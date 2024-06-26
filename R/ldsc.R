ldscore <- function(
    y, x, w, N, M, n_blocks=200
) {

  rlang::check_required(y)
  rlang::check_required(x)
  rlang::check_required(w)
  rlang::check_required(N)
  rlang::check_required(M)
  stopifnot(rlang::is_double(y))
  stopifnot(rlang::is_double(x))
  stopifnot(rlang::is_double(w))
  stopifnot(rlang::is_double(N))
  stopifnot(rlang::is_double(M))
  stopifnot(rlang::is_integerish(n_blocks))


  n_snp <- dim(x)[1]
  n_annot <- dim(x)[2]
  M_tot <- sum(M)
  stopifnot(length(M) == n_annot)
  x_tot <- rowSums(x)

  # provide a starting estimate of heritability
  hsq <- M_tot * (mean(y) - 1) / mean((x_tot * N))

  # first update of weights
  initial_w <- get_weights(ld = x_tot, w_ld = w, N = N, M = M_tot, hsq = hsq)

  # Normalise by the mean of N and add intercept
  Nbar <- mean(N)
  x <- (N*x) / Nbar
  x <- cbind(1, x)

  # the two-step estimator is only implemented for univariate LDscore regression
  if(n_annot == 1) {

    univariate_ldsc(
      y = y,
      x = x,
      x_tot = x_tot,
      w = w,
      N = N,
      initial_w = initial_w,
      M = M,
      Nbar = Nbar
    )

  } else if(n_annot > 1) {

    multivariate_ldsc(
      y = y,
      x = x,
      initial_w = initial_w,
      Nbar = Nbar,
      M = M,
      n_blocks = n_blocks
    )
  }


}

univariate_ldsc <- function(y, x, x_tot,w, N, initial_w, M, Nbar,twostep=30) {


  yp <- y
  ii <- y < twostep
  x1 <- x[ii, ]
  n1 <- sum(ii)
  yp1 <- yp[ii]
  w1 <- w[ii]
  N1 <- N[ii]
  initial_w1 <- initial_w[ii]

  # run two rounds of iterated weighted linear regression to improve weight
  update_func1 <- function(x) update_func(x, ld_tot = x1, w_ld = w1, N = N1, M = M, Nbar = Nbar)
  tmp <- iterated_weights(x = x1, y = yp1, w = initial_w1, update_function = update_func1)
  step1_res <- lstq_jackknife(tmp$x, tmp$y)
  step1_int <- step1_res$full_est[1]

  # estimate h2, with pre-estimated intercept --------------------------------

  updated_s <- update_separators(step1_res$separators, ii)
  yp <- y - step1_int
  x <- x[, 2]
  c <- sum(initial_w * x) / sum(initial_w * x^2)

  update_function2 <- function(x) update_func(x, x_tot, w, N, M, Nbar, step1_int)
  tmp <- iterated_weights(x = x, y = yp, w = initial_w, update_function = update_function2)
  step2_res <- lstq_jackknife(as.matrix(tmp$x), tmp$y, separators = updated_s)
  jknife <- combine_twostep_jknives(step1_res, step2_res, M, c, Nbar)

  res <- extract_jackknife(jknife, Nbar, M)

# done --------------------------------------------------------------------


  res




}

multivariate_ldsc <- function(y, x, initial_w, Nbar, M, n_blocks = 200) {
  initial_w = sqrt(initial_w)
  initial_w = initial_w / sum(initial_w)
  x_weighted <- x*initial_w
  y_weighted <- y*initial_w
  jackknife = lstq_jackknife(x = x_weighted, y = y_weighted, n_blocks = n_blocks)

  extract_jackknife(jackknife, M = M, Nbar = Nbar)

}

combine_twostep_jknives <- function(step1_res, step2_res, M, c, Nbar) {

  est <- c(step1_res$full_est[1], step2_res$full_est)
  delete_values <- matrix(nrow=200, ncol=2)

  delete_values[, 1] <- step1_res$delete_values[,1]
  delete_values[, 2] <- step2_res$delete_values - c*(step1_res$delete_values[,1] - est[1])

  pseudovalues <- delete_values_to_pseudovalues(delete_values, est)
  tmp <- jackknife(pseudovalues)
  list(
    full_est = est,
    cov = tmp$cov,
    var = tmp$var,
    se = tmp$se,
    est = tmp$est,
    delete_values = delete_values
  )

}

update_func <- function(wls, ld_tot, w_ld, N, M, Nbar, intercept = NULL) {
  # ld_tot should be sum of ldscores if partitioned by annotations
  # wls: first entry is intercept, second entry is coef

  if(is.null(intercept)) {

    hsq = M * wls[2] / Nbar
    intercept = max(wls[1])
    ld <- drop(ld_tot[, 2]) # drop intercept from x

  } else {

    hsq = M * wls[1] / Nbar
    ld <- ld_tot

  }
  get_weights(ld, w_ld, N, M, hsq, intercept)

}




# -------------------------------------------------------------------------

get_weights <- function(ld, w_ld, N, M, hsq, intercept=1) {
  stopifnot(!"matrix" %in% class(ld))
  stopifnot(!"w_ld" %in% class(ld))
  stopifnot(is.numeric(N))
  stopifnot(is.numeric(M))
  stopifnot(is.numeric(hsq))

  # hsq is bounded between 0 and 1
  hsq = max(hsq, 0.0)
  hsq = min(hsq, 1.0)
  ld = pmax(ld, 1)
  w_ld = pmax(w_ld, 1)

  c = hsq * N / M
  het_w = 1.0 / (2*  ((intercept + (c*ld))^2) )
  oc_w = 1.0 / w_ld
  w = het_w * oc_w
  w


}

wls <- function(x,y,w) {

  # normalize weights, and weight x and y by w
  w_norm <- (w / sum(w))
  x = w_norm * x
  y = w_norm * y
  XtX = crossprod(x, x)
  Xy = crossprod(x, y)
  res <- solve(XtX, Xy)

  # format the results to return a vector
  drop(t(res))

}

iterated_weights <- function(x,y,w, update_function) {

  w <- sqrt(w)
  for(i in 1:2) {
    ests <- wls(x,y,w)
    updated <- update_function(ests)
    new_w <- sqrt(updated)
    w <- new_w
  }

  w_norm <- (w / sum(w))
  list(
    x = w_norm * x,
    y = w_norm * y
  )

}

update_separators <- function(s, ii) {
  maplist <- seq(1:length(ii))[ii]
  end <- length(s)-1
  update <- s[2:end]
  out <- vector("numeric", length(update))
  for(i in 1:length(update)) {
    out[i] <- maplist[update[i]]
  }

  c(0, out, length(ii))

}

