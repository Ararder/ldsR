gencov <- function(z1,z2, w,x,N1,N2, M, hsq1_tot, hsq2_tot, intercept_hsq1,
                   intercept_hsq2, n_blocks=200, twostep=30) {






  # provide a starting estimate of rho_g
  y <- z1*z2
  N <- sqrt(N1 * N2)
  rho_g <- M * (mean(y) - 0) / mean((x * N))
  x_tot <- x

  # get initial weights
  initial_w <- rg_weights(ld = x, w_ld = w, N1 = N1, N2 = N2, M = M, h1 = hsq1_tot, h2 = hsq2_tot, rho_g = rho_g, intercept_hsq1 = intercept_hsq1, intercept_hsq2 = intercept_hsq2)

  # Normalise by the mean of N and add intercept
  Nbar <- mean(N)
  x <- (N*x) / Nbar
  x <- cbind(1, x)
  yp <- y


  # (1) apply twostep filter for estimation of intercept ---------------------
  ii <- z1^2 < twostep & z2^2 < twostep

  x1 <- x[ii, ]
  yp1 <- yp[ii]
  initial_w1 <- initial_w[ii]

  update_func1 <- function(x) rg_update_func(
    wls = x,
    ld = x1,
    w_ld = w[ii],
    N1 = N1[ii],
    N2 = N2[ii],
    M = M,
    hsq1 = hsq1_tot,
    hsq2 = hsq2_tot,
    intercept_hsq1 = intercept_hsq1,
    intercept_hsq2 = intercept_hsq2,
    Nbar = Nbar
    )


  tmp <- iterated_weights(x = x1, y = yp1, w = initial_w1, update_function = update_func1)
  step1_res <- lstq_jackknife(tmp$x, tmp$y)
  step1_int <- step1_res$full_est[1]

  # (2) estimate h2, with pre-estimated intercept ----------------------------

  updated_s <- update_separators(step1_res$separators, ii)
  yp <- y - step1_int
  x <- x[, 2]
  c <- sum(initial_w * x) / sum(initial_w * x^2)

  update_function2 <- function(wls) { rg_update_func(
    wls = wls,
    ld = x_tot,
    w_ld = w,
    N1 = N1,
    N2 = N2,
    M = M,
    intercept = step1_int,
    hsq1 = hsq1_tot,
    hsq2 = hsq2_tot,
    intercept_hsq1 = intercept_hsq1,
    intercept_hsq2 = intercept_hsq2,
    Nbar = Nbar)
  }


  tmp <- iterated_weights(x = x, y = yp, w = initial_w, update_function = update_function2)
  step2_res <- lstq_jackknife(as.matrix(tmp$x), tmp$y, separators = updated_s)


  # (3) extract and prepare results ------------------------------------------


  jknife <- combine_twostep_jknives(step1_res, step2_res, M, c, Nbar)
  results <- extract_jackknife(jknife, Nbar, M)

  results[["int"]] <- step1_int
  results[["int_se"]] <- step1_res$se[1]
  results[["delete_values"]] <- jknife$delete_values
  results


}
tot_delete_values <- function(jknife, M, Nbar) {
  jknife$delete_values[,2] * M / Nbar

}
rg_weights <- function(ld, w_ld, N1, N2, M, h1, h2, rho_g, intercept_gencov=0, intercept_hsq1, intercept_hsq2) {

  h1 = max(h1, 0) |> min(1)
  h2 = max(h2,0) |> min(1)
  rho_g <- min(rho_g, 1) |> max(-1)
  ld <- pmax(ld, 1)
  w_ld <- pmax(w_ld, 1)

  a <- (N1 * h1 * ld) / M + intercept_hsq1
  b <- (N2 * h2 * ld) / M + intercept_hsq2
  sqrt_n1n2 <- sqrt(N1 * N2)
  c = (sqrt_n1n2 * rho_g * ld) / M + intercept_gencov
  oc_w <- 1 / w_ld
  het_w = 1 / (a*b + c^2)

  # return weights
  het_w * oc_w


}


rg_update_func <-  function(
    wls, ld, w_ld, N1, N2, M, Nbar, intercept=NULL,
    hsq1, hsq2, intercept_hsq1, intercept_hsq2
) {



  if(is.null(intercept)) {
    stopifnot("matrix" %in% class(ld))
    rho_g <-  M * wls[2] / Nbar
    intercept <- wls[1]
    ld <- drop(ld[, 2]) # drop intercept from x
  } else {
    rho_g <-  M * wls[1] / Nbar

  }

  rg_weights(
    ld = ld,
    w_ld = w_ld,
    N1 = N1,
    N2 = N2,
    M = M,
    h1 = hsq1,
    h2 = hsq2,
    rho_g = rho_g,
    intercept_gencov = intercept,
    intercept_hsq1 = intercept_hsq1,
    intercept_hsq2 = intercept_hsq2
  )


}
