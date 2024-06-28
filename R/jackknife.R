lstq_jackknife <- function(x, y, n_blocks = 200,separators = NULL) {

  stopifnot("matrix" %in% class(x))
  stopifnot(is.numeric(y))
  rlang::check_required(x)
  rlang::check_required(y)
  n_snps <- dim(x)[1]
  n_annot <- dim(x)[2]
  stopifnot(length(y) == n_snps)

  if(is.null(separators)) {
    separators <- floor(seq(0, n_snps, length.out = n_blocks + 1))
  }


  blocks <- get_blocks(n_snps, n_blocks, separators)
  XtX <- purrr::map(blocks, \(b) t(x[b,]) %*% x[b,])
  Xy <- purrr::map(blocks, \(b) crossprod(x[b,], y[b]))
  full_XtX <- purrr::reduce(XtX, `+`)
  full_Xy <- purrr::reduce(Xy, `+`)
  est <- t(solve(full_XtX, full_Xy))


  # iterate over the block_values, removing XtX.block from XtX.tot: equivalent to leave-a-block-out jackknife
  delete_values <-
    purrr::map2(XtX, Xy, \(XtX_sub, Xy_sub) {
      XtX <-  full_XtX - XtX_sub
      Xy <-  full_Xy - Xy_sub
      solve(XtX, Xy)
    })


  # clean up formatting and extract estimates
  delete_values <- purrr::reduce(delete_values, cbind) |> t()
  pseudo <- delete_values_to_pseudovalues(delete_values, est)
  res <- jackknife(pseudo)

  # return both pseudo and delete values
  res$pseudo <- pseudo
  res$delete_values <- delete_values
  res$est
  res$full_est <- est
  res$separators <- separators
  res



}


get_blocks <- function(n_snps, n_blocks, s) {
  # shift the separators to the right by 1, python is 0 based, R is 1 based
  s <- s + 1

  blocks <- vector("list", n_blocks)
  for (i in 1:n_blocks) {
    start <- s[i]
    end <- s[i+1] - 1 # python ranges are NOT inclusive, R's is.
    blocks[[i]] <- start:end
  }

  blocks
}


delete_values_to_pseudovalues <- function(delete_values, est) {
  stopifnot("matrix" %in% class(delete_values))
  n_blocks <- nrow(delete_values)
  p <- ncol(delete_values)
  stopifnot(length(est) == p)

  # a little bit of extra code to match dimensions
  est_matrix <- matrix(ncol = p, nrow = n_blocks)
  for(i in 1:n_blocks) est_matrix[i, ] <- est

  (n_blocks * est_matrix) - (n_blocks - 1) * delete_values


}




jackknife <- function(pseudovalues) {
  n_blocks <- nrow(pseudovalues)

  cov <- cov(pseudovalues) / n_blocks
  var <- diag(cov)
  se <- sqrt(var)
  est <- colMeans(pseudovalues)

  list(
    cov = cov,
    var = var,
    se = se,
    est = est
  )


}

extract_jackknife <- function(jknife, Nbar, M) {

  coef <-  jknife$full_est[-1] / Nbar
  coef_cov <- jknife$cov[-1,-1] / Nbar^2
  if("matrix" %in% class(coef_cov)) {
    coef_se <- sqrt(diag(coef_cov))
  } else {
    coef_se <- sqrt(coef_cov)
  }


  cat <- M * coef
  cat_cov <- coef_cov * (M %*% t(M))
  cat_se = sqrt(diag(cat_cov))


  tot <- sum(cat)
  tot_cov <-  sum(cat_cov)
  tot_se <-  sqrt(tot_cov)

  M_prop = M / sum(M)
  enrichment = (cat / M) / ( tot / sum(M))



  list(
    coef = coef,
    coef_se = coef_se,
    cat = cat,
    cat_cov = cat_cov,
    cat_se = cat_se,
    tot = tot,
    tot_se = tot_se,
    tot_cov = tot_cov,
    enrichment = enrichment,
    M_prop = M_prop
  )

}


# -------------------------------------------------------------------------
# efforts to memoise cell-type analysis
# compute_M_T_M <- function(y, x, w, N, M, n_blocks=200) {
#   rlang::check_required(y)
#   rlang::check_required(x)
#   rlang::check_required(w)
#   rlang::check_required(N)
#   rlang::check_required(M)
#   stopifnot(rlang::is_double(y))
#   stopifnot(rlang::is_double(x))
#   stopifnot(rlang::is_double(w))
#   stopifnot(rlang::is_double(N))
#   stopifnot(rlang::is_double(M))
#   stopifnot(rlang::is_integerish(n_blocks))
#
#
#   n_snps <- dim(x)[1]
#   n_annot <- dim(x)[2]
#   M_tot <- sum(M)
#   stopifnot(length(M) == n_annot)
#   x_tot <- rowSums(x)
#
#   # provide a starting estimate of heritability
#   hsq <- M_tot * (mean(y) - 1) / mean((x_tot * N))
#   initial_w <- get_weights(ld = x_tot, w_ld = w, N = N, M = M_tot, hsq = hsq)
#   Nbar <- mean(N)
#   x <- (N*x) / Nbar
#   x <- cbind(1, x)
#   initial_w = sqrt(initial_w)
#   initial_w = initial_w / sum(initial_w)
#
#   x_weighted <- x*initial_w
#   y_weighted <- y*initial_w
#   separators <- floor(seq(0, n_snps, length.out = n_blocks + 1))
#   blocks <- get_blocks(n_snps, n_blocks, separators)
#   compute_xtx(x, blocks)
#
#
# }


# compute_xtx <- function(x, blocks, M_T_M = NULL) {
#   if(is.null(M_T_M)) {
#     purrr::map(blocks, \(b) t(x[b,]) %*% x[b,])
#   } else {
#     j <- x[,ncol(x)]
#     M <- x[, -ncol(x)]
#     stopifnot(length(blocks) == length(M_T_M))
#     purrr::map2(blocks, M_T_M, \(b, m) memoised_xtx(M = M[b,], j = j[b], M_T_M = m))
#   }
# }
#
# memoised_xtx <- function(M, j, M_T_M) {
#
#   MTj <- t(M) %*% j
#   jTj <- t(j) %*% j
#
#   rbind(
#     cbind(M_T_M, MTj),
#     cbind(t(MTj), jTj)
#   )
#
# }


