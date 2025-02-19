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

prop <- function(jknife,cat,tot, Nbar, M) {
  # compute proportional estimates
  n_annot <- length(jknife$est) -1
  n_blocks <- nrow(jknife$delete_values)

  # no intercept
  numer_delete_vals <- jknife$delete_values[,-1]
  for(i in 1:nrow(numer_delete_vals)) numer_delete_vals[i, ] <- (numer_delete_vals[i, ] * M) / Nbar
  denom_delete_vals <- matrix(rowSums(numer_delete_vals))


  denom_delete_vals <- matrix(rep(denom_delete_vals, 53), nrow = 200, ncol = 53, byrow = FALSE)
  ratio <- cat / tot

  pseudovalues <- matrix(nrow = n_blocks, ncol = n_annot)
  for(j in 1:n_blocks) {
    first <- n_blocks * ratio
    second <- (n_blocks - 1) * numer_delete_vals[j,] / denom_delete_vals[j,]
    pseudovalues[j,] <-  first - second
  }

  jackknife <- jackknife(pseudovalues)

  list(
    prop = ratio,
    prop_cov = jackknife$cov,
    prop_se = jackknife$se
  )

}

overlapping_annotations <- function(ldscore_dirs, M, jknife) {
  vals <- read_overlap_matrix(ldscore_dirs)
  overlap_matrix <- vals[["overlap_matrix"]]
  M_tot <- vals[["M_tot"]]
  overlap_matrix_prop <- matrix(nrow= nrow(overlap_matrix), ncol = ncol(overlap_matrix))
  for(i in 1:nrow(overlap_matrix)) {
    overlap_matrix_prop[i,] <- overlap_matrix[i,] / M
  }



  prop_hsq_overlap <- overlap_matrix_prop %*% matrix(jknife$prop$prop, ncol = 1)
  step1 = overlap_matrix_prop %*% jknife$prop$prop_cov
  step2 = step1 %*% t(overlap_matrix_prop)
  prop_hsq_overlap_var <- diag(step2)

  prop_hsq_overlap_se <- sqrt(pmax(prop_hsq_overlap_var,0))
  prop_M_overlap <- M / M_tot
  enrichment = prop_hsq_overlap / prop_M_overlap
  enrichment_se = prop_hsq_overlap_se / prop_M_overlap



  dplyr::tibble(
    prop_snps = prop_M_overlap,
    prop_h2 = drop(prop_hsq_overlap),
    prop_h2_std_error = prop_hsq_overlap_se,
    enrich = drop(enrichment),
    enrich_se = enrichment_se,
  )

}
