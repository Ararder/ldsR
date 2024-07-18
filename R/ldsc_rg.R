
utils::globalVariables(c("annot", "m50", "SNP", "coef", "coef_se", "z", "L2", "rg", "rg_se",
"A1.x","A1.y", "A2.x", "A2.y", "Z.y"))
#' Compute the genetic correlation between two traits
#'
#' @param sumstats1 a [dplyr::tibble()] with atleast the columns SNP, A1, A2, Z, N
#' @param sumstats2 a [dplyr::tibble()] with atleast the columns SNP, A1, A2, Z, N.
#' A list of summary statistics can be passed to estimate rG between sumstats1 and all traits in sumstats2
#' @param weights NOT YET IMPLEMENTED - supply custom weights vector
#' @param M NOT YET IMPLEMENTED - number of SNPs used to calculate weights
#' @param n_blocks Number of blocks to use for jackknife
#'
#' @return a [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#' ldsc_rg(s1, s2)
#' }
ldsc_rg <- function(sumstats1, sumstats2, weights=NULL, M=NULL, n_blocks=200) {
  req_cols <- c("SNP", "Z", "N", "A1", "A2")
  check_columns(req_cols, sumstats1)
  sumstats1 <- dplyr::select(sumstats1, dplyr::all_of(req_cols))
  weights <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"), col_select = c("SNP", "L2"))
  M <- 1173569
  
  
  
  # -------------------------------------------------------------------------
  if("data.frame" %in% class(sumstats2)) sumstats2 <- list(sumstats2)
  

  purrr::map(sumstats2, \(x) .rg(sumstats1, x, M=M, weights=weights, n_blocks = n_blocks),
    .progress = list(type = "tasks", name = "Computing genetic correlations...")
  ) |> 
    purrr::list_rbind()


}



.rg <- function(sumstats1, sumstats2, M,weights, n_blocks = 200,trait1 = NULL, trait2=NULL) {
  req_cols <- c("SNP", "Z", "N", "A1", "A2")

  # check that all columns exist for sumstats2
  sumstats2 <- dplyr::select(sumstats2, dplyr::all_of(req_cols))
  check_columns(req_cols, sumstats2)
  
  before <- nrow(sumstats1)
  m <- dplyr::inner_join(
      weights,
      merge_sumstats(sumstats1, sumstats2),
      by = "SNP"
    )
  cli::cli_alert_warning("{before - nrow(m)} SNPs were removed when merging summary statistics")

  N <- sqrt(m$N.x * m$N.y)
  x <- as.matrix(m$L2)
  res1 <- ldscore(y = m$Z.x^2, x = x, w = m$L2, N = m$N.x, M = M)
  res2 <- ldscore(y = m$Z.y^2, x = x, w = m$L2, N = m$N.y, M = M)
  res3 <- gencov(
    z1 = m$Z.x,z2 = m$Z.y,w = m$L2, x = x, N1 = m$N.x, N2 = m$N.y, M = M,
    hsq1_tot = res1$tot, hsq2_tot = res2$tot, intercept_hsq1 = res1$int,
    intercept_hsq2 = res2$int
  )


  gencov_int <- res3$int
  gencov_int_se <- res3$int_se
  rg_ratio <- res3$tot / sqrt(res1$tot * res2$tot)

  numer <- tot_delete_values(res3, M, mean(N))
  denom <- sqrt(tot_delete_values(res1, M, mean(m$N.x)) * tot_delete_values(res2, M, mean(m$N.y)))

  pseudovalues <- vector("numeric", n_blocks)
  for(j in 1:n_blocks) pseudovalues[j] <- n_blocks * rg_ratio - (n_blocks - 1) * numer[j] / denom[j]
  final <- jackknife(matrix(pseudovalues))


  # report main results ----------------------------------------------------


  out <- dplyr::tibble(
    rg = rg_ratio, 
    rg_se = final$se, 
    h2_trait1 = res1$tot, 
    h2_trait_se = res1$tot_se, 
    h2_trait2 = res2$tot, 
    h2_trait2_se = res2$tot_se,
    gcov = res3$tot,
    gcov_se = res3$tot_se,
    gcov_int = res3$int,
    gcov_int_se = res3$int_se,
    mean_z1z2 = mean(m$Z.x*m$Z.y),
    z = rg /rg_se,
    p = stats::pchisq(z, df = 1, lower.tail = FALSE)
  )
  out$trait1 <- trait1
  out$trait2 <- trait2

  out
}

merge_sumstats <- function(sumstats1, sumstats2) {
  
  m <- dplyr::inner_join(sumstats1, sumstats2, by = "SNP")
  
  
  dplyr::filter(m, A1.x == A1.y & A2.x == A2.y | A1.x == A2.y & A2.x == A1.y) |> 
    # flip Z if alleles are flipped
    dplyr::mutate(Z.y = dplyr::if_else(A1.x == A1.y, Z.y, Z.y*-1)) |> 
    dplyr::select(dplyr::all_of(c("SNP", "Z.x", "Z.y", "N.x", "N.y")))
  
  

  # 
}

