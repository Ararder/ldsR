utils::globalVariables(c("annot", "m50", "SNP", "coef", "coef_se", "z", "L2"))

#' Estimate SNP heritability using LDscore regression
#'
#' @description
#' An R implementation of the LD score regression method to estimate SNP heritability,
#' mimicking `ldsc --h2` from the [ldsc](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation) package. 
#' 
#' LDscores for the European subset of 1000g (2015 release) for 1,290,028 HapMap3 SNPs (with 1,173,569 SNPs with freq > 5%)
#' are bundled within the ldsR package and are used by default. This corresponds to the `eur_w_ld_chr` folder 
#' previously shared at the LDSC [github](https://github.com/bulik/ldsc). 
#' 
#' You can inspect the LDscores used by default: 
#' `arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"))`
#'
#' ldsc_h2 does not perform any quality control on the input summary statistics, except to merge with the
#' 1,290,028 HapMap3 SNPs in the reference panel. See the [munge()] function to mimic the `munge_sumstats.py` function. 
#'
#' 
#'
#'
#' @param sumstat A data.frame or tbl with columns `SNP`, `Z` and `N`
#' @param weights Optional, a data.frame or tbl with columns `SNP`, `L2`
#' @param M Optional, the number of SNPs in the reference panel
#' @param n_blocks Number of blocks to use for the jackknife estimator
#'
#' @return a [dplyr::tibble()] with columns `h2` and `h2_se`
#' @export
#'
#' @examples 
#' p <- system.file("extdata", "eur_w_ld.parquet", package = "ldsR")
#' snps <- arrow::read_parquet(p, col_select = c("SNP"))
#' snps$N <- 130000
#' snps$Z <- rnorm(nrow(snps))
#' ldsc_h2(snps)
#' 
#' 
ldsc_h2 <- function(sumstat, weights=NULL, M=NULL, n_blocks = 200) {
  req_cols <- c("SNP", "Z", "N")
  stopifnot("sumstat has to be a data.frame or tbl" = "data.frame" %in% class(sumstat))
  check_columns(req_cols, sumstat)

  if(is.null(weights)) {
    weights <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"), col_select = c("SNP", "L2"))
    M <- 1173569
  } else {
    stopifnot("`weights` must be a data.frane with columns `SNP` and `L2`" = "data.frame" %in% class(weights))
    check_columns(c("SNP", "L2"), weights)
    stopifnot("To use custom weights, you must also pass `M`" = !is.null(M))
  }


  # merge and run ldscore regression ---------------------------------------

  m <- dplyr::inner_join(weights, dplyr::select(sumstat, dplyr::all_of(req_cols)), by = "SNP")

  results <- ldscore(y = m$Z^2, x = as.matrix(m$L2), w = m$L2, N = m$N, M = M, n_blocks=n_blocks)
  mean_chi2 <- mean(m$Z^2)

  if(mean_chi2 > 1) {
    ratio_se = results$int_se / (mean_chi2 - 1)
    ratio = (results$int - 1) / (mean_chi2 - 1)

  } else {
    ratio_se <- NA_real_
    ratio <- NA_real_
  }

  # -------------------------------------------------------------------------


  dplyr::tibble(
    h2 = results$tot,
    h2_se = results$tot_se,
    int = results$int,
    int_se = results$int_se,
    mean_chi2 = mean_chi2,
    lambda_gc = stats::median(m$Z^2) / 0.4549,

    )


}