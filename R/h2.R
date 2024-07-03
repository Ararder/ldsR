utils::globalVariables(c("annot", "m50", "SNP", "coef", "coef_se", "z", "L2"))

#' Estimate SNP heritability using LDscore regression
#'
#' @description
#' An R implementation of the LD score regression method to estimate SNP heritability,
#' which should produce identical results to ldsc.py --h2 command.
#' The package uses by default the eur_w_ld scores from the LDSC tutorial.
#'
#' The function does not apply any quality control checks, and assumes everything is in order.
#'
#'
#' @param sumstat A data.frame or tbl with columns `SNP`, `Z` and `N`
#' @param weights Optional, a data.frame or tbl with columns `SNP`, `L2`
#' @param M Optional, the number of SNPs in the reference panel
#'
#' @return a [dplyr::tibble()] with columns `h2` and `h2_se`
#' @export
#'
#' @examples \dontrun{
#' ldsc_h2(my_gwas)
#' }
ldsc_h2 <- function(sumstat, weights=NULL, M=NULL) {
  req_cols <- c("SNP", "Z", "N")
  stopifnot("sumstat has to be a data.frame or tbl" = "data.frame" %in% class(sumstat))
  stopifnot("SNP, Z and N are required in `sumstat`" = all(req_cols %in% colnames(sumstat)))

  if(is.null(weights)) {
    weights <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"), col_select = c("SNP", "L2"))
    M <- 1173569
  } else {
    weights <- weights
    stopifnot(!is.null(M))
  }
  sumstat <- dplyr::select(sumstat, dplyr::all_of(req_cols))
  m <- dplyr::inner_join(weights, sumstat, by = "SNP")


  results <- ldscore(y = m$Z^2, x = as.matrix(m$L2), w = m$L2, N = m$N, M = M)

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