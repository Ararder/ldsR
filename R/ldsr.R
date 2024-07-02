utils::globalVariables(c("annot", "m50", "SNP", "coef", "coef_se", "z", "L2"))
req_cols <- c("SNP", "Z", "N")

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




# GENETIC CORRELATIONS ----------------------------------------------------






#' Compute the genetic correlation between two traits
#'
#' @param sumstats1 a [dplyr::tibble()] with columns h2 and h2_se
#' @param sumstats2 a [dplyr::tibble()] with columns h2 and h2_se
#' @param weights Optional, a data.frame or tbl with columns SNP and L2
#' @param M number of snps used to calculate LD in weights
#' @param n_blocks blocks for jackknife
#'
#' @return a [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#' ldsc_rg(s1, s2)
#' }
ldsc_rg <- function(sumstats1, sumstats2, weights=NULL, M=NULL, n_blocks=200) {
  stopifnot("SNP, Z and N are required in `sumstats1`" = all(req_cols %in% colnames(sumstats1)))
  stopifnot("SNP, Z and N are required in `sumstats2`" = all(req_cols %in% colnames(sumstats2)))

  weights <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"))
  M <- 1173569
  sumstats1 <- dplyr::select(sumstats1, dplyr::all_of(req_cols))
  sumstats2 <- dplyr::select(sumstats2, dplyr::all_of(req_cols))

  m <- dplyr::inner_join(weights, sumstats1, by = "SNP") |>
    dplyr::inner_join(sumstats2, by = "SNP")


  # -------------------------------------------------------------------------

  x <- as.matrix(m$L2)
  res1 <- ldscore(y = m$Z.x^2, x = x, w = m$L2, N = m$N.x, M = M)
  res2 <- ldscore(y = m$Z.y^2, x = x, w = m$L2, N = m$N.y, M = M)
  N <- sqrt(m$N.x * m$N.y)
  res3 <- gencov(
    z1 = m$Z.x,z2 = m$Z.y,w = m$L2, x = x, N1 = m$N.x, N2 = m$N.y, M = M,
    hsq1_tot = res1$tot, hsq2_tot = res2$tot, intercept_hsq1 = res1$int,
    intercept_hsq2 = res2$int
  )



  rg_ratio <- res3$tot / sqrt(res1$tot * res2$tot)
  numer <- tot_delete_values(res3, M, mean(N))

  d1_tot <- tot_delete_values(res1, M, mean(m$N.x))
  d2_tot <- tot_delete_values(res2, M, mean(m$N.y))
  denom <- sqrt(d1_tot*d2_tot)
  pseudovalues <- vector("numeric", n_blocks)
  for(j in 1:n_blocks){
    pseudovalues[j] <- n_blocks * rg_ratio - (n_blocks - 1) * numer[j] / denom[j]
  }

  final <- jackknife(matrix(pseudovalues))
  dplyr::tibble(rg = final$est, rg_se = final$se, h2_s1 = res1$tot, h2_1_se = res1$tot_se, h2_s2 = res2$tot, h2_2_se = res2$tot_se)


}




# CELLTYPE H2 -------------------------------------------------------------




#' Run partitioned heritability across many annotations
#'
#' @param sumstat A data.frame or tbl with columns `SNP`, `Z` and `N`
#' @param covariate_dir a directory containing the files `ld.parquet` and `annot.parquet`
#' @param celltype_dir a directory containing the files `ld.parquet` and `annot.parquet`
#' @param weights A numeric vector of regression weights. Comes precomputed within the ldsR package for european LDscores
#'
#' @return a [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#' celltype_analysis(ss_tbl, "/path/to_baseline", "path/to/celltype_ldscores")
#' }
celltype_analysis <- function(sumstat, covariate_dir, celltype_dir, weights = NULL) {

  stopifnot("SNP, Z and N are required in `sumstat`" = all(req_cols %in% colnames(sumstat)))
  sumstat <- dplyr::select(sumstat, dplyr::all_of(req_cols))


  if(is.null(weights)) {
    weights <- arrow::read_parquet(system.file("extdata/eur_w_ld.parquet", package = "ldsR"), col_select = c("SNP", "L2_celltype")) |>
      dplyr::rename(L2 = "L2_celltype") |>
      dplyr::filter(!is.na(L2))
  }

  covars <- parse_parquet_dir(covariate_dir)
  covar_ld <- covars[["ld"]]
  covar_M <- covars[["annot"]][["m50"]]

  celltypes <- parse_parquet_dir(celltype_dir)
  celltypes_ld <- celltypes[["ld"]]
  celltypes_M <- celltypes[["annot"]][["m50"]]
  rm(covars, celltypes)


  #### Merge sumstats with LDscores
  n_snps_before <- nrow(sumstat)
  merged <- dplyr::inner_join(weights, sumstat, by = "SNP") |>
    dplyr::inner_join(covar_ld, by = "SNP") |>
    dplyr::inner_join(celltypes_ld, by = "SNP")

  diff <- nrow(sumstat) - n_snps_before
  cli::cli_alert_info("Removed {.bold {diff}} rows after merging with ldscores")


  ##
  non_ldscores <- unique(c(colnames(sumstat),colnames(weights)))
  celltype_ldscore_names <-  colnames(dplyr::select(celltypes_ld, -c("SNP")))



  purrr::imap(celltype_ldscore_names, \(celltype, i)
             .celltype(
               celltype = celltype,
               merged = merged,
               M = c(covar_M, celltypes_M[i]),
               non_ldscores = non_ldscores,
               celltype_ldscore_names = celltype_ldscore_names
               ),
              .progress = list(
                "name" = "Computing celltype associations...",
                "type" = "tasks"
              )) |>
    purrr::list_rbind()


}


.celltype <- function(celltype, merged, M, non_ldscores, celltype_ldscore_names) {

    cols_to_remove <- unique(c(non_ldscores, celltype_ldscore_names[!celltype_ldscore_names %in% celltype]))
    x <- dplyr::select(merged,-dplyr::all_of(cols_to_remove)) |> as.matrix()
    stopifnot("wrong dimensions between celltype ldscores and M"= length(M) == dim(x)[2])


    res <- ldscore(
      x = x,
      y = merged$Z^2,
      w = merged$L2,
      N = merged$N,
      M = M
    )

    tidy_results(res) |>
      dplyr::filter(annot == celltype)

}





