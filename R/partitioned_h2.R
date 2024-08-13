#' Estimate partitioned SNP heritability
#'
#' @description
#' An R implementation of the LD score regression method to estimate SNP heritability,
#' mimicking `ldsc --h2` from the [ldsc](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation) package.
#'
#' For partitioned heritablity, the noMHC weights are used instead of the eur_w_ld_chr weights.
#'
#' You can inspect the LDscores used by default, in the column `L2_celltype`
#' `arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"))`
#'
#' partitioned_heritability does not perform any quality control on the input summary statistics, except to merge with the
#' SNPs in the reference panel. See the [munge()] function to mimic the `munge_sumstats.py` function.
#'
#' @inheritParams ldsc_h2
#' @param ldscore_dir filepath to a directory with the `annot.parquet` and `ldscores.parquet` files
#'
#' @return a [dplyr::tibble()]
#' @export
#'
#' @examples
#' sumstat1 <- arrow::read_parquet(system.file("extdata", "sumstats.parquet", package = "ldsR")) |>
#' dplyr::select(SNP, Z = Z.x, N=N.x)
#' ldscore_dir <-  system.file("extdata", "baseline1.1_test", package = "ldsR")
#' fs::dir_tree(ldscore_dir)
#' partitioned_h2(sumstat1, ldscore_dir = ldscore_dir)
#'
#'
#'
#'
partitioned_h2 <- function(sumstat, ldscore_dir, weights = NULL, n_blocks=200) {
  stopifnot("sumstat has to be a data.frame or tbl" = "data.frame" %in% class(sumstat))
  check_columns(c("SNP", "Z", "N"), sumstat)


  if(is.null(weights)) {
    weights <- arrow::read_parquet(system.file("extdata/eur_w_ld.parquet", package = "ldsR"), col_select = c("SNP", "L2_celltype")) |>
      dplyr::rename(L2 = "L2_celltype") |>
      dplyr::filter(!is.na(L2))
  }

  covars <- parse_parquet_dir(ldscore_dir)
  covar_ld <- covars[["ld"]]
  covar_M <- covars[["annot"]][["m50"]]


  n_before <- nrow(sumstat)
  merged <- dplyr::inner_join(weights, sumstat, by = "SNP")
  cli::cli_alert_info("Removed {.bold {n_before - nrow(merged)}} rows after merging with weights")
  n_before <- nrow(merged)
  merged <- dplyr::inner_join(merged, covar_ld, by = "SNP")
  cli::cli_alert_info("Removed {.bold {n_before - nrow(merged)}} rows after merging with ldscores")
  remove_cols <- unique(c(colnames(sumstat), colnames(weights)))
  x <- dplyr::select(merged,-dplyr::all_of(remove_cols)) |> as.matrix()


  # run ldscore regression -------------------------------------------------


  res <- ldscore(y = merged$Z^2, x = x, w = merged$L2, N = merged$N, M = as.double(covar_M), n_blocks=n_blocks)



  dplyr::tibble(
    annot = names(res$coef_se),
    coef = res$coef,
    coef_se = res$coef_se,
    enrich = res$enrichment,
    prop = res$M_prop,
    z = coef/coef_se,
    tot = res$tot,
    tot_se = res$tot_se
  ) |>
    dplyr::arrange(dplyr::desc(z))

}
