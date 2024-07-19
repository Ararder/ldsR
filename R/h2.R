utils::globalVariables(c("annot", "m50", "SNP", "coef", "coef_se", "z", "L2"))

#' Estimate SNP heritability using LDscore regression for a single annotation
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
#' @param sumstat A [dplyr::tibble()] with columns `SNP`, `Z` and `N`
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
#' @examples \dontrun{
#' partitioned_heritability(sumstats, "path/to_dir/")
#' }
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












#' Perform cell-type analysis using partitioned heritability
#' 
#' @description
#' A common usage of partitioned heritability is to estimate enrichment and significance for
#' an annotation typical for a specific cell-type.
#' For this analysis, it is useful to adjust for a baseline set of annotations, and then estimate the association for a large set of annotations ("cell-types")
#' 
#'
#' @inheritParams partitioned_heritability
#' @param covariate_dir a directory containing the files `ld.parquet` and `annot.parquet`
#'
#' @return a [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#' celltype_analysis(ss_tbl, "/path/to_baseline", "path/to/celltype_ldscores")
#' }
celltype_analysis <- function(sumstat, covariate_dir, ldscore_dir, weights = NULL) {
  req_cols <- c("SNP", "Z", "N")
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

  celltypes <- parse_parquet_dir(ldscore_dir)
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
      dplyr::arrange(dplyr::desc(abs(z))) |> 
      dplyr::filter(annot == celltype)

}



