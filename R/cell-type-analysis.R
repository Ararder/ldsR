utils::globalVariables(c("annot", "m50", "SNP", "coef", "coef_se", "z", "L2"))
req_cols <- c("SNP", "Z", "N")
#' Run partitioned heritability across many annotations
#'
#' @inheritParams ldsc_h2
#' @param covariate_dir a directory containing the files `ld.parquet` and `annot.parquet`
#' @param celltype_dir a directory containing the files `ld.parquet` and `annot.parquet`
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

tidy_results <- function(res) {
  res
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
    dplyr::arrange(dplyr::desc(abs(z)))

}
