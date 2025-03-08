
#' Perform cell-type analysis using partitioned heritability
#'
#' @description
#' A common usage of partitioned heritability is to estimate enrichment and significance for
#' an annotation typical for a specific cell-type.
#' For this analysis, it is useful to adjust for a baseline set of annotations, and then estimate the association for a large set of annotations ("cell-types")
#'
#'
#' @inheritParams ldsc_h2
#' @param covariate_dir a directory containing the files `ld.parquet` and `annot.parquet`
#' @param ldscore_dir a directory containing the files `ld.parquet` and `annot.parquet`. Each annotation will be executed, while adjusting for all annotations in covariate_dir
#'
#' @return a [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#' sumstat1 <- arrow::read_parquet(system.file("extdata", "sumstats.parquet", package = "ldsR"))
#' sumstat1 <- dplyr::select(sumstat1, SNP, Z = Z.x, N=N.x)
#'
#'
#' celltype_analysis(
#'  sumstat1,
#'  covariate_dir = system.file("extdata", "baseline1.1_test", package = "ldsR"),
#'  ldscore_dir = system.file("extdata", "superclusters", package = "ldsR")
#' )
#' }
#'
celltype_analysis <- function(sumstat, covariate_dir, ldscore_dir, weights = NULL) {
  req_cols <- c("SNP", "Z", "N")
  stopifnot("SNP, Z and N are required in `sumstat`" = all(req_cols %in% colnames(sumstat)))
  sumstat <- dplyr::select(sumstat, dplyr::all_of(req_cols))


  if(is.null(weights)) {
    weights <- arrow::read_parquet(system.file("extdata/1000G_Phase3_weights_hm3_no_MHC.parquet", package = "ldsR"))
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
