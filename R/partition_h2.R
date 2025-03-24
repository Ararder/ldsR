utils::globalVariables(c("common", "enrich"))

#' Flexible partitioning of heritability
#' @description
#' An R implementation of the LD score regression method to estimate SNP heritability, focusing on partioning heritability across multiple annotations.
#'
#' @inheritParams ldsc_h2
#' @param ldscore_dirs a list of directories containing ldscore files
#' @param subset_annots a character vector of annotations to subset the ldscore files. Use [get_annot_names()]
#' @param overlapping_annotations are the annotations overlapping? In such a case, the estimate of the enrichment estimate needs to be adjusted.
#' @returns a tibble with results
#' @export
#'
#' @examples \dontrun{
#' partition_h2(sumstat, ldscore_dirs, subset_annots = c("L2", "L3"))
#' }
partition_h2 <- function(
    sumstat,
    ldscore_dirs,
    subset_annots = NULL,
    overlapping_annotations=FALSE,
    weights = NULL,
    n_blocks = 200
    ) {


  # -------------------------------------------------------------------------
  # basic checks
  rlang::check_required(sumstat)
  check_columns(c("SNP", "Z", "N"), sumstat)
  rlang::check_required(ldscore_dirs)
  stopifnot("sumstat has to be a data.frame or tbl" = "data.frame" %in% class(sumstat))
  purrr::walk(ldscore_dirs, check_is_path)
  stopifnot(is.null(subset_annots) | rlang::is_character(subset_annots))
  stopifnot(rlang::is_scalar_logical(overlapping_annotations))

  if(is.null(weights)) {
    weights <- arrow::read_parquet(system.file("extdata/1000G_Phase3_weights_hm3_no_MHC.parquet", package = "ldsR"))
  }



  # -------------------------------------------------------------------------


  data <- purrr::map(ldscore_dirs, \(x) parse_parquet_dir(x, read_ref = overlapping_annotations, subset_annots = subset_annots)) |>
    purrr::reduce(combine_ld_data, ref = overlapping_annotations, outdir = NULL)

  # Inform users if any annotation was not found
  if(!all(subset_annots %in% data$annot$annot)) {
    cli::cli_alert_warning("Some annotations were not found; {subset_annots[!subset_annots %in% data$annot$annot]}")
  }


  # check that ordering is the same across the data
  data[["annot"]]  <- dplyr::inner_join(dplyr::tibble(annot = colnames(data[["ld"]])[-1]), data[["annot"]], by = "annot")
  stopifnot(all(colnames(data[["ld"]])[-1] == data[["annot"]][["annot"]]))
  if(overlapping_annotations) {
    data[["annot_ref"]] <- dplyr::select(data[["annot_ref"]], dplyr::all_of(colnames(data[["ld"]])))
    stopifnot(all(colnames(data[["annot_ref"]]) == colnames(data[["ld"]])))
  }


  covar_ld <- data[["ld"]]
  covar_M <- data[["annot"]][["m50"]]


  n_before <- nrow(sumstat)
  merged <- dplyr::inner_join(weights, sumstat, by = "SNP")
  cli::cli_alert_info("Removed {.bold {n_before - nrow(merged)}} rows after merging with weights")
  n_before <- nrow(merged)
  merged <- dplyr::inner_join(merged, covar_ld, by = "SNP")
  cli::cli_alert_info("Removed {.bold {n_before - nrow(merged)}} rows after merging with ldscores")
  cli::cli_alert_success("A total of {nrow(merged)} SNPs remain for the regression model")


  remove_cols <- unique(c(colnames(sumstat), colnames(weights)))
  x <- dplyr::select(merged,-dplyr::all_of(remove_cols)) |> as.matrix()


  # run ldscore regression -------------------------------------------------


  res <- ldscore(y = merged$Z^2, x = x, w = merged$L2, N = merged$N, M = as.double(covar_M), n_blocks=n_blocks)

  base_results <- dplyr::tibble(
    annot = names(res$coef_se), coef = res$coef, coef_se = res$coef_se,
    z = coef/coef_se, tot = res$tot, tot_se = res$tot_se, n_snps = covar_M
  ) |>
    dplyr::arrange(dplyr::desc(z))


  if(!overlapping_annotations) {

    return(base_results)

  } else {

    cli::cli_inform("Adjusting enrichment estimates by calculating overlap in annotations")
    freq <- arrow::read_parquet(system.file("extdata/common_snps.parquet", package = "ldsR"))
    m <-  data[["annot_ref"]][freq$common, -1]
    overlap_matrix <- crossprod(as.matrix(m))
    M_tot <- nrow(m)
    overlap_res <- overlapping_annotations(overlap_matrix = overlap_matrix, M_tot = M_tot, M = covar_M, jknife = res)


    dplyr::inner_join(base_results, overlap_res)

  }




}

#' Print names of annotations in a ldscore fileset
#' @param ldscore_dir A filepath to a ldscore directory
#'
#' @returns a character vector of column names
#' @export
#'
#' @examples \dontrun{
#' get_annot_names("ldscore/dir/celltypes")
#' }
get_annot_names <- function(ldscore_dir) {
  arrow::read_parquet(fs::path(ldscore_dir, "annot.parquet")) |> dplyr::pull(annot)
}
