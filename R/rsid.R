utils::globalVariables(c("CHR", "POS", "POS37", "POS38", "REF", "ALT", "A1", "A2", "SNP", ".row"))

#' Map CHR:POS to RSID
#'
#' @description
#' `to_rsid()` adds a `SNP` column with RSIDs to summary statistics that only
#' have genomic coordinates, so that they can be passed to [munge()] and the
#' downstream `ldsc_*()` functions.
#'
#' Variants are matched on `CHR`, `POS` and the unordered allele pair
#' (`A1`/`A2` vs `REF`/`ALT`) against a bundled reference covering all SNPs in the
#' default LD score and weight files (HapMap3, dbSNP155 coordinates for GRCh37 and GRCh38).
#' Rows that do not match the reference are removed.
#'
#' @param dset a [dplyr::tibble()] with columns `CHR`, `POS`, `A1` and `A2`.
#' `CHR` can be given as `1` or `chr1`. Only autosomes are matched.
#' @param build genome build of `POS`: `"guess"` (default) picks the build with
#' the most matching variants, or pass `"37"` or `"38"`.
#'
#' @return `dset` with a `SNP` column, restricted to variants found in the reference
#' @export
#'
#' @examples \dontrun{
#' sumstats |>
#'   to_rsid() |>
#'   munge() |>
#'   ldsc_h2()
#' }
to_rsid <- function(dset, build = c("guess", "37", "38")) {
  stopifnot("data.frame" %in% class(dset))
  build <- rlang::arg_match(build)
  check_columns(c("CHR", "POS", "A1", "A2"), dset)

  if("SNP" %in% colnames(dset)) {
    cli::cli_alert_info("Overwriting existing {.field SNP} column")
    dset <- dplyr::select(dset, -SNP)
  }

  ref <- arrow::read_parquet(system.file("extdata", "rsid_map.parquet", package = "ldsR"))

  keys <- dplyr::mutate(
    dset,
    .row = dplyr::row_number(),
    CHR = suppressWarnings(as.integer(sub("^chr", "", as.character(CHR), ignore.case = TRUE))),
    POS = as.integer(POS),
    A1 = toupper(A1),
    A2 = toupper(A2)
  ) |>
    dplyr::select(".row", "CHR", "POS", "A1", "A2")

  builds <- if(build == "guess") c("37", "38") else build
  matched <- purrr::map(builds, \(b) match_build(keys, ref, b)) |>
    rlang::set_names(builds)
  n_matched <- purrr::map_int(matched, nrow)

  if(build == "guess") {
    cli::cli_alert_info(
      "Guessing genome build: {n_matched[['37']]} variants matched GRCh37, {n_matched[['38']]} matched GRCh38"
    )
    build <- names(which.max(n_matched))
  }

  hits <- matched[[build]]
  if(nrow(hits) == 0) {
    cli::cli_abort("No variants in {.arg dset} matched the reference on GRCh{build}")
  }
  cli::cli_alert_success("Using GRCh{build}")
  cli::cli_alert_warning("Removed {nrow(dset) - nrow(hits)} rows that could not be mapped to an RSID")

  dset[hits$.row, ] |>
    dplyr::mutate(SNP = hits$SNP, .before = 1)
}


match_build <- function(keys, ref, build) {
  pos_col <- paste0("POS", build)

  ref <- dplyr::select(ref, "CHR", POS = dplyr::all_of(pos_col), "SNP", "REF", "ALT")

  dplyr::inner_join(keys, ref, by = c("CHR", "POS"), relationship = "many-to-many") |>
    dplyr::filter((A1 == REF & A2 == ALT) | (A1 == ALT & A2 == REF)) |>
    dplyr::distinct(.row, .keep_all = TRUE) |>
    dplyr::arrange(.row) |>
    dplyr::select(".row", "SNP")
}
