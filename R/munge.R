#' Munge GWAS summary statistics
#'
#' @description
#' `munge()` applies 5 filters (if possible):
#' 1. Deduplication of RSID
#' 2. Filter variants on INFO in INFO column is present
#' 3. Filter variant on effect allele frequency
#' 4. Removes strand ambigious variants
#' 5. Removes variants with `N < round(stats::quantile(step1$N, 0.9) / 1.5)`
#'
#'
#' @param dset a [dplyr::tibble()] with columns `SNP`, `A1` `A2` `Z` `N` and
#' optional columns `EAF` and `INFO`
#' @param info_filter INFO score filter threshold at which to remove rows
#' @param maf_filter  Minor allele frequenc filter at which to remove rows
#'
#' @return a data.frame
#' @export
#'
#' @examples \dontrun{
#' parse_gwas(tbl)
#' }
munge <- function(dset, info_filter = 0.9, maf_filter = 0.01) {
  stopifnot("data.frame" %in% class(dset))
  req_columns <- c("SNP", "Z", "N", "A1", "A2")
  check_columns(req_columns, dset)

  snp_ref <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"), col_select = "SNP")

  before <- nrow(dset)
  step1 <- dplyr::filter(dset, SNP %in% snp_ref$SNP)
  cli::cli_alert_warning("Removed {before - nrow(step1)} rows after filtering on reference panel")


  if("INFO" %in% colnames(dset)) {
    before <- nrow(step1)
    step1 <- dplyr::filter(step1, INFO > info_filter)
    cli::cli_alert_warning("Removed {before - nrow(step1)} rows with INFO below {info_filter}")
  }

  if("EAF" %in% colnames(dset)) {
    before <- nrow(step1)
    step1 <- dplyr::filter(step1, EAF > maf_filter & EAF < (1-maf_filter))
    cli::cli_alert_warning("Removed {before - nrow(step1)} rows with MAF below {maf_filter}")
  }


  # Filter on Unique RSID
  before <- nrow(step1)
  step1 <- dplyr::distinct(step1, SNP, .keep_all = TRUE)
  cli::cli_alert_warning("Removed {before - nrow(step1)} rows with duplicated RSIDs")


  # Strand ambigius SNPS
  before <- nrow(step1)
  step1 <- remove_strand_ambig(step1)
  cli::cli_alert_warning("Removed {before - nrow(step1)} rows due to strand ambigious alleles")

  # Large deviation in sample size
  before <- nrow(step1)
  N_filter <- round(stats::quantile(step1$N, 0.9) / 1.5)
  step1 <- dplyr::filter(step1, N > N_filter)
  cli::cli_alert_warning("Removed {before - nrow(step1)} rows with a sample size smaller than {N_filter}")

  dplyr::select(step1, dplyr::all_of(req_columns))
}


remove_strand_ambig <- function(tbl) {

  dplyr::mutate(tbl, strand_ambig = dplyr::case_when(
    A1 == 'G' & A2 == 'C' ~ TRUE,
    A1 == 'C' & A2 == 'G' ~ TRUE,
    A1 == 'A' & A2 == 'T' ~ TRUE,
    A1 == 'T' & A2 == 'A' ~ TRUE,
    .default = FALSE
  )) |>
    dplyr::filter(!strand_ambig) |>
    dplyr::select(-strand_ambig)
}
