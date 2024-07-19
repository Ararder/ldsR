utils::globalVariables(c("strand_ambig", "INFO", "EAF", "N", "RSID", "EffectAllele", "OtherAllele"))
#' Parse a GWAS summary statistics
#'
#' @param df a in memory [base::data.frame()] or a filepath to ldsc-munged
#' summary statistics file or a filepath to a tidyGWAS folder.
#' @param format format of summary statistics. Either a in-memory dataframe, and LDSC munged file, or tidyGWAS filepath
#'
#' @return a data.frame
#' @export
#'
#' @examples \dontrun{
#' parse_gwas(tbl)
#' }
parse_gwas <- function(df, format = c("dataframe", "ldsc", "tidyGWAS")) {
  format <- rlang::arg_match(format)
  req_columns <- c("SNP", "Z", "N", "A1", "A2")
  ref <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"), col_select = c("SNP"))

  if(rlang::is_scalar_character(path)) {
    cli::cli_inform("Assuming GWAS to be a on-disk file...")
    check_is_path(path)
    df <- arrow::read_tsv_arrow(path)
    check_columns(c("SNP", "A1","A2","Z","N"), df)
    
    tidyr::drop_na(df)
    
  } else if("data.frame" %in% class(df)) {
    cli::cli_inform("In-memory data.frame passed...")
    check_columns(c("SNP", "A1","A2","Z","N"), df)
  
    final <- dplyr::tibble(df) |> 
      tidyr::drop_na() |> 
      dplyr::semi_join(ref, by = "SNP")
  
    munge(final)
    
  } else if("Dataset" %in% class(df) | "arrow_dplyr_query" %in% class (df)) {
    
    df |> 
      dplyr::rename(SNP = RSID, A1 = EffectAllele, A2 = OtherAllele) |>
      dplyr::select(dplyr::any_of(c("SNP", "A1","A2", "Z","N", "INFO", "EAF"))) |>
      dplyr::filter(SNP %in% ref$SNP) |>
      dplyr::collect() |> 
      munge()
  }

  munge(final)

}






#' Munge GWAS summary statistics
#'
#' @param dset a [dplyr::tibble()] with columns `SNP`, `A1` `A2` `Z` `N` and possibly `EAF` and `INFO`
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
  req_columns <- c("SNP", "Z", "N", "EffectAllele", "OtherAllele")

  before <- nrow(dset)
  step1 <- dplyr::distinct(dset, SNP, .keep_all = TRUE)
  cli::cli_alert_warning("Removed {before - nrow(step1)} rows with duplicated RSIDs")


  if("INFO" %in% colnames(dset)) {
    before <- nrow(dset)
    step1 <- dplyr::filter(step1, INFO > info_filter)
    cli::cli_alert_warning("Removed {before - nrow(step1)} rows with INFO below {info_filter}")
  }

  if("EAF" %in% colnames(dset)) {
    before <- nrow(step1)
    step1 <- dplyr::filter(step1, EAF > maf_filter & EAF < (1-maf_filter))
    cli::cli_alert_warning("Removed {before - nrow(step1)} rows with MAF below {maf_filter}")
  }

  before <- nrow(step1)
  step1 <- remove_strand_ambig(step1)
  cli::cli_alert_warning("Removed {before - nrow(step1)} rows due to strand ambigious alleles")

  before <- nrow(step1)
  N_filter <- round(stats::quantile(step1$N, 0.9) / 1.5)
  step1 <- dplyr::filter(step1, N > N_filter)
  cli::cli_alert_warning("Removed {before - nrow(step1)} rows with a sample size smaller than {N_filter}")
  step1
}




parse_parquet_dir <- function(dir) {
  ld_path <- paste0(dir, "/ld.parquet")
  annot_path <- paste0(dir, "/annot.parquet")
  check_is_path(ld_path)
  check_is_path(annot_path)

  ld <- arrow::read_parquet(ld_path)
  annot <- arrow::read_parquet(annot_path)

  if(ncol(ld) != nrow(annot)+1) {
    stop(cli::format_error(
      "The ldscore files provided do not match up.
      The number of columns ({.bold {ncol(ld)}} - 1) in {.path {ld_path}} should
      match the number of rows ({.bold {nrow(annot)}}) in {.path {annot_path}}"
    )
    )
  }

  ld <- check_numeric_columns(ld)
  if(!"SNP" %in% colnames(ld)) {
    stop(cli::format_error(
      "column {.code SNP} is missing from {.path {ld_path}}"
    ))
  }
  if(!"m50" %in% colnames(annot)) {
    stop(cli::format_error(
      "column {.code m50} is missing from {.path {annot_path}}"
    ))
  }

  list(
    "ld" = ld,
    "annot" = annot
  )
}


ldsc_to_parquet <- function(dir, annot_name) {
  ld <- fs::dir_ls(dir, glob = "*ldscore.gz") |>
    purrr::map(arrow::read_tsv_arrow, col_select = c("SNP", "L2")) |>
    purrr::list_rbind() |>
    purrr::set_names(c("SNP", annot_name))

  # annot <- fs::dir_ls(dir, glob = "*annot.gz") |>
  #   purrr::map(arrow::read_tsv_arrow, col_select = c(3,5)) |>
  #   purrr::list_rbind() |>
  #   dplyr::rename({{ annot_name }} := 2)
  m50 <- fs::dir_ls(dir, glob = "*M_5_50") |>
    purrr::map_dbl(\(x) readLines(x) |> as.numeric()) |>
    sum()
  m <- fs::dir_ls(dir, glob = "*M") |>
    purrr::map_dbl(\(x) readLines(x) |> as.numeric()) |>
    sum()

  list(ld, dplyr::tibble(annot_name = annot_name, m50 = m50, m = m))

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
