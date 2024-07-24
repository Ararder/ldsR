utils::globalVariables(c("strand_ambig", "INFO", "EAF", "N", "RSID", "EffectAllele", "OtherAllele"))

parse_gwas <- function(df) {
  req_columns <- c("SNP", "Z", "N", "A1", "A2")
  ref <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"), col_select = c("SNP"))

  if(rlang::is_scalar_character(df)) {
    cli::cli_inform("Assuming GWAS to be a file to be read")
    check_is_path(df)
    df <- arrow::read_tsv_arrow(df)
    check_columns(c("SNP", "A1","A2","Z","N"), df)


  } else if("data.frame" %in% class(df)) {
    cli::cli_inform("In-memory data.frame passed...")
    check_columns(c("SNP", "A1","A2","Z","N"), df)

    df <- dplyr::tibble(df) |>
      tidyr::drop_na() |>
      dplyr::semi_join(ref, by = "SNP")



  } else if("Dataset" %in% class(df) | "arrow_dplyr_query" %in% class (df)) {

    df <- df |>
      dplyr::rename(SNP = RSID, A1 = EffectAllele, A2 = OtherAllele) |>
      dplyr::select(dplyr::any_of(c("SNP", "A1","A2", "Z","N", "INFO", "EAF"))) |>
      dplyr::filter(SNP %in% ref$SNP) |>
      dplyr::collect()

  }

  munge(df)

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
  req_columns <- c("SNP", "Z", "N", "A1", "A2")
  check_columns(req_columns, dset)

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

  dplyr::select(step1, dplyr::all_of(req_columns))
}

read_celltype_parquet <- function(path) {
  check_is_path(path)
  dset <- arrow::open_dataset(path) |> dplyr::collect()

  kk <- purrr::map_dbl(colnames(dset)[-1], \(x) dset[["metadata"]][["r"]][["columns"]][[x]][["attributes"]][["m50"]])


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

  m50 <- fs::dir_ls(dir, glob = "*M_5_50") |>
    purrr::map_dbl(\(x) readLines(x) |> as.numeric()) |>
    sum()
  m <- fs::dir_ls(dir, glob = "*M") |>
    purrr::map_dbl(\(x) readLines(x) |> as.numeric()) |>
    sum()

  list(ld, m50, m)

}


ldsc_to_parquet2 <- function(dir, outdir) {
  ld <- fs::dir_ls(dir, glob = "*ldscore.gz") |>
    purrr::map(arrow::read_tsv_arrow, col_select = -c("CHR", "BP")) |>
    purrr::list_rbind()

  annot_names <- colnames(ld)[-1]

  m50 <-
    fs::dir_ls(dir, glob = "*M_5_50") |>
    purrr::map(\(x) arrow::read_tsv_arrow(x, col_names = FALSE)) |>
    purrr::list_rbind() |>
    dplyr::summarise(dplyr::across(dplyr::everything(), sum)) |>
    purrr::set_names(annot_names) |>
    tidyr::pivot_longer(dplyr::everything(), names_to = "annot", values_to = "m50")

  m <- fs::dir_ls(dir, glob = "*M") |>
    purrr::map(\(x) arrow::read_tsv_arrow(x, col_names = FALSE)) |>
    purrr::list_rbind() |>
    dplyr::summarise(dplyr::across(dplyr::everything(), sum)) |>
    purrr::set_names(annot_names) |>
    tidyr::pivot_longer(dplyr::everything(), names_to = "annot", values_to = "m")

  annot_ref <-
    fs::dir_ls(dir, glob = "*annot.gz") |>
    purrr::map(arrow::read_tsv_arrow, col_select = -c("CHR", "BP")) |>
    purrr::list_rbind()

  annot <- dplyr::inner_join(m50, m, by = "annot")

  arrow::write_parquet(ld, fs::path(outdir, "ld.parquet"))
  arrow::write_parquet(annot_ref, fs::path(outdir, "annot_ref.parquet"))
  arrow::write_parquet(annot, fs::path(outdir, "annot.parquet"))


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
