check_is_path <- function(path) {
  
  if(!rlang::is_scalar_character(path)) {
    stop(cli::format_error(
      "The provided argument is not a string: {path}"
    ))
  }
  # OBS: add more informative path
  if(!file.exists(path)) {
    stop(cli::format_error(
      "The filepath ({.path {path}}) does not exist on the host system"
    ))
  }

}



check_numeric_columns <- function(ld) {
  before <- colnames(ld)
  ld <- dplyr::select(ld, "SNP", dplyr::where(is.numeric))
  after <- colnames(ld)
  
  if(!all(before %in% after)) {
    cli::cli_alert_danger(
      "Columns were removed from ld.parquet because they were non-numeric:
        {.code {before[!before %in% after]}}
        "
    )
  }

  ld
}


parse_gwas <- function(path, ldsc_munge = TRUE) {

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


#
# CONVERT LDSC file structure to parquet
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

