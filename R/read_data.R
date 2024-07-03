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


align_to_ref <- function(dset) {
  # EffectAllele is harmonized to always be the reference allele
  dset |>
    dplyr::mutate(
      EA_is_ref = dplyr::if_else(EffectAllele == REF, TRUE,FALSE),
      tmp = EffectAllele,
      EffectAllele = dplyr::if_else(EA_is_ref, EffectAllele, OtherAllele),
      OtherAllele = dplyr::if_else(EA_is_ref, OtherAllele, tmp),
      B = dplyr::if_else(EA_is_ref, B, B*-1),
      EAF = dplyr::if_else(EA_is_ref, EAF, 1-EAF)
    ) |>
    dplyr::select(-dplyr::all_of(c("EA_is_ref", "tmp")))

}
munge_tidygwas <- function(path, freq=0.9, info=0.9) {
  arrow::open_dataset(path) |> 
    align_to_ref() |> 
    dplyr::filter(!multi-allelic) |> 
    dplyr::filter(!indel) |> 
    # remove strand ambigious SNPs
    # remove MAF
    # filter INFO
    # Filter sample sizes outside 2-3x sd
    # from ldsc:    n_min = args.n_min if args.n_min else dat.N.quantile(0.9) / 1.5
    # merge-alleles

}

ldsc_munge <- function(sumstats) {
  stopifnot("data.frame" %in% class(sumstats) | rlang::is_scalar_character(sumstats))
  
  
}

parse_gwas <- function(tbl, ldsc_munge = TRUE) {

  # check if tbl is a filepath or in-memory dataframe or arrow::dataset connection
  
  # should ldsc_munge steps be applied?

  # should munging be done?

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

