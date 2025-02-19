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





# note to self:
# should be able to take any set of ldscores
read_overlap_matrix <- function(ldscore_dir, ordered=TRUE) {
  stopifnot("`ordered` should be either TRUE or FALSE" = rlang::is_bool(ordered))

  if(ordered) {
    annots <- arrow::read_parquet(fs::path(ldscore_dir, "annot_ref.parquet")) |>
      dplyr::select(-dplyr::any_of(c("SNP","CM")))
    freq <- arrow::read_parquet(fs::path(ldscore_dir, "snp_freq.parquet")) |>
      dplyr::select(-dplyr::any_of(c("SNP")))

    annots <- dplyr::bind_cols(annots, freq) |>
      dplyr::filter(dplyr::.data[["MAF"]] > 0.05) |>
      dplyr::select(-dplyr::any_of(c("MAF")))

  } else {
    annots <- arrow::read_parquet(fs::path(ldscore_dir, "annot_ref.parquet"))
    freq <- arrow::read_parquet(fs::path(ldscore_dir, "snp_freq.parquet"))
    annots <- dplyr::semi_join(annots, dplyr::filter(freq, dplyr::.data[["MAF"]] > 0.05), by = "SNP")
    annots <- dplyr::select(annots,-dplyr::any_of(c("SNP", "CM")))

  }


  annots

}

create_overlap_matrix <- function(ldscore_dirs) {

  m <- purrr::map(ldscore_dirs, read_overlap_matrix) |>
    purrr::list_cbind()


  overlap <- crossprod(as.matrix(m))
  M_tot <- nrow(m)

  list("overlap_matrix" = overlap, M_tot = M_tot)

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
      The number of columns (ncol = {.bold {ncol(ld)} - 1} ) in {.path {ld_path}} should
      match the number of rows (nrow = {.bold {nrow(annot)}}) in {.path {annot_path}}"
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




#' Transform LDSC formatted annotation ldscores to ldsR format
#'
#' @param dir directory with ldscores in LDSC format
#' @param outdir directory to save the parquet files
#'
#' @return NULL
#' @export
#'
#' @examples \dontrun{
#' ldsc_to_parquet2("path/to/ldsc", "path/to/outdir")
#' }
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



ldsc_to_parquet <- function(dir) {
  annot_name <- fs::path_file(dir)

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

  annot_ref <-
    fs::dir_ls(dir, glob = "*annot.gz") |>
    purrr::map(\(x) arrow::read_tsv_arrow(x, col_select = c(5))) |>
    purrr::list_rbind() |>
    purrr::set_names(annot_name)

  list("ld" = ld, "m50" = m50, "m" = m, "annot" = annot_ref)

}

get_snps <- function(dir) {


  fs::dir_ls(dir, glob = "*annot.gz") |>
    purrr::map(\(x) arrow::read_tsv_arrow(x, col_select = c(3) )) |>
    purrr::list_rbind()


}

to_celltype_dataset <- function(parent_dir, outdir) {
  fs::dir_create(outdir)
  stopifnot(fs::dir_exists(parent_dir))

  # read in list data
  list_data <- purrr::map(fs::dir_ls(parent_dir), ldsc_to_parquet, .progress = list(type = "tasks", name = "reading in raw ldscore data"))

  # SNPs should the same in all directories, can get from first directory
  snps_in_ref <- get_snps(fs::dir_ls(parent_dir)[1])

  # merge LDscore columns ---------------------------------------------------

  all_ld <- purrr::map(list_data, "ld") |>
    purrr::map(\(x) dplyr::select(x, -1)) |>
    unname() |>
    purrr::list_cbind() |>
    janitor::clean_names()

  snp <- dplyr::select(list_data[[1]][["ld"]], "SNP")
  all_ld <- dplyr::bind_cols(snp, all_ld)

  # Merge m50 ---------------------------------------------------------------

  all_m50 <- purrr::map(list_data, "m50") |>
    purrr::imap(\(val, name) dplyr::tibble(annot = fs::path_file(name) |> janitor::make_clean_names(), m50 = val))  |>
    purrr::list_rbind()


  # merge M -----------------------------------------------------------------

  all_m <- purrr::map(list_data, "m") |>
    purrr::imap(\(val, name) dplyr::tibble(annot = fs::path_file(name) |> janitor::make_clean_names(), m = val))  |>
    purrr::list_rbind()

  # merge annot_ref ---------------------------------------------------------
  # bind_cols, and then
  annot_ref <- purrr::map(list_data, "annot") |>
    unname() |>
    purrr::list_cbind()
  annot_ref <- dplyr::bind_cols(snps_in_ref, annot_ref)



  annot <- dplyr::inner_join(all_m50, all_m, by = "annot")
  arrow::write_parquet(annot, fs::path(outdir, "annot.parquet"))
  arrow::write_parquet(all_ld, fs::path(outdir, "ld.parquet"))
  arrow::write_parquet(annot_ref, fs::path(outdir, "annot_ref.parquet"))

}

