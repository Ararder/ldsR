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






#' Read LDscore format from ldsRs parquet format
#'
#' @param dir directory containing the files 'annot.parquet', 'annot_ref.parquet' and 'ld.parquet'
#' @param read_ref logical, whether to read the reference file 'annot_ref.parquet'
#'
#' @return a [list()]
#' @export
#'
#' @examples \dontrun{
#' files <- parse_parquet_dir("ldscores/atac")
#' }
#'
parse_parquet_dir <- function(dir, read_ref=FALSE) {
  ld_path <- paste0(dir, "/ld.parquet")
  annot_path <- paste0(dir, "/annot.parquet")
  check_is_path(ld_path)
  check_is_path(annot_path)


  ld <- arrow::read_parquet(ld_path)
  annot <- arrow::read_parquet(annot_path)

  if(isTRUE(read_ref)) {
    annot_ref_path <- paste0(dir, "/annot_ref.parquet")
    check_is_path(annot_ref_path)
    annot_ref <- arrow::read_parquet(annot_ref_path)

  }

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



  if(isTRUE(read_ref)) {
    list(
      "ld" = ld,
      "annot" = annot,
      "annot_ref" = annot_ref
    )

  } else{
    list(
      "ld" = ld,
      "annot" = annot
    )
  }


}




#' Transform a directory of LDscores to parquet format
#'
#' @param dir directory containing the python LDSC ldscore format
#' @param thin If the --thin flag has been used, provide a character vector of RSIDs for the full dataset used to calculate LDscores
#'
#' @return a [list()]
#' @export
#'
#' @examples \dontrun{
#' ldsc_to_parquet("/directory/ldsc")
#' }
ldsc_to_parquet <- function(dir, thin=FALSE) {
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

  annot <- dplyr::tibble(annot = annot_name, m50 = m50, m = m)


  if(!isTRUE(thin)) {
    annot_ref <-
      fs::dir_ls(dir, glob = "*annot.gz") |>
      purrr::map(\(x) arrow::read_tsv_arrow(x, col_select = c(5))) |>
      purrr::list_rbind() |>
      purrr::set_names(annot_name)

  } else {
    annot_ref <-
      fs::dir_ls(dir, glob = "*annot.gz") |>
      purrr::map(\(x) arrow::read_tsv_arrow(x)) |>
      purrr::list_rbind() |>
      purrr::set_names(annot_name)

  }



  list("ld" = ld, "annot" = annot, "annot_ref" = annot_ref)

}

get_snps <- function(dir) {


  fs::dir_ls(dir, glob = "*annot.gz") |>
    purrr::map(\(x) arrow::read_tsv_arrow(x, col_select = c(3) )) |>
    purrr::list_rbind()


}

#' Convert the output of many LDscores to a single dataframe suitable for ldsR
#'
#' @param parent_dir a directory with subdirectories containing LDscore data
#' @param outdir directory to save the parquet files
#' @param thin a character vector of RSIDs corresponding to SNPs used in the
#' full dataset used to generate LDscores
#'
#' @return NULL
#' @export
#'
#' @examples \dontrun{
#' to_celltype_dataset("files/ldsc", "files/ldsc_parquet")
#' }
to_celltype_dataset <- function(parent_dir, outdir, thin = NULL) {
  fs::dir_create(outdir)
  stopifnot(fs::dir_exists(parent_dir))

  # read in list data
  list_data <- purrr::map(fs::dir_ls(parent_dir, type = "dir"),\(x) ldsc_to_parquet(x, thin = thin), .progress = list(type = "tasks", name = "reading in raw ldscore data"))

  # SNPs should the same in all directories, can get from first directory
  if(is.null(thin)) {
    snps_in_ref <- get_snps(fs::dir_ls(parent_dir)[1])
  } else {
    snps_in_ref <- dplyr::tibble(SNP = thin)
  }

  # merge LDscore columns ---------------------------------------------------

  all_ld <- purrr::map(list_data, "ld") |>
    purrr::map(\(x) dplyr::select(x, -1)) |>
    unname() |>
    purrr::list_cbind() |>
    janitor::clean_names()

  snp <- dplyr::select(list_data[[1]][["ld"]], "SNP")
  all_ld <- dplyr::bind_cols(snp, all_ld)

  # Merge annot ---------------------------------------------------------------

  annot <- purrr::map(list_data, "annot") |>
    purrr::list_rbind()

  # merge annot_ref ---------------------------------------------------------
  # bind_cols, and then
  annot_ref <- purrr::map(list_data, "annot") |>
    unname() |>
    purrr::list_cbind()
  annot_ref <- dplyr::bind_cols(snps_in_ref, annot_ref)


  arrow::write_parquet(annot, fs::path(outdir, "annot.parquet"))
  arrow::write_parquet(all_ld, fs::path(outdir, "ld.parquet"))
  arrow::write_parquet(annot_ref, fs::path(outdir, "annot_ref.parquet"))

}




#' Merge two sets of LD data, and save to parquet
#'
#' @param list1 output of [ldsc_to_parquet()] or [parse_parquet_dir()]
#' @param list2 output of [ldsc_to_parquet()] or [parse_parquet_dir()]
#' @param outdir directory to store merged LD data
#'
#' @return NULL
#' @export
#'
#' @examples \dontrun{
#' combe_ld_data(l1, l2, "path/to/storage")
#' }
combine_ld_data <- function(list1, list2, outdir) {
  stopifnot(all(names(list1) == names(list2)))

  ld <-  dplyr::bind_cols(list1$ld, dplyr::select(list2$ld, -dplyr::any_of(c("SNP"))))
  annot <- dplyr::bind_rows(list1$annot, list2$annot)
  annot_ref <- dplyr::bind_cols(list1$annot_ref, dplyr::select(list2$annot_ref, -dplyr::any_of(c("SNP"))))


  ld_path <- paste0(outdir, "/ld.parquet")
  annot_path <- paste0(outdir, "/annot.parquet")
  ref_path <- paste0(outdir, "/annot_ref.parquet")


  arrow::write_parquet(ld, ld_path)
  arrow::write_parquet(annot, annot_path)
  arrow::write_parquet(annot_ref, ref_path)


}


