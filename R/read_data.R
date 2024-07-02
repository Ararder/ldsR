# check_is_path() <- function(path) {
#   rlang::is_scalar_character(path)
#   # OBS: add more informative path
#   stopifnot(file.exists(path))

# }



}




parse_gwas <- function(path, ldsc_munge = TRUE) {

}





parse_parquet_dir <- function(dir) {
  ld <- arrow::read_parquet(paste0(dir, "/ld.parquet"))
  annot <- arrow::read_parquet(paste0(dir, "/annot.parquet"))

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

