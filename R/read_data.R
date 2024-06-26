parse_parquet <- function(gwas, ldscores, weights) {
  stopifnot(fs::dir_exists(ldscores))
  stopifnot(fs::file_exists(gwas))
  stopifnot(fs::file_exists(weights))

  g <- arrow::read_parquet(gwas, col_select = c("SNP", "Z", "N"))
  w <- arrow::read_parquet(weights, col_select = c("SNP", "L2"))

  unpack <- parse_parquet_dir(dir)
  ld <- unpack[[1]]
  annot <- unpack[[2]]
  rm(unpack)


  merged <- dplyr::inner_join(g, w, by = "SNP") |>
    dplyr::inner_join(ld, by = "SNP")

  remove <- c(colnames(g), colnames(w))

  y <- merged$Z^2
  x <- dplyr::select(merged,-dplyr::all_of(remove)) |> as.matrix()
  w <- merged$L2
  N <- merged$N
  M <- annot$m50

  list(
    y = y,
    x = x,
    w = w,
    N = N,
    M = M
  )



}

parse_parquet_dir <- function(dir) {
  ld <- arrow::read_parquet(paste0(dir, "/ld.parquet"))
  annot <- arrow::read_parquet(paste0(dir, "/annot.parquet"))

  list(ld, annot)
}

# dir = "~/Downloads/superclusters/amygdala_excitatory/"
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
