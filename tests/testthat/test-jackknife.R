
test_that("celltype_ldsc works", {
  skip()
  sumstat <- arrow::read_parquet(test_path("fixtures/move/speed_scz2022_eur.parquet"), col_select = c("SNP", "Z","N"))
  weights <- arrow::read_parquet(test_path("fixtures/move/weights2.parquet"), col_select = c("SNP", "L2"))
  baseline_ldscores <- arrow::read_parquet(test_path("fixtures/move/ld.parquet"))
  baseline_M <- arrow::read_parquet(test_path("fixtures/move/annot.parquet"))
  celltype_ldscores <- arrow::read_parquet(test_path("fixtures/superclusters/ld.parquet")) |>
    dplyr::mutate(SNP = stringr::str_remove(SNP, "rs") |> as.numeric())
  celltype_M <- arrow::read_parquet(test_path("fixtures/superclusters/annot.parquet"))

  c <- c(colnames(weights), colnames(sumstat))
  merged <- dplyr::inner_join(weights, sumstat, by = "SNP") |>
    dplyr::inner_join(baseline_ldscores, by = "SNP")


  test1 <- dplyr::inner_join(merged, dplyr::select(celltype_ldscores, 1,2), by = "SNP")
  test2 <- dplyr::inner_join(merged, dplyr::select(celltype_ldscores, 1,3), by = "SNP")


})


