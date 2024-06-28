
test_that("celltype_ldsc works", {
  skip()
  sumstat <- arrow::read_parquet(test_path("fixtures/move/speed_scz2022_eur.parquet"), col_select = c("SNP", "Z","N"))
  weights <- arrow::read_parquet(test_path("fixtures/move/weights2.parquet"), col_select = c("SNP", "L2"))
  baseline_ldscores <- arrow::read_parquet(test_path("fixtures/move/ld.parquet"))
  baseline_M <- arrow::read_parquet(test_path("fixtures/move/annot.parquet"))
  celltype_ldscores <- arrow::read_parquet(test_path("fixtures/superclusters/ld.parquet")) |>
    dplyr::mutate(SNP = stringr::str_remove(SNP, "rs") |> as.numeric())
  celltype_M <- arrow::read_parquet(test_path("fixtures/superclusters/annot.parquet"))


})





