
test_that("celltype_ldsc works", {
  skip()


  sumstat <- arrow::read_parquet(test_path("fixtures/test_data.parquet"), col_select = c("SNP", "Z" ="Z.x", "N" = "N.x")) |>
    dplyr::select(dplyr::all_of(c("SNP", "Z" ="Z.x", "N" = "N.x")))
  covariate_dir <- "~/projects/move"
  celltype_dir <- "~/projects/move/superclusters/"
  weights = NULL


})

