
test_that("celltype_ldsc works", {
  skip()
  
  
  sumstat <- arrow::read_parquet(test_path("fixtures/test_data.parquet"), col_select = c("SNP", "Z" ="Z.x", "N" = "N.x")) |>
    dplyr::select(dplyr::all_of(c("SNP", "Z" ="Z.x", "N" = "N.x")))
  covariate_dir <- "~/projects/move"
  celltype_dir <- "~/projects/move/superclusters/"  
  weights = NULL
  

})


test_that("Genetic correlatios works", {
  skip()
  
  sumstat <- arrow::read_parquet(test_path("fixtures/test_data.parquet"))
  
  s1 <- dplyr::select(sumstat, SNP, Z = Z.x, N = N.x)
  s2 <- dplyr::select(sumstat, SNP, Z = Z.y, N = N.y)
  
  expect_no_error(test <- ldsc_rg(s1,s2))
  # profvis::profvis(ldsc_rg(s1, s2))
  # microbenchmark::microbenchmark(
  #   ldsc_rg(s1, s2),
  #   times = 10L
  # )
  

})








