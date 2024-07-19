
test_that("partitioned_heratbility", {
  skip()
  testdata <- arrow::read_parquet(test_path("fixtures/test_data.parquet"))
  sumstat <- dplyr::select(testdata, SNP, Z = Z.x, N = N.x)
  covars <- parse_parquet_dir("~/projects/move/superclusters/")

  


})

