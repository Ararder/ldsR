test_that("parse_parquet_dir works", {
  l <- parse_parquet_dir("~/projects/move")
  l <- parse_parquet_dir("~/projects/move/superclusters")
})
