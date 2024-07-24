



test_that("parse_gwas", {
  skip()
  expect_no_error(
    parse_gwas(
      test_path("fixtures/test_ldsc.sumstats.gz"),
      format = "ldsc"
    )
  )


  expect_error(
    parse_gwas(
      test_path("fixtures/test_ldsc.sumstats.gz"),
      format = "dataframe"
    )
  )





})

test_that("munge works", {
  df <- arrow::read_tsv_arrow(test_path("fixtures/test_ldsc.sumstats.gz"))


  expect_no_error(munge(tidyr::drop_na(df)))


})





test_that("ldsc_to_parquet works", {
  skip()
  snps <- arrow::read_parquet("tests/testthat/fixtures/baseline_v1.1/ld.parquet") |>
    dplyr::slice_head(n = 100000)
  write_parquet("")


  ld <- arrow::read_parquet(fs::path(outdir, "ld.parquet"))
  arrow::read_parquet(fs::path("annot_ref.parquet"))
  annot <- arrow::read_parquet(fs::path(outdir, "annot.parquet"))







})
