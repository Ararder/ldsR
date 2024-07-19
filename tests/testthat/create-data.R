snps <- arrow::read_parquet("tests/testthat/fixtures/baseline_v1.1/ld.parquet") |> 
    dplyr::slice_head(n = 100000)

snps |> 
  arrow::write_parquet(test_path("testdata/baseline/ld.parquet"))
fs::file_copy("tests/testthat/fixtures/baseline_v1.1/annot.parquet", test_path("testdata/baseline/annot.parquet"))




# celll-types ------------------------------------------------------------
snps <- arrow::read_parquet("~/projects/move/superclusters/ld.parquet") |> 
    dplyr::slice_head(n = 100000)

arrow::write_parquet(snps, test_path("testdata/superclusters/ld.parquet"))
fs::file_copy("~/projects/move/superclusters/annot.parquet", test_path("testdata/superclusters/annot.parquet"))