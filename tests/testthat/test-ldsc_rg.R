testdata <- arrow::read_parquet(test_path("fixtures/test_data.parquet"))
weights <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"))
testdata <- dplyr::inner_join(weights, testdata, by = "SNP")
tmp <- testdata |>
    dplyr::slice_head(n = 1000000) |> 
    dplyr::mutate(A1 = "A", A2 = "T")
sumstats1 <- dplyr::select(tmp, SNP, Z = Z.x, N = N.x, A1, A2)
sumstats2 <- dplyr::select(tmp, SNP, Z = Z.y, N = N.y, A1, A2)


test_that("Test inner implementation",{

  
  expect_no_error(.rg(sumstats1, sumstats2, M = 1173569, weights=weights, n_blocks=200))

  expect_no_error(.rg(sumstats1, sumstats2, M = 1173569, weights=weights, n_blocks=200))

})

test_that("ldsc_rg works", {

  res <- ldsc_rg(sumstats1, list("scz" = sumstats2, "bip" = sumstats2))

})

test_that("check_columns detects missing columns", {

  expect_error(
    check_columns(
      c("SNP", "A1", "A2"),
      dplyr::tibble("SNP" = "rs", "A1" = "A")
    )
  )

  expect_no_error(
    check_columns(
      c("SNP", "A1", "A2"),
      dplyr::tibble("SNP" = "rs", "A1" = "A", "A2" = "T")
    )
  )

  s1 <- dplyr::select(sumstats1, -A1, -A2)
  expect_error(
    ldsc_rg(
      s1,
      sumstats2
    )
  )



})
