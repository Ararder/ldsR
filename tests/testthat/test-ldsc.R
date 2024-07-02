testdata <- arrow::read_parquet(test_path("fixtures/test_data.parquet"))
weights <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"))
testdata <- dplyr::inner_join(weights, testdata, by = "SNP")

test_that("ldsc_h2 runs and reproduces LDSC for scz", {

  expect_no_error(res <- ldsc_h2(dplyr::select(testdata, SNP, Z = Z.x, N = N.x)))
  expect_equal(res$h2, 0.35893491)
  expect_equal(res$int, 1.08191455)


})

test_that("ldsc_h2 runs and reproduces LDSC for bip", {

  expect_no_error(res <- ldsc_h2(dplyr::select(testdata, SNP, Z = Z.y, N = N.y)))
  expect_equal(res$h2, 0.07079196, tolerance = 1e-05)
  expect_equal(res$int, 1.02470645)


})







