testdata <- arrow::read_parquet(test_path("fixtures/test_data.parquet"))
weights <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"))
testdata <- dplyr::inner_join(weights, testdata, by = "SNP")

test_that("ldsc_h2 per chrom", {
  # data <- fs::dir_ls("~/Downloads/sldsc_ref/1000G_Phase3_frq/", glob = "*frq") |>
  #   purrr::map(readr::read_table) |>
  #   purrr::list_rbind()
  #
  # x <- testdata |>
  #   dplyr::left_join(dplyr::select(data, CHR, SNP), by = c("SNP" = "SNP")) |>
  #   dplyr::rename(Z = Z.x, N = N.x)
  #
  #
  # res <- split(x, x$CHR) |>
  #   purrr::map(ldsc_h2)
  #
  expect_no_error(res <- ldsc_h2(dplyr::select(testdata, SNP, Z = Z.x, N = N.x)))
  expect_equal(res$h2, 0.35893491)
  expect_equal(res$int, 1.08191455)


})


test_that("pop and samp prev works", {

  res <- ldsc_h2(dplyr::select(testdata, SNP, Z = Z.x, N = N.x), pop_prev = 0.01, sample_prev = 53000/(73000 + 53000))
  expect_equal(res$lia_h2, 0.203, tolerance = 1e-02)
  expect_equal(res$int, 1.08191455)


})

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

# test_that("benchmark and profile h2", {
#   skip()
#   s1 <- dplyr::select(testdata, SNP, Z = Z.y, N = N.y)
#
#   bench <- microbenchmark::microbenchmark(
#     ldsc_h2(s1),
#     times = 10L
#   )
#
#   time <- profvis::profvis(ldsc_h2(s1))
#
#
#
# })

test_that("partitioned heritability runs and reproduces results", {

  s1 <- dplyr::select(testdata, SNP, Z = Z.x, N = N.x)
  res <- partitioned_h2(
    sumstat = s1,
    ldscore_dir = test_path("fixtures/baseline_v1.1")
  )
  expect_equal(res$tot[1], 0.3347495, tolerance = 1e-06)



})

test_that("cell-type analysis runs and reproduces results", {

  s1 <- dplyr::select(testdata, SNP, Z = Z.x, N = N.x)
  expect_no_error(
    res <- celltype_analysis(
      sumstat = s1,
      covariate_dir = test_path("testdata/baseline"),
     ldscore_dir = test_path("testdata/superclusters")
    )
  )




})

