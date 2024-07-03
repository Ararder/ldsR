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

#

test_that("benchmarking", {
  skip()
  weights <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"), col_select = c("SNP", "L2"))
  M <- 1173569
  df <- dplyr::select(testdata, SNP, Z = Z.y, N = N.y)
  m <- dplyr::inner_join(weights, df)
  # profvis::profvis(res <- ldsc_h2(dplyr::select(testdata, SNP, Z = Z.y, N = N.y)))
  # microbenchmark::microbenchmark(
  #   ldscore(y = m$Z^2, x = m$L2 |> as.matrix(), w = m$L2, N = m$N,M = M),
  #   times = 10L
  # )
  # profvis::profvis(ldscore(y = m$Z^2, x = m$L2 |> as.matrix(), w = m$L2, N = m$N,M = M))
  


})





test_that("speeding up", {
  skip()
  weights <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"), col_select = c("SNP", "L2"))
  df <- dplyr::select(testdata, SNP, Z = Z.y, N = N.y)
  m <- dplyr::inner_join(weights, df)
  
  y = m$Z^2
  x = m$L2 |> as.matrix()
  w = m$L2
  N = m$N
  M <- 1173569
  
  n_snp <- dim(x)[1]
  n_annot <- dim(x)[2]
  M_tot <- sum(M)
  x_tot <- rowSums(x)

  # provide a starting estimate of heritability
  hsq <- M_tot * (mean(y) - 1) / mean((x_tot * N))

  # first update of weights
  initial_w <- get_weights(ld = x_tot, w_ld = w, N = N, M = M_tot, hsq = hsq)

  # Normalise by the mean of N and add intercept
  Nbar <- mean(N)
  x <- (N*x) / Nbar
  x <- cbind(1, x)
})
