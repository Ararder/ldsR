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


test_that("partitioned heritability runs and reproduces results", {

  s1 <- dplyr::select(testdata, SNP, Z = Z.x, N = N.x)
  res <- partitioned_h2(
    sumstat = s1,
    ldscore_dir = test_path("fixtures/baseline_v1.1")
  )
  expect_equal(res$tot[1], 0.3347495, tolerance = 1e-06)



})

test_that("cell-type analysis runs and reproduces results", {
  skip()
  s1 <- dplyr::select(testdata, SNP, Z = Z.x, N = N.x)
  expect_no_error(
    res <- celltype_analysis(
      sumstat = s1,
      covariate_dir = test_path("testdata/baseline"),
      ldscore_dir = test_path("testdata/superclusters")
    )
  )

})



test_that("partitioned_h2 can adjust for overlapping annotations", {
  skip("Requires external data")
  ldscore_dir = "~/Desktop/baseline_v1.1/"
  testdata <- arrow::read_parquet(test_path("fixtures/test_data.parquet"))
  sumstat <- dplyr::select(testdata, SNP, Z = Z.x, N = N.x)



  weights <- arrow::read_parquet(system.file("extdata/eur_w_ld.parquet", package = "ldsR"), col_select = c("SNP", "L2_celltype")) |>
    dplyr::rename(L2 = "L2_celltype") |>
    dplyr::filter(!is.na(L2))


  covars <- parse_parquet_dir(ldscore_dir)
  covar_ld <- covars[["ld"]]
  covar_M <- covars[["annot"]][["m50"]]

  n_before <- nrow(sumstat)
  merged <- dplyr::inner_join(weights, sumstat, by = "SNP")
  cli::cli_alert_info("Removed {.bold {n_before - nrow(merged)}} rows after merging with weights")
  n_before <- nrow(merged)
  merged <- dplyr::inner_join(merged, covar_ld, by = "SNP")
  cli::cli_alert_info("Removed {.bold {n_before - nrow(merged)}} rows after merging with ldscores")
  remove_cols <- unique(c(colnames(sumstat), colnames(weights)))
  x <- dplyr::select(merged,-dplyr::all_of(remove_cols)) |> as.matrix()

  y = merged$Z^2
  x = x
  w = merged$L2
  N = merged$N
  M = as.double(covar_M)
  n_blocks=200

  # -------------------------------------------------------------------------
  res <- ldscore(y = y, x = x, w = w, N = N, M = M, n_blocks = 200)



  # -------------------------------------------------------------------------


  # n_snp <- dim(x)[1]
  # n_annot <- dim(x)[2]
  # M_tot <- sum(M)
  # stopifnot("The number of annotations in x should be equal to the of M" = length(M) == n_annot)
  # x_tot <- rowSums(x)
  #
  # # provide a starting estimate of heritability
  # hsq <- M_tot * (mean(y) - 1) / mean((x_tot * N))
  #
  # # first update of weights
  # initial_w <- get_weights(ld = x_tot, w_ld = w, N = N, M = M_tot, hsq = hsq)
  #
  # # Normalise by the mean of N and add intercept
  # Nbar <- mean(N)
  # x <- (N*x) / Nbar
  # x <- cbind(1, x)
  #
  # initial_w = sqrt(initial_w)
  # initial_w = initial_w / sum(initial_w)
  # x_weighted <- x*initial_w
  # y_weighted <- y*initial_w
  # jknife = lstq_jackknife(x = x_weighted, y = y_weighted, n_blocks = n_blocks)
  #
  # ex_jk <- extract_jackknife(jknife, M = M, Nbar = Nbar)
  # cat <- ex_jk$cat
  # tot <- ex_jk$tot
  # kk <- prop(jknife = jknife, cat = cat, tot =tot, Nbar=Nbar)
  #
  #



  # -------------------------------------------------------------------------





})

