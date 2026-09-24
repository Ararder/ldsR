rsid_map <- arrow::read_parquet(system.file("extdata", "rsid_map.parquet", package = "ldsR"))
sumstats <- arrow::read_parquet(system.file("extdata", "sumstats.parquet", package = "ldsR")) |>
  dplyr::select(SNP, Z = Z.x, N = N.x)

# rsid-based sumstats with alleles, plus the same rows keyed only on coordinates
with_alleles <- dplyr::inner_join(
  sumstats,
  dplyr::distinct(rsid_map, SNP, .keep_all = TRUE),
  by = "SNP"
) |>
  dplyr::rename(A1 = ALT, A2 = REF)

coords <- function(build) {
  dplyr::select(with_alleles, CHR, POS = dplyr::all_of(paste0("POS", build)), A1, A2, Z, N)
}


test_that("to_rsid recovers RSIDs on both builds", {
  for(b in c("37", "38")) {
    res <- to_rsid(coords(b), build = b)
    expect_equal(nrow(res), nrow(with_alleles))
    expect_equal(res$SNP, with_alleles$SNP)
    expect_equal(colnames(res), c("SNP", "CHR", "POS", "A1", "A2", "Z", "N"))
  }
})


test_that("to_rsid guesses the build", {
  expect_message(res37 <- to_rsid(coords("37")), "Using GRCh37")
  expect_message(res38 <- to_rsid(coords("38")), "Using GRCh38")
  expect_equal(res37$SNP, with_alleles$SNP)
  expect_equal(res38$SNP, with_alleles$SNP)
})


test_that("wrong build finds few matches", {
  res <- to_rsid(coords("37"), build = "38")
  expect_lt(nrow(res), 0.01 * nrow(with_alleles))
})


test_that("to_rsid matches alleles regardless of order and case", {
  dset <- coords("38") |>
    dplyr::slice_head(n = 1000) |>
    dplyr::mutate(tmp = A1, A1 = tolower(A2), A2 = tmp, .keep = "unused")
  res <- to_rsid(dset, build = "38")
  expect_equal(res$SNP, with_alleles$SNP[1:1000])
})


test_that("to_rsid drops rows with non-matching alleles", {
  dset <- dplyr::slice_head(coords("38"), n = 100)
  dset$A1[1:10] <- "I"
  dset$A2[1:10] <- "D"
  res <- to_rsid(dset, build = "38")
  expect_equal(res$SNP, with_alleles$SNP[11:100])
})


test_that("to_rsid handles multi-allelic reference sites", {
  multi <- rsid_map |>
    dplyr::add_count(SNP) |>
    dplyr::filter(n > 1) |>
    dplyr::slice_head(n = 2)
  dset <- dplyr::select(multi, CHR, POS = POS38, A1 = ALT, A2 = REF)
  res <- to_rsid(dset, build = "38")
  expect_equal(nrow(res), 2)
  expect_equal(res$SNP, multi$SNP)
})


test_that("to_rsid accepts chr-prefixed and character chromosomes", {
  dset <- dplyr::slice_head(coords("38"), n = 100) |>
    dplyr::mutate(CHR = paste0("chr", CHR))
  dset$CHR[1] <- "chrX"
  res <- to_rsid(dset, build = "38")
  expect_equal(res$SNP, with_alleles$SNP[2:100])
})


test_that("to_rsid overwrites an existing SNP column", {
  dset <- dplyr::slice_head(coords("38"), n = 10) |>
    dplyr::mutate(SNP = "1:1")
  expect_message(res <- to_rsid(dset, build = "38"), "Overwriting")
  expect_equal(res$SNP, with_alleles$SNP[1:10])
})


test_that("to_rsid errors on missing columns, bad build and no matches", {
  expect_error(to_rsid(dplyr::select(coords("38"), -POS)), "POS")
  expect_error(to_rsid(coords("38"), build = "19"))
  expect_error(
    to_rsid(dplyr::tibble(CHR = 1L, POS = 1L, A1 = "A", A2 = "G")),
    "No variants"
  )
})


test_that("to_rsid |> munge |> ldsc_h2 reproduces the RSID-based estimate", {
  expected <- with_alleles |>
    dplyr::select(SNP, A1, A2, Z, N) |>
    munge() |>
    ldsc_h2()

  res <- coords("37") |>
    to_rsid() |>
    munge() |>
    ldsc_h2()

  expect_equal(res, expected)
})
