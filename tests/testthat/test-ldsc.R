test_that("ldsc_h2 runs", {

  skip()

  gwas <- arrow::read_parquet(test_path("fixtures/move/speed_scz2022_eur.parquet"), col_select = c("SNP", "Z","N")) |>
    dplyr::mutate(SNP = stringr::str_c("rs", as.integer(SNP)))

  tictoc::tic()
  ldsc_h2(gwas)
  tictoc::toc()


})



test_that("partitioned ldscores", {
  skip()
  gwas <- arrow::read_parquet("tests/testthat/fixtures/move/speed_scz2022_eur.parquet")
  weights <- arrow::read_parquet("tests/testthat/fixtures/move/weights2.parquet") |>
    dplyr::arrange(CHR, BP)
  ld <- arrow::read_parquet("tests/testthat/fixtures/move/ld.parquet")
  annot <- arrow::read_parquet("tests/testthat/fixtures/move/annot.parquet")


  merged <- dplyr::inner_join(weights, ld, by = "SNP") |>
    dplyr::inner_join(gwas, by = "SNP")
  remove <- c(colnames(weights), colnames(gwas))

  expect_no_error(res <- ldscore(
    y = merged$Z^2,
    x = dplyr::select(merged,-dplyr::all_of(remove)) |> as.matrix(),
    w = merged$L2,
    N = merged$N,
    M = annot$m50
    ))
  tidy_results(res)



})


