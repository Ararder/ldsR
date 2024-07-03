test_that("parse_parquet_dir works", {
  l <- parse_parquet_dir("~/projects/move")
  l <- parse_parquet_dir("~/projects/move/superclusters")
})

path <- "/work/users/a/r/arvhar/tidyGWAS_stuff/output2/scz2022_eur_noSWE/tidyGWAS_hivestyle/"
path2 <- "/work/users/a/r/arvhar/tidyGWAS_stuff/output2/scz2022_eur/tidyGWAS_hivestyle/"
# path <- "~/shared/gwas_sumstats/updated_sumstats/scz2022_core/tidyGWAS_hivestyle/"


test_that("munge works", {
  dset <- arrow::open_dataset(path)
  tictoc::tic()
  munged <- ldsc_munge(dset)
  tictoc::toc()
  munged <- ldsc_munge(dset)
  ldsc_h2(dplyr::mutate(munged, N = as.double(N)))


})
