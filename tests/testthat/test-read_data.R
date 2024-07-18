# test_that("parse_parquet_dir works", {
#   l <- parse_parquet_dir("~/projects/move")
#   l <- parse_parquet_dir("~/projects/move/superclusters")
# })
#
# path <- "/work/users/a/r/arvhar/tidyGWAS_stuff/output2/scz2022_eur_noSWE/tidyGWAS_hivestyle/"
# path2 <- "/work/users/a/r/arvhar/tidyGWAS_stuff/output2/scz2022_eur/tidyGWAS_hivestyle/"
# # path <- "~/shared/gwas_sumstats/updated_sumstats/scz2022_core/tidyGWAS_hivestyle/"
#
#
# test_that("munge works", {
#   dset <- arrow::open_dataset(path)
#   tictoc::tic()
#   munged <- ldsc_munge(dset)
#   tictoc::toc()
#   munged <- ldsc_munge(dset)
#   ldsc_h2(dplyr::mutate(munged, N = as.double(N)))
#
#
# })
#
#
#
# # Parse GWAS handles three distinct situations
# # LDSC format: sumstats have already been munged with LDSC
# # in-memory data-frame: Req columns are SNP, Z, N. Optional are : A1, A2, INFO, MaF
#
# test_that("parse_gwas", {
#   path <- "~/Downloads/ldsc.sumstats.gz"
#   path2 <- test_path("fixtures/test_data.parquet")
#   dset <- arrow::read_tsv_arrow(path) |>
#     tidyr::drop_na() |>
#     # drop multi-allelic SNPs
#     dplyr::distinct(RSID, .keep_all = TRUE)
#     step1 <- dplyr::filter(dset, .data[["SNP"]] %in% ref$SNP)
#
#
#     ref <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"), col_select = c("SNP"))
#
#
#
#
#
# })
