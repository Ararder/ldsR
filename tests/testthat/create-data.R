



df <- arrow::read_parquet(system.file("extdata", "sumstats.parquet", package = "ldsR"))

ld <- arrow::read_parquet("~/Downloads/ldsR_ldscores/siletti_superclusters/ld.parquet") |>
  semi_join(df, by = "SNP")

annot_ref <- arrow::read_parquet("~/Downloads/ldsR_ldscores/siletti_superclusters/annot_ref.parquet") |>
  semi_join(df, by = "SNP")

annot <- arrow::read_parquet("~/Downloads/ldsR_ldscores/siletti_superclusters/annot.parquet")



arrow::write_parquet(
  ld,
  "inst/extdata/superclusters/ld.parquet"
)
arrow::write_parquet(
  annot_ref,
  "inst/extdata/superclusters/annot_ref.parquet"
)
arrow::write_parquet(
  annot,
  "inst/extdata/superclusters/annot.parquet"
)

# now for baseline

ld <- arrow::read_parquet("~/Downloads/ldsR_ldscores/baseline_model-v1.2/ld.parquet") |>
  dplyr::semi_join(df, by = "SNP")

annot_ref <- arrow::read_parquet("~/Downloads/ldsR_ldscores/baseline_model-v1.2/annot_ref.parquet") |>
  dplyr::semi_join(df, by = "SNP")

annot <- arrow::read_parquet("~/Downloads/ldsR_ldscores/baseline_model-v1.2/annot.parquet")


dir.create("inst/extdata/baseline1.2")
arrow::write_parquet(
  ld,
  "inst/extdata/baseline1.2/ld.parquet"
)
arrow::write_parquet(
  annot_ref,
  "inst/extdata/baseline1.2/annot_ref.parquet"
)
arrow::write_parquet(
  annot,
  "inst/extdata/baseline1.2/annot.parquet"
)










snps |>
  arrow::write_parquet(test_path("testdata/baseline/ld.parquet"))
fs::file_copy("tests/testthat/fixtures/baseline_v1.1/annot.parquet", test_path("testdata/baseline/annot.parquet"))




# celll-types ------------------------------------------------------------
snps <- arrow::read_parquet("~/projects/move/superclusters/ld.parquet") |>
    dplyr::slice_head(n = 100000)

arrow::write_parquet(snps, test_path("testdata/superclusters/ld.parquet"))




ld_baseline <- arrow::read_parquet("tests/testthat/fixtures/baseline_v1.1/ld.parquet") |>
  dplyr::select(1,2,3,5, "H3K4me1_peaks_Trynka.bedL2") |>
  dplyr::slice_head(n = 100000)
arrow::write_parquet(ld_baseline, "inst/extdata/baseline1.1_test/ld.parquet")


annot <- arrow::read_parquet("tests/testthat/fixtures/baseline_v1.1/annot.parquet") |>
  dplyr::filter(annot %in% colnames(ld_baseline))

arrow::write_parquet(annot, "inst/extdata/baseline1.1_test/annot.parquet")

s1 <- arrow::read_parquet("inst/extdata/sumstats.parquet") |>
  dplyr::select(SNP, Z = Z.x, N = N.x)


partitioned_h2(
  s1,
  ldscore_dir = "inst/extdata/baseline1.1_test"
)



# use siletti ldscores to inst/extdata
ld_baseline <- arrow::read_parquet("tests/testthat/testdata/superclusters/ld.parquet") |>
  dplyr::select(1:5)
annot <- arrow::read_parquet("tests/testthat/testdata/superclusters/annot.parquet") |>
  dplyr::filter(annot %in% colnames(ld_baseline))

arrow::write_parquet(ld_baseline, "inst/extdata/superclusters/ld.parquet")
arrow::write_parquet(annot, "inst/extdata/superclusters/annot.parquet")
