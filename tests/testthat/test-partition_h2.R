
test_that("multiplication works", {
  skip("only works on local machine")

  annots <- c("baseL2", "Conserved_LindbladToh.bedL2")





  ldscore_dirs <- c("~/Downloads/ldsR_ldscores/baseline1.1","/Users/arvhar/Downloads/ldsR_ldscores/siletti_superclusters", "/Users/arvhar/Downloads/ldsR_ldscores/fetal_cCRE_union")
  subset_annots <- c("amygdala_excitatory", "baseL2", "Coding_UCSC.bedL2", "fetal_ATAC", "pooop")






  bl <- get_annot_names("~/Downloads/ldsR_ldscores/baseline1.1")
  subset_annots <- c(bl[-c(4,5)], "amygdala_excitatory", "phyloPam")
  sumstat <- arrow::read_tsv_arrow("~/Downloads/ldsR_ldscores/scz_example.tsv.gz")
  k <- partition_h2(
    sumstat = sumstat,
    ldscore_dirs = c("~/Downloads/ldsR_ldscores/baseline1.1", "/Users/arvhar/Downloads/ldsR_ldscores/siletti_superclusters", "/Users/arvhar/Downloads/ldsR_ldscores/zoonomia_sullivan2023/"),
    subset_annots = bl,
    overlapping_annotations = TRUE
  )


    k <- partition_h2(
    sumstat = sumstat,
    ldscore_dirs = c("~/Downloads/ldsR_ldscores/baseline1.1", "/Users/arvhar/Downloads/ldsR_ldscores/siletti_superclusters", "/Users/arvhar/Downloads/ldsR_ldscores/zoonomia_sullivan2023/"),
    subset_annots = c(bl, "amygdala_excitatory"),
    overlapping_annotations = FALSE
  )

  k <- partition_h2(
    sumstat = sumstat,
    ldscore_dirs = c("~/Downloads/ldsR_ldscores/baseline1.1", "/Users/arvhar/Downloads/ldsR_ldscores/siletti_superclusters", "/Users/arvhar/Downloads/ldsR_ldscores/zoonomia_sullivan2023/","/Users/arvhar/Downloads/ldsR_ldscores/siletti_clusters"),
    subset_annots = c(bl, "phastCons", "amygdala_excitatory_405", "hippocampal_ca1_3_188"),
    overlapping_annotations = FALSE
  )



})

#
test_that("multiplication works", {
  skip("only works on local machine")

  bdir <- "~/Library/CloudStorage/OneDrive-KarolinskaInstitutet/ldsR_ldscores/"

  bl_names <- get_annot_names(paste0(bdir, "/baseline_model-v1.2"))
  yl_names <- get_annot_names(paste0(bdir, "/yangli2023_subclass"))

  sumstat <- arrow::read_tsv_arrow(paste0(bdir, "/scz_example.tsv.gz"))
  ldscore_dirs <- c(
    paste0(bdir, "/baseline_model-v1.2"),
    paste0(bdir, "/yangli2023_subclass"),
    paste0(bdir, "/zoonomia_sullivan2023")
  )
  subset_annots <- c(bl_names[-c(4:5)], "phastCons","phyloPam")
  weights = NULL
  n_blocks = 200
  overlapping_annotations=TRUE

  res <- purrr::map(yl_names, \(x) {
    partition_h2(
      sumstat,
      ldscore_dirs,
      c(subset_annots,x),
      overlapping_annotations = TRUE
    )

  }, .progress = list(type = "tasks"))

})


