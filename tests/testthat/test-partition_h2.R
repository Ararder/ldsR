
test_that("multiplication works", {
  skip("only works on local machine")

  annots <- c("baseL2", "Conserved_LindbladToh.bedL2")

  



  ldscore_dirs <- c("~/Downloads/ldsR_ldscores/baseline1.1","/Users/arvhar/Downloads/ldsR_ldscores/siletti_superclusters", "/Users/arvhar/Downloads/ldsR_ldscores/fetal_cCRE_union")
  subset_annots <- c("amygdala_excitatory", "baseL2", "Coding_UCSC.bedL2", "fetal_ATAC", "pooop")





  get_annot_names("/Users/arvhar/Downloads/ldsR_ldscores/siletti_clusters/")
  bl <- get_annot_names("~/Downloads/ldsR_ldscores/baseline1.1")
  subset_annots <- c(bl[-c(4,5)], "amygdala_excitatory", "phyloPam")
  sumstat <- arrow::read_tsv_arrow("~/Downloads/ldsR_ldscores/scz_example.tsv.gz")
  weights = NULL
  n_blocks=200
  k <- partition_h2(
    sumstat = sumstat,
    ldscore_dirs = c("~/Downloads/ldsR_ldscores/baseline1.1", "/Users/arvhar/Downloads/ldsR_ldscores/siletti_superclusters", "/Users/arvhar/Downloads/ldsR_ldscores/zoonomia_sullivan2023/"),
    subset_annots = subset_annots,
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
  df <- readr::read_tsv("~/Downloads/eo.sumstats.gz") |> 
    dplyr::filter(!is.na(Z))

  bl <- "~/Downloads/ldsR_ldscores/baseline1.1"
  bl_names <- get_annot_names(bl)
  roadmap <- "~/Downloads/ldsR_ldscores/roadmap2018/"
  kk <- get_annot_names(roadmap)



  
  

  res <- partition_h2(
    sumstat = df,
    ldscore_dirs = c(bl, roadmap),
    subset_annots = c(bl_names, "fetal_brain_female_d_nase")
  )

})

