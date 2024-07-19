



test_that("parse_gwas", {
  skip()
  expect_no_error(
    parse_gwas(
      test_path("fixtures/test_ldsc.sumstats.gz"),
      format = "ldsc"
    )
  )
    
    
  expect_error(
    parse_gwas(
      test_path("fixtures/test_ldsc.sumstats.gz"),
      format = "dataframe"
    )
  )





})

test_that("munge works", {
  df <- arrow::read_tsv_arrow(test_path("fixtures/test_ldsc.sumstats.gz"))

  
  expect_no_error(munge(tidyr::drop_na(df)))
    

})


# # 

# test <- dplyr::tibble(SNP = "rs100", A1 = "T", A2 = "C", B = 0.05)
# path <- "~/Downloads/ldsc.sumstats.gz"

# else if("data.frame" %in% class(df)) {
#   check_columns(c("SNP", "A1","A2","Z","N"), df)

#   final <- dplyr::tibble(df) |> 
#     tidyr::drop_na() |> 
#     dplyr::semi_join(ref, by = "SNP")

#   munge(final)
  
# }




test_that("ldsc_to_parquet works", {
  skip()
  snps <- arrow::read_parquet("tests/testthat/fixtures/baseline_v1.1/ld.parquet") |> 
    dplyr::slice_head(n = 100000)
  write_parquet("")


  ld <- arrow::read_parquet(fs::path(outdir, "ld.parquet"))
  arrow::read_parquet(fs::path("annot_ref.parquet"))
  annot <- arrow::read_parquet(fs::path(outdir, "annot.parquet"))
  

  # ds1 <- arrow::read_parquet("/Users/arvhar/projects/move/clusters/annot.parquet")
  # df <- arrow::open_dataset("/Users/arvhar/projects/move/clusters/ld.parquet") |> dplyr::collect()
  # m50 <- ds1$m50
  
  # for(idx in 2:ncol(df)) {
  #   name <- colnames(df)[idx]
  #   attr(df[[name]], "m50") <- m50[idx-1]
  # }
  
  # arrow::write_parquet(df, "/Users/arvhar/projects/move/superclusters/clusters.parquet")
  # l <- arrow::open_dataset("test.parquet")

  # purrr::map_dbl(colnames(l)[-1], \(x) l[["metadata"]][["r"]][["columns"]][[x]][["attributes"]][["m50"]])



  
  


  
})