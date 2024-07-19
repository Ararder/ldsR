



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


# 

test <- dplyr::tibble(SNP = "rs100", A1 = "T", A2 = "C", B = 0.05)
path <- "~/Downloads/ldsc.sumstats.gz"

else if("data.frame" %in% class(df)) {
  check_columns(c("SNP", "A1","A2","Z","N"), df)

  final <- dplyr::tibble(df) |> 
    tidyr::drop_na() |> 
    dplyr::semi_join(ref, by = "SNP")

  munge(final)
  
}


