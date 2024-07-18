



test_that("parse_gwas", {
  
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
