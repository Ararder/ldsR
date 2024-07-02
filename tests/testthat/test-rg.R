test_that("rg_update_function runs", {
  expect_no_error(
    rg_update_func(
      c(0.05, 0.02),
      ld = matrix(rnorm(100), ncol = 2),
      w_ld = rnorm(100),
      N1 = rnorm(100, 100),
      N2 = rnorm(100, 100),
      M = 100,
      Nbar = 100,
      hsq1 = 0.1,
      hsq2 = 0.2,
      intercept_hsq1 = 0.1,
      intercept_hsq2 = 0.2
    )
  )

  expect_no_error(
    rg_update_func(
      c(0.05, 0.02),
      ld = matrix(rnorm(100), ncol = 2),
      w_ld = rnorm(100),
      N1 = rnorm(100, 100),
      N2 = rnorm(100, 100),
      M = 100,
      intercept =0.1,
      Nbar = 100,
      hsq1 = 0.1,
      hsq2 = 0.2,
      intercept_hsq1 = 0.1,
      intercept_hsq2 = 0.2
    )
  )
})



# test_that("Ratio jackknife works", {
#   list <- list(
#     "rg_ratio" = rg_ratio,
#     "gencov_tot_delete_values" = gencov_tot_delete_values,
#     "denom" = denom
#     )
#   ll <- readr::read_rds("~/ratio_jknife_tdata.rds")
#   est <- ll$rg_ratio
#   denom <- ll$denom
#   numer <- ll$gencov_tot_delete_values
#   n_blocks <- 200
#   pseudovalues <- vector("numeric", n_blocks)
#   for(j in 1:n_blocks){
#     pseudovalues[j] <- n_blocks * est - (n_blocks - 1) * numer[j] / denom[j]
#     # print(pseudovalues[j])
#   }
#
#   final <- jackknife(matrix(pseudovalues))
#
#
#
# })
