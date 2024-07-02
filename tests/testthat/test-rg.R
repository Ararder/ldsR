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
