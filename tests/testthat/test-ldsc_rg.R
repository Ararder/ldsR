test_that("multiplication works", {

    # read in data -----------------------------------------------------------



    sumstat <- arrow::read_parquet(test_path("fixtures/test_data.parquet"))
    sumstats1 <- dplyr::select(sumstat, SNP, Z = Z.x, N = N.x)
    sumstats2 <- dplyr::select(sumstat, SNP, Z = Z.y, N = N.y)

    # other test-suite
    sumstats1 <- arrow::read_tsv_arrow("~/Desktop/iPSYCH2012.sumstats.gz")
    sumstats2 <- arrow::read_tsv_arrow("~/Desktop/recur.sumstats.gz")
    res <- ldsc_rg(sumstats1, sumstats2)


    weights <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"))
    M <- 1173569
    sumstats1 <- dplyr::select(sumstats1, dplyr::all_of(req_cols))
    sumstats2 <- dplyr::select(sumstats2, dplyr::all_of(req_cols))
  
    m <- dplyr::inner_join(weights, sumstats1, by = "SNP") |>
      dplyr::inner_join(sumstats2, by = "SNP")
    n_blocks <- 200
    x <- as.matrix(m$L2)
  
    # -------------------------------------------------------------------------
  
    res1 <- ldscore(y = m$Z.x^2, x = x, w = m$L2, N = m$N.x, M = M)
    res2 <- ldscore(y = m$Z.y^2, x = x, w = m$L2, N = m$N.y, M = M)
    N <- sqrt(m$N.x * m$N.y)
    res3 <- gencov(
      z1 = m$Z.x,z2 = m$Z.y,w = m$L2, x = x, N1 = m$N.x, N2 = m$N.y, M = M,
      hsq1_tot = res1$tot, hsq2_tot = res2$tot, intercept_hsq1 = res1$int,
      intercept_hsq2 = res2$int
    )
  
  
    gencov_int <- res3$int
    gencov_int_se <- res3$int_se
    rg_ratio <- res3$tot / sqrt(res1$tot * res2$tot)
    numer <- tot_delete_values(res3, M, mean(N))
  
    d1_tot <- tot_delete_values(res1, M, mean(m$N.x))
    d2_tot <- tot_delete_values(res2, M, mean(m$N.y))
    denom <- sqrt(d1_tot*d2_tot)
    pseudovalues <- vector("numeric", n_blocks)
    for(j in 1:n_blocks){
      pseudovalues[j] <- n_blocks * rg_ratio - (n_blocks - 1) * numer[j] / denom[j]
    }
  

    final <- jackknife(matrix(pseudovalues))
    dplyr::tibble(
      rg = final$est, 
      rg_se = final$se, 
      h2_trait1 = res1$tot, 
      h2_trait_se = res1$tot_se, 
      h2_trait2 = res2$tot, 
      h2_trait2_se = res2$tot_se,
      gcov = res3$tot,
      gcov_se = res3$tot_se,
      gcov_int = res3$int,
      gcov_int_se = res3$int_se,
      mean_z1z2 = mean(m$Z.x*m$Z.y),
      z = rg /rg_se,
      p = stats::pchisq(z, df = 1, lower.tail = FALSE)
    )
  
  
  
})
