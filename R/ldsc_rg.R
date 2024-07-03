
utils::globalVariables(c("annot", "m50", "SNP", "coef", "coef_se", "z", "L2", "rg", "rg_se"))
req_cols <- c("SNP", "Z", "N")
#' Compute the genetic correlation between two traits
#'
#' @param sumstats1 a [dplyr::tibble()] with columns h2 and h2_se
#' @param sumstats2 a [dplyr::tibble()] with columns h2 and h2_se
#' @param weights Optional, a data.frame or tbl with columns SNP and L2
#' @param M number of snps used to calculate LD in weights
#' @param n_blocks blocks for jackknife
#'
#' @return a [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#' ldsc_rg(s1, s2)
#' }
ldsc_rg <- function(sumstats1, sumstats2, weights=NULL, M=NULL, n_blocks=200) {
  stopifnot("SNP, Z and N are required in `sumstats1`" = all(req_cols %in% colnames(sumstats1)))
  stopifnot("SNP, Z and N are required in `sumstats2`" = all(req_cols %in% colnames(sumstats2)))

  weights <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"))
  M <- 1173569
  sumstats1 <- dplyr::select(sumstats1, dplyr::all_of(req_cols))
  sumstats2 <- dplyr::select(sumstats2, dplyr::all_of(req_cols))

  m <- dplyr::inner_join(weights, sumstats1, by = "SNP") |>
    dplyr::inner_join(sumstats2, by = "SNP")


  # -------------------------------------------------------------------------

  x <- as.matrix(m$L2)
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


}




# CELLTYPE H2 -------------------------------------------------------------





