
utils::globalVariables(c("annot", "m50", "SNP", "coef", "coef_se", "z", "L2", "rg", "rg_se",
"A1.x","A1.y", "A2.x", "A2.y", "Z.y", "N.x", "N.y"))
#' Compute the genetic correlation between two traits
#'
#' @param sumstats1 a [dplyr::tibble()] with atleast the columns SNP, A1, A2, Z, N.
#' To perform quality checks, use [munge()] before running `ldsc_rg()`
#' @param sumstats2 a [dplyr::tibble()] with atleast the columns SNP, A1, A2, Z, N.
#' To estimate rG with between sumstats1 and several traits, wrap the data.frames in a list.
#' If the names are lost, the column `trait2` will correspond to the index of the list.
#' Use a named list the retain the names in the `trait2` column
#' `list(dataframe1, dataframe2, ...)`
#' @inheritParams ldsc_h2
#' 
#' 
#'
#' @return a [dplyr::tibble()]
#' @export
#'
#' @examples
#' path <- system.file("extdata", "eur_w_ld.parquet", package = "ldsR")
#' snps <- arrow::read_parquet(path, col_select = c("SNP"))
#' # to make example faster
#' snps <- dplyr::slice_head(snps, n = 100000)
#' snps$A1 <- "A"
#' snps$A2 <- 
#' snps$N <- 130000
#' snps$Z <- rnorm(nrow(snps))
#' snps2 <- snps
#' snps2$N <- 75000
#' snps2$Z <- rnorm(nrow(snps))
#' ldsc_rg(snps, snps2)
#' # to run estimate the genetic correlations for many traits, wrap s2 in a list
#' # ldsc_rg(snps, list(snps2, snps2))
#' # use a named list to create the `trait2` column in the output
#' # ldsc_rg(snps, list("trait2" = s2, "trait3" = s3, "trait4" = s4))
ldsc_rg <- function(sumstats1, sumstats2, weights=NULL, M=NULL, n_blocks=200) {
  req_cols <- c("SNP", "Z", "N", "A1", "A2")
  check_columns(req_cols, sumstats1)
  
  if(is.null(weights)) {
    weights <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"), col_select = c("SNP", "L2"))
    M <- 1173569
  } else {
    stopifnot("`weights` must be a data.frane with columns `SNP` and `L2`" = "data.frame" %in% class(weights))
    check_columns(c("SNP", "L2"), weights)
    stopifnot("To use custom weights, you must also pass `M`" = !is.null(M))
  }
  
  sumstats1 <- dplyr::select(sumstats1, dplyr::all_of(req_cols))
  
  
  
  # -------------------------------------------------------------------------
  if("data.frame" %in% class(sumstats2)) sumstats2 <- list(sumstats2)
  

  purrr::map(sumstats2, \(x) .rg(sumstats1, sumstats2=x, M=M, weights=weights, n_blocks = n_blocks),
    .progress = list(type = "tasks", name = "Computing genetic correlations...")
  ) |> 
    purrr::list_rbind(names_to = "trait2")


}



.rg <- function(sumstats1, sumstats2, M, weights, n_blocks) {
  req_cols <- c("SNP", "Z", "N", "A1", "A2")

  # check that all columns exist for sumstats2
  check_columns(req_cols, sumstats2)
  sumstats2 <- dplyr::select(sumstats2, dplyr::all_of(req_cols)) |> 
    dplyr::filter(SNP %in% weights$SNP)
  
  before <- nrow(sumstats1)
  # sumstats1 and 2 are merged and Z-score is standardized on alleles in merge_sumstats
  m <- dplyr::inner_join(weights, merge_sumstats(sumstats1, sumstats2), by = "SNP") |> 
    dplyr::mutate(
      N.x = as.double(N.x),
      N.y = as.double(N.y),
    )
  cli::cli_alert_warning("{before - nrow(m)} SNPs were removed when merging summary statistics")


  # start rg calculation ---------------------------------------------------
  
  N <- sqrt(m$N.x * m$N.y)
  x <- as.matrix(m$L2)

  # first get ldscore regression results for each trait
  res1 <- ldscore(y = m$Z.x^2, x = x, w = m$L2, N = m$N.x, M = M)
  res2 <- ldscore(y = m$Z.y^2, x = x, w = m$L2, N = m$N.y, M = M)
  
  # then for product of z-scores
  res3 <- gencov(
    z1 = m$Z.x,z2 = m$Z.y,w = m$L2, x = x, N1 = m$N.x, N2 = m$N.y, M = M,
    hsq1_tot = res1$tot, hsq2_tot = res2$tot, intercept_hsq1 = res1$int,
    intercept_hsq2 = res2$int
  )

  # get intercept and rg
  gencov_int <- res3$int
  gencov_int <- res3$int
  gencov_int_se <- res3$int_se
  rg_ratio <- res3$tot / sqrt(res1$tot * res2$tot)

  # get standard error by computing ratio jackknifes
  numer <- tot_delete_values(res3, M, mean(N))
  denom <- sqrt(tot_delete_values(res1, M, mean(m$N.x)) * tot_delete_values(res2, M, mean(m$N.y)))
  pseudovalues <- vector("numeric", n_blocks)
  for(j in 1:n_blocks) pseudovalues[j] <- n_blocks * rg_ratio - (n_blocks - 1) * numer[j] / denom[j]
  final <- jackknife(matrix(pseudovalues))


  # report main results ----------------------------------------------------


  out <- dplyr::tibble(
    rg = rg_ratio, 
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

  out
}

merge_sumstats <- function(sumstats1, sumstats2) {
  
  m <- dplyr::inner_join(sumstats1, sumstats2, by = "SNP")
  
  
  dplyr::filter(m, A1.x == A1.y & A2.x == A2.y | A1.x == A2.y & A2.x == A1.y) |> 
    # flip Z if alleles are flipped
    dplyr::mutate(Z.y = dplyr::if_else(A1.x == A1.y, Z.y, Z.y*-1)) |> 
    dplyr::select(dplyr::all_of(c("SNP", "Z.x", "Z.y", "N.x", "N.y")))
  
  

  # 
}

