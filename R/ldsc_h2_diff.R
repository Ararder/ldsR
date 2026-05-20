utils::globalVariables(c("L2", "Z_1", "Z_2", "N_1", "N_2"))


#' Test whether two SNP heritabilities differ
#'
#' @description
#' Block jackknife test of the null hypothesis that h²(trait 1) equals
#' h²(trait 2).
#'
#' Both summary statistics are merged with the weights to a single common SNP
#' set before either h² is computed, so the jackknife blocks correspond
#' one-to-one across the two traits. The per-block h² delete values from each
#' trait are differenced; pseudovalues
#' `P_i = n * d - (n - 1) * D_i` with `d = h²_1 - h²_2` give a jackknife
#' estimate of the difference and its variance, and `z = mean(P) / sqrt(var(P)/n)`.
#'
#' If `pop_prev1` (and optionally `pop_prev2`) is supplied, the corresponding
#' h² and per-block delete values are converted to liability scale via
#' [liability_h2()] before differencing. Comparing observed-scale h² for one
#' trait against liability-scale h² for the other is allowed but rarely makes
#' sense.
#'
#' @param sumstat1,sumstat2 Summary statistics tibbles with at least
#'   `SNP`, `Z`, `N`.
#' @param pop_prev1,sample_prev1 Optional liability-scale parameters for
#'   `sumstat1`. Leave `pop_prev1 = NULL` to test on the observed scale.
#' @param pop_prev2,sample_prev2 Same, for `sumstat2`.
#' @param return_blocks If `TRUE`, return the per-block h² delete values for both traits
#' @inheritParams ldsc_h2
#'
#' @return a [dplyr::tibble()] with the two h² estimates (on the chosen
#'   scale), the global difference, the jackknife-corrected difference, its
#'   standard error, and the z statistic with chi-squared p-value.
#' @export
#'
#' @examples
#' \dontrun{
#' # Observed-scale comparison
#' ldsc_h2_diff(trait_a, trait_b)
#'
#' # Liability-scale comparison of two binary traits
#' ldsc_h2_diff(
#'   trait_a, trait_b,
#'   pop_prev1 = 0.01, sample_prev1 = 0.5,
#'   pop_prev2 = 0.05, sample_prev2 = 0.5
#' )
#' }
ldsc_h2_diff <- function(sumstat1, sumstat2,
                         return_blocks = FALSE,
                         pop_prev1 = NULL, sample_prev1 = 0.5,
                         pop_prev2 = NULL, sample_prev2 = 0.5,
                         weights = NULL, M = NULL, n_blocks = 200) {
  req_cols <- c("SNP", "Z", "N")
  stopifnot("sumstat1 has to be a data.frame or tbl" = "data.frame" %in% class(sumstat1))
  stopifnot("sumstat2 has to be a data.frame or tbl" = "data.frame" %in% class(sumstat2))
  check_columns(req_cols, sumstat1)
  check_columns(req_cols, sumstat2)

  if(is.null(weights)) {
    weights <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"), col_select = c("SNP", "L2"))
    M <- 1173569
  } else {
    stopifnot("`weights` must be a data.frane with columns `SNP` and `L2`" = "data.frame" %in% class(weights))
    check_columns(c("SNP", "L2"), weights)
    stopifnot("To use custom weights, you must also pass `M`" = !is.null(M))
  }


  # merge both sumstats + weights onto a common SNP set ---------------------
  prep <- function(df, suffix) {
    df |>
      dplyr::select(dplyr::all_of(req_cols)) |>
      dplyr::rename_with(\(x) paste0(x, "_", suffix), -SNP)
  }
  before <- nrow(weights)
  m <- weights |>
    dplyr::inner_join(prep(sumstat1, "1"), by = "SNP") |>
    dplyr::inner_join(prep(sumstat2, "2"), by = "SNP")
  cli::cli_alert_warning("{before - nrow(m)} SNPs were removed when merging summary statistics")


  # run ldscore regression for each trait on the common set -----------------
  x <- as.matrix(m$L2)
  res1 <- ldscore(y = m$Z_1^2, x = x, w = m$L2, N = m$N_1, M = M, n_blocks = n_blocks)
  res2 <- ldscore(y = m$Z_2^2, x = x, w = m$L2, N = m$N_2, M = M, n_blocks = n_blocks)


  # per-block h² delete values (observed scale) -----------------------------
  H1 <- tot_delete_values(res1, M, mean(m$N_1))
  H2 <- tot_delete_values(res2, M, mean(m$N_2))
  h2_1 <- res1$tot
  h2_2 <- res2$tot


  # optional liability-scale conversion (linear, so applies to delete values too) ----
  if(!is.null(pop_prev1)) {
    H1   <- liability_h2(H1,   pop_prev1, sample_prev1)
    h2_1 <- liability_h2(h2_1, pop_prev1, sample_prev1)
  }
  if(!is.null(pop_prev2)) {
    H2   <- liability_h2(H2,   pop_prev2, sample_prev2)
    h2_2 <- liability_h2(h2_2, pop_prev2, sample_prev2)
  }


  # block jackknife test for the difference ---------------------------------
  d <- h2_1 - h2_2
  D <- H1 - H2
  P <- n_blocks * d - (n_blocks - 1) * D
  m_pseudo <- mean(P)
  v_pseudo <- sum((P - m_pseudo)^2) / (n_blocks - 1)
  diff_se <- sqrt(v_pseudo / n_blocks)
  z <- m_pseudo / diff_se

  if(isTRUE(return_blocks)) {
    return(list(
      H1 = H1,
      H2 = H2,
      h2_1 = h2_1,
      h2_2 = h2_2,
      n_blocks = n_blocks
    ))

  }

  dplyr::tibble(
    h2_1 = h2_1,
    h2_2 = h2_2,
    diff = d,
    diff_jackknife = m_pseudo,
    diff_se = diff_se,
    z = z,
    p = stats::pchisq(z^2, df = 1, lower.tail = FALSE)
  )
}
