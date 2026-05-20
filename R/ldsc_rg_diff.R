utils::globalVariables(c("L2",
                         "A1_A", "A2_A", "A1_B", "A2_B", "A1_C", "A2_C", "A1_D", "A2_D",
                         "Z_A", "Z_B", "Z_C", "Z_D",
                         "N_A", "N_B", "N_C", "N_D"))


#' Test whether two genetic correlations differ
#'
#' @description
#' Block jackknife test of the null hypothesis that the genetic correlation
#' between traits A and B equals the genetic correlation between traits C and D.
#' A trait may appear in both pairs (e.g. rG(A,B) vs rG(C,B)) by being passed
#' twice.
#'
#' All summary statistics are merged with the weights to a single common SNP
#' set before either rG is computed, so the jackknife blocks correspond
#' one-to-one across the two correlations. This captures the covariance
#' between the two rGs, which is non-zero whenever a trait (or just SNPs)
#' is shared.
#'
#' For each block `i`, delete values `R(A,B)_i` and `R(C,D)_i` for the rG
#' are differenced to form `D_i`. Pseudovalues are
#' `P_i = n * d - (n - 1) * D_i`, where `d = rG(A,B) - rG(C,D)` is the
#' global difference and `n` is the number of blocks. The test statistic is
#' `z = mean(P) / sqrt(var(P) / n)`.
#'
#' @param pair1 A list of two summary statistics tibbles defining the first
#'   genetic correlation. Each tibble must contain `SNP`, `A1`, `A2`, `Z`, `N`.
#' @param pair2 Same shape as `pair1`. May share a tibble with `pair1` (pass
#'   the same data frame in both pairs to share a trait).
#' @inheritParams ldsc_rg
#'
#' @return a [dplyr::tibble()] with the two global rGs, the global difference,
#'   the jackknife-corrected difference, its standard error, and the z
#'   statistic with chi-squared p-value.
#' @export
#'
#' @examples
#' \dontrun{
#' # Four distinct traits
#' ldsc_rg_diff(
#'   pair1 = list(trait_a, trait_b),
#'   pair2 = list(trait_c, trait_d)
#' )
#'
#' # Three traits: rG(A,B) vs rG(C,B), with B shared
#' ldsc_rg_diff(
#'   pair1 = list(trait_a, trait_b),
#'   pair2 = list(trait_c, trait_b)
#' )
#' }
ldsc_rg_diff <- function(pair1, pair2, weights = NULL, M = NULL, n_blocks = 200) {
  stopifnot("`pair1` must be a list of two summary statistics data frames" =
              is.list(pair1) && length(pair1) == 2)
  stopifnot("`pair2` must be a list of two summary statistics data frames" =
              is.list(pair2) && length(pair2) == 2)

  req_cols <- c("SNP", "A1", "A2", "Z", "N")
  purrr::walk(c(pair1, pair2), \(s) check_columns(req_cols, s))

  if(is.null(weights)) {
    weights <- arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"), col_select = c("SNP", "L2"))
    M <- 1173569
  } else {
    stopifnot("`weights` must be a data.frane with columns `SNP` and `L2`" = "data.frame" %in% class(weights))
    check_columns(c("SNP", "L2"), weights)
    stopifnot("To use custom weights, you must also pass `M`" = !is.null(M))
  }


  # merge all four sumstats + weights onto a common SNP set ------------------
  before <- nrow(weights)
  m <- merge_four_sumstats(pair1[[1]], pair1[[2]], pair2[[1]], pair2[[2]], weights)
  cli::cli_alert_warning("{before - nrow(m)} SNPs were removed when merging summary statistics")


  # compute rG and per-block delete values for each pair --------------------
  ab <- rg_delete_values(m, "A", "B", M)
  cd <- rg_delete_values(m, "C", "D", M)


  # block jackknife test for the difference ---------------------------------
  d <- ab$rg - cd$rg
  D <- ab$R - cd$R
  P <- n_blocks * d - (n_blocks - 1) * D
  m_pseudo <- mean(P)
  v_pseudo <- sum((P - m_pseudo)^2) / (n_blocks - 1)
  diff_se <- sqrt(v_pseudo / n_blocks)
  z <- m_pseudo / diff_se

  dplyr::tibble(
    rg1 = ab$rg,
    rg2 = cd$rg,
    diff = d,
    diff_jackknife = m_pseudo,
    diff_se = diff_se,
    z = z,
    p = stats::pchisq(z^2, df = 1, lower.tail = FALSE)
  )
}



rg_delete_values <- function(m, t1, t2, M) {
  N1 <- as.double(m[[paste0("N_", t1)]])
  N2 <- as.double(m[[paste0("N_", t2)]])
  Z1 <- m[[paste0("Z_", t1)]]
  Z2 <- m[[paste0("Z_", t2)]]
  x <- as.matrix(m$L2)

  res1 <- ldscore(y = Z1^2, x = x, w = m$L2, N = N1, M = M)
  res2 <- ldscore(y = Z2^2, x = x, w = m$L2, N = N2, M = M)
  res3 <- gencov(
    z1 = Z1, z2 = Z2, w = m$L2, x = x, N1 = N1, N2 = N2, M = M,
    hsq1_tot = res1$tot, hsq2_tot = res2$tot,
    intercept_hsq1 = res1$int, intercept_hsq2 = res2$int
  )

  N <- sqrt(N1 * N2)
  numer <- tot_delete_values(res3, M, mean(N))
  denom <- sqrt(tot_delete_values(res1, M, mean(N1)) * tot_delete_values(res2, M, mean(N2)))

  list(
    rg = res3$tot / sqrt(res1$tot * res2$tot),
    R = numer / denom
  )
}



merge_four_sumstats <- function(sA, sB, sC, sD, weights) {
  req <- c("SNP", "A1", "A2", "Z", "N")
  prep <- function(df, suffix) {
    df |>
      dplyr::select(dplyr::all_of(req)) |>
      dplyr::rename_with(\(x) paste0(x, "_", suffix), -SNP)
  }

  m <- weights |>
    dplyr::inner_join(prep(sA, "A"), by = "SNP") |>
    dplyr::inner_join(prep(sB, "B"), by = "SNP") |>
    dplyr::inner_join(prep(sC, "C"), by = "SNP") |>
    dplyr::inner_join(prep(sD, "D"), by = "SNP")

  match_or_flip <- function(a1_ref, a2_ref, a1, a2) {
    (a1_ref == a1 & a2_ref == a2) | (a1_ref == a2 & a2_ref == a1)
  }
  keep <- match_or_flip(m$A1_A, m$A2_A, m$A1_B, m$A2_B) &
          match_or_flip(m$A1_A, m$A2_A, m$A1_C, m$A2_C) &
          match_or_flip(m$A1_A, m$A2_A, m$A1_D, m$A2_D)
  m <- m[keep, ]

  for(nm in c("B", "C", "D")) {
    flip <- m$A1_A == m[[paste0("A2_", nm)]] & m$A2_A == m[[paste0("A1_", nm)]]
    m[[paste0("Z_", nm)]] <- dplyr::if_else(flip, -m[[paste0("Z_", nm)]], m[[paste0("Z_", nm)]])
  }

  dplyr::select(m,
                SNP, L2,
                Z_A, Z_B, Z_C, Z_D,
                N_A, N_B, N_C, N_D)
}
