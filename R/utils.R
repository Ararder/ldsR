#' Convert an estimate of observed-scale heritability to liability scale heritability
#'
#' @param obs_h2 observed-scale heritability
#' @param pop_prev prevalence of the disorder in the general population
#' @param sample_prev the prevalence of the disorder in the sample.
#'  Default value is 0.5, reflecting a case-control study using effective N as sample size
#'
#' @return a double
#' @export
#'
#' @examples
#' liability_h2(0.25, 0.02)
liability_h2 <- function(obs_h2, pop_prev, sample_prev = 0.5) {
  K <- pop_prev
  P <- sample_prev
  zv <- stats::dnorm(stats::qnorm(K))

  obs_h2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2

}
