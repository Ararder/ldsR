utils::globalVariables(c("B","SE", "Z"))

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


#' Calculate the effective sample size of a case-control GWAS
#'
#' @param N_case number of cases
#' @param N_control Number of controls
#'
#' @returns a double
#' @export
#'
#' @examples
#' calc_effective_n(500, 1000)
calc_effective_n <- function(N_case, N_control) {

  4 / (((1 / N_case) + (1 / N_control)))
}


#' Parse GWAS format of `tidyGWAS::tidyGWAS()`
#'
#' @param tbl a [dplyr::tibble()]
#' @param n Column name of sample size, default is "N".
#'
#' @return a munged [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#' munged <- from_tidyGWAS("path/tidyGWAS/cleaned/tidyGWAS_hivestyle")
#' }
from_tidyGWAS <- function(tbl, n = c("N", "EffectiveN")) {
  n <- rlang::arg_match(n)

  if("data.frame" %in% class(tbl)) {

    tbl |>
      dplyr::mutate(Z = B/SE) |>
      dplyr::select(
        SNP = RSID,
        A1 = EffectAllele,
        A2 = OtherAllele,
        Z,
        N = {{ n }},
        dplyr::any_of(c("INFO", "EAF"))
      ) |>
      munge()

  } else if(rlang::is_scalar_vector(tbl)) {

    if(fs::path_file(tbl) != "tidyGWAS_hivestyle") {
      tbl <- fs::path(tbl, "tidyGWAS_hivestyle")
    }


    arrow::open_dataset(tbl)  |>
      dplyr::select(
        "SNP" = "RSID",
        "A1" = "EffectAllele",
        "A2" = "OtherAllele",
        "Z",
        "N" = {{ n }},
        dplyr::any_of(c("INFO", "EAF"))
      ) |>
      dplyr::collect() |>
      munge()

  }
}



#' Print names of annotations in a ldscore fileset
#' @param ldscore_dir A filepath to a ldscore directory
#'
#' @returns a character vector of column names
#' @export
#'
#' @examples \dontrun{
#' get_annot_names("ldscore/dir/celltypes")
#' }
get_annot_names <- function(ldscore_dir) {
  arrow::read_parquet(fs::path(ldscore_dir, "annot.parquet")) |> dplyr::pull(annot)
}
