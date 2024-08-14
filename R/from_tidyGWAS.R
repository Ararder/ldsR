utils::globalVariables(c("B","SE", "Z"))
#' Parse GWAS format of [tidyGWAS::tidyGWAS()]
#'
#' @param tbl a [dplyr::tibble()]
#' @param n Column name of sample s
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
