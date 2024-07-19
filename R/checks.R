
check_columns <- function(columns, df) {
  arg <- deparse(substitute(df))
  missing <- columns[!columns %in% colnames(df)]

  if(!rlang::is_empty(missing)) {
    cli::cli_abort(
      "Required columns were missing from the data.frame {.arg {arg}}:
        {missing}",
      
    ) 
  }
}

check_is_path <- function(path) {

  if(!rlang::is_scalar_character(path)) {
    stop(cli::format_error(
      "The provided argument is not a string"
    ))
  }
  # OBS: add more informative path
  if(!file.exists(path)) {
    stop(cli::format_error(
      "The filepath ({.path {path}}) does not exist on the host system"
    ))
  }

}



check_numeric_columns <- function(ld) {
  before <- colnames(ld)
  ld <- dplyr::select(ld, "SNP", dplyr::where(is.numeric))
  after <- colnames(ld)

  if(!all(before %in% after)) {
    cli::cli_alert_danger(
      "Columns were removed from ld.parquet because they were non-numeric:
        {.code {before[!before %in% after]}}
        "
    )
  }

  ld
}

