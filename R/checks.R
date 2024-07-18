
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
