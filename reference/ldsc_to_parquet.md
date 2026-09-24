# Transform a directory of LDscores to parquet format

Transform a directory of LDscores to parquet format

## Usage

``` r
ldsc_to_parquet(dir, thin = FALSE)
```

## Arguments

- dir:

  directory containing the python LDSC ldscore format

- thin:

  If the –thin flag has been used, provide a character vector of RSIDs
  for the full dataset used to calculate LDscores

## Value

a [`list()`](https://rdrr.io/r/base/list.html)

## Examples

``` r
if (FALSE) { # \dontrun{
ldsc_to_parquet("/directory/ldsc")
} # }
```
