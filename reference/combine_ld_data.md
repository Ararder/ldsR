# Combine two ldsc_to_parquet outputs

Combine two ldsc_to_parquet outputs

## Usage

``` r
combine_ld_data(list1, list2, ref = TRUE, outdir = NULL)
```

## Arguments

- list1, list2:

  output of ldsc_to_parquet

- ref:

  if TRUE, combine the reference annotations

- outdir:

  if NULL, return the combined data, otherwise save to outdir

## Value

a list or NULL

## Examples

``` r
if (FALSE) { # \dontrun{
combine_ld_data("files/ldsc_parquet/CD4_T_cells", "files/ldsc_parquet/CD8_T_cells")
} # }
```
