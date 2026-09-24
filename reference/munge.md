# Munge GWAS summary statistics

`munge()` applies 5 filters (if possible):

1.  Deduplication of RSID

2.  Filter variants on INFO in INFO column is present

3.  Filter variant on effect allele frequency

4.  Removes strand ambigious variants

5.  Removes variants with `N < round(stats::quantile(N, 0.9) / 1.5)`

## Usage

``` r
munge(dset, info_filter = 0.9, eaf_filter = 0.01)
```

## Arguments

- dset:

  a
  [`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
  with columns `SNP`, `A1` `A2` `Z` `N` and optional columns `EAF` and
  `INFO`

- info_filter:

  INFO score filter threshold at which to remove rows

- eaf_filter:

  effective allele filter at which to remove rows. eaf_filter=0.01 would
  filter to eaf \> 0.01 & eaf \< 0.99

## Value

a data.frame

## Examples

``` r
if (FALSE) { # \dontrun{
munge(tbl)
} # }
```
