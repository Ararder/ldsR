# Compute the genetic correlation between two traits

Compute the genetic correlation between two traits

## Usage

``` r
ldsc_rg(sumstats1, sumstats2, weights = NULL, M = NULL, n_blocks = 200)
```

## Arguments

- sumstats1:

  a
  [`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
  with atleast the columns SNP, A1, A2, Z, N. To perform quality checks,
  use [`munge()`](https://ararder.github.io/ldsR/reference/munge.md)
  before running `ldsc_rg()`

- sumstats2:

  a
  [`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
  with atleast the columns SNP, A1, A2, Z, N. To estimate rG with
  between sumstats1 and several traits, wrap the data.frames in a list.
  If the names are lost, the column `trait2` will correspond to the
  index of the list. Use a named list the retain the names in the
  `trait2` column `list(dataframe1, dataframe2, ...)`

- weights:

  Optional, a data.frame or tbl with columns `SNP`, `L2`

- M:

  Optional, the number of SNPs in the reference panel

- n_blocks:

  Number of blocks to use for the jackknife estimator

## Value

a
[`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)

## Examples

``` r
path <- system.file("extdata", "eur_w_ld.parquet", package = "ldsR")
snps <- arrow::read_parquet(path, col_select = c("SNP"))
# to make example faster
snps <- dplyr::slice_head(snps, n = 100000)
snps$A1 <- "A"
snps$A2 <-
snps$N <- 130000
snps$Z <- rnorm(nrow(snps))
snps2 <- snps
snps2$N <- 75000
snps2$Z <- rnorm(nrow(snps))
ldsc_rg(snps, snps2)
#> ! 0 SNPs were removed when merging summary statistics
#> Warning: NaNs produced
#> Warning: NaNs produced
#> # A tibble: 1 × 13
#>      rg rg_se h2_trait1 h2_trait_se h2_trait2 h2_trait2_se     gcov gcov_se
#>   <dbl> <dbl>     <dbl>       <dbl>     <dbl>        <dbl>    <dbl>   <dbl>
#> 1   NaN    NA  0.000682     0.00392  -0.00572      0.00633 0.000787 0.00392
#> # ℹ 5 more variables: gcov_int <dbl>, gcov_int_se <dbl>, mean_z1z2 <dbl>,
#> #   z <dbl>, p <dbl>
# to run estimate the genetic correlations for many traits, wrap s2 in a list
# ldsc_rg(snps, list(snps2, snps2))
# use a named list to create the `trait2` column in the output
# ldsc_rg(snps, list("trait2" = s2, "trait3" = s3, "trait4" = s4))
```
