# Estimate SNP heritability using LDscore regression for a single annotation

An R implementation of the LD score regression method to estimate SNP
heritability, mimicking `ldsc --h2` from the
[ldsc](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation)
package.

LDscores for the European subset of 1000g (2015 release) for 1,290,028
HapMap3 SNPs (with 1,173,569 SNPs with freq \> 5%) are bundled within
the ldsR package and are used by default. This corresponds to the
`eur_w_ld_chr` folder previously shared at the LDSC
[github](https://github.com/bulik/ldsc).

You can inspect the LDscores used by default:
`arrow::read_parquet(system.file("extdata", "eur_w_ld.parquet", package = "ldsR"))`

ldsc_h2 does not perform any quality control on the input summary
statistics, except to merge with the 1,290,028 HapMap3 SNPs in the
reference panel. See the
[`munge()`](https://ararder.github.io/ldsR/reference/munge.md) function
to mimic the `munge_sumstats.py` function.

## Usage

``` r
ldsc_h2(
  sumstat,
  pop_prev = NULL,
  sample_prev = 0.5,
  weights = NULL,
  M = NULL,
  n_blocks = 200
)
```

## Arguments

- sumstat:

  A
  [`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
  with columns `SNP`, `Z` and `N`

- pop_prev:

  prevalence of the disorder in the general population

- sample_prev:

  the prevalence of the disorder in the sample. Default value is 0.5,
  reflecting a case-control study using effective N as sample size

- weights:

  Optional, a data.frame or tbl with columns `SNP`, `L2`

- M:

  Optional, the number of SNPs in the reference panel

- n_blocks:

  Number of blocks to use for the jackknife estimator

## Value

a
[`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
with columns `h2` and `h2_se`

## Examples

``` r
p <- system.file("extdata", "eur_w_ld.parquet", package = "ldsR")
snps <- arrow::read_parquet(p, col_select = c("SNP"))
snps$N <- 130000
snps$Z <- rnorm(nrow(snps))
ldsc_h2(snps)
#> # A tibble: 1 × 6
#>          h2   h2_se   int  int_se mean_chi2 lambda_gc
#>       <dbl>   <dbl> <dbl>   <dbl>     <dbl>     <dbl>
#> 1 -0.000113 0.00110 1.000 0.00312     0.999     0.996

```
