# Test whether two SNP heritabilities differ

Block jackknife test of the null hypothesis that h²(trait 1) equals
h²(trait 2).

Both summary statistics are merged with the weights to a single common
SNP set before either h² is computed, so the jackknife blocks correspond
one-to-one across the two traits. The per-block h² delete values from
each trait are differenced; pseudovalues `P_i = n * d - (n - 1) * D_i`
with `d = h²_1 - h²_2` give a jackknife estimate of the difference and
its variance, and `z = mean(P) / sqrt(var(P)/n)`.

If `pop_prev1` (and optionally `pop_prev2`) is supplied, the
corresponding h² and per-block delete values are converted to liability
scale via
[`liability_h2()`](https://ararder.github.io/ldsR/reference/liability_h2.md)
before differencing. Comparing observed-scale h² for one trait against
liability-scale h² for the other is allowed but rarely makes sense.

## Usage

``` r
ldsc_h2_diff(
  sumstat1,
  sumstat2,
  return_blocks = FALSE,
  pop_prev1 = NULL,
  sample_prev1 = 0.5,
  pop_prev2 = NULL,
  sample_prev2 = 0.5,
  weights = NULL,
  M = NULL,
  n_blocks = 200
)
```

## Arguments

- sumstat1, sumstat2:

  Summary statistics tibbles with at least `SNP`, `Z`, `N`.

- return_blocks:

  If `TRUE`, return the per-block h² delete values for both traits

- pop_prev1, sample_prev1:

  Optional liability-scale parameters for `sumstat1`. Leave
  `pop_prev1 = NULL` to test on the observed scale.

- pop_prev2, sample_prev2:

  Same, for `sumstat2`.

- weights:

  Optional, a data.frame or tbl with columns `SNP`, `L2`

- M:

  Optional, the number of SNPs in the reference panel

- n_blocks:

  Number of blocks to use for the jackknife estimator

## Value

a
[`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
with the two h² estimates (on the chosen scale), the global difference,
the jackknife-corrected difference, its standard error, and the z
statistic with chi-squared p-value.

## Examples

``` r
if (FALSE) { # \dontrun{
# Observed-scale comparison
ldsc_h2_diff(trait_a, trait_b)

# Liability-scale comparison of two binary traits
ldsc_h2_diff(
  trait_a, trait_b,
  pop_prev1 = 0.01, sample_prev1 = 0.5,
  pop_prev2 = 0.05, sample_prev2 = 0.5
)
} # }
```
