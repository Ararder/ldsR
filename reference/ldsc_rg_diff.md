# Test whether two genetic correlations differ

Block jackknife test of the null hypothesis that the genetic correlation
between traits A and B equals the genetic correlation between traits C
and D. A trait may appear in both pairs (e.g. rG(A,B) vs rG(C,B)) by
being passed twice.

All summary statistics are merged with the weights to a single common
SNP set before either rG is computed, so the jackknife blocks correspond
one-to-one across the two correlations. This captures the covariance
between the two rGs, which is non-zero whenever a trait (or just SNPs)
is shared.

For each block `i`, delete values `R(A,B)_i` and `R(C,D)_i` for the rG
are differenced to form `D_i`. Pseudovalues are
`P_i = n * d - (n - 1) * D_i`, where `d = rG(A,B) - rG(C,D)` is the
global difference and `n` is the number of blocks. The test statistic is
`z = mean(P) / sqrt(var(P) / n)`.

## Usage

``` r
ldsc_rg_diff(pair1, pair2, weights = NULL, M = NULL, n_blocks = 200)
```

## Arguments

- pair1:

  A list of two summary statistics tibbles defining the first genetic
  correlation. Each tibble must contain `SNP`, `A1`, `A2`, `Z`, `N`.

- pair2:

  Same shape as `pair1`. May share a tibble with `pair1` (pass the same
  data frame in both pairs to share a trait).

- weights:

  Optional, a data.frame or tbl with columns `SNP`, `L2`

- M:

  Optional, the number of SNPs in the reference panel

- n_blocks:

  Number of blocks to use for the jackknife estimator

## Value

a
[`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
with the two global rGs, the global difference, the jackknife-corrected
difference, its standard error, and the z statistic with chi-squared
p-value.

## Examples

``` r
if (FALSE) { # \dontrun{
# Four distinct traits
ldsc_rg_diff(
  pair1 = list(trait_a, trait_b),
  pair2 = list(trait_c, trait_d)
)

# Three traits: rG(A,B) vs rG(C,B), with B shared
ldsc_rg_diff(
  pair1 = list(trait_a, trait_b),
  pair2 = list(trait_c, trait_b)
)
} # }
```
