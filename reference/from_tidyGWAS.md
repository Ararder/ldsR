# Parse GWAS format of `tidyGWAS::tidyGWAS()`

Parse GWAS format of `tidyGWAS::tidyGWAS()`

## Usage

``` r
from_tidyGWAS(tbl, n = c("N", "EffectiveN"))
```

## Arguments

- tbl:

  a
  [`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)

- n:

  Column name of sample size, default is "N".

## Value

a munged
[`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)

## Examples

``` r
if (FALSE) { # \dontrun{
munged <- from_tidyGWAS("path/tidyGWAS/cleaned/tidyGWAS_hivestyle")
} # }
```
