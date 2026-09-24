# Map CHR:POS to RSID

`to_rsid()` adds a `SNP` column with RSIDs to summary statistics that
only have genomic coordinates, so that they can be passed to
[`munge()`](https://ararder.github.io/ldsR/reference/munge.md) and the
downstream `ldsc_*()` functions.

Variants are matched on `CHR`, `POS` and the unordered allele pair
(`A1`/`A2` vs `REF`/`ALT`) against a bundled reference covering all SNPs
in the default LD score and weight files (HapMap3, dbSNP155 coordinates
for GRCh37 and GRCh38). Rows that do not match the reference are
removed.

## Usage

``` r
to_rsid(dset, build = c("guess", "37", "38"))
```

## Arguments

- dset:

  a
  [`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
  with columns `CHR`, `POS`, `A1` and `A2`. `CHR` can be given as `1` or
  `chr1`. Only autosomes are matched.

- build:

  genome build of `POS`: `"guess"` (default) picks the build with the
  most matching variants, or pass `"37"` or `"38"`.

## Value

`dset` with a `SNP` column, restricted to variants found in the
reference

## Examples

``` r
if (FALSE) { # \dontrun{
sumstats |>
  to_rsid() |>
  munge() |>
  ldsc_h2()
} # }
```
