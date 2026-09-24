# Flexible partitioning of heritability

An R implementation of the LD score regression method to estimate SNP
heritability, focusing on partioning heritability across multiple
annotations.

## Usage

``` r
partition_h2(
  sumstat,
  ldscore_dirs,
  subset_annots = NULL,
  overlapping_annotations = FALSE,
  weights = NULL,
  n_blocks = 200
)
```

## Arguments

- sumstat:

  A
  [`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
  with columns `SNP`, `Z` and `N`

- ldscore_dirs:

  a list of directories containing ldscore files

- subset_annots:

  a character vector of annotations to subset the ldscore files. Use
  [`get_annot_names()`](https://ararder.github.io/ldsR/reference/get_annot_names.md)

- overlapping_annotations:

  are the annotations overlapping? In such a case, the estimate of the
  enrichment estimate needs to be adjusted.

- weights:

  Optional, a data.frame or tbl with columns `SNP`, `L2`

- n_blocks:

  Number of blocks to use for the jackknife estimator

## Value

a tibble with results

## Examples

``` r
if (FALSE) { # \dontrun{
partition_h2(sumstat, ldscore_dirs, subset_annots = c("L2", "L3"))
} # }
```
