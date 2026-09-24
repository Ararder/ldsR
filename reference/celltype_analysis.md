# Perform cell-type analysis using partitioned heritability

A common usage of partitioned heritability is to estimate enrichment and
significance for an annotation typical for a specific cell-type. For
this analysis, it is useful to adjust for a baseline set of annotations,
and then estimate the association for a large set of annotations
("cell-types")

## Usage

``` r
celltype_analysis(sumstat, covariate_dir, ldscore_dir, weights = NULL)
```

## Arguments

- sumstat:

  A
  [`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)
  with columns `SNP`, `Z` and `N`

- covariate_dir:

  a directory containing the files `ld.parquet` and `annot.parquet`

- ldscore_dir:

  a directory containing the files `ld.parquet` and `annot.parquet`.
  Each annotation will be executed, while adjusting for all annotations
  in covariate_dir

- weights:

  Optional, a data.frame or tbl with columns `SNP`, `L2`

## Value

a
[`dplyr::tibble()`](https://dplyr.tidyverse.org/reference/reexports.html)

## Examples

``` r
if (FALSE) { # \dontrun{
sumstat1 <- arrow::read_parquet(system.file("extdata", "sumstats.parquet", package = "ldsR"))
sumstat1 <- dplyr::select(sumstat1, SNP, Z = Z.x, N=N.x)


celltype_analysis(
 sumstat1,
 covariate_dir = system.file("extdata", "baseline1.1_test", package = "ldsR"),
 ldscore_dir = system.file("extdata", "superclusters", package = "ldsR")
)
} # }
```
