# Read LDscore format from ldsRs parquet format

Read LDscore format from ldsRs parquet format

## Usage

``` r
parse_parquet_dir(dir, read_ref = FALSE, subset_annots = NULL)
```

## Arguments

- dir:

  directory containing the files 'annot.parquet', 'annot_ref.parquet'
  and 'ld.parquet'

- read_ref:

  logical, whether to read the reference file 'annot_ref.parquet'

- subset_annots:

  Specify a character vector of annotations to read in, if NULL, all
  annotations are read in

## Value

a [`list()`](https://rdrr.io/r/base/list.html)

## Examples

``` r
if (FALSE) { # \dontrun{
files <- parse_parquet_dir("ldscores/atac")
} # }
```
