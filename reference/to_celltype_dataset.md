# Convert the output of many LDscores to a single dataframe suitable for ldsR

Convert the output of many LDscores to a single dataframe suitable for
ldsR

## Usage

``` r
to_celltype_dataset(parent_dir, outdir, ref_snps = NULL)
```

## Arguments

- parent_dir:

  a directory with subdirectories containing LDscore data

- outdir:

  directory to save the parquet files

- ref_snps:

  a character vector of RSIDs corresponding to SNPs used in the full
  dataset used to generate LDscores. Provide if –thin-annot has been
  used in the LDSC command to generate ldscores

## Examples

``` r
if (FALSE) { # \dontrun{
to_celltype_dataset("files/ldsc", "files/ldsc_parquet")
} # }
```
