# Package index

## ldsR

### SNP heritability and genetic correlations

- [`munge()`](https://ararder.github.io/ldsR/reference/munge.md) : Munge
  GWAS summary statistics
- [`ldsc_h2()`](https://ararder.github.io/ldsR/reference/ldsc_h2.md) :
  Estimate SNP heritability using LDscore regression for a single
  annotation
- [`ldsc_rg()`](https://ararder.github.io/ldsR/reference/ldsc_rg.md) :
  Compute the genetic correlation between two traits
- [`ldsc_h2_diff()`](https://ararder.github.io/ldsR/reference/ldsc_h2_diff.md)
  : Test whether two SNP heritabilities differ
- [`ldsc_rg_diff()`](https://ararder.github.io/ldsR/reference/ldsc_rg_diff.md)
  : Test whether two genetic correlations differ

### Partition Heritability

- [`partition_h2()`](https://ararder.github.io/ldsR/reference/partition_h2.md)
  : Flexible partitioning of heritability
- [`celltype_analysis()`](https://ararder.github.io/ldsR/reference/celltype_analysis.md)
  : Perform cell-type analysis using partitioned heritability

### Helpful functions

- [`liability_h2()`](https://ararder.github.io/ldsR/reference/liability_h2.md)
  : Convert an estimate of observed-scale heritability to liability
  scale heritability

- [`calc_effective_n()`](https://ararder.github.io/ldsR/reference/calc_effective_n.md)
  : Calculate the effective sample size of a case-control GWAS

- [`from_tidyGWAS()`](https://ararder.github.io/ldsR/reference/from_tidyGWAS.md)
  :

  Parse GWAS format of `tidyGWAS::tidyGWAS()`

- [`to_rsid()`](https://ararder.github.io/ldsR/reference/to_rsid.md) :
  Map CHR:POS to RSID

### Read and transform LDscore data

- [`parse_parquet_dir()`](https://ararder.github.io/ldsR/reference/parse_parquet_dir.md)
  : Read LDscore format from ldsRs parquet format
- [`ldsc_to_parquet()`](https://ararder.github.io/ldsR/reference/ldsc_to_parquet.md)
  : Transform a directory of LDscores to parquet format
- [`to_celltype_dataset()`](https://ararder.github.io/ldsR/reference/to_celltype_dataset.md)
  : Convert the output of many LDscores to a single dataframe suitable
  for ldsR
- [`get_annot_names()`](https://ararder.github.io/ldsR/reference/get_annot_names.md)
  : Print names of annotations in a ldscore fileset
- [`combine_ld_data()`](https://ararder.github.io/ldsR/reference/combine_ld_data.md)
  : Combine two ldsc_to_parquet outputs
