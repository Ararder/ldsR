# 

``` r

library(ldsR)
library(testthat)
library(fs)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

## Introduction

ldsR is a R rewrite of the [ldsc](https://github.com/bulik/ldsc)
software. The aim is to provide a significantly simplified user
experience.

Currently, ldsR is able to reproduce all features of ldsc, except:

1\. LDscore calculation

2\. Quantitative annotations in partitioned heritability.

The showcase here is based on a small example file of 100,000 rows, with
effect sizes from schizophrenia and bipolar disorder.

``` r

df <- arrow::read_parquet(system.file("extdata", "sumstats.parquet", package = "ldsR"))
head(df)
#> # A tibble: 6 × 5
#>   SNP          Z.x    N.x    Z.y    N.y
#>   <chr>      <dbl>  <dbl>  <dbl>  <dbl>
#> 1 rs3094315  0.289 126849 -0.101 405254
#> 2 rs3131972  0.269 126849 -0.299 405254
#> 3 rs3131969 -0.074 126849 -0.296 405254
#> 4 rs1048488  0.273 126849 -0.051 405254
#> 5 rs3115850  0.263 126849 -0.045 405254
#> 6 rs2286139  0.034 126849 -0.68  387638
```

## Munge

To perfectly reproduce LDSC results, the input summary statistics needs
to filtered in an identical manner.
[`munge()`](https://ararder.github.io/ldsR/reference/munge.md) mirrors
the functionality of `munge_sumstats.py`. Munge requires the following
columns:

`SNP` `A1` `A2` `Z` `N`

``` r

# create alleles so that we can showcase how munge works
test_df <- mutate(df, A1 = "G", A2 = "T", Z = Z.x, N=N.x)
test_df <- munge(test_df)
#> ! Removed 0 rows after filtering on reference panel
#> ! Removed 0 rows with non-RSIDs in SNP column
#> ! Removed 0 rows with duplicated RSIDs
#> ! Removed 0 rows due to strand ambigious alleles
#> ! Removed 0 rows with a sample size smaller than 87096
```

## SNP heritability

[`ldsc_h2()`](https://ararder.github.io/ldsR/reference/ldsc_h2.md)
estimates the SNP-heritability for a single trait, one at a time. By
default, regression weights and LDscores are from the European 1000G
reference panel.

The only required argument is the summary statistics. If you want to use
other regression and LDscore weights (for example from another ancestry,
you can use the `weights` argument and the `M` argument).

``` r

ldsc_h2(test_df)
#> # A tibble: 1 × 6
#>      h2  h2_se   int int_se mean_chi2 lambda_gc
#>   <dbl>  <dbl> <dbl>  <dbl>     <dbl>     <dbl>
#> 1 0.371 0.0358  1.09 0.0328      2.14      1.80
```

For binary traits, it’s often useful to convert the heritability
estimate to the liability scale. This can be done using the `pop_prev`
and `sample_prev` arguments.

In this particular case, the N column in the example file reflects the
total number of cases and controls. In such a case, we need both the
population prevalence and sample prevalence

``` r

n_cont <- 77344
n_cas <- 53300
s_prev <- n_cas/(n_cas + n_cont)
ldsc_h2(test_df, pop_prev = 0.01, sample_prev = s_prev)
#> # A tibble: 1 × 8
#>   lia_h2 lia_se    h2  h2_se   int int_se mean_chi2 lambda_gc
#>    <dbl>  <dbl> <dbl>  <dbl> <dbl>  <dbl>     <dbl>     <dbl>
#> 1  0.212 0.0205 0.371 0.0358  1.09 0.0328      2.14      1.80
```

Ideally, for estimation of heritability, the per SNP sample size should
be calculated using the effective sample size summed across cohorts, see
[here](https://pmc.ncbi.nlm.nih.gov/articles/PMC10066905/). If the N
column corresponds to the effective sample size, only the `pop_prev`
argument is required (a 50% sample prevalence is assumed).

``` r

eff_n <- calc_effective_n(n_cas, n_cont)
test_df$N <- eff_n
ldsc_h2(test_df, pop_prev = 0.01, sample_prev = s_prev)
#> # A tibble: 1 × 8
#>   lia_h2 lia_se    h2  h2_se   int int_se mean_chi2 lambda_gc
#>    <dbl>  <dbl> <dbl>  <dbl> <dbl>  <dbl>     <dbl>     <dbl>
#> 1  0.219 0.0212 0.384 0.0371  1.09 0.0329      2.14      1.80
```

Lastly, you can estimate the observed scale heritability first, and then
calculate the liability scale heritability at your convenience using
[`liability_h2()`](https://ararder.github.io/ldsR/reference/liability_h2.md).

``` r

liability_h2(0.37, 0.01, 0.4)
#> [1] 0.2127143
```

## Genetic correlations

ldsR also implements the estimation of genetic correlations in
[`ldsc_rg()`](https://ararder.github.io/ldsR/reference/ldsc_rg.md).
Regression weights and LDscores are from the European 1000G reference
panel, same as with
[`ldsc_h2()`](https://ararder.github.io/ldsR/reference/ldsc_h2.md).

Technically,
[`ldsc_h2()`](https://ararder.github.io/ldsR/reference/ldsc_h2.md) only
requires SNP, Z and N.

[`ldsc_rg()`](https://ararder.github.io/ldsR/reference/ldsc_rg.md) on
the other hand requires both `A1` and `A2`, to be able to align Z scores
to the same reference allele.

Using the example data:

``` r

sumstat1 <- dplyr::rename(df, Z = Z.x, N = N.x) |> 
  dplyr::mutate(A1 = "A", A2 = "G")
sumstat2 <- dplyr::rename(df, Z = Z.y, N = N.y) |> 
  dplyr::mutate(A1 = "A", A2 = "G")

ldsc_rg(sumstat1, sumstat2)
#> ! 0 SNPs were removed when merging summary statistics
#> # A tibble: 1 × 13
#>      rg  rg_se h2_trait1 h2_trait_se h2_trait2 h2_trait2_se  gcov gcov_se
#>   <dbl>  <dbl>     <dbl>       <dbl>     <dbl>        <dbl> <dbl>   <dbl>
#> 1 0.734 0.0496     0.371      0.0358    0.0732      0.00685 0.121  0.0113
#> # ℹ 5 more variables: gcov_int <dbl>, gcov_int_se <dbl>, mean_z1z2 <dbl>,
#> #   z <dbl>, p <dbl>
```

Additionally,
[`ldsc_rg()`](https://ararder.github.io/ldsR/reference/ldsc_rg.md)
attemps to make it easier to calculate genetic correlations among many
traits.

sumstats2 can be both a
[`dplyr::tibble()`](https://tibble.tidyverse.org/reference/tibble.html),
and a list containing multiple
[`dplyr::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)’s.
If a list is passed to the `sumstats2` argument, the genetic correlation
will be estimated between `sumstats1` and all entries in the list.

``` r

sumstats <- list(
  "trait1" = sumstat2,
  "trait2" = sumstat2,
  "trait4" = sumstat2
)
          
ldsc_rg(sumstat1, sumstats)
#> ! 0 SNPs were removed when merging summary statistics
#> ⠙ 1/3 ETA:  3s | Computing genetic correlations...
#> ! 0 SNPs were removed when merging summary statistics
#> ⠙ 1/3 ETA:  3s | Computing genetic correlations...  ! 0 SNPs were removed when merging summary statistics
#> ⠙ 1/3 ETA:  3s | Computing genetic correlations...  
#> # A tibble: 3 × 14
#>   trait2    rg  rg_se h2_trait1 h2_trait_se h2_trait2 h2_trait2_se  gcov gcov_se
#>   <chr>  <dbl>  <dbl>     <dbl>       <dbl>     <dbl>        <dbl> <dbl>   <dbl>
#> 1 trait1 0.734 0.0496     0.371      0.0358    0.0732      0.00685 0.121  0.0113
#> 2 trait2 0.734 0.0496     0.371      0.0358    0.0732      0.00685 0.121  0.0113
#> 3 trait4 0.734 0.0496     0.371      0.0358    0.0732      0.00685 0.121  0.0113
#> # ℹ 5 more variables: gcov_int <dbl>, gcov_int_se <dbl>, mean_z1z2 <dbl>,
#> #   z <dbl>, p <dbl>
```

If the list is not named, the names will correspond to the index in the
list.

``` r

sumstats <- list(
  sumstat2,
  sumstat2,
  sumstat2
)
          
ldsc_rg(sumstat1, sumstats)
#> ! 0 SNPs were removed when merging summary statistics
#> ⠙ 1/3 ETA:  3s | Computing genetic correlations...
#> ! 0 SNPs were removed when merging summary statistics
#> ⠙ 1/3 ETA:  3s | Computing genetic correlations...  ! 0 SNPs were removed when merging summary statistics
#> ⠙ 1/3 ETA:  3s | Computing genetic correlations...  
#> # A tibble: 3 × 14
#>   trait2    rg  rg_se h2_trait1 h2_trait_se h2_trait2 h2_trait2_se  gcov gcov_se
#>    <int> <dbl>  <dbl>     <dbl>       <dbl>     <dbl>        <dbl> <dbl>   <dbl>
#> 1      1 0.734 0.0496     0.371      0.0358    0.0732      0.00685 0.121  0.0113
#> 2      2 0.734 0.0496     0.371      0.0358    0.0732      0.00685 0.121  0.0113
#> 3      3 0.734 0.0496     0.371      0.0358    0.0732      0.00685 0.121  0.0113
#> # ℹ 5 more variables: gcov_int <dbl>, gcov_int_se <dbl>, mean_z1z2 <dbl>,
#> #   z <dbl>, p <dbl>
```

## Partitioning heritability

`ldsR` deviates from `ldsc` in the file format used for the partitioned
LD scores.

The example files used below are subset of the full LDscores, which are
available [here](https://zenodo.org/records/15194124).

The original baseline annotations in ldsc LDscore format are available
[here](https://zenodo.org/records/10515792).

The primary difference between ldsR and ldsc partitioned ldscores is the
number of files. LDSC splits out the files across chromosomes and
data-type, resulting in 88 files per data set (22 chromosomes, .ldscore,
.annot, .m50, .m5_50).

ldsR collapses this to three files, stored in parquet files.

``` r

dir_tree(system.file("extdata", "superclusters/", package = "ldsR"))
#> /home/runner/work/_temp/Library/ldsR/extdata/superclusters/
#> ├── annot.parquet
#> ├── annot_ref.parquet
#> └── ld.parquet
```

`annot.parquet` contains information on the number of SNPs, ie the M50
and M5_50 files per annotation

``` r

arrow::read_parquet(system.file("extdata", "superclusters/annot.parquet", package = "ldsR")) |>
  head()
#> # A tibble: 4 × 3
#>   annot                    m50       m
#>   <chr>                  <dbl>   <dbl>
#> 1 amygdala_excitatory   911324 1540052
#> 2 astrocyte             713230 1207523
#> 3 bergmann_glia         696875 1183544
#> 4 cerebellar_inhibitory 690546 1168161
```

`ld.parquet` contains the LDscores per annotation

``` r

arrow::read_parquet(system.file("extdata", "superclusters/ld.parquet", package = "ldsR")) |> 
  head()
#> # A tibble: 6 × 5
#>   SNP       amygdala_excitatory astrocyte bergmann_glia cerebellar_inhibitory
#>   <chr>                   <dbl>     <dbl>         <dbl>                 <dbl>
#> 1 rs3094315                2.65      42.3         0.041                 0.034
#> 2 rs3131972                2.74      42.4         0.061                 0.175
#> 3 rs3131969                3.62      54.8         0.298                 0.625
#> 4 rs1048488                2.77      42.5         0.024                -0.024
#> 5 rs3115850                2.61      42.0        -0.025                -0.105
#> 6 rs2286139                4.08      55.0         0.25                  0.86
```

`annot_ref.parquet` contains the SNP to annotation mapping for each
annotation.

A value of 1 means that the SNP is mapped to that annotation.

``` r

arrow::read_parquet(system.file("extdata", "superclusters/annot_ref.parquet", package = "ldsR")) |> 
  head()
#> # A tibble: 6 × 32
#>   SNP       amygdala_excitatory astrocyte bergmann_glia cerebellar_inhibitory
#>   <chr>                   <int>     <int>         <int>                 <int>
#> 1 rs3094315                   0         0             0                     0
#> 2 rs3131972                   0         0             0                     0
#> 3 rs3131969                   0         0             0                     0
#> 4 rs1048488                   0         1             0                     0
#> 5 rs3115850                   0         1             0                     0
#> 6 rs2286139                   0         1             0                     0
#> # ℹ 27 more variables: cge_interneuron <int>, choroid_plexus <int>,
#> #   committed_oligodendrocyte_precursor <int>,
#> #   deep_layer_corticothalamic_and_6b <int>,
#> #   deep_layer_intratelencephalic <int>, deep_layer_near_projecting <int>,
#> #   eccentric_medium_spiny_neuron <int>, ependymal <int>, fibroblast <int>,
#> #   hippocampal_ca1_3 <int>, hippocampal_ca4 <int>,
#> #   hippocampal_dentate_gyrus <int>, lamp5_lhx6_and_chandelier <int>, …
```

The main function to partition heritability is
[`partition_h2()`](https://ararder.github.io/ldsR/reference/partition_h2.md).

1.  `sumstats` is the first argument, corresponding to a dataframe with
    SNP, Z and N. Use
    [`munge()`](https://ararder.github.io/ldsR/reference/munge.md)
    first.
2.  `ldscore_dirs` takes a vector of filepaths, to directories with
    partitioned ldscores
3.  `subset_annots` allows for flexible subsetting of annotations to use
    in the model. Default is NULL; which means that all of the annotions
    will be used in the model
4.  `overlapping_annotations` which will result in the enrichment for
    each annotation also being estimated.

In this first example we estimate the partitioned heritability across
the 4 annotations present in this example data.

``` r

# use the dataframe we munged earlier

partition_h2(
  test_df,
  ldscore_dirs = system.file("extdata", "baseline1.1_test/", package = "ldsR"),
)
#> ℹ Removed 449 rows after merging with weights
#> ℹ Removed 5961 rows after merging with ldscores
#> ✔ A total of 93590 SNPs remain for the regression model
#> # A tibble: 4 × 7
#>   annot                                coef    coef_se     z   tot tot_se n_snps
#>   <chr>                               <dbl>      <dbl> <dbl> <dbl>  <dbl>  <int>
#> 1 Conserved_LindbladToh.bedL2  0.000000880     2.63e-7  3.35 0.350 0.0393 1.53e5
#> 2 baseL2                       0.0000000274    8.79e-9  3.12 0.350 0.0393 5.96e6
#> 3 H3K4me1_peaks_Trynka.bedL2   0.0000000717    4.81e-8  1.49 0.350 0.0393 1.01e6
#> 4 Coding_UCSC.bedL2           -0.000000236     1.66e-7 -1.42 0.350 0.0393 8.51e4
```

A common analysis is to estimate the enrichment of a particular
annotation

Perhaps we are interested in doing “cell-type analysis”, testing the
enrichment of a specific annotation while adjusting for the 53 baseline
annotations. With
[`partition_h2()`](https://ararder.github.io/ldsR/reference/partition_h2.md),
we simply add the filepath to the directory with additional ldscores.

Here we use a subsetted example of LDscores derived from human brain
cell-types.

To be more specific in which LDscores we want to include in the model,
we need pass a list of annotation names that will be included in the
model. For this we need to know what the annotation names are per
ldscore dataset. Here we use
[`get_annot_names()`](https://ararder.github.io/ldsR/reference/get_annot_names.md),
which will simply return the names of all annotations in a ldscore
dataset.

``` r

all_ldscores <- c(
  system.file("extdata", "baseline1.1_test/", package = "ldsR"), 
  system.file("extdata", "superclusters/", package = "ldsR")
  )

baseline_names <- get_annot_names(all_ldscores[1])
supercluster_names <- get_annot_names(all_ldscores[2])
```

``` r

baseline_names
#> [1] "baseL2"                      "Coding_UCSC.bedL2"          
#> [3] "Conserved_LindbladToh.bedL2" "H3K4me1_peaks_Trynka.bedL2"
```

``` r

supercluster_names
#> [1] "amygdala_excitatory"   "astrocyte"             "bergmann_glia"        
#> [4] "cerebellar_inhibitory"
```

With this, we specify which annotations to include in the model through
`subset_annots`.

``` r

partition_h2(
  test_df,
  ldscore_dirs = all_ldscores,
  # here we test the association with hippocampal neurons
  subset_annot = c(baseline_names,"amygdala_excitatory")
) |> 
  filter(annot == "amygdala_excitatory")
#> ℹ Removed 449 rows after merging with weights
#> ℹ Removed 5961 rows after merging with ldscores
#> ✔ A total of 93590 SNPs remain for the regression model
#> # A tibble: 1 × 7
#>   annot                       coef      coef_se     z   tot tot_se n_snps
#>   <chr>                      <dbl>        <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1 amygdala_excitatory 0.0000000312 0.0000000180  1.73 0.350 0.0373 911324
```

With this setup we can fit any model we might like, for example by
adding additional cell-types to run conditional analysis.

``` r

partition_h2(
  test_df,
  ldscore_dirs = all_ldscores,
  # here we test the association with hippocampal neurons
  subset_annot = c(baseline_names,"amygdala_excitatory","astrocyte")
) |> 
  filter(annot %in% c("amygdala_excitatory","astrocyte"))
#> ℹ Removed 449 rows after merging with weights
#> ℹ Removed 5961 rows after merging with ldscores
#> ✔ A total of 93590 SNPs remain for the regression model
#> # A tibble: 2 × 7
#>   annot                        coef      coef_se     z   tot tot_se n_snps
#>   <chr>                       <dbl>        <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1 amygdala_excitatory  0.0000000272 0.0000000184  1.48 0.349 0.0379 911324
#> 2 astrocyte           -0.0000000186 0.0000000157 -1.18 0.349 0.0379 713230
```

Lastly, we can use `overlapping_annotations=TRUE` to estimate the
enrichment of each annotation.

``` r

partition_h2(
  test_df,
  ldscore_dirs = all_ldscores,
  subset_annot = c(baseline_names,"amygdala_excitatory","astrocyte"),
  overlapping_annotations=TRUE
) |> 
  filter(annot %in% c("amygdala_excitatory","astrocyte"))
```

## Cell-type analysis

Lastly,
[`celltype_analysis()`](https://ararder.github.io/ldsR/reference/celltype_analysis.md)
is a convenience function for estimating the association across several
annotations, while adjusting for a set of “baseline” annotations.

1.  `sumstat`
2.  `covariate_dir` filepath to annotations that are used as controls
3.  `ldscore_dir` filepath to annotations that are to be tested, one at
    a time while adjusting for all annotations in covariate_dir

[`celltype_analysis()`](https://ararder.github.io/ldsR/reference/celltype_analysis.md)
will loop over each annotation in ldscore_dir, testing it’s association
while adjusting for all annotations in covariate_dir.

``` r

celltype_analysis(
  sumstat = test_df,
  covariate_dir = all_ldscores[1],
  ldscore_dir = all_ldscores[2]
  
)
#> ℹ Removed 0 rows after merging with ldscores
#> # A tibble: 4 × 8
#>   annot                     coef      coef_se enrich   prop      z   tot tot_se
#>   <chr>                    <dbl>        <dbl>  <dbl>  <dbl>  <dbl> <dbl>  <dbl>
#> 1 amygdala_excitatory    3.12e-8 0.0000000180  0.724 0.112   1.73  0.350 0.0373
#> 2 astrocyte             -1.49e-8 0.0000000165 -0.336 0.0900 -0.905 0.352 0.0400
#> 3 bergmann_glia          4.55e-9 0.0000000159  0.101 0.0881  0.286 0.357 0.0396
#> 4 cerebellar_inhibitory -1.84e-8 0.0000000143 -0.411 0.0874 -1.29  0.354 0.0404
```

## 
