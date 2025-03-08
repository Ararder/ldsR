---
title: "reference-files"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{reference-files}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(ldsR)
library(purrr)
library(arrow)
library(fs)
library(dplyr)
```

# Files shipped with the ldsR package
Three set of files are made available through the ldsR package.
1.  `eur_w_ld_chr` - used for univariate ldscore regression (estimating heritability and genetic correlations)
2.  `w_hm3.snplist.parquet` used for the munge function. This file exists to recreate the output of the munge.py function
3.  `1000G_Phase3_weights_hm3_no_MHC.parquet` - a second set of regression weights - used for partitioned heritability


# How were the reference files used by ldsR construced?
There are a lot of different versions of the original LDscore files floating around the internet. 
for example, the tutorial from the ldsc [github](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation)
references the `eur_w_ld_chr` folder, which as far as i can tell is longer available online.

The readme from the previous location https://alkesgroup.broadinstitute.org/LDSCORE/README_new_data_location.txt references
the S-LDSC directory from Steven Gaza at [Zenodo](https://doi.org/10.5281/zenodo.7768714)

However, the files corresponding to univariate LDscore regression, 
`1000G_Phase3_ldscores.tgz` are not identical to the previously standard, `eur_w_ld_chr`
I have made available the version of `eur_w_ld_chr` that our group had available at [zenodo](https://zenodo.org/records/14993076)
This file is shipped with the ldsR package. 

```{r, eval = FALSE}
setwd(tempdir())
# Define the URL and destination file
url <- "https://zenodo.org/records/14993076/files/eur_w_ld_chr.tgz"
destfile <- "eur_w_ld_chr.tgz"

# Download the file
download.file(url, destfile, mode = "wb")
# Extract the contents
untar(destfile, exdir = "eur_w_ld_chr")

# Remove the tar file after extraction
file.remove(destfile)

# contains LD split by chr. 
list.files("eur_w_ld_chr/eur_w_ld_chr")

# ordered matters here for the jackknife. Reorder as intended
ld_files <- paste0("eur_w_ld_chr/eur_w_ld_chr/", c(1:22), ".l2.ldscore.gz")

ld <- map(ld_files, \(x) read_delim_arrow(x, delim = "\t")) |> 
  list_rbind()

# Verify that they are ordered by physical location.
all(dplyr::arrange(ld, CHR, BP)$SNP == ld$SNP)

# get M
m_files <- paste0("eur_w_ld_chr/eur_w_ld_chr/", c(1:22), ".l2.M_5_50")
m <- map_dbl(m_files, \(x) readLines(x) |> as.numeric()) |> sum()
# m : 1173569

dplyr::select(ld, SNP, L2) |> 
  arrow::write_parquet("inst/extdata/eur_w_ld.parquet")


# The ldsc munge function actual restrices the regression SNPs even further, with the set of SNPs in w_hm3.snplist" 
df <- arrow::read_tsv_arrow("eur_w_ld_chr/eur_w_ld_chr/w_hm3.snplist") |> dplyr::tibble()

arrow::write_parquet(df, "inst/extdata/w_hm3.snplist.parquet")

```

# partitioned heritability (multivariate LDscore regression)

For partitioned heritability, cell-type analysis and other applications with multiple annotations, another set of weights is the standard use case.
Here we turn the files on zenodo from Steven Gazal at [Zenodo](https://doi.org/10.5281/zenodo.7768714).


```{r, eval = FALSE}
setwd(tempdir())
# Define the URL and destination file
url <- "https://zenodo.org/records/7768714/files/1000G_Phase3_weights_hm3_no_MHC.tgz"
destfile <- "1000G_Phase3_weights_hm3_no_MHC.tgz"

# Download the file
download.file(url, destfile, mode = "wb")
# Extract the contents
untar(destfile, exdir = "1000G_Phase3_weights_hm3_no_MHC")

# Remove the tar file after extraction
file.remove(destfile)

# contains LD split by chr. 
list.files("1000G_Phase3_weights_hm3_no_MHC/1000G_Phase3_weights_hm3_no_MHC")

# again, order matters
ld_files <- paste0("1000G_Phase3_weights_hm3_no_MHC/1000G_Phase3_weights_hm3_no_MHC/", "weights.hm3_noMHC.", c(1:22),".l2.ldscore.gz")
ld <- map(ld_files, \(x) read_delim_arrow(x, delim = "\t")) |> 
  list_rbind()

all(dplyr::arrange(ld, CHR, BP)$SNP == ld$SNP)
arrow::write_parquet(dplyr::select(ld, SNP, L2), "inst/extdata/1000G_Phase3_weights_hm3_no_MHC.parquet")
```



# partitioned heritability, baseline annotationss
```{r, eval = FALSE}
setwd("~/")
# Define the URL and destination file
options(timeout = max(3000, getOption("timeout")))
downloaded <- TRUE
if(!downloaded) {
  url <- "https://zenodo.org/records/7768714/files/1000G_Phase3_baseline_v1.2_ldscores.tgz"
  destfile <- "1000G_Phase3_baseline_v1.2_ldscores.tgz"
  
  # Download the file
  download.file(url, destfile, mode = "wb")
  # Extract the contents
  untar(destfile, exdir = "baseline_v1.2")
  # Remove the tar file after extraction
  file.remove(destfile)
}

dir <- "~/baseline_v1.2/"
naming <- fs::dir_ls(dir, glob = "*ldscore.gz")[1] |> 
  fs::path_file() |> 
  stringr::str_remove(".1.l2.ldscore.gz")

paste0(naming, ".1.l2.ldscore.gz")
stopifnot(fs::file_exists(fs::path(dir, paste0(naming, ".1.l2.ldscore.gz"))))

ld_paths <- fs::path(dir, paste0(naming,".", c(1:22), ".l2.ldscore.gz"))

m50_paths <- fs::path(dir, paste0(naming,".", c(1:22), ".l2.M_5_50"))
m_paths <- fs::path(dir, paste0(naming,".", c(1:22), ".l2.M"))
annot_path <- fs::path(dir, paste0(naming,".", c(1:22), ".annot.gz"))

ld <- ld_paths |> 
  purrr::map(arrow::read_tsv_arrow) |>
  purrr::list_rbind()
ld <- select(ld, -CHR,-BP)

m50 <- m50_paths |>
  purrr::map(\(x) read_tsv_arrow(x, col_names = FALSE)) |> 
  map(\(x) tidyr::pivot_longer(x, everything()) |> pull(value)) |> 
  reduce(`+`)

m <- m_paths |>
  purrr::map(\(x) read_tsv_arrow(x, col_names = FALSE)) |> 
  map(\(x) tidyr::pivot_longer(x, everything()) |> pull(value)) |> 
  reduce(`+`)

annot <- dplyr::tibble(annot = colnames(ld)[-1], m50 = m50, m = m)



annot_ref <-
  annot_path |>
  purrr::map(\(x) arrow::read_tsv_arrow(x)) |>
  purrr::list_rbind()
annot_ref <- select(annot_ref, -CHR, -CM, -BP)




outdir <- "~/baseline_model-v1.2" |> dir_create()
arrow::write_parquet(annot, fs::path(outdir, "annot.parquet"))
arrow::write_parquet(ld, fs::path(outdir, "ld.parquet"))
arrow::write_parquet(annot_ref, fs::path(outdir, "annot_ref.parquet"))

```
