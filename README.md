
# ldsR

<!-- badges: start -->
<!-- badges: end -->

The core functions of the [ldsc python package](https://github.com/bulik/ldsc) such as 
estimating heritability, intercept, genetic correlations and partitioned heritability are extremely common in the analysis of genetic data. 

However, the interface to these functions are no longer user friendly with difficulties in installation of python 2.7 dependencies and getting the nessecary reference data.

In addition, R is commonly used for the everyday tasks of inspecting, munging and working with data, making it time costly to swap to a command line interface for analysis that takes seconds. 

Here we introduce the core ldsc algorithms rewritten in R and make the most common reference data available within the R package making the estimation of genetic correlations, heritability, intercept and partitioned heritabiliy easier than ever.

Compare the code to estimate heritability with ldsc using the command line interface with ldsR:

``` bash
# NOTE that these links are no longer valid 
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
tar -jxvf eur_w_ld_chr.tar.bz2
unzip -o pgc.cross.scz.zip
bunzip2 w_hm3.snplist.bz2

munge_sumstats.py \
--sumstats pgc.cross.SCZ17.2013-05.txt \
--N 17115 \
--out scz \
--merge-alleles w_hm3.snplist

ldsc.py \
--h2 scz.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out scz
```

In ldsR, the ldscore data contained in the eur_w_ld_chr folder comes with the R package, and the core functions take only a few seconds to run on a modern computer.
``` r
library(ldsR)
sumstats <- readr::read_tsv("my_sumstats.tsv") |>
    munge()
h2_res <- ldsc_h2(sumstats)

# estimate rg with itself
rg_est = ldsc_rg(sumstats, sumstats)
```

## Installation

You can install the development version of ldsR from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Ararder/ldsR")
```

