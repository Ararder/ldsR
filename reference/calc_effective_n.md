# Calculate the effective sample size of a case-control GWAS

Calculate the effective sample size of a case-control GWAS

## Usage

``` r
calc_effective_n(N_case, N_control)
```

## Arguments

- N_case:

  number of cases

- N_control:

  Number of controls

## Value

a double

## Examples

``` r
calc_effective_n(500, 1000)
#> [1] 1333.333
```
