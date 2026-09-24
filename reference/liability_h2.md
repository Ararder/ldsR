# Convert an estimate of observed-scale heritability to liability scale heritability

Convert an estimate of observed-scale heritability to liability scale
heritability

## Usage

``` r
liability_h2(obs_h2, pop_prev, sample_prev = 0.5)
```

## Arguments

- obs_h2:

  observed-scale heritability

- pop_prev:

  prevalence of the disorder in the general population

- sample_prev:

  the prevalence of the disorder in the sample. Default value is 0.5,
  reflecting a case-control study using effective N as sample size

## Value

a double

## Examples

``` r
liability_h2(0.25, 0.02)
#> [1] 0.1638687
```
