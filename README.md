
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ritest

<!-- badges: start -->

[![R-CMD-check](https://github.com/grantmcdermott/ritest/workflows/R-CMD-check/badge.svg)](https://github.com/grantmcdermott/ritest/actions)
[![Codecov test
coverage](https://codecov.io/gh/grantmcdermott/ritest/branch/master/graph/badge.svg)](https://app.codecov.io/gh/grantmcdermott/ritest?branch=master)
<!-- badges: end -->

Conduct [**randomization
inference**](https://dimewiki.worldbank.org/Randomization_Inference) on
R model objects.

This R package is a port of the
[`-ritest-`](https://github.com/simonheb) Stata routine (by [Simon
Heß](https://github.com/simonheb)). It doesn’t (yet) try to support all
of the features in the Stata version and is currently limited to `lm()`
and `fixest::feols()` models. But it does appear to be significantly
faster, and aims to support a variety of model classes once it is fully
baked.

[`Installation`](#Installation) \| [`Useage`](#Useage) \|
[`Benchmarks`](#Benchmarks) \| [`License`](#)

## Installation

``` r
# install.packages("remotes")
remotes::install_github("grantmcdermott/ritest")
```

## Useage

Here follows a basic example using the `colombia` dataset that comes
bundled with the package. Note that I’ll also use the `fixest::feols()`
function to estimate the parametric model first.

``` r
library(ritest)       ## This package
library(fixest)       ## For fast (high-dimensional) fixed-effect regressions

data("colombia")

## Parametric model using fixest::feols()
co_est = 
  feols(
    dayscorab ~ b_treat + b_dayscorab + miss_b_dayscorab | b_pair + round2 + round3, 
    vcov = ~b_block, data = colombia
    )
```

Let’s conduct RI on the treatment variable (`b_treat`), whilst taking
account of the stratified and clustered experimental design.

``` r
co_ri = ritest(co_est, 'b_treat', strata='b_pair', cluster='b_block', reps=1e3, seed=546L)
co_ri
#> 
#>           Call: feols(fml = dayscorab ~ b_treat + b_dayscorab + miss_b_dayscorab | b_pair + round2 + round3, data = colombia, vcov = ~b_block)
#>    Res. var(s): b_treat
#>             H0: b_treat=0
#>  Strata var(s): b_pair
#>         Strata: 31
#> Cluster var(s): b_block
#>       Clusters: 63
#>      Num. reps: 1000
#> ──────────────────────────────────────────────────────────────────────────────── 
#>   T(obs)         c         n     p=c/n     SE(p)   CI 2.5%  CI 97.5%  
#>  -0.1807       102      1000     0.102   0.01575   0.07609    0.1279  
#> ──────────────────────────────────────────────────────────────────────────────── 
#> Note: Confidence interval is with respect to p=c/n. 
#> Note: c = #{|T| >= |T(obs)|}
```

A default plotting method is provided.

``` r
plot(co_ri)
```

<img src="man/figures/README-plot_co_ri-1.png" width="100%" />

For more examples and demonstration of additional features, please see
the introductory vignette.

## Benchmarks

I generally observe a speed increase of between 25x–50x compared to the
Stata version. Again, see the introductory vignette for timed examples.

## License

The software code contained within this repository is made available
under the [MIT license](http://opensource.org/licenses/mit-license.php).
