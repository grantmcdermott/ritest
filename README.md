
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

This R package is a port of the excellent
[`-ritest-`](https://github.com/simonheb/ritest) Stata routine by Simon
Heß. It doesn’t (yet) try to support all of the features in the Stata
version and is currently limited to `lm()` and `fixest::feols()` models.
But it does appear to be significantly faster, and aims to support a
variety of model classes once it is fully baked.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("grantmcdermott/ritest")
```

## Usage

A detailed walkthrough of the package is provided in the introductory
vignette. See `vignette("ritest")`.

As a quickstart, here follows a basic example using data from a
randomized control trial (RCT) that was conducted in Colombia. The
dataset is provided with this package.

First, we use the `fixest::feols()` function to estimate the parametric
model.

``` r
library(ritest)  ## This package
library(fixest)  ## For fast (high-dimensional) fixed-effect regressions

data("colombia")

## Parametric model using fixest::feols()
co_est = 
  feols(
    dayscorab ~ b_treat + b_dayscorab | b_pair + miss_b_dayscorab + round2 + round3, 
    vcov = ~b_block, data = colombia
    )
co_est
#> OLS estimation, Dep. Var.: dayscorab
#> Observations: 2,346 
#> Fixed-effects: b_pair: 31,  miss_b_dayscorab: 2,  round2: 2,  round3: 2
#> Standard-errors: Clustered (b_block) 
#>              Estimate Std. Error  t value  Pr(>|t|)    
#> b_treat     -0.180738   0.078174 -2.31201  0.024113 *  
#> b_dayscorab  0.524761   0.029423 17.83478 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> RMSE: 1.91167     Adj. R2: 0.282038
#>                 Within R2: 0.21367
```

Our key treatment variable (`b_treat`) is deemed to be statistically
significant (p-value of 0.024), even as we cluster the standard errors.

But let’s see if this result is robust to randomization inference (RI).
We’ll perform 1,000 RI permutations on `b_treat`, whilst taking into
account the stratified and clustered experimental design of the
underlying RCT.

``` r
ritest(co_est, 'b_treat', strata='b_pair', cluster='b_block', reps=1e3, seed=1234)
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
#>  -0.1807       107      1000     0.107   0.01609   0.08054    0.1335  
#> ──────────────────────────────────────────────────────────────────────────────── 
#> Note: Confidence interval is with respect to p=c/n. 
#> Note: c = #{|T| >= |T(obs)|}
```

The simulated p-value is noticeably larger than the parametric one
(0.107 vs 0.024), suggesting that our parametric model is overstating
the effectiveness of treatment.

For more examples and additional features — plotting, regression tables,
etc. — please see the [**introductory
vignette**](http://grantmcdermott.com/ritest/articles/ritest.html).

## Benchmarks

I generally observe a speed increase of between 25x–50x compared to the
Stata version. Again, see the introductory vignette for timed examples.

## Other software

Apart from the Stata [`-ritest-`](https://github.com/simonheb/ritest)
routine, there are several other packages for conducting randomization
inference in R. For example, the
[**ri**](https://cran.r-project.org/web/packages/ri/index.html) has been
available for nearly a decade. More recently, the successor
[**ri2**](https://cran.r-project.org/web/packages/ri2/index.html)
package extends upon the original, with updated syntax and
functionality. An advantage of **ri2** is that it integrates with the
wider [DeclareDesign](https://declaredesign.org/) suite of R packages
for experimental and empirical research design. This enables researchers
to build RI and other experimental considerations into the incipient
design process. On the other hand, this places some restrictions on
conducting RI *ex post* or in quasi-experimental settings (e.g. a study
that leverages a natural experiment). For example, you can’t pass an
existing regression model object to `ri2::conduct_ri()`, which is what
`ritest::ritest()` was designed for. Your use case will likely determine
which software is optimal for you.

## License

The software code contained within this repository is made available
under the [MIT license](http://opensource.org/licenses/mit-license.php).
