# Glance at an \`ritest\` object

Glance at an \`ritest\` object

## Usage

``` r
# S3 method for class 'ritest'
glance(x, ...)
```

## Arguments

- x:

  An object produced by the \`ritest\` function.

- ...:

  Additional arguments. Currently ignored.

## Value

A data frame of goodness-of-fit statistics that conforms to the
\`broom\` package specifications.

## Examples

``` r
est = lm(yield ~ N + P + K, data = npk)
est_ri = ritest(est, 'N', reps = 1e3, seed = 1234L, pcores = 2L)
tidy(est_ri)
#>   term estimate   std.error p.value    conf.low  conf.high
#> 1   N1 5.616667 0.007461833   0.021 0.008726377 0.03327362
glance(est_ri)
#>     H0 Num.Reps
#> 1 N1=0     1000
```
