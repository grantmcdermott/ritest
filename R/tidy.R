#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance


#' Tidy an `ritest` object
#'
#' @param x An object produced by the `ritest` function.
#' @param conf.int Logical indicating whether or not to include a confidence
#'   interval.
#' @param conf.level The confidence level to use for the confidence interval if
#'   `conf.int=TRUE`. Must be strictly greater than 0 and less than 1. Defaults
#'   to 0.95, which corresponds to a 95 percent confidence interval.
#' @inheritParams ritest
#' @return A "tidy" `data.frame` of summary statistics that conforms to the
#'   `broom` package specification.
#' @export
#' @examples
#' est = lm(yield ~ N + P + K, data = npk)
#' est_ri = ritest(est, 'N', reps = 1e3, seed = 1234L, pcores = 2L)
#' tidy(est_ri)
#' glance(est_ri)
tidy.ritest <- function(x,
                        conf.int = TRUE,
                        conf.level = 0.95,
                        ...) {

  ret = data.frame(term = x$resampvar,
                   estimate = x$beta_parm,
                   std.error = x$se,
                   p.value = x$pval)

  # confidence intervals
  if (conf.int) {
    ci = x$ci
    if (attr(x$pval, 'side')=='<=') ci[2] = Inf
    if (attr(x$pval, 'side')=='>=') ci[1] = -Inf
    ret$conf.low = ci[1]
    ret$conf.high = ci[2]
  }

  return(ret)
}


#' Glance at an `ritest` object
#'
#' @param x An object produced by the `ritest` function.
#' @return A "tidy" `data.frame` of goodness-of-fit statistics that conforms to
#'   the `broom` package specification.
#' @inherit tidy.ritest params examples
#' @export
glance.ritest <- function(x, ...) {
  ret = data.frame(H0 = x$h0, Num.Reps = as.integer(x$reps))
  return(ret)
}
