#' @title Perform randomization inference on a model object
#' @aliases ritest
#'
#' @description Perform randomization inference (RI) testing on a model object, e.g. a
#'   coefficient from a linear regression model. It tries to mimic the `-ritest-`
#'   Stata routine (Heß, 2017) in its design and functionality. The package is
#'   quite experimental and only a subset of this functionality is currently
#'   supported. However, it does appear to be significantly faster.
#'
#' @param object Model object containing the `resampvar` variable. At present,
#'   only `stats::lm` and `fixest::feols` models are supported.
#' @param resampvar Character or one-sided formula. The variable (coefficient)
#'   that you want to perform RI on.
#'   By default, the RI procedure will conduct a standard two-sided test against
#'   a sharp null hypothesis of zero (i.e. H0: resampvar = 0). Other null
#'   hypotheses may be specified as part of a character string. These must
#'   take the form of the resampled variable, followed by an (in)equality sign,
#'   and then a number. For example 'x>=0' or 'x<=0' would yield respective
#'   one-sided tests against the zero null. Similarly, you can test against
#'   values other than zero (e.g. 'x=1' or 'x>=1'). However, note that multiple
#'   multiple comparison tests (e.g. 'x1-x2=0') are not yet supported. See
#'   Examples.
#' @param reps Integer. The number of repetitions (permutation draws) in the RI
#'   simulation. Default is 100, but you probably want more that that. Young
#'   (2019) finds that rejection rates stabilise at around 2,000 draws.
#' @param strata Character or one-sided formula. Permute `resampvar` within
#'   strata (AKA blocks)? See Details and Examples below.
#' @param cluster Character or one-sided formula. Keep `resampvar` constant
#'   within clusters? See Details and Examples below.
#' @param level Numeric. The desired confidence level. Default if 0.95.
#' @param stack Logical. Should the permuted data be stacked in memory all at
#'   once, rather than being recalculated during each iteration? Stacking takes
#'   advantage of vectorisation and is thus more efficient. (It also helps to
#'   ensure reproducibility when the results are being generated in parallel ---
#'   see the note on random-number generation below.) But it does
#'   require additional memory. If no explicit choice is provided, then the
#'   function will automatically stack as long this implies an additional memory
#'   overhead less than the `stack_lim` argument. Note that stacking is only
#'   relevant if at least one of `strata` or `cluster` are defined.
#' @param stack_lim Numeric. What is the memory limit (in gigabytes) for
#'   determining whether the permuted data should be stacked in memory ahead of
#'   time? Default is 1 GB. You probably want to increase this if you are
#'   working with a large dataset and have lots of RAM.
#' @param parallel Logical. Should the permuted fits be executed in parallel?
#'   Default is TRUE, with additional options being passed to the `ptype` and
#'   `pcores` arguments.
#' @param ptype Character. What type of parallel strategy should be used? The
#'   default behaviour on Linux and Mac is parallel forking ("fork"), while on
#'   Windows it will revert to parallel sockets ("psock"). Note that forking
#'   is more efficient, but unavailable on Windows.
#' @param pcores Integer. How many parallel cores should be used? If none is
#'   provided, then it will default to half of the total available CPU cores
#'   on the user's machine.
#' @param seed Integer. Random seed for reproducible results. Note that the
#'   choice of parallel behaviour can alter results even when using the same
#'   seed. See the Note on random number generation below.
#' @param pb Logical. Display a progress bar? Default is FALSE. Progress bars
#'   can add a surprising amount of computational overhead to iterative
#'   functions. I've therefore set the number of updating steps --- i.e. via the
#'   `nout` argument of [pbapply::pboptions()] --- to 5, which should limit this
#'   kind of overhead in the event that a user invokes the `pb = TRUE` argument.
#' @param verbose Logical. Display the underlying model `object` summary and
#'   `ritest` return value? Default is `FALSE`.
#' @param ... Additional arguments. Currently ignored.
#'
#' @details This function is experimental and functionality is still quite
#'   limited. Albeit, that it does support the most likely use case for RI on a
#'   regression model, i.e. permutation testing of a coefficient value. Present
#'   limitations include: only `lm` and `fixest::feols` model objects are
#'   supported; only one permutation (RI) test is allowed; and only one strata
#'   and/or cluster variable, respectively, can be supplied. I hope to resolve
#'   these limitations as time permits.
#'
#' @note The use of parallelism in computation introduces some well-known
#'   complications with respect to random number generation (RNG). See the
#'   vignette accompanying the `parallel` package for more discussion. The
#'   `ritest()` function adopts various best practices to facilitate reliably
#'   reproducible results, irrespective of the choice of parallel strategy.
#'   These include the use of "L'Ecuyer-CMRG" RNG kind and the default behaviour
#'   of "stacking" the permuted results in memory before passing them on to the
#'   fitting stage of the randomization inference routine (where the fitting would
#'   still be done in parallel). To briefly dwell on the latter, note that this
#'   stacking behaviour is very similar to --- I nearly wrote "parallels" ---
#'   the approach adopted by `boot` and other packages that trade off RNG with
#'   parallelism and reproducibility. In the phrasing
#'   of the `parallel` vignette: "One way to avoid any difficulties is (where
#'   possible) to do all the randomization in the master process" (p. 5).
#'
#'   The upshot is that running `ritest()` under a given set of arguments should
#'   generally yield the same, reproducible results regardless of when or where
#'   they were run. This should be true even if users turn the default parallel
#'   behaviour of the function off or back on, or whether they change the
#'   `ptype` argument to "psock" or "fork" and back again. However, it cannot be
#'   guaranteed under every scenario. Perhaps the most obvious case is when
#'   `ritest()` is run without any declared strata or clusters. (Reason:
#'   Stacking in this simple case is turned off because it yields no efficiency
#'   gains.) In this case, the RNG stream will be sensitive to the _number_ of
#'   available cores on a user's computer. Two parallel cores will yield a
#'   slightly different result than four cores, or six cores, etc. Of course,
#'   users can always ensure perfect reproducibility by explicitly defining the
#'   number of required cores in the `ritest()` call via the `pcores` argument.
#'   Stepping back, the focus on reproducible exactitude rather misses the
#'   point in an exercise like randomization inference testing. Much more important is
#'   running enough permutation trials so that your rejection rates are stable
#'   and any minor differences due to different RNG seeds are moot. Remember,
#'   the ultimate goal of inference (and research) isn't simply to generate
#'   reproducible results under a very specific set of circumstances. Rather, it
#'   is to generate consistent insights that hold even under varying
#'   circumstances.
#'
#'
#' @references Simon Heß (2017). \cite{Robust Randomization inference with
#'   Stata: A guide and software}, The Stata Journal, 17, Number 3, pp. 630--651
#' @references Alwyn Young (2019). \cite{Channeling Fisher: Randomization
#'   Tests and the Statistical Insignificance of Seemingly Significant
#'   Experimental Results}, The Quarterly Journal of Economics, 134, Issue 2,
#'   pp. 557--598
#' @return An list object of class `ritest`. Default print and plotting methods
#'   are supported.
#' @seealso [print.ritest()], [plot.ritest()]
#' @import data.table
#' @importFrom pbapply pbsapply
#' @importFrom stats .lm.fit coefficients complete.cases confint model.matrix qnorm sd
#' @importFrom utils capture.output head tail
#' @export
#' @examples
#' #
#' ## Example 1: Basic functionality
#' #
#'
#' # First estimate a simple linaer regression on the base 'npk' dataset. For
#' # this first example, we won't worry about strata or clusters, or other
#' # experimental design complications.
#' est = lm(yield ~ N + P + K, data = npk)
#'
#' # Conduct RI on the 'N' (i.e. nitrogen) coefficient. We'll do 1,000
#' # simulations and, just for illustration, limit the number of parallel cores
#' # to 2 (default is half of the available cores). The 'verbose = TRUE'
#' # argument simply prints the results upon completion, including the original
#' # regression model summary.
#' est_ri = ritest(est, 'N', reps = 1e3, seed = 1234L, verbose = TRUE)
#'
#' # Result: The RI rejection rate (0.021) is very similar to the parametric
#' # p-value (0.019).
#'
#' # We can plot the results and various options are available to customise the appearance.
#' plot(est_ri)
#' plot(est_ri, type = 'hist')
#' # etc
#'
#' # Aside: By default, ritest() conducts a standard two-sided test against a
#' # sharp null hypothesis of zero. You can can specify other null hypotheses as
#' # part of the 'resampvar' string argument. For example, a (left) one-sided
#' # test...
#' plot(ritest(est, 'N<=0', reps = 1e3, seed = 1234L, pcores = 2L))
#' # ... or, null values different from zero.
#' plot(ritest(est, 'N=2', reps = 1e3, seed = 1234L, pcores = 2L))
#'
#'
#' #
#' ## Example 2: Real-life example
#' #
#'
#' # Now that we've seen the basic functionality, here's a more realistic RI
#' # example using data from a randomized control trial conducted in Colombia.
#' # More details on the dataset -- kindly provided by the study authors -- can
#' # be found in the accompanying helpfile ("?colombia"). The most important
#' # thing to note is that we need to control for the stratified (aka "blocked")
#' # and clustered experimental design.
#'
#' data("colombia")
#'
#' # We'll use the fixest package to estimate our parametric regression model,
#' # specifying the strata (here: treatment-control pairs) as fixed-effects and
#' # clustering the standard errors by location (here: city blocks).
#' library(fixest)
#' co_est = feols(dayscorab ~ b_treat + b_dayscorab + miss_b_dayscorab |
#'                  b_pair + round2 + round3,
#'                vcov = ~b_block, data = colombia)
#' co_est
#'
#' # Run RI on the 'b_treat' variable, specifying the strata and clusters.
#' co_ri = ritest(co_est, 'b_treat', strata='b_pair', cluster='b_block',
#'                reps=1e3, seed=123L)
#' co_ri
#' plot(co_ri, type = 'hist', highlight = 'fill')
#'
#' # This time, the RI rejection rate (0.11) is noticeably higher than the
#' # parametric p-value (0.024) from the regression model.
#'
ritest = function(object,
                  resampvar,
                  reps = 100,
                  strata = NULL,
                  cluster = NULL,
                  # fixlevels = NULL,
                  level = 0.95,
                  parallel = TRUE,
                  ptype = c('auto', 'fork', 'psock'),
                  pcores = NULL,
                  stack = NULL,
                  stack_lim = 1L,
                  seed = NULL,
                  pb = FALSE,
                  verbose = FALSE,
                  ...) {

  ## Silence NSE notes in R CMD check. See:
  ## https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html#globals
  orig_order = .ii = treat = nn = treat_samp = NULL

  ## Temp vars that we'll be using later
  DT = DATA = strata_split = cluster_split = split_list = cl = NULL

  if (!pb) {
    opb = pbapply::pboptions(type = if (interactive()) "timer" else "none")
    pbapply::pboptions(type = 'none')
    on.exit(pbapply::pboptions(opb))
  }

  if (inherits(resampvar, "formula")) {
    resampvar = paste0(resampvar)[2]
    h0_symbol = '='
    h0_value = 0L
  } else {
    resampvar = gsub('[[:space:]]', '', resampvar)
    h0 = resampvar
    resampvar = gsub('=.*|>.*|<.*', '', resampvar)
    if (resampvar==h0) {
      h0_symbol = '='
      h0_value = 0L
    } else {
      h0_value = as.numeric(gsub('.*=|.*>|.*<', '', gsub(resampvar, '', h0)))
      h0_symbol = sub(h0_value, '', gsub(resampvar, '', h0, fixed = TRUE))
      if (!grepl('=$', h0_symbol)) h0_symbol = paste0(h0_symbol, '=')
    }
  }

  if(!is.null(seed)) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
  }

  fixest_obj = inherits(object, c('fixest', 'fixest_multi'))

  if (inherits(object, c('lm'))) {
    Ymat = object$model[, 1]
  } else if (fixest_obj) {
    Ymat = model.matrix(object, type = 'lhs', as.matrix = TRUE)
  } else {
    stop('\nModel or object class not currently supported. See help documentation.\n')
  }

  object_summ = summary(object, lean = TRUE)
  object_summ_string = paste(capture.output(object_summ), collapse = '\n')
  call_string = paste(deparse(object$call), collapse = '')
  call_string = gsub("\"", "\'", trimws(gsub("\\s+", " ", call_string)))

  Xmat = model.matrix(object)

  fmat = Xmat_dm = Ymat_dm = NULL
  if (fixest_obj && !is.null(object$fixef_vars)) {
    fmat = model.matrix(object, type = 'fixef')
    Xmat_dm = fixest::demean(Xmat, fmat)
    Ymat_dm = fixest::demean(Ymat, fmat)
  }

  Xnames = colnames(Xmat)
  resampvar_pos = match(resampvar, Xnames)
  resampvar_orig = NULL
  if (is.na(resampvar_pos)) {
    resampvar_orig = resampvar
    resampvar = paste0(resampvar, 1) ## Catch for single-level numeric factors.
    resampvar_pos = match(resampvar, Xnames)
    if (is.na(resampvar_pos)) {
      if (verbose) cat(object_summ_string)
      stop('\n\nThe specified resampling variable "', resampvar,
           '" was not found among the model explanatory variables:\n\n',
           call_string,
           '\n\nIf you transformed this variable directly in the model call, ',
           'e.g. `I(x^2)`, then please use the transformed string as a verbatim ',
           'argument instead, e.g. `ritest(resampvar = "I(x^2)", ...)`.',
           '\n\nSimilarly, if you want to resample a factor variable with more ',
           'than two levels then you should specify the level that you want, ',
           'e.g. `ritest(resampvar = "cyl6", ...)`. To see a list of all ',
           'possible coefficients and their verbatim names, consider re-running ',
           'the `ritest` function with the `verbose = TRUE` argument.')
      }
  }
  onames = setdiff(Xnames, resampvar) ## other (non-treatment vars)
  Xtreat = Xmat[,resampvar_pos]

  prep_split_var = function(x) {
    if (inherits(x, "formula")) {
      x_string = strsplit(paste0(x)[2], split = ' \\+ ')[[1]]
    } else {
      x_string = paste(x)
    }
    if (!is.null(fmat) && x_string %in% colnames(fmat)) {
      ret = fmat[, x_string]
    } else if (x_string %in% colnames(Xmat)) {
      ret = Xmat[, x_string]
    } else {
      if (is.null(DATA)) {
        DATA = eval(object$call$data)
        ## Assign to parent env to avoid redoing if possible
        assign('DATA', DATA)
      }
      if (x_string %in% names(DATA)) {
        all_vars = sapply(list(Ymat, Xmat, fmat), colnames)
        if (inherits(all_vars, 'list')) all_vars = do.call('c', all_vars)
        all_vars = union(union(all_vars, resampvar_orig), x_string)
        DATA = data.frame(DATA)[, intersect(colnames(DATA), all_vars)]
        DATA = DATA[complete.cases(DATA), ]
        ret = model.matrix(~0+., DATA[, x_string, drop=FALSE])
      } else {
        stop(paste0('Could not find ', x, '. Please provide a valid input.\n'))
      }
    }
    attr(ret, 'string') = x_string
    return(ret)
  }

  if (!is.null(strata)) {
    strata_split = prep_split_var(strata)
    strata = attr(strata_split, 'string')
    attr(strata, 'levels') = length(unique(strata_split))
    }
  if (!is.null(cluster)) {
    cluster_split = prep_split_var(cluster)
    cluster = attr(cluster_split, 'string')
    attr(cluster, 'levels') = length(unique(cluster_split))
    }

  ## Next line mostly for indexing
  if (!is.null(strata_split) || !is.null(cluster_split)) {
    split_list = list(strata = strata_split, cluster = cluster_split)
  }

  if (!is.null(split_list)) {
    DT = data.table(cbind(Xtreat, strata_split, cluster_split))
    colnames(DT) = c('treat', c('strata', 'cluster')[!sapply(split_list, is.null)])
    DT[, orig_order := seq_len(.N)]
    # Add dummy column in cases where we are clustering without strata
    if (is.null(strata_split)) {
      DT[, strata := TRUE]
    }

    DT_prep = function(DT) {
      data.table::setkey(DT, .ii, strata)
      if(is.null(cluster)) {
        DT_prepped = DT[,
                        list(orig_order, treat_samp = treat[sample.int(.N)]),
                        by = list(.ii, strata)]
        data.table::setorder(DT_prepped, .ii, orig_order) ## Back to original order for fitting
      } else {
        DT2 = DT[data.table::rowid(.ii, strata, cluster)==1L]
        DT2[, treat_samp := treat[sample.int(.N)], by = list(.ii, strata)]
        DT_prepped = DT[DT2[, list(.ii, strata, cluster, treat_samp)],
                        on = list(.ii, strata, cluster)]
        data.table::setorder(DT_prepped, .ii, orig_order) ## Back to original order for fitting
        DT_prepped = DT_prepped[, list(.ii, treat_samp)]#; gc()
        data.table::setkey(DT_prepped, .ii)
      }
      DT_prepped = DT_prepped[,list(Xtreat_samp = list(as.matrix(treat_samp))), by=.ii]
      return(DT_prepped)
    }

    if (is.null(stack) || stack) {
      sdict = c('logical' = 4, 'integer' = 4, 'numeric' = 8)
      ## Need to add an extra 4 byte (integer) variable for .ii index column
      stacked_DT_size = mean(c(sdict[sapply(DT, class)], 4))*NROW(DT)*NCOL(DT) * reps / 1e9
      stack = stacked_DT_size < stack_lim
    }

    ## stacked?
    if (stack) {
      if (verbose) cat('\nPermuted data of approximate size',
                       sprintf("%.2g", stacked_DT_size),
                       'GB being stacked in memory.')
      DT = rbindlist(lapply(1:reps, function(ii) {DT; DT$.ii = ii; DT}))
      DT = DT_prep(DT)
    }

  } else {
    ## No stacking if neither strata nor cluster vars are defined
    stack = FALSE
  }

  ## Simulation function
  ri_sims =
    function(i) {

      ## First get our appropriately permuted treatment variable
      if (is.null(split_list)) { ## Case 1. No strata or clustering vars. Also means cannot be stacked.
        ## Sample the treatment variable
        Xtreat_samp = sample(Xtreat)
      } else { ## Case 2. Strata and/or cluster vars defined. Can also be stacked.
        if (!stack) {
          ## For non-stacked data we'll prep on the fly during each iteration
          DT$.ii = i
          DT = DT_prep(DT)
        }
        ## Sample the treatment variable
        Xtreat_samp = DT[.ii==i, Xtreat_samp[[1]]]
      }

      ## Then fit on the permuted data
      if (!is.null(fmat)) {
        ## We need to FE demean the sampled treatment vector
        Xtreat_samp_dm = fixest::demean(Xtreat_samp, fmat, nthreads = 1)
        beta_samp = coefficients(.lm.fit(cbind(Xtreat_samp_dm, Xmat_dm[, onames]), Ymat_dm))[1]
      } else {
        beta_samp = coefficients(.lm.fit(cbind(Xtreat_samp, Xmat[, onames]), Ymat))[1]
      }

      return(beta_samp)
    }


  ## Parallel strategy
  if (parallel) {

    ptype = match.arg(ptype)
    if (ptype=='auto') {
      if (.Platform$OS.type == "windows") {
        ptype = 'psock'
      } else {
        ptype = 'fork'
      }
    }

    if (is.null(pcores)) {
      pcores = ceiling(parallelly::availableCores()/2L)
    } else if (!is.integer(pcores)) {
      pcores = as.integer(pcores)
    }
    if (is.na(pcores) || pcores > parallelly::availableCores()) {
      message("Invalid number of parallel cores selected. Reverting to default behaviour.")
      pcores = ceiling(parallelly::availableCores()/2L)
    }

    if (ptype=='fork') {

      cl = pcores

      if (verbose) cat("\nRunning", reps, "parallel RI simulations as forked",
                       "processes across", pcores, "CPU cores.")

    } else {

      if (!stack) {
        ## Experimenting, I find a sharp decrease in psock performance with
        ## pcores > 4 when the data aren't stacked and we require that sampling
        ## remain constant within clusters. I'm not entirely sure why (probably
        ## (something to do with the overhead), so we'll set this as the maximum
        ## level to safeguard.
        if (!is.null(cluster)) pcores = min(pcores, 4L)
        setDTthreads(1L)        ## avoid nested parallelism
        on.exit(setDTthreads(), add = TRUE) ## reset after function exists
      }

      cl = parallel::makeCluster(pcores)
      on.exit(parallel::stopCluster(cl), add = TRUE)

      parallel::clusterEvalQ(cl, c(library(data.table)))
      parallel::clusterExport(
        cl, c("DT", "split_list", "Xtreat", "stack", "Xmat", "Xmat_dm",
              "Ymat", "Ymat_dm", "fmat", "onames"),
        envir=environment()
      )
      parallel::clusterSetRNGStream(cl, seed)

      if (verbose) cat("\nRunning", reps, "parallel RI simulations in a PSOCK",
                       "cluster across", pcores, "CPU cores.")

      }

  } else {

    if (verbose) cat("\nRunning", reps, "RI simulations sequentially.")

  }

  ## Run the sims
  betas = pbapply::pbsapply(1:reps, ri_sims, cl = cl)

  ## Parametric values
  h0 = paste0(resampvar, h0_symbol, h0_value)
  ci_sides = 1:2
  if (h0_value==0 && h0_symbol=='=') {
    if (fixest_obj) {
      coeftab = fixest::coeftable(object)
    } else {
      coeftab = summary(object)$coefficients
    }
    ci_parm = confint(object, level = level)
    ci_parm = ci_parm[rownames(ci_parm)==resampvar,]
  } else {
    glht_summ =
      if (fixest_obj) {
        ## Need two minor adjustments for fixest objects, since there's
        ## currently no dedicated glht.fixest method:
        ## 1) Suppress harmless (but annoying) warning message about
        ## unrecognised complete = FALSE argument (passed via ...)
        ## 2) Need manual DoF override to get right SEs / CIs.
        summary(suppressWarnings(
          multcomp::glht(object, h0,
                         df=fixest::degrees_freedom(object, type = 'resid'), ...)
          ))
        } else {
          summary(multcomp::glht(object, h0, ...))
        }
    coeftab = data.frame(glht_summ$test[c("coefficients", "sigma", "tstat", "pvalues")])
    ci_parm = confint(glht_summ)$confint[, c('lwr', 'upr')[ci_sides]]
  }
  # beta_parm = coeftab[resampvar, 'Estimate']
  # pval_parm = coeftab[resampvar, 'Pr(>|t|)']
  beta_parm = coeftab[resampvar, 1]
  pval_parm = coeftab[resampvar, 4]


  if (h0_symbol=='=') {
    probs = abs(betas) >= abs(beta_parm - h0_value)
  # } else {
  #   probs = eval(parse(text = paste(betas, h0_symbol, beta_parm, '-', h0_value)))
  # }
  } else if (h0_symbol=='<=') {
    ## Note: we have to reverse the inequality
    probs = betas >= beta_parm - h0_value
  } else {
    probs = betas <= beta_parm - h0_value
  }
  count = sum(probs)
  pval = count/reps
  attr(pval, 'side') = h0_symbol
  se = qnorm(level) * sd(probs) / sqrt(reps)
  ci = qnorm(level) * se
  ci = pval + c(-ci, ci)[ci_sides]
  alpha = 1-level
  names(ci) = paste0('CI ', c(alpha/2, 1-alpha/2)*100, '%')[ci_sides]

  out = list(call = call_string,
             resampvar = resampvar,
             h0 = h0,
             reps = reps,
             strata = strata,
             cluster = cluster,
             count = count,
             pval = pval,
             se = se,
             ci = ci,
             betas = betas,
             object_summ_string = object_summ_string,
             pval_parm = pval_parm,
             beta_parm = beta_parm,
             ci_parm = ci_parm)

  class(out) = "ritest"

  if (verbose) print(out, verbose = TRUE)

  return(out)
}

