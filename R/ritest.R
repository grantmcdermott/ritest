#' @title Perform random inference on a model object
#' @aliases ritest
#'
#' @description Perform random inference (RI) testing on a model object, e.g. a
#'   coefficient from a linear regression model. It tries to mimic the `ritest`
#'   Stata routine (Heß, 2017) in its design and functionality. The package is
#'   quite experimental and only a subset of this functionality is currently
#'   supported. However, it does appear to be significantly faster.
#'
#' @param resampvar Character. The name of the variable (coefficient) that you want
#'   to test. At present, only a single variable is permitted.
#' @param object Model object containing the `resampvar` variable. At present,
#'   only `stats::lm` and `fixest::feols` models are supported.
#' @param reps Integer. The number of repetitions (permutations) in the RI
#'   simulation. Default is 100, but you probably want more that that. (Alwyn
#'   Young has suggested at least 2000 in a research setting.)
#' @param strata Character or one-sided formula. Permute `resampvar` within
#'   strata? See Details and Examples below.
#' @param cluster Character or one-sided formula. Keep `resampvar` constant
#'   within clusters? See Details and Examples below.#'
#' @param pvals Character. How should the test values should be computed? The
#'   default is "both", which means that two-sided p-values will be computed.
#'   Alternatively, users can specify one-sided p-values with either "left" or
#'   "right".
#' @param level Numeric. The desired confidence level. Default if 0.95.
#' @param stack Logical. Should the permuted data be stacked in memory all at
#'   once, rather than being recalculated during each iteration? Stacking takes
#'   advantage of vectorisation and is thus more efficient. But it does
#'   require additional memory. If no explicit choice is provided, then the
#'   function will automatically stack as long this implies an additional memory
#'   overhead less than the `stack_lim` argument. (Note that stacking is only
#'   relevant if at least one of `strata` or `cluster` are defined.)
#' @param stack_lim Numeric. What is the memory limit (in gigabytes) for
#'   determining whether the permuted data should be stacked in memory ahead of
#'   time? Default is 1 GB.
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
#' @param seed Integer. Random seed for reproducible results.
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
#' @references Simon Heß (2017)
#' \cite{Robust Randomization inference with Stata: A guide and software}, The
#' Stata Journal, 17, Number 3, pp. 630--651
#' @return An list object of class `ritest`. Default print and plotting methods
#'   are supported.
#' @seealso [print.ritest()], [plot.ritest()]
#' @import data.table
#' @export
#' @examples
#' library(fixest)
#' mod = feols(yield ~ N + P + K | block, data = npk)
#'
#' # Conduct RI on the 'N' (nitrogen) coefficient. We'll do 1,000 simulations
#' # and permute within the 'block' variable strata. The 'verbose = TRUE'
#' # argument simply prints the results upon completion.
#' mod_ri = ritest('N', mod, reps = 1e3, strata = ~block, verbose = TRUE)
#'
#' # We can also plot the results and various options are available to
#' # customise the appearance.
#' plot(mod_ri)
#' plot(mod_ri, highlight_par = TRUE) ## Add the parametric CI lines
ritest = function(resampvar,
                  object,
                  reps = 100,
                  pvals = c("both", "left", "right"),
                  strata = NULL,
                  cluster = NULL,
                  # fixlevels = NULL,
                  # h0 = NULL,
                  level = 0.95,
                  parallel = TRUE,
                  ptype = c('auto', 'fork', 'psock'),
                  pcores = NULL,
                  stack = NULL,
                  stack_lim = 1L,
                  seed = NULL,
                  verbose = FALSE,
                  ...) {

  pvals = match.arg(pvals)

  if(!is.null(seed)) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    }

  if (inherits(object, c('lm'))) {
    # Ymat = object$model[, 1, drop = FALSE]
    Ymat = object$model[, 1]
  } else if (inherits(object, c('fixest', 'fixest_multi'))) {
    Ymat = model.matrix(object, type = 'lhs', as.matrix = TRUE)
  } else {
    stop('\nModel or object class not currently supported. See help documentation.\n')
  }

  object_summ = summary(object, lean = TRUE)
  object_summ_string = paste(capture.output(object_summ), collapse = '\n')
  call_string = paste(deparse(object$call), collapse = '')
  call_string = gsub("\"", "\'", trimws(gsub("\\s+", " ", call_string)))

  Xmat = model.matrix(object)

  fmat = NULL
  if (inherits(object, c('fixest', 'fixest_multi')) && !is.null(object$fixef_vars)) {
    fmat = model.matrix(object, type = 'fixef')
    Xmat_dm = fixest::demean(Xmat, fmat)
    Ymat_dm = fixest::demean(Ymat, fmat)
  }

  Xnames = colnames(Xmat)
  resampvar_pos = match(resampvar, Xnames)
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

  DATA = NULL
  strata_split = NULL
  cluster_split = NULL
  split_list = NULL
  # Xtreat_split = NULL

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
        all_vars = do.call('c', sapply(list(Ymat, Xmat, fmat), colnames))
        all_vars = union(all_vars, x_string) ## probably redundant
        DATA = eval(object$call$data)
        DATA = data.frame(DATA)[, all_vars]
        DATA = DATA[complete.cases(DATA), ]
        ## Assign to parent env to avoid redoing if possible
        assign('DATA', DATA)
      }
      if (x_string %in% names(DATA)) {
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
    }
  if (!is.null(cluster)) {
    cluster_split = prep_split_var(cluster)
    cluster = attr(cluster_split, 'string')
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
        DT_prepped = DT[DT[ , .I[sample(.N,.N)] , by = .(.ii, strata)]$V1,
                        .(treat_samp = treat),
                        keyby = .ii]
      } else {
        DT2 = DT[data.table::rowid(.ii, strata, cluster)==1L]
        DT2[, nn := data.table::rowid(.ii, strata)]
        DT2[, treat_samp := treat[sample(nn)], by = .(.ii, strata)]
        DT_prepped = DT[DT2[, .(.ii, strata, cluster, treat_samp)],
                        on = .(.ii, strata, cluster)]
        data.table::setorder(DT_prepped, .ii, orig_order) ## Back to original order for fitting
        DT_prepped = DT_prepped[, .(.ii, treat_samp)]#; gc()
        data.table::setkey(DT_prepped, .ii)
      }
      DT_prepped = DT_prepped[,.(Xtreat_samp = list(as.matrix(treat_samp))), by=.ii]
      return(DT_prepped)
    }

    if (is.null(stack) || stack) {
      sdict = c('integer' = 4, 'numeric' = 8)
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

  ## Run the simulations
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

    if (ptype=='psock') {
      if (!stack) {
        ## Experimenting, I find a sharp decrease in psock performance with
        ## pcores > 4 when the data aren't stacked and we require that sampling
        ## remain constant within clusters. I'm not entirely sure why (probably
        ## (something to do with the overhead), so we'll set this as the maximum
        ## level to safeguard.
        if (!is.null(cluster)) pcores = min(pcores, 4L)
        setDTthreads(1L) ## avoid nested parallelism
        }
      cl = parallel::makeCluster(pcores)
      parallel::clusterEvalQ(cl, c(library(data.table)))
      parallel::clusterExport(
        cl, c("DT", "split_list", "Xtreat", "stack", "Xmat", "Xmat_dm",
              "Ymat", "Ymat_dm", "fmat", "onames"),
        envir=environment()
        )
      if (verbose) cat("\nRunning", reps, "parallel RI simulations in a PSOCK",
                       "cluster across", pcores, "CPU cores.")
      betas = parallel::parLapply(cl = cl, 1:reps, ri_sims)
      parallel::stopCluster(cl)
      if (!stack) setDTthreads()
    } else {
      if (verbose) cat("\nRunning", reps, "parallel RI simulations as forked",
                       "processes across", pcores, "CPU cores.")
      betas = parallel::mclapply(1:reps, ri_sims, mc.cores = pcores)
    }
    betas = simplify2array(betas)

  } else {
    if (verbose) cat("\nRunning", reps, "RI simulations sequentially.")
    betas = sapply(1:reps, ri_sims)
  }

  ## Parametric values
  beta_par = coeftable(object)[resampvar, 'Estimate']
  pval_par = coeftable(object)[resampvar, 'Pr(>|t|)']

  if (pvals=='both') {
    probs = abs(betas) >= abs(beta_par)
  } else if (pvals=='left') {
    probs = betas <= beta_par
  } else {
    probs = betas >= beta_par
  }
  count = sum(probs)
  pval = count/reps
  attr(pval, 'side') = pvals
  # se = sd(probs)/sqrt(reps)
  # ci = quantile(probs, abs(c(0, 1) - alpha))
  se = qnorm(level) * sd(probs) / sqrt(reps)
  ci = qnorm(level) * se
  ci = pval + c(-ci, ci)
  alpha = 1-level
  names(ci) = paste0('CI ', c(alpha/2, 1-alpha/2)*100, '%')
  ci_par = confint(object)
  ci_par = ci_par[rownames(ci_par)==resampvar,]

  out = list(call = call_string,
             resampvar = resampvar,
             reps = reps,
             strata = strata,
             cluster = cluster,
             count = count,
             pval = pval,
             se = se,
             ci = ci,
             betas = betas,
             object_summ_string = object_summ_string,
             pval_par = pval_par,
             beta_par = beta_par,
             ci_par = ci_par)

  class(out) = "ritest"

  if (verbose) print(out, verbose = TRUE)

  return(out)
}


#' @title A print method for ritest objects
#' @name print.ritest
#' @description Printed display of ritest objects. Tries to mimic the display
#'   of the equivalent Stata routine.
#' @param x An ritest object.
#' @param verbose Logical. Should we display the the original model summary too?
#'   Default is FALSE.
#' @param ... Currently ignored.
#' @export
print.ritest = function(x, verbose = FALSE, ...) {

  ri_mat = c(`T(obs)` = x$beta_par, `c` = x$count, `n` = x$reps,
             `p=c/n` = x$pval, `SE(p)` = x$se)
  ri_mat = c(ri_mat, x$ci)
  ri_mat_nms = names(ri_mat)
  ri_mat = sprintf(ri_mat, fmt = '%.4g')
  names(ri_mat) = ri_mat_nms
  pval_string = if (attr(x$pval, 'side')=='both') {
    "Note: c = #{|T| >= |T(obs)|}"
  } else if (attr(x$pval, 'side')=='left') {
    "Note: c = #{T <= T(obs)}"
  } else {
    "Note: c = #{T >= T(obs)}"
  }

  if (verbose) {
    cat("\n",
        "******************",
        "* ORIGINAL MODEL *",
        "******************",
        sep = "\n")
    cat("\n", x$object_summ_string, "\n\n", sep = "")
    cat("******************",
        "* RITEST RESULTS *",
        "******************",
        sep = "\n")
  }


  cat("\nCall: ", x$call, "\n", sep = "")
  cat("Res. var(s): ", x$resampvar, "\n", sep = "")
  cat("Strata var(s): ", x$strata, "\n", sep = "")
  cat("Cluster var(s): ", x$cluster, "\n", sep = "")
  cat("Num. reps: ", x$reps, "\n", sep = "")
  cat("---", "\n")
  print(ri_mat, quote = FALSE, print.gap = 2L)
  cat("---", "\n")
  cat("Note: Confidence interval is with respect to p=c/n", "\n")
  cat(pval_string, "\n")
  cat("\n")
  # invisible(x)
}


#' @title A plot method for ritest objects
#' @name plot.ritest
#' @description Nice plots of your ritest objects.
#' @param x An ritest object.
#' @param type Character. What type of plot do you want?
#' @param highlight Character. How do you want to highlight the H0 rejection
#'   regions in the distribution tails?
#' @param highlight_par Logical. Should we highlight the parametric H0 rejection
#'   regions too?
#' @param breaks Character. Histogram plot only. What type of breaks do you
#'   want? The default method creates more breaks than the standard R behaviour.
#'   You can revert to the latter by selecting NULL.
#' @param family Character. The font family. Defaults to 'HersheySans' instead
#'   of R's normal Arial plotting font.
#' @param ... Other plot arguments. Currently ignored.
#' @export
plot.ritest = function(x, type = c('density', 'hist'),
                       highlight = c('lines', 'fill', 'both', 'none'),
                       highlight_par = FALSE,  breaks = 'auto',
                       family = NULL, ...) {
  type = match.arg(type)
  highlight = match.arg(highlight)
  if (is.null(family)) family = 'HersheySans'
  if (type=='density') {
    dens = density(x$betas)
    plot(dens,
         main = '', xlab = '',
         xlim = range(x$betas, x$beta_par),
         family = family)
    if (highlight %in% c('both', 'fill')) {
      x1 = 1
      x2 = tail(which(dens$x <= -abs(x$beta_par)), 1)
      x3 = head(which(dens$x >= abs(x$beta_par)), 1)
      x4 = length(dens$x)
      fcol = rgb(1,0,0,0.5)
      with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=fcol, border = FALSE))
      with(dens, polygon(x=c(x[c(x3,x3:x4,x4)]), y= c(0, y[x3:x4], 0), col=fcol, border = FALSE))
    }
  } else {
    if (is.null(breaks)) {
      breaks = "Sturges"
    } else if (breaks=='auto') {
      nbreaks = min(length(unique(x$betas)), 101)
      breaks = seq(min(x$betas), max(x$betas), l=nbreaks)
    }

    hist_df = hist(x$betas, breaks = breaks, plot = FALSE)
    if (highlight %in% c('both', 'fill')) {
      hist_col = ifelse(abs(hist_df$breaks) > abs(x$beta_par), rgb(1,0,0,0.5), rgb(0.2,0.2,0.2,0.2))
    } else {
      hist_col = rgb(0.2,0.2,0.2,0.2)
    }
    plot(hist_df,
         col = hist_col,
         border = FALSE,
         main = NULL, xlab = NULL,
         xlim = range(x$betas, x$beta_par),
         family = family)
  }
  title(main = paste('Random Inference:', x$resampvar),
        sub = paste0('Simulated p-val: ',  sprintf('%.3g', x$pval),
                     '. Parametric p-val: ', sprintf('%.3g', x$pval_par),'.'),
        cex.sub = 0.75,
        xlab = 'Simulated values',
        family = family)
  if (highlight %in% c('both', 'lines')) {
    abline(v = x$beta_par, col = "red", lty = 1)
    abline(v = -x$beta_par, col = "red", lty = 1)
  }
  if (highlight_par) {
    abline(v = sweep(x$ci_par, 2, rowMeans(x$ci_par)), lty = 2, col = 'red')
  }
}

