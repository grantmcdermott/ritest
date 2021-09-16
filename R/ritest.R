#' @title Perform random inference on a model object
#' @aliases ritest
#'
#' @description Perform random inference (RI) testing on a model object, e.g. a
#'   coefficient from a linear regression model. It tries to mimic the `ritest`
#'   Stata routine (Hess, 2017) in its design and functionality. Only a subset
#'   of this functionality is currently supported. However, it does appear to be
#'   a lot faster.
#'
#' @param resampvar Character. The name of the variable (coefficient) that you want
#'   to test. At present, only a single variable is permitted.
#' @param object Model object containing the `resampvar` variable. At present,
#'   only `stats::lm` and `fixest::feols` models are supported.
#' @param reps Integer. The number of repetitions (permutations) in the RI
#'   simulation. Default is 100, but you probably want more that that. (Alwyn
#'   Young has suggested at least 2000 in a research setting.)
#' @param pvals Character. How should the p-values should be computed? The
#'   default is "both", which means that two-sided p-values will be computed.
#'   Alternatively, users can specify one-sided p-values with either "left" or
#'   "right".
#' @param strata Character or one-sided formula. Permute `resampvar` with strata
#'   strata? See Details and Examples below.
#' @param cluster Character or one-sided formula. Keep `resampvar` constant with
#'   clusters? See Details and Examples below.
#' @param level Numeric. The desired confidence level. Default if 0.95.
#' @param seed Integer. Random seed for reproducible results.
#' @param ... Additional arguments. Currently ignored.
#'
#' @details This function can be used with a regression or classification tree
#'   containing one or (at most) two numeric predictors.
#' @seealso [geom_parttree()], [rpart::rpart()], [partykit::ctree()].
#' @return A data frame comprising seven columns: the leaf node, its path, a set
#'   of coordinates understandable to `ggplot2` (i.e., xmin, xmax, ymin, ymax),
#'   and a final column corresponding to the predicted value for that leaf.
#' @importFrom data.table :=
#' @export
#' @examples
#' library(fixest)
#' (mod = feols(yield ~ N + P + K | block, data = npk))
#' ritest('N1', mod, strata = ~block) ## We need 'N1' b/c that's the factored name in the model
ritest = function(resampvar,
                  object,
                  reps = 100,
                  pvals = c("both", "left", "right"),
                  strata = NULL,
                  cluster = NULL,
                  # fixlevels = NULL,
                  # h0 = NULL,
                  level = 0.95,
                  seed = NULL,
                  ...) {

  pvals = match.arg(pvals)
  if(!is.null(seed)) set.seed(seed)

  if (inherits(object, c('lm'))) {
    # Ymat = object$model[, 1, drop = FALSE]
    Ymat = object$model[, 1]
  } else if (inherits(object, c('fixest', 'fixest_multi'))) {
    Ymat = model.matrix(object, type = 'lhs', as.matrix = TRUE)
  }

  Xmat = model.matrix(object)

  fmat = NULL
  if (inherits(object, c('fixest', 'fixest_multi')) && !is.null(object$fixef_vars)) {
    fmat = model.matrix(object, type = 'fixef')
    Xmat_dm = demean(Xmat, fmat)
    Ymat_dm = demean(Ymat, fmat)
  }

  Xnames = colnames(Xmat)
  # resampvar_orig = resampvar
  resampvar_pos = grep(resampvar, Xnames)
  onames = setdiff(Xnames, resampvar) ## other (non-treatment vars)
  Xtreat = Xmat[,resampvar_pos]

  DATA = NULL
  strata_split = NULL
  cluster_split = NULL
  split_list = NULL
  # Xtreat_split = NULL

  if (!is.null(strata)) {
    strata_split = prep_split_var(strata, Ymat=Ymat, Xmat=Xmat, fmat=fmat, DATA=DATA, object = object)
    strata = attr(strata_split, 'string')
    }
  if (!is.null(cluster)) {
    cluster_split = prep_split_var(cluster, Ymat=Ymat, Xmat=Xmat, fmat=fmat, DATA=DATA, object = object)
    cluster = attr(cluster_split, 'string')
    }

  if (!is.null(strata_split) || !is.null(cluster_split)) {
    split_list = list(strata = strata_split, cluster = cluster_split)
    split_list[sapply(split_list, is.null)] = NULL
  }

  # ## Split the treatment variable by the appropriate strata and/or cluster vars
  # if (!is.null(split_list)) {
  #   # Xtreat_split = split(Xmat, split_list)
  #   Xtreat_split = split(Xtreat, split_list)
  #   # Xtreat_split[sapply(Xtreat_split, function(x) length(x)==0)] = NULL
  # }

  if (!is.null(split_list)) {
  # if (!is.null(cluster)) {
    DT = data.table::data.table(cbind(Xtreat, strata_split, cluster_split))
    colnames(DT) = c(resampvar, strata, cluster)
    DT[, orig_order := seq_len(.N)]
  }

  betas =
    sapply(
      1:reps,
      function(i) {
        if (!is.null(strata) || !is.null(cluster)) {
        # if (!is.null(Xtreat_split)) {
          if(is.null(cluster)) {
            ## Base
            # Xtreat_samp = unsplit(lapply(Xtreat_split, sample), split_list)
            # Might as well use data.table
            Xtreat_samp = DT[DT[ , .I[sample(.N,.N)] , by = strata]$V1, ..resampvar]
            Xtreat_samp = as.matrix(Xtreat_samp)
          } else {

            DT$rorder = runif(nrow(DT))
            DT[data.table::rowidv(DT, c(strata, cluster))==1, ind := 1] ## Get first obs by strata+cluster
            data.table::setorderv(DT, c(strata, 'ind'), na.last = TRUE) ## Put these first obs at top w/in each strata
            DT[!is.na(ind), nn := data.table::rowidv(DT[!is.na(ind)], c(strata, 'ind'))]
            data.table::setorderv(DT, c(strata, 'ind', 'rorder'), na.last=TRUE) ## shuffle
            DT[!is.na(ind), newt := .SD[nn], by = c(strata, 'ind'), .SDcols = resampvar]
            data.table::setorderv(DT, c(strata, cluster, 'ind'), na.last = TRUE)
            DT[, newt := data.table::nafill(newt, type = "locf"), by = c(strata, cluster)]
            data.table::setorder(DT, orig_order) ## Back to original order for fitting
            Xtreat_samp = as.matrix(DT[, 'newt'])

            ## Base
            # Xtreat_samp =
            #   unsplit(lapply(Xtreat_split, function(x) {
            #     ifelse(length(x)<=1, x, replace(x, seq_along(x), sample(x, 1)))
            #   }),
            #   split_list)
          }
        } else {
          Xtreat_samp = sample(Xtreat)
        }
        if (is.null(fmat)) {
          coefficients(.lm.fit(cbind(Xtreat_samp, Xmat[, onames]), Ymat))[1]
        } else {
          Xtreat_samp_dm = demean(Xtreat_samp, fmat)
          coefficients(.lm.fit(cbind(Xtreat_samp_dm, Xmat_dm[, onames]), Ymat_dm))[1]
        }

      }
    )

  ## Parametric values
  # par_vals = coeftable(object)[rownames(coeftable(object))==resampvar, ]
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

  call_string = paste(deparse(object$call), collapse = '')
  call_string = gsub("\"", "\'", trimws(gsub("\\s+", " ", call_string)))

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
             pval_par = pval_par,
             beta_par = beta_par,
             ci_par = ci_par)

  class(out) = "ritest"

  return(out)
}

#' @title Prep split variables
#' @name prep_split_var
#' @description Internal function for prepping the split vars, i.e. strata
#'   and/or clusters. Will extract and assign the underlying data (used in the
#'   model object) to the parent environment if required.
#' @param x String or one-sided formula.
#' @param ... Additional helpers objects, e.g. from a model.matrix.
#' @NoRd
prep_split_var = function(x, ...) {
  dots = list(...)
  Ymat = dots$Ymat; Xmat = dots$Xmat; fmat = dots$fmat; DATA = dots$DATA; object = dots$object
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
