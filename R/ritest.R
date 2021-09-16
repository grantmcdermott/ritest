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
#' @import data.table
#' @export
#' @examples
#' library(fixest)
#' (mod = feols(yield ~ N + P + K | block, data = npk))
#' mod_ri = ritest('N1', mod, strata = ~block) ## We need 'N1' b/c that's the factored name in the model
#' mod_ri
#' plot(mod_ri)
#' plot(mod_ri, type = 'hist', highlight = 'fill', highlight_par = TRUE)
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

  # ## Split the treatment variable by the appropriate strata and/or cluster vars
  # if (!is.null(split_list)) {
  #   # Xtreat_split = split(Xmat, split_list)
  #   Xtreat_split = split(Xtreat, split_list)
  #   # Xtreat_split[sapply(Xtreat_split, function(x) length(x)==0)] = NULL
  # }

  if (!is.null(split_list)) {
    DT = data.table(cbind(Xtreat, strata_split, cluster_split))
    colnames(DT) = c('treat', c('strata', 'cluster')[!sapply(split_list, is.null)])
    DT[, orig_order := seq_len(.N)]
    # Add dummy column in case where clustering without strata
    if (is.null(strata_split)) {
      DT[, strata := TRUE]
    }
  }

  betas =
    sapply(
      1:reps,
      function(i) {
        if (!is.null(split_list)) {
          if(is.null(cluster)) {
            ## Base
            # Xtreat_samp = unsplit(lapply(Xtreat_split, sample), split_list)
            # Might as well use data.table
            Xtreat_samp = DT[DT[ , .I[sample(.N,.N)] , by = strata]$V1, treat]
            Xtreat_samp = as.matrix(Xtreat_samp)
          } else {

            DT$rorder = runif(nrow(DT))
            DT[, ind := rowid(strata, cluster)==1L] ## Get first obs by strata+cluster
            DT[(ind), nn := rowid(strata, ind)]
            setorder(DT, strata, ind, rorder) ## shuffle
            DT[(ind), newt := treat[nn], by = .(strata, ind)]
            setorder(DT, strata, cluster, -ind, na.last = TRUE)
            DT[, newt := nafill(newt, type = "locf"), by=.(strata, cluster)]
            setorder(DT, orig_order) ## Back to original order for fitting
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


#' @title A print method for ritest objects
#' @name print.ritest
#' @description Printed display of ritest objects. Tries to mimic the display
#'   of the equivalent Stata routine.
#' @param x An ritest object.
#' @param ... Currently ignored.
#' @export
print.ritest = function(x, ...) {

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
#' @param break Character. Histogram plot only. What type of breaks do you want?
#'   The default method creates more breaks than the standard R behaviour. You
#'   can revert to the latter by selecting NULL.
#' @param highlight Character. How do you want to highlight the H0 rejection
#'   regions in the distribution tails?
#' @param highlight_par Logical. Should we highlight the parametric H0 rejection
#'   regions too?
#' @param family Character. The font family. Defaults to HersheySans instead of
#'   R's normal Arial plotting font.
#' @param ... Other plot arguments. Currently ignored.
#' @export
plot.ritest = function(x, type = c('density', 'hist'), breaks = 'auto',
                       highlight = c('lines', 'fill', 'both', 'none'),
                       highlight_par = FALSE,
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

