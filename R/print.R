#' @title A print method for ritest objects
#' @name print.ritest
#' @description Printed display of ritest objects. Tries to mimic the display
#'   of the equivalent Stata routine.
#' @param x An ritest object.
#' @param verbose Logical. Should we display the original model summary too?
#'   Default is FALSE.
#' @param ... Currently ignored.
#' @inherit ritest examples
#' @export
print.ritest = function(x, verbose = FALSE, ...) {

  ri_mat = c(`T(obs)` = x$beta_parm, `c` = x$count, `n` = x$reps,
             `p=c/n` = x$pval, `SE(p)` = x$se)
  ri_mat = c(ri_mat, x$ci)
  ri_mat_nms = names(ri_mat)
  ri_mat = sprintf(ri_mat, fmt = '%.4g')
  names(ri_mat) = ri_mat_nms
  pval_string = if (attr(x$pval, 'side') %in% c('=', '==')) {
    "Note: c = #{|T| >= |T(obs)|}"
  } else if (attr(x$pval, 'side')=='<=') {
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
    cat(x$object_summ_string, "\n\n", sep = "")
    cat("******************",
        "* RITEST RESULTS *",
        "******************",
        sep = "\n")
  }


  cat(cli::style_bold("\n          Call: "), x$call, "\n", sep = "")
  cat(cli::style_bold("   Res. var(s): "), x$resampvar, "\n", sep = "")
  cat(cli::style_bold("            H0: "), x$h0, "\n", sep = "")
  if (!is.null(x$strata)) {
    cat(cli::style_bold(" Strata var(s): "), x$strata, "\n", sep = "")
    cat(cli::style_bold("        Strata: "), attr(x$strata, 'levels'), "\n", sep = "")
  }
  if (!is.null(x$cluster)) {
    cat(cli::style_bold("Cluster var(s): "), x$cluster, "\n", sep = "")
    cat(cli::style_bold("      Clusters: "), attr(x$cluster, 'levels'), "\n", sep = "")
  }
  cat(cli::style_bold("     Num. reps: "), x$reps, "\n", sep = "")
  cat(cli::rule(), "\n")
  print(ri_mat, quote = FALSE, print.gap = 2L)
  cat(cli::rule(), "\n")
  cat("Note: Confidence interval is with respect to p=c/n.", "\n")
  cat(pval_string, "\n")
  cat("\n")
  # invisible(x)
}
