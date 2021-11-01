#' @title A plot method for ritest objects
#' @name plot.ritest
#' @description Nice plots of your ritest objects.
#' @param x An ritest object.
#' @param type Character. What type of plot do you want?
#' @param highlight Character. How do you want to highlight the H0 rejection
#'   regions in the distribution tails?
#' @param show_parm Logical. Should we highlight the parametric H0 rejection
#'   regions too?
#' @param breaks Character. Histogram plot only. What type of breaks do you
#'   want? The default method creates more breaks than the standard R behaviour.
#'   You can revert to the latter by selecting NULL.
#' @param family Character. The font family. Defaults to 'HersheySans' instead
#'   of R's normal Arial plotting font.
#' @param ... Other plot arguments. Currently ignored.
#' @importFrom grDevices rgb
#' @importFrom graphics abline hist title
#' @importFrom stats density
#' @inherit ritest examples
#' @export
plot.ritest = function(x, type = c('density', 'hist'),
                       highlight = c('lines', 'fill', 'both', 'none'),
                       show_parm = FALSE,  breaks = 'auto',
                       family = NULL, ...) {
  type = match.arg(type)
  highlight = match.arg(highlight)
  if (is.null(family)) family = 'HersheySans'

  h0_value = as.numeric(gsub('.*=|<|>', '', x$h0))

  beta_parm = x$beta_parm - h0_value
  side = attr(x$pval, 'side')

  ci_sides =
    if (side=='<=') {
      1
    } else if (side=='>=') {
      2
    } else {
      1:2
    }

  hi_lines = beta_parm * (c(-1, 1)[ci_sides])
  if (beta_parm < 0) hi_lines = -hi_lines

  xlim =
    if (show_parm) {
      # range(x$betas, x$beta_parm, x$ci_parm - abs(x$beta_parm), finite = TRUE)
      range(x$betas, beta_parm, hi_lines, finite = TRUE)
    } else {
      range(x$betas, x$beta_parm, finite = TRUE)
    }

  if (type=='density') {
    dens = density(x$betas)
    plot(dens,
         main = '', xlab = '',
         xlim = xlim,
         family = family)
    if (highlight %in% c('both', 'fill')) {
      fcol = rgb(1,0,0,0.5)
      x1 = x2 = x3 = x4 = NA
      if (1 %in% ci_sides) {
        x1 = if (beta_parm < 0 && length(ci_sides)==1) length(dens$x) else 1
        x2 = tail(which(dens$x <= -abs(beta_parm)), 1)
        with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=fcol, border = FALSE))
      }
      if (2 %in% ci_sides) {
        x3 = head(which(dens$x >= abs(beta_parm)), 1)
        x4 = length(dens$x)
        with(dens, polygon(x=c(x[c(x3,x3:x4,x4)]), y= c(0, y[x3:x4], 0), col=fcol, border = FALSE))
      }
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
      hist_cond =
        if (side=='<=' && beta_parm >= 0) {
          hist_df$breaks <= hi_lines#-abs(beta_parm)
        } else if (side=='<=' && beta_parm < 0) {
          hist_df$breaks >= hi_lines
        } else if (side=='>=') {
          hist_df$breaks >= hi_lines#abs(beta_parm)
        } else {
          abs(hist_df$breaks) >= abs(beta_parm)
        }
      hist_col = ifelse(hist_cond, rgb(1,0,0,0.5), rgb(0.2,0.2,0.2,0.2))
    } else {
      hist_col = rgb(0.2,0.2,0.2,0.2)
    }
    plot(hist_df,
         col = hist_col,
         border = FALSE,
         main = NULL, xlab = NULL,
         xlim = xlim,
         family = family)
  }
  title(main = paste('Randomization Inference:', x$h0),
        sub = paste0('Simulated p-val: ',  sprintf('%.3f', x$pval),
                     '. Parametric p-val: ', sprintf('%.3f', x$pval_parm),'.'),
        cex.sub = 0.75,
        xlab = 'Simulated values',
        family = family)
  if (highlight %in% c('both', 'lines')) {
    abline(v = hi_lines, col = "red", lty = 1)
  }
  if (show_parm) {
    abline(v = x$ci_parm - beta_parm, lty = 2, col = 'blue')
  }
}
