#' @title S3 Class for PS Inference Test Results
#'
#' @description
#' The \code{ps_test} class is the unified output object returned by all
#' inferential functions in \pkg{PSinference}. It stores the test result,
#' the simulated null distribution, and all relevant metadata, and provides
#' \code{print}, \code{summary}, and \code{plot} methods for convenient
#' inspection and reporting.
#'
#' @section Slots:
#' \describe{
#'   \item{statistic}{Numeric. Observed value of the test statistic.}
#'   \item{p.value}{Numeric. Monte Carlo p-value.}
#'   \item{alpha}{Numeric. Significance level used.}
#'   \item{decision}{Character. \code{"Reject H0"} or
#'     \code{"Fail to Reject H0"}.}
#'   \item{null.dist}{Numeric vector. Simulated null distribution.}
#'   \item{test}{Character. One of \code{"gv"}, \code{"sphericity"},
#'     \code{"independence"}, \code{"regression"}.}
#'   \item{n}{Integer. Original sample size.}
#'   \item{M}{Integer. Number of synthetic releases.}
#'   \item{N}{Integer. Effective sample size \eqn{N = Mn}.}
#'   \item{p}{Integer. Number of variables.}
#'   \item{conf.int}{Numeric vector of length 2 or \code{NULL}.
#'     Confidence interval (generalised variance only).}
#'   \item{sigma2.hat}{Numeric or \code{NULL}. Plug-in estimator of
#'     \eqn{\sigma^2} (sphericity only).}
#'   \item{Delta.hat}{Matrix or \code{NULL}. Plug-in estimator of
#'     \eqn{\Delta} (regression only).}
#'   \item{iterations}{Integer. Number of Monte Carlo iterations used
#'     to calibrate the null distribution.}
#' }
#'
#' @name ps_test-class
NULL

## ------------------------------------------------------------------
## Internal constructor
## ------------------------------------------------------------------

#' @keywords internal
new_ps_test <- function(statistic,
                        p.value,
                        alpha,
                        decision,
                        null.dist,
                        test,
                        n,
                        M,
                        p,
                        conf.int   = NULL,
                        sigma2.hat = NULL,
                        Delta.hat  = NULL,
                        lbl1       = NULL,
                        lbl2       = NULL) {

  stopifnot(
    is.numeric(statistic) || is.na(statistic),
    is.numeric(p.value)   || is.na(p.value),
    is.numeric(alpha),
    is.character(decision),
    is.numeric(null.dist),
    is.character(test),
    is.numeric(n),
    is.numeric(M),
    is.numeric(p)
  )

  structure(
    list(
      statistic  = statistic,
      p.value    = p.value,
      alpha      = alpha,
      decision   = decision,
      null.dist  = null.dist,
      test       = test,
      n          = as.integer(n),
      M          = as.integer(M),
      N          = as.integer(M * n),
      p          = as.integer(p),
      conf.int   = conf.int,
      sigma2.hat = sigma2.hat,
      Delta.hat  = Delta.hat,
      lbl1       = lbl1,
      lbl2       = lbl2,
      iterations = length(null.dist)
    ),
    class = "ps_test"
  )
}

## ------------------------------------------------------------------
## print method
## ------------------------------------------------------------------

#' @title Print a \code{ps_test} Object
#' @description Prints a concise, human-readable summary of the test
#'   result stored in a \code{ps_test} object.
#' @param x An object of class \code{ps_test}.
#' @param ... Further arguments (currently ignored).
#' @return Invisibly returns \code{x}.
#' @exportS3Method print ps_test
#' @examples
#' data(ps_attitude)
#' V <- simSynthData(ps_attitude, M = 3)
#' res <- sphericity_test(V, M = 3)
#' print(res)
print.ps_test <- function(x, ...) {

  test_label <- .test_label(x$test)
  bar        <- strrep("-", 56)

  cat("\n")
  cat("PSInference:", test_label, "\n")
  cat(bar, "\n")
  cat(sprintf("  Original sample size  n = %d\n",   x$n))
  cat(sprintf("  Number of variables   p = %d\n",   x$p))
  cat(sprintf("  Number of releases    M = %d\n",   x$M))
  cat(sprintf("  Effective sample size N = Mn = %d\n", x$N))
  cat(bar, "\n")
  cat(sprintf("  Test statistic  : %.6g\n", x$statistic))
  cat(sprintf("  p-value         : %s\n",  .ps_fmt_pvalue(x$p.value, x$iterations)))
  cat(sprintf("  alpha           : %.2f\n",   x$alpha))
  cat(sprintf("  Decision        : %s\n",     x$decision))

  # Test-specific extras
  if (x$test == "gv") {
    cat(sprintf("  H0 : |Sigma| = |Sigma_0|  (two-sided test)\n"))
    if (!is.null(x$conf.int)) {
      cat(sprintf(
        "  %.0f%% CI for |Sigma|: (%.4e, %.4e)\n",
        100 * (1 - x$alpha), x$conf.int[1], x$conf.int[2]))
    }
  }
  if (x$test == "sphericity" && !is.null(x$sigma2.hat)) {
    cat(sprintf("  H0 : Sigma = sigma^2 * I_p\n"))
    cat(sprintf("  sigma^2_hat = %.4f  (under H0)\n", x$sigma2.hat))
  }
  if (x$test == "independence") {
    if (!is.null(x$lbl1) && !is.null(x$lbl2)) {
      cat(sprintf("  H0 : {%s}  independent of  {%s}\n",
                  x$lbl1, x$lbl2))
    } else {
      cat("  H0 : Sigma_12 = 0\n")
    }
  }
  if (x$test == "regression") {
    if (!is.null(x$lbl1) && !is.null(x$lbl2)) {
      cat(sprintf("  H0 : Delta = Delta_0  (regress {%s} on {%s})\n",
                  x$lbl1, x$lbl2))
    } else {
      cat("  H0 : Delta = Delta_0\n")
    }
    if (!is.null(x$Delta.hat)) {
      cat("  Delta_hat (plug-in slope):\n")
      print(round(x$Delta.hat, 4))
    }
  }
  cat(bar, "\n")
  cat(sprintf("  Monte Carlo iterations: %d\n", x$iterations))
  cat("\n")
  invisible(x)
}

## ------------------------------------------------------------------
## summary method
## ------------------------------------------------------------------

#' @title Summarise a \code{ps_test} Object
#' @description Prints a detailed summary including the null
#'   distribution quantiles and a comparison with the observed
#'   statistic.
#' @param object An object of class \code{ps_test}.
#' @param ... Further arguments (currently ignored).
#' @return Invisibly returns \code{object}.
#' @exportS3Method summary ps_test
#' @examples
#' data(ps_attitude)
#' V <- simSynthData(ps_attitude, M = 3)
#' res <- sphericity_test(V, M = 3)
#' summary(res)
summary.ps_test <- function(object, ...) {

  print(object)

  nd <- object$null.dist
  cat("Null distribution summary (Monte Carlo):\n")
  q <- quantile(nd, probs = c(0.01, 0.025, 0.05,
                              0.25, 0.50, 0.75,
                              0.95, 0.975, 0.99))
  print(round(q, 6))
  cat(sprintf("\nObserved statistic : %.6g\n", object$statistic))

  # Position of observed statistic in null distribution
  pct <- round(100 * mean(nd <= object$statistic), 1)
  cat(sprintf("Percentile in null : %.1f%%\n\n", pct))

  invisible(object)
}

## ------------------------------------------------------------------
## plot method
## ------------------------------------------------------------------

#' @title Plot a \code{ps_test} Object
#'
#' @description
#' Produces a density plot of the simulated null distribution with the
#' observed test statistic and critical value(s) marked, and the
#' rejection region shaded.
#'
#' The x-axis always shows both the null distribution and the observed
#' statistic. For the generalised variance and regression tests a log10
#' scale is used automatically because these statistics span many orders
#' of magnitude. Key information is placed below the plot as text so it
#' never overlaps the density curve. The function is multi-panel aware:
#' inside \code{par(mfrow = ...)} it uses compact in-plot annotations
#' and does not modify the outer margins.
#'
#' @param x An object of class \code{ps_test}.
#' @param main Optional title string. Auto-generated if \code{NULL}.
#' @param shade_col Colour for the rejection-region shading.
#' @param dist_col  Colour for the null-distribution density fill.
#' @param stat_col  Colour for the observed-statistic line.
#' @param crit_col  Colour for the critical-value line(s).
#' @param ...       Further arguments passed to \code{plot()}.
#' @return Invisibly returns \code{x}.
#' @exportS3Method plot ps_test
#' @examples
#' data(ps_attitude)
#' V <- simSynthData(ps_attitude, M = 3)
#' plot(sphericity_test(V, M = 3))
plot.ps_test <- function(x,
                         main      = NULL,
                         shade_col = grDevices::adjustcolor("tomato",    0.45),
                         dist_col  = grDevices::adjustcolor("steelblue", 0.22),
                         stat_col  = "firebrick",
                         crit_col  = "steelblue4",
                         ...) {

  nd  <- x$null.dist
  obs <- x$statistic

  ## ---- log10 scale for GV and regression --------------------------------
  use_log <- x$test %in% c("gv", "regression")

  if (use_log) {
    nd_pos   <- nd[nd > 0]
    nd_plot  <- log10(nd_pos)
    obs_plot <- if (!is.na(obs) && obs > 0) log10(obs) else NA_real_
    xlab_str <- expression(log[10](italic(T)^"*"))
  } else {
    nd_plot  <- nd
    obs_plot <- obs
    xlab_str <- expression(italic(T)^"*")
  }

  ## ---- tail type --------------------------------------------------------
  two_sided <- x$test == "gv"
  left_tail <- x$test %in% c("sphericity", "independence")

  ## ---- critical values on plot scale ------------------------------------
  .cv <- function(p) {
    v <- as.numeric(stats::quantile(nd, p))
    if (use_log) log10(max(v, .Machine$double.eps)) else v
  }
  .cv_raw <- function(p) as.numeric(stats::quantile(nd, p))

  if (two_sided) {
    cl_plot <- .cv(x$alpha / 2);        ch_plot <- .cv(1 - x$alpha / 2)
    cl_raw  <- .cv_raw(x$alpha / 2);    ch_raw  <- .cv_raw(1 - x$alpha / 2)
  } else if (left_tail) {
    cl_plot <- .cv(x$alpha);  ch_plot <- NA_real_
    cl_raw  <- .cv_raw(x$alpha); ch_raw <- NA_real_
  } else {
    cl_plot <- NA_real_; ch_plot <- .cv(1 - x$alpha)
    cl_raw  <- NA_real_; ch_raw  <- .cv_raw(1 - x$alpha)
  }

  ## ---- density on plot scale --------------------------------------------
  d <- stats::density(nd_plot, bw = "SJ")

  ## ---- x-axis: always includes null dist AND observed stat --------------
  nd_range <- range(d$x)
  pad      <- diff(nd_range) * 0.08
  if (!is.na(obs_plot) && is.finite(obs_plot)) {
    x_lo <- min(nd_range[1L] - pad, obs_plot - pad * 0.3)
    x_hi <- max(nd_range[2L] + pad, obs_plot + pad * 0.3)
  } else {
    x_lo <- nd_range[1L] - pad
    x_hi <- nd_range[2L] + pad
  }
  ylim <- c(0, max(d$y) * 1.08)

  ## ---- title and annotation strings ------------------------------------
  main_ <- if (!is.null(main)) main else
    sprintf("%s  (M = %d,  N = %d)", .test_label(x$test), x$M, x$N)

  .fmt <- function(v_plot) {
    if (is.na(v_plot) || !is.finite(v_plot)) return("NA")
    if (use_log) .ps_fmt_sci(10^v_plot) else sprintf("%.4g", v_plot)
  }

  ann_obs  <- sprintf("Observed T* = %s", .fmt(obs_plot))
  ann_pval <- sprintf("p-value = %s", .ps_fmt_pvalue(x$p.value, x$iterations))
  ann_dec  <- x$decision
  if (two_sided) {
    ann_crit <- sprintf("Reject if T* < %s  or  T* > %s",
                        .fmt(cl_plot), .fmt(ch_plot))
  } else if (left_tail) {
    ann_crit <- sprintf("Reject if T* < %s", .fmt(cl_plot))
  } else {
    ann_crit <- sprintf("Reject if T* > %s", .fmt(ch_plot))
  }

  ## ---- multi-panel awareness -------------------------------------------
  multi_panel <- prod(graphics::par("mfrow")) > 1L

  if (multi_panel) {
    op <- graphics::par(
      mar = c(4.8, 4.0, 2.5, 0.8),
      mgp = c(2.4, 0.6, 0)
    )
  } else {
    op <- graphics::par(
      oma = c(6.0, 0, 0, 0),   # taller outer margin
      mar = c(4.0, 4.5, 2.8, 1.2),  # extra bottom mar creates gap
      mgp = c(2.8, 0.7, 0)     # x-label sits at line 2.8
    )
  }
  on.exit(graphics::par(op), add = TRUE)

  ## ---- base plot -------------------------------------------------------
  plot(d,
       main = main_,
       xlab = xlab_str,
       ylab = "Density",
       xlim = c(x_lo, x_hi),
       ylim = ylim,
       lwd  = 2.5,
       col  = "steelblue",
       axes = FALSE,
       zero.line = FALSE,
       ...)
  graphics::axis(1L, las = 1L)
  graphics::axis(2L, las = 1L)
  graphics::box()

  ## null distribution fill
  graphics::polygon(c(d$x, rev(d$x)), c(d$y, rep(0L, length(d$y))),
                    col = dist_col, border = NA)
  graphics::lines(d, col = "steelblue", lwd = 2.5)

  ## ---- shade rejection region(s) ----------------------------------------
  .shade <- function(thresh, left) {
    if (is.na(thresh) || !is.finite(thresh)) return(invisible(NULL))
    idx <- if (left) d$x <= thresh else d$x >= thresh
    if (sum(idx) < 2L) return(invisible(NULL))
    graphics::polygon(c(d$x[idx], rev(d$x[idx])),
                      c(d$y[idx], rep(0, sum(idx))),
                      col = shade_col, border = NA)
  }
  if (two_sided)      { .shade(cl_plot, TRUE);  .shade(ch_plot, FALSE)
  } else if (left_tail) { .shade(cl_plot, TRUE)
  } else                { .shade(ch_plot, FALSE) }

  graphics::lines(d, col = "steelblue", lwd = 2.5)  # redraw on top

  ## ---- observed statistic line -----------------------------------------
  if (!is.na(obs_plot) && is.finite(obs_plot))
    graphics::abline(v = obs_plot, col = stat_col, lwd = 2.2, lty = 2L)

  ## ---- critical value line(s) ------------------------------------------
  if (!is.na(cl_plot) && is.finite(cl_plot))
    graphics::abline(v = cl_plot, col = crit_col, lwd = 1.6, lty = 3L)
  if (!is.na(ch_plot) && is.finite(ch_plot))
    graphics::abline(v = ch_plot, col = crit_col, lwd = 1.6, lty = 3L)

  ## ---- annotation below the plot ---------------------------------------
  dec_col <- if (grepl("^Reject", x$decision)) "firebrick" else "darkgreen"

  if (multi_panel) {
    ## Two compact lines inside the plot's bottom margin
    sub1 <- sprintf("%s   |   %s", ann_obs, ann_pval)
    sub2 <- sprintf("%s   |   %s", ann_dec, ann_crit)
    graphics::mtext(sub1, side = 1L, line = 3.0, cex = 0.65,
                    col = "black", font = 1L)
    graphics::mtext(sub2, side = 1L, line = 4.0, cex = 0.65,
                    col = dec_col, font = 1L)
  } else {
    ## Four lines in the outer margin (single-panel)
    ann  <- list(
      list(ann_obs,  stat_col, 2L),
      list(ann_pval, "black",  1L),
      list(ann_dec,  dec_col,  2L),
      list(ann_crit, crit_col, 1L)
    )
    graphics::par(xpd = NA)
    for (k in seq_along(ann))
      graphics::mtext(ann[[k]][[1L]],
                      side = 1L, outer = TRUE,
                      line = (k - 1L) * 1.20 + 0.8,  # start at 0.8, step 1.20
                      cex = 0.90, adj = 0.5,
                      font = ann[[k]][[3L]],
                      col  = ann[[k]][[2L]])
  }

  invisible(x)
}


## ------------------------------------------------------------------
## is / as helpers
## ------------------------------------------------------------------

#' @title Test if an Object is of Class \code{ps_test}
#' @param x Any R object.
#' @return Logical.
#' @export
is.ps_test <- function(x) inherits(x, "ps_test")

## ------------------------------------------------------------------
## Internal helpers
## ------------------------------------------------------------------

#' @keywords internal
.test_label <- function(test) {
  switch(test,
         gv           = "Generalised Variance Test",
         sphericity   = "Sphericity Test",
         independence = "Independence Test",
         regression   = "Regression Test",
         test)
}

#' @keywords internal
.ps_fmt_sci <- function(v) {
  if (abs(v) >= 0.001 && abs(v) < 1e5) return(sprintf("%.4g", v))
  sprintf("%.3e", v)
}

#' @keywords internal
#' Format a Monte Carlo p-value for display.
#' Never returns "0.0000"; uses "< 1/B" notation when the MC count is zero.
#' @param pval  Numeric p-value in [0, 1].
#' @param B     Number of Monte Carlo iterations used (determines resolution).
.ps_fmt_pvalue <- function(pval, B = NULL) {
  if (is.null(B)) B <- 10000L          # conservative default resolution
  min_detectable <- 1.0 / B           # smallest non-zero MC p-value

  if (is.na(pval) || !is.finite(pval)) return("NA")

  if (pval == 0 || pval < min_detectable) {
    # Express as "< 0.0001" (4 sig fig of min_detectable)
    thresh <- signif(min_detectable, 1)
    # Format threshold without trailing zeros
    if (thresh >= 0.001) {
      return(sprintf("< %.4f", thresh))
    } else {
      return(sprintf("< %.2e", thresh))
    }
  }

  if (pval >= 0.001) return(sprintf("%.4f", pval))
  sprintf("%.2e", pval)
}
