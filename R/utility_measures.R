#' @title Utility Measures for Plug-in Sampling Synthetic Data
#'
#' @description
#' Computes a battery of utility measures comparing the statistical
#' properties of PS synthetic data to the original confidential data.
#' Five complementary measures are reported: Frobenius distance between
#' covariance matrices, standardised mean differences (SMD),
#' variance ratios, propensity score MSE (pMSE), and confidence interval
#' overlap. For \eqn{M > 1} releases, per-release statistics are
#' computed and averaged.
#'
#' @param X A numeric matrix (\eqn{n \times p}), the original
#'   confidential data.
#' @param V A numeric matrix (\eqn{Mn \times p}), the stacked PS
#'   synthetic data returned by \code{\link{simSynthData}}.
#' @param M A positive integer. Number of synthetic releases
#'   (default \code{1L}).
#' @param alpha Significance level for confidence interval overlap
#'   (default \code{0.05}).
#' @param verbose Logical. If \code{TRUE} (default), prints a formatted
#'   summary.
#'
#' @return A list of class \code{ps_utility} (invisibly) with components:
#' \describe{
#'   \item{frobenius}{Numeric. Frobenius distance
#'     \eqn{\|\hat{\bm{\Sigma}}_{\mathrm{orig}} -
#'     \bar{\hat{\bm{\Sigma}}}_{\mathrm{synth}}\|_F}.}
#'   \item{smd}{Named numeric vector. Per-variable standardised mean
#'     difference \eqn{|\bar{v}_j - \bar{x}_j| / s_{x_j}}.}
#'   \item{smd_mean}{Numeric. Mean SMD across all \eqn{p} variables.}
#'   \item{var_ratio}{Named numeric vector. Per-variable variance ratio
#'     \eqn{\bar{s}^2_{v_j} / s^2_{x_j}}.}
#'   \item{pmse}{Numeric. Propensity score MSE
#'     \eqn{n^{-1}\sum_i(\hat{p}_i - c)^2} where \eqn{c = n/(n + Mn)}.}
#'   \item{pmse_null}{Numeric. Expected pMSE under a correctly
#'     specified model: \eqn{c(1-c)(p+1)/(n + Mn)}.}
#'   \item{pmse_ratio}{Numeric. \code{pmse / pmse_null}; values near 1
#'     indicate good utility, large values indicate distinguishability.}
#'   \item{ci_overlap}{Named numeric vector. Per-variable confidence
#'     interval overlap (Karr et al., 2006).}
#'   \item{ci_overlap_mean}{Numeric. Mean CI overlap.}
#'   \item{M}{Integer. Number of releases.}
#'   \item{n}{Integer. Original sample size.}
#'   \item{p}{Integer. Number of variables.}
#' }
#'
#' @details
#' \subsection{Frobenius distance}{
#' Measures preservation of the multivariate covariance structure:
#' \deqn{
#'   d_F = \|\hat{\bm{\Sigma}}_X -
#'   \bar{\hat{\bm{\Sigma}}}_V\|_F
#'   = \sqrt{\sum_{j,k}(\hat\sigma_{jk} -
#'   \bar{\hat\sigma}^{\star}_{jk})^2},
#' }
#' where \eqn{\bar{\hat{\bm{\Sigma}}}_V = M^{-1}\sum_{m=1}^M
#' \hat{\bm{\Sigma}}_{V_m}} is the average per-release sample
#' covariance. Smaller values indicate better covariance utility.
#' }
#'
#' \subsection{Standardised Mean Difference}{
#' For each variable \eqn{j}:
#' \deqn{
#'   \mathrm{SMD}_j = \frac{|\bar{v}_j - \bar{x}_j|}{s_{x_j}},
#' }
#' where \eqn{\bar{v}_j} is the mean of the pooled synthetic data and
#' \eqn{s_{x_j}} is the original standard deviation. Values below 0.10
#' indicate negligible difference; 0.10--0.20 is small; above 0.20
#' warrants investigation.
#' }
#'
#' \subsection{Variance ratio}{
#' Per-variable ratio of average synthetic variance to original variance:
#' \deqn{
#'   \mathrm{VR}_j = \bar{s}^2_{v_j} / s^2_{x_j}.
#' }
#' Values near 1 indicate good variance preservation.
#' }
#'
#' \subsection{Propensity Score MSE (pMSE)}{
#' Combines original (\eqn{n} rows, label 0) and synthetic (\eqn{Mn}
#' rows, label 1) data, fits a logistic regression of the label on all
#' \eqn{p} variables, and computes:
#' \deqn{
#'   \mathrm{pMSE} = \frac{1}{n + Mn}\sum_{i=1}^{n+Mn}
#'   (\hat{p}_i - c)^2,
#'   \quad c = \frac{Mn}{n + Mn}.
#' }
#' The null value (expected under a correctly specified model) is
#' \eqn{c(1-c)(p+1)/(n+Mn)}. The ratio
#' \eqn{\mathrm{pMSE}/\mathrm{pMSE}_{\mathrm{null}}} should be near 1
#' for high-utility synthetic data; large values indicate the synthetic
#' data is easily distinguished from the original.
#' }
#'
#' \subsection{Confidence interval overlap}{
#' For each variable \eqn{j}, compute the \eqn{(1-\alpha)} CI from
#' the original and from the synthetic data, and measure their overlap
#' (Karr et al., 2006):
#' \deqn{
#'   J_j = \frac{\max(0,\; u^{V}_j - l^{X}_j,\;
#'   u^{X}_j - l^{V}_j)}{\max(u^X_j - l^X_j,\; u^V_j - l^V_j)},
#' }
#' where \eqn{[l^X_j, u^X_j]} and \eqn{[l^V_j, u^V_j]} are the CIs
#' from original and synthetic data. Values near 1 indicate high
#' inferential utility.
#' }
#'
#' @references
#' Karr, A. F., Kohnen, C. N., Oganian, A., Reiter, J. P. and
#' Sanil, A. P. (2006). A framework for evaluating the utility of
#' data altered to protect confidentiality.
#' \emph{The American Statistician}, 60, 224--232.
#'
#' Snoke, J., Raab, G. M., Nowok, B., Dibben, C. and Slavkovic, A.
#' (2018). General and specific utility measures for synthetic data.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' 181, 663--688.
#'
#' Woo, M.-J., Reiter, J. P., Oganian, A. and Karr, A. F. (2009).
#' Global measures of data utility for microdata masked for disclosure
#' limitation. \emph{Journal of Privacy and Confidentiality}, 1,
#' 111--124.
#'
#' @seealso \code{\link{simSynthData}}, \code{\link{ps_test}}
#'
#' @export
#'
#' @examples
#' data(brittany_soil_ps)
#' set.seed(1)
#' V3 <- simSynthData(brittany_soil_ps, M = 3)
#' utility_measures(brittany_soil_ps, V3, M = 3)
utility_measures <- function(X, V, M = 1L,
                             alpha   = 0.05,
                             verbose = TRUE) {

  # ------------------------------------------------------------------
  # Input validation
  # ------------------------------------------------------------------
  X <- .validate_X(X)
  V <- .validate_X(V)
  M <- as.integer(M)

  n <- nrow(X); p <- ncol(X)
  N <- nrow(V)  # = M * n

  if (N != M * n)
    stop(sprintf(
      "nrow(V) = %d is not equal to M * nrow(X) = %d * %d = %d.",
      N, M, n, M * n), call. = FALSE)

  if (ncol(V) != p)
    stop(sprintf(
      "ncol(V) = %d must equal ncol(X) = %d.", ncol(V), p),
      call. = FALSE)

  vnames <- if (!is.null(colnames(X))) colnames(X) else
    paste0("V", seq_len(p))

  # ------------------------------------------------------------------
  # Per-release statistics
  # ------------------------------------------------------------------
  # Split V into M releases of n rows each
  rel_idx <- lapply(seq_len(M), function(m) (m - 1L) * n + seq_len(n))

  # Per-release sample covariance matrices
  S_list  <- lapply(rel_idx, function(idx) stats::cov(V[idx, ,drop=FALSE]))
  # Per-release column means
  mu_list <- lapply(rel_idx, function(idx) colMeans(V[idx, ,drop=FALSE]))

  # Averages across releases
  S_avg   <- Reduce("+", S_list) / M
  mu_avg  <- Reduce("+", mu_list) / M

  # Original statistics
  S_orig  <- stats::cov(X)
  mu_orig <- colMeans(X)
  sd_orig <- apply(X, 2L, stats::sd)

  # ------------------------------------------------------------------
  # 1. Frobenius distance
  # ------------------------------------------------------------------
  frob <- norm(S_orig - S_avg, "F")

  # ------------------------------------------------------------------
  # 2. Standardised Mean Difference (SMD)
  # ------------------------------------------------------------------
  smd      <- abs(mu_avg - mu_orig) / sd_orig
  names(smd) <- vnames
  smd_mean <- mean(smd)

  # ------------------------------------------------------------------
  # 3. Variance ratio
  # ------------------------------------------------------------------
  var_syn   <- diag(S_avg)
  var_orig  <- diag(S_orig)
  var_ratio <- var_syn / var_orig
  names(var_ratio) <- vnames

  # ------------------------------------------------------------------
  # 4. pMSE (propensity score MSE)
  # ------------------------------------------------------------------
  c_prop <- M * n / (n + M * n)   # proportion of synthetic rows
  df_pool <- data.frame(
    rbind(X, V),
    .label = c(rep(0L, n), rep(1L, M * n))
  )
  fit_pmse <- tryCatch(
    stats::glm(.label ~ ., data = df_pool,
               family = stats::binomial(link = "logit")),
    error = function(e) NULL,
    warning = function(w) suppressWarnings(
      stats::glm(.label ~ ., data = df_pool,
                 family = stats::binomial(link = "logit")))
  )

  if (!is.null(fit_pmse)) {
    phat      <- stats::fitted(fit_pmse)
    pmse      <- mean((phat - c_prop)^2)
    pmse_null <- c_prop * (1 - c_prop) * (p + 1L) / (n + M * n)
    pmse_ratio <- pmse / pmse_null
  } else {
    pmse <- pmse_null <- pmse_ratio <- NA_real_
  }

  # ------------------------------------------------------------------
  # 5. Confidence interval overlap (Karr et al. 2006)
  # ------------------------------------------------------------------
  t_crit <- stats::qt(1 - alpha / 2, df = n - 1L)
  t_crit_s <- stats::qt(1 - alpha / 2, df = M * n - 1L)

  ci_overlap <- vapply(seq_len(p), function(j) {
    mx <- mu_orig[j];       sx <- sd_orig[j]
    mv <- colMeans(V)[j];   sv <- sqrt(diag(stats::cov(V))[j])
    lx <- mx - t_crit   * sx / sqrt(n)
    ux <- mx + t_crit   * sx / sqrt(n)
    lv <- mv - t_crit_s * sv / sqrt(M * n)
    uv <- mv + t_crit_s * sv / sqrt(M * n)
    overlap <- max(0, min(ux, uv) - max(lx, lv))
    width   <- max(ux - lx, uv - lv)
    if (width < .Machine$double.eps) 1 else overlap / width
  }, numeric(1L))
  names(ci_overlap) <- vnames
  ci_overlap_mean <- mean(ci_overlap)

  # ------------------------------------------------------------------
  # Print
  # ------------------------------------------------------------------
  if (verbose) {
    bar <- strrep("-", 60)
    cat("\nPS Synthetic Data Utility Assessment\n")
    cat(bar, "\n")
    cat(sprintf("  n = %d | p = %d | M = %d releases (N = %d)\n",
                n, p, M, M * n))
    cat(bar, "\n\n")

    cat("1. Frobenius distance  ||cov(X) - avg cov(V)||_F\n")
    cat(sprintf("   %.4f\n\n", frob))

    cat("2. Standardised Mean Differences  |mean(V_j) - mean(X_j)| / sd(X_j)\n")
    smd_df <- data.frame(
      variable = vnames,
      SMD      = round(smd, 4),
      flag     = ifelse(smd > 0.20, "**", ifelse(smd > 0.10, "*", ""))
    )
    print(smd_df, row.names = FALSE)
    cat(sprintf("   Mean SMD = %.4f", smd_mean))
    cat("  (* > 0.10  ** > 0.20)\n\n")

    cat("3. Variance Ratios  var(V_j) / var(X_j)\n")
    vr_df <- data.frame(
      variable  = vnames,
      var_ratio = round(var_ratio, 4),
      flag      = ifelse(abs(var_ratio - 1) > 0.20, "**",
                         ifelse(abs(var_ratio - 1) > 0.10, "*", ""))
    )
    print(vr_df, row.names = FALSE)
    cat("  (* |VR-1| > 0.10  ** |VR-1| > 0.20)\n\n")

    cat("4. Propensity Score MSE (pMSE)\n")
    if (!is.na(pmse)) {
      cat(sprintf("   pMSE       = %.4f\n", pmse))
      cat(sprintf("   pMSE_null  = %.4f  (expected under good utility)\n",
                  pmse_null))
      cat(sprintf("   pMSE ratio = %.4f  (near 1 = good)\n\n",
                  pmse_ratio))
    } else {
      cat("   Could not fit propensity model (possible separation).\n\n")
    }

    cat("5. Confidence Interval Overlap  (alpha =", alpha, ")\n")
    ci_df <- data.frame(
      variable   = vnames,
      CI_overlap = round(ci_overlap, 4),
      flag       = ifelse(ci_overlap < 0.70, "**",
                          ifelse(ci_overlap < 0.90, "*", ""))
    )
    print(ci_df, row.names = FALSE)
    cat(sprintf("   Mean CI overlap = %.4f", ci_overlap_mean))
    cat("  (* < 0.90  ** < 0.70)\n\n")
    cat(bar, "\n")
  }

  # ------------------------------------------------------------------
  # Return
  # ------------------------------------------------------------------
  result <- list(
    frobenius       = round(frob, 4),
    smd             = round(smd, 4),
    smd_mean        = round(smd_mean, 4),
    var_ratio       = round(var_ratio, 4),
    pmse            = round(pmse, 4),
    pmse_null       = round(pmse_null, 4),
    pmse_ratio      = round(pmse_ratio, 4),
    ci_overlap      = round(ci_overlap, 4),
    ci_overlap_mean = round(ci_overlap_mean, 4),
    M               = M,
    n               = n,
    p               = p,
    alpha           = alpha
  )
  class(result) <- "ps_utility"
  invisible(result)
}

#' @title Print Method for \code{ps_utility} Objects
#' @param x An object of class \code{ps_utility}.
#' @param ... Further arguments (ignored).
#' @exportS3Method print ps_utility
print.ps_utility <- function(x, ...) {
  cat(sprintf(
    "PS Utility: Frob=%.4f | Mean SMD=%.4f | Mean CI overlap=%.4f | pMSE ratio=%.4f\n",
    x$frobenius, x$smd_mean, x$ci_overlap_mean, x$pmse_ratio))
  invisible(x)
}
