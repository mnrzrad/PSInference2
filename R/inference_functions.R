#' @title Generalised Variance Test and Confidence Interval
#'
#' @description
#' Tests \eqn{H_0 : |\Sigma| = |\Sigma_0|} and computes an exact
#' \eqn{(1-\alpha)}-level confidence interval for the generalised
#' variance \eqn{|\Sigma|}, based on \eqn{M} released plug-in
#' sampling synthetic datasets stacked into \code{V}.
#' Setting \code{M = 1} recovers the single-release procedure of
#' Klein et al. (2021).
#'
#' @param V Stacked synthetic dataset (\eqn{Mn \times p} matrix),
#'   as returned by \code{\link{simSynthData}}.
#' @param M Positive integer. Number of synthetic releases
#'   (default \code{1L}).
#' @param Sigma A \eqn{p \times p} positive-definite matrix specifying
#'   the null value \eqn{\Sigma_0}. Typically \code{cov(X)}.
#' @param alpha Significance level (default \code{0.05}).
#' @param iterations Monte Carlo sample size (default \code{10000L}).
#'
#' @return An object of class \code{\link{ps_test}} with component
#'   \code{conf.int} giving the exact \eqn{(1-\alpha)} CI for
#'   \eqn{|\Sigma|}. All S3 methods (\code{print}, \code{summary},
#'   \code{plot}) are available.
#'
#' @seealso \code{\link{GVdist}}, \code{\link{ps_test}},
#'   \code{\link{simSynthData}}
#'
#' @references
#' Klein, M., Moura, R., and Sinha, B. (2021). Multivariate normal
#' inference based on singly imputed synthetic data under plug-in
#' sampling. \emph{Sankhya B}, 83, 273--287.
#'
#' @export
#'
#' @examples
#' data(ps_attitude)
#' set.seed(1)
#' V1 <- simSynthData(ps_attitude)
#' res <- gv_test(V1, M = 1, Sigma = cov(ps_attitude))
#' print(res)
#' plot(res)
#'
#' V5 <- simSynthData(ps_attitude, M = 5)
#' gv_test(V5, M = 5, Sigma = cov(ps_attitude))
gv_test <- function(V, M = 1L, Sigma,
                    alpha = 0.05, iterations = 10000L,
                    null_dist = NULL) {

  V   <- .validate_X(V)
  M   <- .validate_M(M)
  N   <- nrow(V)
  n   <- N / M
  p   <- ncol(V)

  .check_N(N, p)

  if (missing(Sigma))
    stop(paste0(
      "'Sigma' must be supplied: the null covariance matrix Sigma_0 ",
      "for H0: |Sigma| = |Sigma_0|. Typically use cov(X)."),
      call. = FALSE)

  Sigma <- as.matrix(Sigma)
  if (!isTRUE(all.equal(dim(Sigma), c(p, p))))
    stop(sprintf("'Sigma' must be a %d x %d matrix.", p, p),
         call. = FALSE)
  .check_pd(Sigma, "'Sigma'")

  S_star <- .compute_S_star(V)

  # GVdist now takes n (original sample size) and M separately.
  # The null distribution is prod(chi2_{Mn-j}) * prod(chi2_{n-j}).
  # The test statistic uses (n-1)^p, matching the B_j ~ chi2_{n-j} factor.
  null_d <- if (!is.null(null_dist)) null_dist else
    GVdist(nsample = n, pvariates = p, M = M, iterations = iterations)

  q_lo <- as.numeric(quantile(null_d, probs = alpha / 2))
  q_hi <- as.numeric(quantile(null_d, probs = 1 - alpha / 2))

  # Scale factor: (n-1)^p matches the B_j degrees of freedom in GVdist
  scale <- (n - 1L)^p * det(S_star)
  T_obs <- scale / det(Sigma)
  pval  <- min(2 * min(mean(null_d <= T_obs),
                       mean(null_d >= T_obs)), 1)
  dec   <- if (T_obs < q_lo || T_obs > q_hi)
    "Reject H0" else "Fail to Reject H0"

  ci <- c(lower = scale / q_hi, upper = scale / q_lo)

  new_ps_test(
    statistic = T_obs,
    p.value   = pval,
    alpha     = alpha,
    decision  = dec,
    null.dist = null_d,
    test      = "gv",
    n         = n,
    M         = M,
    p         = p,
    conf.int  = ci
  )
}

#' @title Generalised Variance (alias for \code{gv_test})
#' @description Backward-compatible alias for \code{\link{gv_test}}.
#' @inheritParams gv_test
#' @export
gv_ci <- function(V, M = 1L, Sigma,
                  alpha = 0.05, iterations = 10000L,
                  null_dist = NULL) {
  gv_test(V = V, M = M, Sigma = Sigma,
          alpha = alpha, iterations = iterations, null_dist = null_dist)
}


#' @title Sphericity Test
#'
#' @description
#' Tests \eqn{H_0 : \Sigma = \sigma^2 I_p} (all variables uncorrelated
#' with equal variance) based on \eqn{M} released plug-in sampling
#' synthetic datasets stacked into \code{V}. The test is left-tailed.
#' Setting \code{M = 1} recovers Klein et al. (2021).
#'
#' @param V Stacked synthetic dataset (\eqn{Mn \times p} matrix).
#' @param M Positive integer. Number of synthetic releases
#'   (default \code{1L}).
#' @param alpha Significance level (default \code{0.05}).
#' @param iterations Monte Carlo sample size (default \code{10000L}).
#'
#' @return An object of class \code{\link{ps_test}}. Component
#'   \code{sigma2.hat} gives the plug-in estimator
#'   \eqn{\hat\sigma^2 = \mathrm{tr}(S^\star)/(p(N-1))} under \eqn{H_0}.
#'
#' @seealso \code{\link{Sphdist}}, \code{\link{ps_test}}
#'
#' @references
#' Klein, M., Moura, R., and Sinha, B. (2021). Multivariate normal
#' inference based on singly imputed synthetic data under plug-in
#' sampling. \emph{Sankhya B}, 83, 273--287.
#'
#' @export
#'
#' @examples
#' data(ps_attitude)
#' set.seed(1)
#' V5 <- simSynthData(ps_attitude, M = 5)
#' res <- sphericity_test(V5, M = 5)
#' print(res)
#' plot(res)
sphericity_test <- function(V, M = 1L, alpha = 0.05,
                            iterations = 10000L,
                            null_dist = NULL) {

  V  <- .validate_X(V)
  M  <- .validate_M(M)
  N  <- nrow(V)
  n  <- N / M
  p  <- ncol(V)

  .check_N(N, p)

  S_star <- .compute_S_star(V)
  T_obs  <- det(S_star)^(1.0 / p) / (sum(diag(S_star)) / p)

  null_d    <- if (!is.null(null_dist)) null_dist else
    Sphdist(nsample = n, pvariates = p, M = M, iterations = iterations)
  crit      <- as.numeric(quantile(null_d, probs = alpha))
  pval      <- mean(null_d <= T_obs)
  dec       <- if (T_obs < crit) "Reject H0" else "Fail to Reject H0"
  sigma2hat <- sum(diag(S_star)) / (p * (N - 1L))

  new_ps_test(
    statistic  = T_obs,
    p.value    = pval,
    alpha      = alpha,
    decision   = dec,
    null.dist  = null_d,
    test       = "sphericity",
    n          = n,
    M          = M,
    p          = p,
    sigma2.hat = sigma2hat
  )
}


#' @title Independence Test
#'
#' @description
#' Tests \eqn{H_0 : \Sigma_{12} = \mathbf{0}} (independence between two
#' subsets of variables) based on \eqn{M} released plug-in sampling
#' synthetic datasets stacked into \code{V}. Setting \code{M = 1}
#' recovers Klein et al. (2021).
#'
#' The two variable blocks can be specified in two ways (exactly one
#' must be used):
#' \describe{
#'   \item{\code{part} (integer scalar)}{First \code{part} columns
#'     form Block 1; remaining columns form Block 2. Original interface;
#'     backward-compatible.}
#'   \item{\code{group_a} and \code{group_b}}{Integer indices \emph{or}
#'     column names identifying each block. Together they must cover all
#'     columns of \code{V} exactly. \code{V} must have \code{colnames}
#'     when names are used.}
#' }
#'
#' @param V Stacked synthetic dataset (\eqn{Mn \times p} matrix).
#' @param M Positive integer. Number of synthetic releases
#'   (default \code{1L}).
#' @param part Integer scalar. Size of Block 1 (first \code{part}
#'   columns). Ignored when \code{group_a}/\code{group_b} are supplied.
#' @param group_a Integer indices or column names identifying Block 1.
#' @param group_b Integer indices or column names identifying Block 2.
#'   Together with \code{group_a} must cover all columns.
#' @param alpha Significance level (default \code{0.05}).
#' @param iterations Monte Carlo sample size (default \code{10000L}).
#'
#' @return An object of class \code{\link{ps_test}}. The null hypothesis
#'   string in \code{$null.value} names the two blocks explicitly.
#'
#' @seealso \code{\link{Inddist}}, \code{\link{ps_test}}
#'
#' @references
#' Klein, M., Moura, R., and Sinha, B. (2021). Multivariate normal
#' inference based on singly imputed synthetic data under plug-in
#' sampling. \emph{Sankhya B}, 83, 273--287.
#'
#' @export
#'
#' @examples
#' data(ps_attitude)
#' set.seed(1)
#' V5 <- simSynthData(ps_attitude, M = 5)
#'
#' ## Integer interface
#' independence_test(V5, M = 5, part = 2L)
#'
#' ## Named interface: shows variable names in output
#' independence_test(V5, M = 5,
#'                   group_a = c("rating", "complaints"),
#'                   group_b = c("privileges", "learning"))
independence_test <- function(V, M = 1L,
                              part    = NULL,
                              group_a = NULL,
                              group_b = NULL,
                              alpha = 0.05,
                              iterations = 10000L,
                              null_dist = NULL) {

  V  <- .validate_X(V)
  M  <- .validate_M(M)
  N  <- nrow(V)
  n  <- N / M
  p  <- ncol(V)

  .check_N(N, p)

  res  <- .resolve_blocks(V, p, part, group_a, group_b,
                          "independence_test")
  V    <- res$V
  p1   <- res$p1
  p2   <- p - p1
  lbl1 <- res$lbl1
  lbl2 <- res$lbl2

  S_star <- .compute_S_star(V)
  idx1   <- seq_len(p1)
  idx2   <- seq_len(p2) + p1
  S_11   <- S_star[idx1, idx1, drop = FALSE]
  S_22   <- S_star[idx2, idx2, drop = FALSE]
  T_obs  <- det(S_star) / (det(S_11) * det(S_22))

  null_d <- if (!is.null(null_dist)) null_dist else
    Inddist(part = p1, nsample = n, pvariates = p, M = M, iterations = iterations)
  crit   <- as.numeric(quantile(null_d, probs = alpha))
  pval   <- mean(null_d <= T_obs)
  dec    <- if (T_obs < crit) "Reject H0" else "Fail to Reject H0"

  new_ps_test(
    statistic  = T_obs,
    p.value    = pval,
    alpha      = alpha,
    decision   = dec,
    null.dist  = null_d,
    test       = "independence",
    n          = n,
    M          = M,
    p          = p,
    lbl1       = lbl1,
    lbl2       = lbl2
  )
}


#' @title Regression Test
#'
#' @description
#' Tests \eqn{H_0 : \Delta = \Delta_0} for the population regression
#' matrix \eqn{\Delta = \Sigma_{12}\Sigma_{22}^{-1}}, based on \eqn{M}
#' released plug-in sampling synthetic datasets stacked into \code{V}.
#' Requires \eqn{p_1 \leq p_2}.
#' Setting \code{M = 1} recovers Klein et al. (2021).
#'
#' The two variable blocks can be specified in two ways (exactly one
#' must be used):
#' \describe{
#'   \item{\code{part} (integer scalar)}{First \code{part} columns are
#'     the response block (Block 1); remaining columns are the predictor
#'     block (Block 2). Original interface; backward-compatible.}
#'   \item{\code{response} and \code{predictors}}{Integer indices
#'     \emph{or} column names identifying the response and predictor
#'     blocks. Together they must cover all columns of \code{V}.}
#' }
#'
#' @param V Stacked synthetic dataset (\eqn{Mn \times p} matrix).
#' @param M Positive integer. Number of synthetic releases
#'   (default \code{1L}).
#' @param part Integer scalar. Size of the response block (Block 1,
#'   first \code{part} columns). Must satisfy \eqn{p_1 \leq p_2}.
#'   Ignored when \code{response}/\code{predictors} are supplied.
#' @param Delta0 A \eqn{p_1 \times p_2} matrix giving the null value
#'   \eqn{\Delta_0}. Default is the zero matrix (test of zero
#'   regression).
#' @param response Integer or character vector identifying the response
#'   block (Block 1).
#' @param predictors Integer or character vector identifying the
#'   predictor block (Block 2). Together with \code{response} must
#'   cover all columns.
#' @param alpha Significance level (default \code{0.05}).
#' @param iterations Monte Carlo sample size (default \code{10000L}).
#'
#' @return An object of class \code{\link{ps_test}}. Component
#'   \code{Delta.hat} gives the plug-in slope estimator
#'   \eqn{\hat\Delta = S_{12}^\star (S_{22}^\star)^{-1}}.
#'   The null hypothesis string in \code{$null.value} names both blocks.
#'
#' @seealso \code{\link{canodist}}, \code{\link{ps_test}}
#'
#' @references
#' Klein, M., Moura, R., and Sinha, B. (2021). Multivariate normal
#' inference based on singly imputed synthetic data under plug-in
#' sampling. \emph{Sankhya B}, 83, 273--287.
#'
#' @export
#'
#' @examples
#' data(ps_attitude)
#' set.seed(1)
#' V5 <- simSynthData(ps_attitude, M = 5)
#'
#' ## Integer interface (zero regression)
#' regression_test(V5, M = 5, part = 2L)
#'
#' ## Named interface with true Delta0 (should fail to reject)
#' b      <- partition(cov(ps_attitude),
#'                     part1 = c("rating", "complaints"))
#' Delta0 <- b$B %*% solve(b$D)
#' regression_test(V5, M = 5,
#'                  response   = c("rating", "complaints"),
#'                  predictors = c("privileges", "learning"),
#'                  Delta0     = Delta0)
regression_test <- function(V, M = 1L,
                            part       = NULL,
                            Delta0     = NULL,
                            response   = NULL,
                            predictors = NULL,
                            alpha = 0.05,
                            iterations = 10000L,
                            null_dist = NULL) {

  V  <- .validate_X(V)
  M  <- .validate_M(M)
  N  <- nrow(V)
  n  <- N / M
  p  <- ncol(V)

  .check_N(N, p)

  res  <- .resolve_blocks(V, p, part, response, predictors,
                          "regression_test")
  V    <- res$V
  p1   <- res$p1
  p2   <- p - p1
  lbl1 <- res$lbl1
  lbl2 <- res$lbl2

  if (p1 > p2)
    stop(sprintf(
      "Regression test requires p1 <= p2. Got p1 = %d, p2 = %d.",
      p1, p2), call. = FALSE)

  if (is.null(Delta0)) {
    Delta0 <- matrix(0.0, nrow = p1, ncol = p2)
  } else {
    Delta0 <- as.matrix(Delta0)
    if (!isTRUE(all.equal(dim(Delta0), c(p1, p2))))
      stop(sprintf("'Delta0' must be a %d x %d matrix.", p1, p2),
           call. = FALSE)
  }

  idx1   <- seq_len(p1)
  idx2   <- seq_len(p2) + p1

  S_star    <- .compute_S_star(V)
  S_11      <- S_star[idx1, idx1, drop = FALSE]
  S_12      <- S_star[idx1, idx2, drop = FALSE]
  S_22      <- S_star[idx2, idx2, drop = FALSE]
  S22_inv   <- solve(S_22)
  Delta_hat <- S_12 %*% S22_inv
  diff      <- Delta_hat - Delta0
  num       <- diff %*% S_22 %*% t(diff)
  den       <- S_11 - S_12 %*% S22_inv %*% t(S_12)
  T_obs     <- det(num) / det(den)

  null_d <- if (!is.null(null_dist)) null_dist else
    canodist(part = p1, nsample = n, pvariates = p, M = M, iterations = iterations)
  crit   <- as.numeric(quantile(null_d, probs = 1 - alpha))
  pval   <- mean(null_d >= T_obs)
  dec    <- if (T_obs > crit) "Reject H0" else "Fail to Reject H0"

  new_ps_test(
    statistic = T_obs,
    p.value   = pval,
    alpha     = alpha,
    decision  = dec,
    null.dist = null_d,
    test      = "regression",
    n         = n,
    M         = M,
    p         = p,
    Delta.hat = Delta_hat,
    lbl1      = lbl1,
    lbl2      = lbl2
  )
}
