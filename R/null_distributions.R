#' @title Simulate the Generalised Variance Null Distribution
#'
#' @description
#' Simulates the null distribution of the generalised variance
#' pivotal statistic \eqn{T_1^\star} under plug-in sampling.
#'
#' Under the multiple-release stacking result,
#' \eqn{\mathbf{S}^\star_{\mathrm{complete}} \mid \mathbf{S}
#' \sim \mathcal{W}_p(Mn-1, \mathbf{S}/(n-1))},
#' the Bartlett decomposition gives:
#' \deqn{
#'   T_1^\star = (n-1)^p \frac{|\mathbf{S}^\star|}{|\Sigma|}
#'   \;\overset{d}{=}\;
#'   \Bigl(\prod_{j=1}^p A_j\Bigr)\Bigl(\prod_{j=1}^p B_j\Bigr),
#' }
#' where \eqn{A_j \sim \chi^2_{Mn-j}} (from the synthetic Wishart,
#' degrees of freedom \eqn{Mn-1}) and \eqn{B_j \sim \chi^2_{n-j}}
#' (from the original sample Wishart, degrees of freedom \eqn{n-1}),
#' all \eqn{2p} variables mutually independent.
#' For \eqn{M = 1} we have \eqn{A_j \overset{d}{=} B_j \sim
#' \chi^2_{n-j}}, recovering the single-release result of
#' Klein et al. (2021).
#'
#' @param nsample Original sample size \eqn{n}.
#' @param pvariates Number of variables \eqn{p}.
#' @param iterations Number of Monte Carlo draws (default 10000).
#' @param M Number of synthetic releases (default \code{1L}).
#'   The effective sample size is \eqn{N = Mn}.
#'
#' @return A numeric vector of length \code{iterations} containing
#'   draws from the null distribution of \eqn{T_1^\star}.
#'
#' @seealso \code{\link{gv_test}}
#'
#' @references
#' Klein, M., Moura, R., and Sinha, B. (2021). Multivariate normal
#' inference based on singly imputed synthetic data under plug-in
#' sampling. \emph{Sankhya B}, 83, 273--287.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' # Single release (M = 1)
#' nd1 <- GVdist(nsample = 50, pvariates = 4)
#' quantile(nd1, probs = c(0.025, 0.975))
#'
#' # Five releases (M = 5): A_j ~ chi2_{250-j}, B_j ~ chi2_{50-j}
#' nd5 <- GVdist(nsample = 50, pvariates = 4, M = 5)
#' quantile(nd5, probs = c(0.025, 0.975))  # narrower
GVdist <- function(nsample, pvariates, iterations = 10000L,
                   M = 1L) {

  n  <- as.integer(nsample)
  p  <- as.integer(pvariates)
  it <- as.integer(iterations)
  M  <- as.integer(M)
  N  <- n * M   # effective sample size

  if (n <= p)
    stop("'nsample' must exceed 'pvariates'.", call. = FALSE)
  if (it < 1L)
    stop("'iterations' must be a positive integer.", call. = FALSE)
  if (M < 1L)
    stop("'M' must be a positive integer.", call. = FALSE)

  # A_j ~ chi2_{N-j}  (synthetic Wishart df = N-1 = Mn-1)
  # B_j ~ chi2_{n-j}  (original sample Wishart df = n-1)
  # For M=1: N=n so A_j = B_j ~ chi2_{n-j}
  df_A <- (N - 1L) - seq_len(p) + 1L   # Mn-1, Mn-2, ..., Mn-p
  df_B <- (n - 1L) - seq_len(p) + 1L   # n-1,  n-2,  ..., n-p

  replicate(it, {
    A <- prod(stats::rchisq(p, df = df_A))
    B <- prod(stats::rchisq(p, df = df_B))
    A * B
  })
}

#' @title Simulate the Sphericity Null Distribution
#'
#' @description
#' Simulates the null distribution of the sphericity pivotal statistic
#' \eqn{T_2^\star} under plug-in sampling for \eqn{M \geq 1} releases.
#'
#' Under the stacking result the compound Wishart representation gives:
#' \deqn{
#'   T_2^\star \overset{d}{=}
#'   \frac{|W_1 W_2|^{1/p}}{\mathrm{tr}(W_1 W_2)/p},
#' }
#' where \eqn{W_1 \sim \mathcal{W}_p(n-1,\,(n-1)^{-1} I_p)} (df \eqn{n-1})
#' and \eqn{W_2 \sim \mathcal{W}_p(Mn-1,\, I_p)} (df \eqn{Mn-1}), independently.
#' For \eqn{M=1} both have df \eqn{n-1}, recovering Klein et al. (2021).
#'
#' @param nsample Original sample size \eqn{n}.
#' @param pvariates Number of variables \eqn{p}.
#' @param iterations Number of Monte Carlo draws (default 10000).
#' @param M Number of synthetic releases (default \code{1L}).
#' @return Numeric vector of length \code{iterations}.
#' @seealso \code{\link{sphericity_test}}
#' @export
#' @examples
#' set.seed(1)
#' Sphdist(nsample = 50, pvariates = 4, M = 1) |> quantile(0.05)
#' Sphdist(nsample = 50, pvariates = 4, M = 5) |> quantile(0.05)
Sphdist <- function(nsample, pvariates, iterations = 10000L, M = 1L) {
  n  <- as.integer(nsample)
  p  <- as.integer(pvariates)
  it <- as.integer(iterations)
  M  <- as.integer(M)
  N  <- n * M
  if (n <= p) stop("'nsample' must exceed 'pvariates'.", call. = FALSE)
  if (it < 1L) stop("'iterations' must be a positive integer.", call. = FALSE)
  Ip <- diag(p)
  replicate(it, {
    W1 <- stats::rWishart(1L, df = n - 1L, Sigma = Ip / (n - 1L))[,,1L]
    W2 <- stats::rWishart(1L, df = N - 1L, Sigma = Ip)[,,1L]
    WW <- W1 %*% W2
    det(WW)^(1.0/p) / (sum(diag(WW))/p)
  })
}


#' @title Simulate the Independence Null Distribution
#'
#' @description
#' Simulates the null distribution of the independence pivotal statistic
#' \eqn{T_3^\star} under plug-in sampling for \eqn{M \geq 1} releases.
#'
#' Under the stacking result the compound Wishart representation gives:
#' \deqn{
#'   T_3^\star \overset{d}{=}
#'   \frac{|\Omega_2|}{|\Omega_{2,11}||\Omega_{2,22}|},
#' }
#' where \eqn{\Omega_1 \sim \mathcal{W}_p(n-1,\,I_p)} (df \eqn{n-1})
#' and \eqn{\Omega_2 \mid \Omega_1 \sim \mathcal{W}_p(Mn-1,\,(n-1)^{-1}\Omega_1)}
#' (df \eqn{Mn-1}).
#'
#' @param part Size of the first variable block \eqn{p_1}.
#' @param nsample Original sample size \eqn{n}.
#' @param pvariates Total number of variables \eqn{p}.
#' @param iterations Number of Monte Carlo draws (default 10000).
#' @param M Number of synthetic releases (default \code{1L}).
#' @return Numeric vector of length \code{iterations}.
#' @seealso \code{\link{independence_test}}
#' @export
#' @examples
#' set.seed(1)
#' Inddist(part = 2, nsample = 50, pvariates = 4, M = 1) |> quantile(0.05)
#' Inddist(part = 2, nsample = 50, pvariates = 4, M = 5) |> quantile(0.05)
Inddist <- function(part, nsample, pvariates, iterations = 10000L, M = 1L) {
  p1 <- as.integer(part)
  n  <- as.integer(nsample)
  p  <- as.integer(pvariates)
  it <- as.integer(iterations)
  M  <- as.integer(M)
  N  <- n * M
  p2 <- p - p1
  if (p1 < 1L || p1 >= p)
    stop("'part' must satisfy 1 <= part < pvariates.", call. = FALSE)
  if (n <= p) stop("'nsample' must exceed 'pvariates'.", call. = FALSE)
  if (it < 1L) stop("'iterations' must be a positive integer.", call. = FALSE)
  idx1 <- seq_len(p1);  idx2 <- seq_len(p2) + p1
  replicate(it, {
    Om1 <- stats::rWishart(1L, df = n - 1L, Sigma = diag(p))[,,1L]
    Om2 <- stats::rWishart(1L, df = N - 1L, Sigma = Om1/(n-1L))[,,1L]
    det(Om2) / (det(Om2[idx1,idx1,drop=FALSE]) * det(Om2[idx2,idx2,drop=FALSE]))
  })
}


#' @title Simulate the Regression Canonical Null Distribution
#'
#' @description
#' Simulates the null distribution of the regression pivotal statistic
#' \eqn{T_4^\star} under plug-in sampling for \eqn{M \geq 1} releases.
#'
#' Same compound Wishart structure as \code{\link{Inddist}}:
#' \eqn{\Omega_1 \sim \mathcal{W}_p(n-1,\,I_p)} and
#' \eqn{\Omega_2 \mid \Omega_1 \sim \mathcal{W}_p(Mn-1,\,(n-1)^{-1}\Omega_1)}.
#'
#' @param part Size of the first variable block \eqn{p_1}
#'   (\eqn{p_1 \leq p_2}).
#' @param nsample Original sample size \eqn{n}.
#' @param pvariates Total number of variables \eqn{p}.
#' @param iterations Number of Monte Carlo draws (default 10000).
#' @param M Number of synthetic releases (default \code{1L}).
#' @return Numeric vector of length \code{iterations}.
#' @seealso \code{\link{regression_test}}
#' @export
#' @examples
#' set.seed(1)
#' canodist(part = 2, nsample = 50, pvariates = 4, M = 1) |> quantile(0.95)
#' canodist(part = 2, nsample = 50, pvariates = 4, M = 5) |> quantile(0.95)
canodist <- function(part, nsample, pvariates, iterations = 10000L, M = 1L) {
  p1 <- as.integer(part)
  n  <- as.integer(nsample)
  p  <- as.integer(pvariates)
  it <- as.integer(iterations)
  M  <- as.integer(M)
  N  <- n * M
  p2 <- p - p1
  if (p1 < 1L || p1 >= p)
    stop("'part' must satisfy 1 <= part < pvariates.", call. = FALSE)
  if (p1 > p2) stop("Regression test requires p1 <= p2.", call. = FALSE)
  if (n <= p) stop("'nsample' must exceed 'pvariates'.", call. = FALSE)
  if (it < 1L) stop("'iterations' must be a positive integer.", call. = FALSE)
  idx1 <- seq_len(p1);  idx2 <- seq_len(p2) + p1
  replicate(it, {
    Om1   <- stats::rWishart(1L, df = n - 1L, Sigma = diag(p))[,,1L]
    Om2   <- stats::rWishart(1L, df = N - 1L, Sigma = Om1/(n-1L))[,,1L]
    Om_11 <- Om2[idx1,idx1,drop=FALSE];  Om_12 <- Om2[idx1,idx2,drop=FALSE]
    Om_21 <- Om2[idx2,idx1,drop=FALSE];  Om_22 <- Om2[idx2,idx2,drop=FALSE]
    Om22i <- solve(Om_22)
    num   <- Om_12 %*% Om22i %*% Om_21
    det(num) / det(Om_11 - num)
  })
}
