#' @title Generate Plug-in Sampling Synthetic Dataset(s)
#'
#' @description
#' Generates \eqn{M \geq 1} independent fully synthetic datasets from an
#' original numeric matrix \code{X} using the plug-in sampling (PS)
#' mechanism under a multivariate normal model, and returns them as a
#' single \emph{stacked} \eqn{Mn \times p} matrix.
#'
#' The unknown population parameters \eqn{\boldsymbol{\mu}} and
#' \eqn{\boldsymbol{\Sigma}} are replaced by the sample mean
#' \eqn{\bar{\mathbf{x}}} and sample covariance matrix
#' \eqn{\hat{\boldsymbol{\Sigma}}}, and \eqn{Mn} synthetic observations
#' are drawn independently from
#' \eqn{\mathcal{N}_p(\bar{\mathbf{x}}, \hat{\boldsymbol{\Sigma}})}.
#'
#' Setting \code{M = 1} (the default) produces a single synthetic dataset
#' of size \eqn{n}, identical to the classical single-release PS procedure
#' of Klein et al. (2021). Setting \code{M > 1} produces the stacked
#' dataset \eqn{\mathbf{V}_{\mathrm{complete}}} used by the
#' multiple-release procedures:
#' \deqn{
#'   \mathbf{V}_{\mathrm{complete}} =
#'   \begin{pmatrix} \mathbf{V}_1 \\ \vdots \\ \mathbf{V}_M
#'   \end{pmatrix} \in \mathbb{R}^{Mn \times p}.
#' }
#'
#' @param X A numeric matrix or data frame of original confidential
#'   observations (\eqn{n \times p}). Rows are observations, columns are
#'   variables. Must satisfy \eqn{n > p}.
#' @param M A positive integer. Number of independent synthetic releases
#'   to generate (default \code{1L}). The returned matrix has \eqn{Mn}
#'   rows.
#'
#' @return An \eqn{Mn \times p} numeric matrix. Column names are preserved
#'   from \code{X}. Row names encode the release index and observation
#'   index in the form \code{"release_j.obs_i"} when \code{M > 1};
#'   for \code{M = 1} the row names of \code{X} are used if available.
#'
#' @details
#' The stacked representation is statistically justified because all
#' \eqn{Mn} rows are conditionally independent and identically distributed
#' given the original data, so the stacked sufficient statistic
#' \eqn{\mathbf{S}^\star_{\mathrm{complete}}} satisfies:
#' \deqn{
#'   \mathbf{S}^\star_{\mathrm{complete}} \mid \mathbf{S}
#'   \sim \mathcal{W}_p\!\left(Mn - 1,\; \tfrac{1}{n-1}\mathbf{S}\right).
#' }
#'
#' @seealso \code{\link{ps_test}}, \code{\link{gv_test}},
#'   \code{\link{sphericity_test}}, \code{\link{independence_test}},
#'   \code{\link{regression_test}}
#'
#' @references
#' Klein, M., Moura, R., and Sinha, B. (2021). Multivariate normal
#' inference based on singly imputed synthetic data under plug-in
#' sampling. \emph{Sankhya B}, 83, 273--287.
#' \doi{10.1007/s13571-019-00215-9}
#'
#' @export
#'
#' @examples
#' data(ps_attitude)
#'
#' # Single release (M = 1, default)
#' set.seed(1)
#' V1 <- simSynthData(ps_attitude)
#' dim(V1)   # 30 x 4  (same dimensions as ps_attitude)
#'
#' # Five releases stacked (M = 5)
#' set.seed(1)
#' V5 <- simSynthData(ps_attitude, M = 5)
#' dim(V5)   # 150 x 4  (5 * 30 rows)
#'
#' # M = 1 is the default: pass to any test function as before
#' sphericity_test(V1, M = 1)
#'
#' # M = 5: stacked data passed with M = 5
#' sphericity_test(V5, M = 5)
simSynthData <- function(X, M = 1L) {

  X <- .validate_X(X)
  M <- .validate_M(M)
  n <- nrow(X)
  p <- ncol(X)

  xbar  <- colMeans(X)
  Sigma <- cov(X)
  .check_pd(Sigma, "sample covariance matrix of X")

  if (M == 1L) {
    # Single release: preserve original row names if present
    V           <- MASS::mvrnorm(n = n, mu = xbar, Sigma = Sigma)
    colnames(V) <- colnames(X)
    if (!is.null(rownames(X))) rownames(V) <- rownames(X)
    return(V)
  }

  # Multiple releases: stack M independent datasets row-wise
  V_list <- lapply(seq_len(M), function(j) {
    Vj           <- MASS::mvrnorm(n = n, mu = xbar, Sigma = Sigma)
    colnames(Vj) <- colnames(X)
    rownames(Vj) <- paste0("release_", j, ".obs_", seq_len(n))
    Vj
  })

  do.call(rbind, V_list)
}

## ------------------------------------------------------------------
## Internal validators (used by all inference functions)
## ------------------------------------------------------------------

#' @keywords internal
.validate_X <- function(X, name = "X") {
  X <- as.matrix(X)
  if (!is.numeric(X))
    stop(sprintf("'%s' must be a numeric matrix or data frame.", name),
         call. = FALSE)
  n <- nrow(X); p <- ncol(X)
  if (n <= 1L)
    stop(sprintf("'%s' must have at least 2 rows.", name),
         call. = FALSE)
  if (p < 1L)
    stop(sprintf("'%s' must have at least 1 column.", name),
         call. = FALSE)
  if (n <= p)
    stop(sprintf(
      "Sample size n = %d must exceed the number of variables p = %d.",
      n, p), call. = FALSE)
  if (any(!is.finite(X)))
    stop(sprintf("'%s' contains non-finite values (NA, NaN, Inf).", name),
         call. = FALSE)
  X
}

#' @keywords internal
.validate_M <- function(M) {
  M <- suppressWarnings(as.integer(M))
  if (is.na(M) || M < 1L)
    stop("'M' must be a positive integer.", call. = FALSE)
  M
}

#' @keywords internal
.validate_part <- function(part, p) {
  p1 <- suppressWarnings(as.integer(part))
  if (is.na(p1) || p1 < 1L || p1 >= p)
    stop(sprintf(
      "'part' must be an integer in {1, ..., p-1}. Got %s with p = %d.",
      part, p), call. = FALSE)
  p1
}

#' @keywords internal
.check_pd <- function(A, name = "matrix") {
  ev <- eigen(A, only.values = TRUE)$values
  if (any(ev <= .Machine$double.eps * max(abs(ev)) * nrow(A)))
    stop(sprintf("The %s is not positive definite.", name),
         call. = FALSE)
}

#' @keywords internal
.check_N <- function(N, p) {
  if (N <= p + 1L)
    stop(sprintf(
      "Effective sample size N = Mn = %d must exceed p + 1 = %d. Increase n or M.",
      N, p + 1L), call. = FALSE)
}

#' @keywords internal
.compute_S_star <- function(V) {
  Vc <- sweep(V, 2L, colMeans(V), "-")
  crossprod(Vc)
}

## ------------------------------------------------------------------
## Internal block-resolution helper (shared by independence & regression)
## ------------------------------------------------------------------

#' @keywords internal
.resolve_blocks <- function(V, p, part, group_a, group_b,
                            fun_name) {
  cn <- colnames(V)

  # ---- integer interface (original) ------------------------------------
  if (!is.null(part) && is.null(group_a) && is.null(group_b)) {
    p1   <- .validate_part(part, p)
    idx1 <- seq_len(p1)
    idx2 <- seq_len(p - p1) + p1
    lbl1 <- if (!is.null(cn)) paste(cn[idx1], collapse = ", ") else
      paste0("cols 1-", p1)
    lbl2 <- if (!is.null(cn)) paste(cn[idx2], collapse = ", ") else
      paste0("cols ", p1 + 1L, "-", p)
    return(list(V    = V,
                p1   = p1,
                idx1 = idx1,
                idx2 = idx2,
                lbl1 = lbl1,
                lbl2 = lbl2))
  }

  # ---- named / indexed interface ---------------------------------------
  if (!is.null(group_a) && !is.null(group_b)) {
    .to_idx <- function(g, which_arg) {
      if (is.character(g)) {
        if (is.null(cn))
          stop(sprintf(
            "'%s' requires V to have column names when using names.",
            fun_name), call. = FALSE)
        bad <- setdiff(g, cn)
        if (length(bad))
          stop(sprintf(
            "'%s' contains names not found in colnames(V): %s.",
            which_arg, paste(bad, collapse = ", ")),
            call. = FALSE)
        match(g, cn)
      } else {
        as.integer(g)
      }
    }
    idx1 <- .to_idx(group_a, "group_a")
    idx2 <- .to_idx(group_b, "group_b")

    if (anyDuplicated(c(idx1, idx2)))
      stop(sprintf(
        "In '%s': 'group_a' and 'group_b' must not overlap.", fun_name),
        call. = FALSE)
    if (!setequal(c(idx1, idx2), seq_len(p)))
      stop(sprintf(
        "In '%s': 'group_a' and 'group_b' must cover all %d columns.",
        fun_name, p), call. = FALSE)

    p1   <- length(idx1)
    V    <- V[, c(idx1, idx2), drop = FALSE]
    lbl1 <- if (!is.null(cn)) paste(cn[idx1], collapse = ", ") else
      paste(idx1, collapse = ", ")
    lbl2 <- if (!is.null(cn)) paste(cn[idx2], collapse = ", ") else
      paste(idx2, collapse = ", ")
    return(list(V    = V,
                p1   = p1,
                idx1 = seq_len(p1),
                idx2 = seq_len(p - p1) + p1,
                lbl1 = lbl1,
                lbl2 = lbl2))
  }

  stop(sprintf(
    "In '%s': supply either 'part' (integer) or both 'group_a' and 'group_b'.",
    fun_name), call. = FALSE)
}
