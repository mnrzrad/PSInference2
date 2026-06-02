#' @title Unified Wrapper for PS Inference Tests
#'
#' @description
#' Dispatches to the appropriate exact inferential procedure based on
#' the \code{test} argument, and optionally produces the diagnostic
#' plot immediately. This is the main entry point for users who prefer
#' a single function over the four individual test functions.
#'
#' @param V Stacked synthetic dataset (\eqn{Mn \times p} matrix), as
#'   returned by \code{\link{simSynthData}}.
#' @param M Positive integer. Number of synthetic releases
#'   (default \code{1L}). Setting \code{M = 1} recovers the
#'   single-release procedures of Klein et al. (2021).
#' @param test Character string specifying the test. One of
#'   \code{"gv"} (generalised variance),
#'   \code{"sphericity"},
#'   \code{"independence"},
#'   \code{"regression"}.
#' @param plot Logical. If \code{TRUE}, calls \code{plot()} on the
#'   result before returning it, so the null distribution diagnostic
#'   is displayed automatically. Default \code{FALSE}.
#' @param ... Additional arguments passed to the corresponding test
#'   function or, when \code{plot = TRUE}, to \code{plot.ps_test()}.
#'   Arguments intended for the test function (e.g. \code{part},
#'   \code{Sigma}, \code{Delta0}, \code{iterations}) and graphical
#'   arguments (e.g. \code{main}, \code{shade_col}) are separated
#'   automatically.
#'
#' @return An object of class \code{\link{ps_test}}, invisibly when
#'   \code{plot = TRUE}.
#'
#' @seealso \code{\link{gv_test}}, \code{\link{sphericity_test}},
#'   \code{\link{independence_test}}, \code{\link{regression_test}},
#'   \code{\link{plot.ps_test}}
#'
#' @export
#'
#' @examples
#' data(ps_attitude)
#' set.seed(1)
#' V <- simSynthData(ps_attitude, M = 3)
#'
#' # Run and print only
#' ps_test(V, M = 3, test = "sphericity")
#'
#' # Run and plot in one call
#' ps_test(V, M = 3, test = "sphericity", plot = TRUE)
#'
#' # Independence with named blocks, plot automatically
#' ps_test(V, M = 3, test = "independence",
#'         group_a = c("rating", "complaints"),
#'         group_b = c("privileges", "learning"),
#'         plot = TRUE)
#'
#' # Generalised variance with reference covariance
#' ps_test(V, M = 3, test = "gv",
#'         Sigma = cov(ps_attitude), plot = TRUE)
#'
#' # M = 1 recovers Klein et al. (2021)
#' ps_test(simSynthData(ps_attitude), M = 1,
#'         test = "sphericity", plot = TRUE)
ps_test <- function(V,
                    M    = 1L,
                    test = c("gv", "sphericity",
                             "independence", "regression"),
                    plot = FALSE,
                    ...) {
  test <- match.arg(test)

  ## ---- separate test arguments from plot arguments --------------------
  ## Plot arguments recognised by plot.ps_test:
  plot_arg_names <- c("main", "shade_col", "dist_col",
                      "stat_col", "crit_col")
  dots      <- list(...)
  plot_args <- dots[intersect(names(dots), plot_arg_names)]
  test_args <- dots[setdiff(names(dots), plot_arg_names)]

  ## ---- run the test ---------------------------------------------------
  res <- switch(test,
                gv           = do.call(gv_test,
                                       c(list(V = V, M = M), test_args)),
                sphericity   = do.call(sphericity_test,
                                       c(list(V = V, M = M), test_args)),
                independence = do.call(independence_test,
                                       c(list(V = V, M = M), test_args)),
                regression   = do.call(regression_test,
                                       c(list(V = V, M = M), test_args))
  )

  ## ---- optional plot --------------------------------------------------
  if (plot) {
    do.call(graphics::plot, c(list(x = res), plot_args))
    return(invisible(res))
  }

  res
}


#' Partition a Matrix into Four Blocks
#'
#' Splits a numeric matrix \eqn{\mathbf{M}} into four sub-matrices according
#' to a two-group partition of its rows and columns:
#' \deqn{\mathbf{M} = \begin{bmatrix} \mathbf{A} & \mathbf{B} \\
#'                                     \mathbf{C} & \mathbf{D} \end{bmatrix}}
#' where
#' \describe{
#'   \item{\eqn{\mathbf{A}}}{part1 \eqn{\times} part1}
#'   \item{\eqn{\mathbf{B}}}{part1 \eqn{\times} part2}
#'   \item{\eqn{\mathbf{C}}}{part2 \eqn{\times} part1}
#'   \item{\eqn{\mathbf{D}}}{part2 \eqn{\times} part2}
#' }
#'
#' \strong{Two interfaces:}
#' \describe{
#'   \item{Integer interface (original)}{Supply \code{nrows} and \code{ncols}.
#'     The first \code{nrows} rows and \code{ncols} columns form part1;
#'     the rest form part2.}
#'   \item{Name interface}{Supply \code{part1} as a character vector of
#'     row/column names. The matrix is reordered so that \code{part1} rows
#'     and columns appear first. \code{part2} is optional; if supplied it
#'     is validated as the complement of \code{part1}.}
#' }
#'
#' @param Matrix  A numeric matrix. Must have dimnames when using the name
#'   interface.
#' @param nrows   Integer; number of part1 rows. Ignored when \code{part1}
#'   is supplied.
#' @param ncols   Integer; number of part1 columns. Ignored when \code{part1}
#'   is supplied.
#' @param part1   Character vector of variable names forming the first group.
#'   Used for both row and column reordering of a square matrix.
#' @param part2   Optional character vector naming the second group. If
#'   supplied it is checked to equal the complement of \code{part1}.
#'
#' @return A named list of class \code{"ps_partition"} with elements
#'   \code{A}, \code{B}, \code{C}, \code{D}.
#'   Numeric indexing \code{[[1]]}--\code{[[4]]} also works.
#'
#' @export
#'
#' @examples
#' ## Integer interface
#' M <- matrix(1:16, 4, 4,
#'             dimnames = list(c("A","B","C","D"), c("A","B","C","D")))
#' b <- partition(M, nrows = 2, ncols = 2)
#' b$A   # rows A,B x cols A,B
#'
#' ## Name interface
#' b2 <- partition(M, part1 = c("A","B"))
#' b2$A  # rows A,B x cols A,B  (part1 x part1)
#' b2$D  # rows C,D x cols C,D  (part2 x part2)
#'
#' ## Covariance matrix of ps_attitude
#' data(ps_attitude)
#' b3 <- partition(cov(ps_attitude),
#'                 part1 = c("rating", "complaints"),
#'                 part2 = c("privileges", "learning"))
#' b3$A  # 2x2: rating, complaints
#' b3$D  # 2x2: privileges, learning
partition <- function(Matrix, nrows = NULL, ncols = NULL,
                      part1 = NULL, part2 = NULL) {

  ## ── input checks ────────────────────────────────────────────────────────
  if (!is.matrix(Matrix) || !is.numeric(Matrix))
    stop("'Matrix' must be a numeric matrix.", call. = FALSE)
  if (anyNA(Matrix))
    stop("'Matrix' contains NA values.", call. = FALSE)

  r  <- nrow(Matrix)
  cc <- ncol(Matrix)
  rn <- rownames(Matrix)
  cn <- colnames(Matrix)

  ## ── name interface ───────────────────────────────────────────────────────
  if (!is.null(part1)) {

    if (!is.character(part1) || length(part1) < 1L)
      stop("'part1' must be a non-empty character vector of names.",
           call. = FALSE)
    if (is.null(rn) || is.null(cn))
      stop("'Matrix' must have row and column names when 'part1' is supplied.",
           call. = FALSE)

    bad <- setdiff(part1, rn)
    if (length(bad))
      stop("'part1' contains names not found in rownames(Matrix): ",
           paste(bad, collapse = ", "), ".", call. = FALSE)
    bad <- setdiff(part1, cn)
    if (length(bad))
      stop("'part1' contains names not found in colnames(Matrix): ",
           paste(bad, collapse = ", "), ".", call. = FALSE)

    part2_rows <- setdiff(rn, part1)
    part2_cols <- setdiff(cn, part1)
    if (length(part2_rows) == 0L || length(part2_cols) == 0L)
      stop("'part1' covers all rows or columns -- part2 would be empty.",
           call. = FALSE)

    if (!is.null(part2)) {
      if (!is.character(part2))
        stop("'part2' must be a character vector of names.", call. = FALSE)
      expected <- setdiff(rn, part1)
      if (!setequal(part2, expected))
        stop("'part2' must be the complement of 'part1'. Got: ",
             paste(sort(part2), collapse = ", "), ". Expected: ",
             paste(sort(expected), collapse = ", "), ".", call. = FALSE)
    }

    Matrix <- Matrix[c(part1, part2_rows), c(part1, part2_cols),
                     drop = FALSE]
    nrows  <- length(part1)
    ncols  <- length(part1)
    r      <- nrow(Matrix)
    cc     <- ncol(Matrix)

  } else {
    ## ── integer interface ──────────────────────────────────────────────────
    if (is.null(nrows) || is.null(ncols))
      stop("Supply either 'part1' (character vector) or both 'nrows' and ",
           "'ncols' (integers).", call. = FALSE)
    if (!is.numeric(nrows) || length(nrows) != 1L || is.na(nrows) ||
        nrows < 1L || nrows != round(nrows) || nrows >= r)
      stop(sprintf("'nrows' must be a single integer in [1, %d). Got: %s.",
                   r, as.character(nrows)), call. = FALSE)
    if (!is.numeric(ncols) || length(ncols) != 1L || is.na(ncols) ||
        ncols < 1L || ncols != round(ncols) || ncols >= cc)
      stop(sprintf("'ncols' must be a single integer in [1, %d). Got: %s.",
                   cc, as.character(ncols)), call. = FALSE)
  }

  ## ── split ────────────────────────────────────────────────────────────────
  idx_r1 <- seq_len(nrows)
  idx_r2 <- seq.int(nrows + 1L, r)
  idx_c1 <- seq_len(ncols)
  idx_c2 <- seq.int(ncols + 1L, cc)

  A <- Matrix[idx_r1, idx_c1, drop = FALSE]
  B <- Matrix[idx_r1, idx_c2, drop = FALSE]
  C <- Matrix[idx_r2, idx_c1, drop = FALSE]
  D <- Matrix[idx_r2, idx_c2, drop = FALSE]

  structure(list(A = A, B = B, C = C, D = D),
            class = "ps_partition")
}

#' @export
`[[.ps_partition` <- function(x, i) unclass(x)[[i]]

#' @export
print.ps_partition <- function(x, ...) {
  cat("ps_partition:",
      " A (", nrow(x$A), "x", ncol(x$A), ")",
      " B (", nrow(x$B), "x", ncol(x$B), ")",
      " C (", nrow(x$C), "x", ncol(x$C), ")",
      " D (", nrow(x$D), "x", ncol(x$D), ")\n",
      sep = "")
  invisible(x)
}

