# Thin wrapper around survey::calibrate(). Does NOT reimplement it: this
# translates the 2D matrix problem into the long-format table calibrate()
# expects, calls it, and translates the result back into a matrix. The
# guarantees, limitations and convergence behaviour are those of the
# survey package, not of RASbound.

#' @keywords internal
#' @noRd
.matrix_to_long <- function(x0) {
  m <- nrow(x0); n <- ncol(x0)
  data.frame(
    row_id = factor(rep(seq_len(m), times = n)),
    col_id = factor(rep(seq_len(n), each = m)),
    w = as.vector(x0)
  )
}

#' Drop rows/columns whose target is exactly zero (internal)
#'
#' Under non-negative bounds a row or column with target 0 forces every
#' cell in it to 0, so it can be removed from the problem and reinserted
#' as zeros afterwards without changing the solution. This is only valid
#' when no cell in it has a positive lower bound (otherwise the problem is
#' infeasible and an error is raised).
#' @keywords internal
#' @noRd
.drop_zero_margins <- function(x0, row_targets, col_targets, lower, upper) {
  zr <- row_targets == 0
  zc <- col_targets == 0
  if (any(zr) && any(lower[zr, , drop = FALSE] > 0)) {
    stop("A row with target 0 contains a cell with a positive `lower` bound; ",
         "the problem is infeasible.", call. = FALSE)
  }
  if (any(zc) && any(lower[, zc, drop = FALSE] > 0)) {
    stop("A column with target 0 contains a cell with a positive `lower` bound; ",
         "the problem is infeasible.", call. = FALSE)
  }
  keep_rows <- which(!zr)
  keep_cols <- which(!zc)
  list(
    x0 = x0[keep_rows, keep_cols, drop = FALSE],
    row_targets = row_targets[keep_rows],
    col_targets = col_targets[keep_cols],
    lower = lower[keep_rows, keep_cols, drop = FALSE],
    upper = upper[keep_rows, keep_cols, drop = FALSE],
    keep_rows = keep_rows, keep_cols = keep_cols,
    m = length(keep_rows), n = length(keep_cols)
  )
}

#' Reinsert zero rows/columns removed by `.drop_zero_margins()` (internal)
#' @keywords internal
#' @noRd
.reinsert_zero_margins <- function(x_sub, drop_info, m, n) {
  out <- matrix(0, m, n)
  out[drop_info$keep_rows, drop_info$keep_cols] <- x_sub
  out
}

#' Closed-form solution when only one row or column (or none) remains (internal)
#'
#' `survey::calibrate()` cannot handle a main-effects formula with a
#' single-level factor, and the problem is trivial in that case: with one
#' remaining row, each cell equals its column target (and symmetrically
#' for one column). Errors if that value violates the cell's bounds.
#' @keywords internal
#' @noRd
.calibrate_trivial <- function(drop, m0, n0, row_targets0, col_targets0,
                               lower, upper, tol, call) {
  m <- drop$m; n <- drop$n
  if (m == 0L || n == 0L) {
    x_sub <- matrix(0, m, n)
  } else {
    x_sub <- if (m == 1L) matrix(drop$col_targets, 1L, n) else matrix(drop$row_targets, m, 1L)
    if (any(x_sub < lower - tol) || any(x_sub > upper + tol)) {
      stop("Only one row or column has a non-zero target, which forces every ",
           "remaining cell to a single value, and that value violates the ",
           "supplied bounds; the problem is infeasible.", call. = FALSE)
    }
  }
  x_hat <- .reinsert_zero_margins(x_sub, drop, m0, n0)
  row_err <- abs(rowSums(x_hat) - row_targets0)
  col_err <- abs(colSums(x_hat) - col_targets0)
  .new_bras_fit(
    matrix = x_hat, converged = max(row_err, col_err) < tol, feasible = NA,
    n_outer_iter = 0L,
    fixed_cells = data.frame(row = integer(), col = integer(),
                              bound = character(), value = numeric(),
                              outer_pass = integer()),
    row_margin_error = row_err, col_margin_error = col_err,
    infeasible_at = NULL, method = "calibrate (closed form, single row/column)",
    call = call
  )
}

#' Bounded RAS via `survey::calibrate(calfun = "logit")` (wrapper)
#'
#' **This function is a thin wrapper around [survey::calibrate()].** It
#' gives the same arguments and the same `bras_fit` result as the rest of
#' RASbound, but it is not an independent implementation. The algorithm,
#' its guarantees and its limits are those of the `survey` package. For
#' example, `survey` documents that the Newton-Raphson iterations of
#' `calfun = "logit"` can fail to converge when the bounds are very tight.
#' If that happens, the error comes from `survey`, not from RASbound. See
#' `citation("survey")`.
#'
#' Internally, each cell of `x0` becomes one row of a long table, with the
#' cell value as its calibration weight and a formula with row and column
#' effects only. `survey` enforces the bounds as real constraints inside
#' its Newton-Raphson solve, not as a correction afterwards. The `"logit"`
#' and `"linear"` versions differ only in the distance function they
#' minimize. This is a different approach from [bras_freeze_sim()] and the
#' other freeze functions. It optimizes the whole matrix in one step
#' instead of solving, finding violations, fixing them and repeating.
#'
#' Rows and columns with a target of exactly 0 are removed before calling
#' `survey` and put back as zeros afterwards. An error is raised if such a
#' row or column has a positive `lower` bound.
#'
#' @section This is not the same solution as classic RAS/IPF:
#' Only [bras_calibrate_raking()] minimizes the same multiplicative
#' distance as classic RAS. Without bounds, it gives the same matrix as
#' [ras()], up to numerical tolerance. `calfun = "logit"` (this function)
#' and `calfun = "linear"` ([bras_calibrate_linear()]) minimize other
#' distance functions (a bounded logistic one and a least-squares one).
#' Without bounds, both still meet the same margins as classic RAS, but
#' they generally return a different matrix. See
#' `vignette("RASbound-diagnostics")`.
#'
#' @section Known fragility on real large tables:
#' On a real 57 x 57 input-output-style table, all three
#' `bras_calibrate_*()` wrappers failed inside `survey::calibrate()`, while
#' the `bras_freeze_*()` functions converged. Removing the zero-target rows
#' and columns did not fix this. The cause is not known. Differences in
#' scale between cells are one possibility that has not been tested. If a
#' wrapper fails on your data, use [bras_freeze_sim()], [bras_freeze_seq()]
#' or [bras_freeze_hybrid()].
#'
#' @section Logit calibration cannot represent an unbounded cell:
#' The logit calibration of `survey` maps each weight onto the interval
#' `(lower, upper)` with a logistic transform. That transform is not
#' defined for an infinite bound, and `survey` stops with the error
#' `"Logit calibration requires finite bounds"`. The other functions of
#' this package use `upper = Inf` for "no upper bound". So this wrapper
#' replaces every infinite upper bound by
#' `unbounded_mult * max(row_targets, col_targets)` (20 times by default).
#' No cell can reach such a value, and it is far enough from the
#' saturation zone of the logistic transform. It also scales with the
#' marginal totals, not with the matrix size.
#'
#' A tighter replacement is not used. The exact maximum of a cell,
#' `min(row_targets[i], col_targets[j])`, is correct but numerically
#' fragile: a bound that close to the true supremum is where the logistic
#' transform saturates and the Newton-Raphson iterations stall (`survey`
#' then reports `"Failed to converge"`). Increase `unbounded_mult` if that
#' still happens on your data. An infinite lower bound (unusual, since the
#' default is `0`) is replaced by a large negative number.
#'
#' @section Failure mode: a cell at a real, tight bound:
#' Suppose a cell has a real, tight `upper` and its seed in `x0` is at or
#' above it. The seed is moved to exactly that bound before the
#' calibration starts (see `vignette("RASbound-diagnostics")`). If most
#' other cells have no bound and get the replacement above, one cell sits
#' on a tight bound next to cells whose bound is tens of times larger.
#' This combination can make the Newton-Raphson solve diverge, with the
#' error `"infinite or missing values"` from an internal routine, instead
#' of converging slowly. `vignette("RASbound-diagnostics")` has a
#' reproducible example.
#'
#' If you know that some cells need a tight bound, give real finite bounds
#' to the other cells too, even generous ones, instead of relying on the
#' replacement of `Inf`. Or use [bras_calibrate_raking()] or
#' [bras_calibrate_linear()]. Both accept `Inf` and need no replacement,
#' so they are more robust when most cells are effectively unbounded.
#'
#' @inheritParams bras_freeze_sim
#' @param tol Tolerance used to decide `converged`. It is checked against
#'   the margin error of the returned matrix. It is looser than the
#'   default of the freeze functions, because the solution comes from
#'   `survey` and is reshaped back to a matrix.
#' @param unbounded_mult Multiplier used to replace an infinite `upper`
#'   (see "Logit calibration cannot represent an unbounded cell"). The
#'   replacement is `unbounded_mult * max(row_targets, col_targets)`.
#'   Increase it if `survey::calibrate()` still fails to converge. Decrease
#'   it only for a specific reason, since a tighter replacement is more
#'   fragile.
#' @param ... Additional arguments passed on to [survey::calibrate()].
#' @return An object of class `bras_fit`; see `?bras_fit`. `converged` is
#'   computed from the margin error of the returned matrix. `feasible` is
#'   `TRUE` with the default `force = FALSE`, because
#'   `survey::calibrate()` stops with an error instead of returning an
#'   unconverged result. If you pass `force = TRUE` through `...`, this
#'   `TRUE` is optimistic, so check `converged` as well.
#' @examples
#' if (requireNamespace("survey", quietly = TRUE) &&
#'     requireNamespace("MASS", quietly = TRUE)) {
#'   x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
#'   # upper = Inf on most cells is fine here: internally substituted with
#'   # unbounded_mult * max(row_target, col_target), since logit
#'   # calibration requires finite bounds (see "Logit calibration cannot
#'   # represent an unbounded cell").
#'   fit <- bras_calibrate_logit(x0, c(20, 20), c(15, 25),
#'                                upper = matrix(c(Inf, Inf, 12, Inf), nrow = 2))
#'   fit$matrix
#' }
#' @export
bras_calibrate_logit <- function(x0, row_targets, col_targets, lower = 0, upper = Inf,
                                  tol = 1e-6, unbounded_mult = 20, ...) {
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package 'survey' is required for bras_calibrate_logit() but is not installed.",
         call. = FALSE)
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required by survey::calibrate() for bounded calibration ",
         "but is not installed. Install it with install.packages(\"MASS\").",
         call. = FALSE)
  }
  .validate_inputs(x0, row_targets, col_targets, lower, upper)
  m <- nrow(x0); n <- ncol(x0)
  lower_mat <- .as_bound_matrix(lower, m, n)
  upper_mat <- .as_bound_matrix(upper, m, n)
  if (any(lower_mat > upper_mat)) stop("`lower` must be <= `upper` in every cell.", call. = FALSE)
  x0 <- .clamp_seed_to_bounds(x0, lower_mat, upper_mat)
  m0 <- m; n0 <- n
  row_targets0 <- row_targets; col_targets0 <- col_targets
  drop <- .drop_zero_margins(x0, row_targets, col_targets, lower_mat, upper_mat)
  x0 <- drop$x0; row_targets <- drop$row_targets; col_targets <- drop$col_targets
  lower_mat <- drop$lower; upper_mat <- drop$upper
  m <- drop$m; n <- drop$n
  if (m <= 1L || n <= 1L) {
    return(.calibrate_trivial(drop, m0, n0, row_targets0, col_targets0,
                              lower_mat, upper_mat, tol, match.call()))
  }

  # Logit calibration requires finite bounds (see the "Logit calibration
  # cannot represent an unbounded cell" section above). Substitute any
  # non-finite upper bound with `unbounded_mult` times the largest single
  # row or column target -- loose enough to stay well clear of the
  # logistic transform's saturation zone (see that section for why the
  # tighter, mathematically exact marginal cap was tried and rejected),
  # and unlike a fixed sentinel or the total problem mass, this scales the
  # same way regardless of matrix size: it is always exactly
  # `unbounded_mult` times the biggest marginal, not `nrow`/`ncol` times it.
  scale_hint <- max(1, sum(abs(row_targets)), sum(abs(col_targets)))
  upper_mat[!is.finite(upper_mat)] <- max(row_targets, col_targets) * unbounded_mult
  lower_mat[!is.finite(lower_mat)] <- -1e9 * scale_hint
  # Defensive floor: a degenerate `row_targets`/`col_targets` pair (very
  # small or negative) could otherwise push the substituted upper below
  # `lower`, producing an invalid bound that `survey::calibrate()` would
  # fail on internally with an opaque error rather than a clear one from us.
  if (any(lower_mat > upper_mat)) {
    stop("`lower` > `upper` at one or more cells after substituting finite ",
         "bounds for calfun = \"logit\" (this can happen with a very small ",
         "or non-positive row/column target). Check `row_targets` and ",
         "`col_targets` are all reasonably positive.", call. = FALSE)
  }

  df <- .matrix_to_long(x0)
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)

  row_lv <- levels(df$row_id); col_lv <- levels(df$col_id)
  pop <- c(
    `(Intercept)` = sum(row_targets),
    stats::setNames(row_targets[-1], paste0("row_id", row_lv[-1])),
    stats::setNames(col_targets[-1], paste0("col_id", col_lv[-1]))
  )

  cal <- tryCatch(
    survey::calibrate(
      des, formula = ~row_id + col_id, population = pop,
      calfun = "logit",
      bounds = list(lower = as.vector(lower_mat), upper = as.vector(upper_mat)),
      bounds.const = TRUE,
      ...
    ),
    error = function(e) {
      stop("survey::calibrate(calfun = \"logit\") failed: ", conditionMessage(e), ". ",
           "This calibration function requires finite bounds and can be numerically ",
           "fragile when most cells are only loosely bounded (see the \"Logit ",
           "calibration cannot represent an unbounded cell\" section of ",
           "?bras_calibrate_logit). Try increasing `unbounded_mult`, switching to ",
           "bras_calibrate_linear() instead -- it does not require finite bounds and ",
           "does not need this substitution -- or supplying a tighter, meaningful ",
           "`upper` yourself.",
           call. = FALSE)
    }
  )

  x_hat <- .reinsert_zero_margins(matrix(stats::weights(cal), nrow = m, ncol = n),
                                   drop, m0, n0)
  row_err <- abs(rowSums(x_hat) - row_targets0)
  col_err <- abs(colSums(x_hat) - col_targets0)

  .new_bras_fit(
    matrix = x_hat, converged = max(row_err, col_err) < tol, feasible = TRUE,
    n_outer_iter = 1L,
    fixed_cells = data.frame(row = integer(), col = integer(),
                              bound = character(), value = numeric(),
                              outer_pass = integer()),
    row_margin_error = row_err, col_margin_error = col_err,
    infeasible_at = NULL, method = "calibrate_logit (survey, wrapper)",
    call = match.call()
  )
}

#' Bounded RAS via `survey::calibrate(calfun = "linear")` (wrapper)
#'
#' **This function is a thin wrapper around [survey::calibrate()]**, for
#' the same reason as [bras_calibrate_logit()]. Like the logit version, it
#' enforces `lower` and `upper` as real constraints inside the same
#' Newton-Raphson solve (Deville, Saerndal and Sautory, 1993). It is not a
#' correction made afterwards. It differs from [bras_calibrate_logit()] in
#' the distance function (linear, that is, least squares) and in one
#' practical point: **its bounds are optional and `Inf` is accepted**, so
#' no replacement is needed. See the section "This is not the same
#' solution as classic RAS/IPF" in [bras_calibrate_logit()]. Even without
#' bounds, linear calibration generally returns a different matrix from
#' classic RAS, with the same margins. See `citation("survey")`.
#'
#' Rows and columns with a target of exactly 0 are removed before calling
#' `survey` and put back as zeros afterwards.
#'
#' @inheritParams bras_calibrate_logit
#' @return An object of class `bras_fit`; see `?bras_fit`. `converged` is
#'   computed from the margin error of the returned matrix. `feasible` is
#'   always `NA`, because the bounds are optional here and a successful
#'   return does not prove that they could be met. Check `converged`
#'   instead.
#' @examples
#' if (requireNamespace("survey", quietly = TRUE) &&
#'     requireNamespace("MASS", quietly = TRUE)) {
#'   x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
#'   fit <- bras_calibrate_linear(x0, c(20, 20), c(15, 25),
#'                                 upper = matrix(c(Inf, Inf, 12, Inf), nrow = 2))
#'   fit$matrix
#'   fit$converged
#' }
#' @export
bras_calibrate_linear <- function(x0, row_targets, col_targets, lower = 0, upper = Inf,
                                   tol = 1e-6, ...) {
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package 'survey' is required for bras_calibrate_linear() but is not installed.",
         call. = FALSE)
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required by survey::calibrate() for bounded calibration ",
         "but is not installed. Install it with install.packages(\"MASS\").",
         call. = FALSE)
  }
  .validate_inputs(x0, row_targets, col_targets, lower, upper)
  m <- nrow(x0); n <- ncol(x0)
  lower_mat <- .as_bound_matrix(lower, m, n)
  upper_mat <- .as_bound_matrix(upper, m, n)
  if (any(lower_mat > upper_mat)) stop("`lower` must be <= `upper` in every cell.", call. = FALSE)
  x0 <- .clamp_seed_to_bounds(x0, lower_mat, upper_mat)
  m0 <- m; n0 <- n
  row_targets0 <- row_targets; col_targets0 <- col_targets
  drop <- .drop_zero_margins(x0, row_targets, col_targets, lower_mat, upper_mat)
  x0 <- drop$x0; row_targets <- drop$row_targets; col_targets <- drop$col_targets
  lower_mat <- drop$lower; upper_mat <- drop$upper
  m <- drop$m; n <- drop$n
  if (m <= 1L || n <= 1L) {
    return(.calibrate_trivial(drop, m0, n0, row_targets0, col_targets0,
                              lower_mat, upper_mat, tol, match.call()))
  }

  df <- .matrix_to_long(x0)
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)

  row_lv <- levels(df$row_id); col_lv <- levels(df$col_id)
  pop <- c(
    `(Intercept)` = sum(row_targets),
    stats::setNames(row_targets[-1], paste0("row_id", row_lv[-1])),
    stats::setNames(col_targets[-1], paste0("col_id", col_lv[-1]))
  )

  cal <- tryCatch(
    survey::calibrate(
      des, formula = ~row_id + col_id, population = pop,
      calfun = "linear",
      bounds = list(lower = as.vector(lower_mat), upper = as.vector(upper_mat)),
      bounds.const = TRUE,
      ...
    ),
    error = function(e) {
      stop("survey::calibrate(calfun = \"linear\") failed: ", conditionMessage(e), ". ",
           "See ?bras_calibrate_linear.", call. = FALSE)
    }
  )

  x_hat <- .reinsert_zero_margins(matrix(stats::weights(cal), nrow = m, ncol = n),
                                   drop, m0, n0)
  row_err <- abs(rowSums(x_hat) - row_targets0)
  col_err <- abs(colSums(x_hat) - col_targets0)

  .new_bras_fit(
    matrix = x_hat, converged = max(row_err, col_err) < tol, feasible = NA,
    n_outer_iter = 1L,
    fixed_cells = data.frame(row = integer(), col = integer(),
                              bound = character(), value = numeric(),
                              outer_pass = integer()),
    row_margin_error = row_err, col_margin_error = col_err,
    infeasible_at = NULL, method = "calibrate_linear (survey, wrapper)",
    call = match.call()
  )
}

#' Bounded RAS via `survey::calibrate(calfun = "raking")` (wrapper)
#'
#' **This function is a thin wrapper around [survey::calibrate()]**, for
#' the same reason as [bras_calibrate_logit()] and
#' [bras_calibrate_linear()]. Unlike those two, `calfun = "raking"`
#' minimizes the same multiplicative distance as classic RAS/IPF. Without
#' bounds, it is **expected to give essentially the same matrix as
#' [ras()]**, up to numerical tolerance. It gets there with a different
#' algorithm: the Newton-Raphson solve of `survey` (Deville, Saerndal and
#' Sautory, 1993) instead of this package's alternating row and column
#' scaling.
#'
#' The package's cross-validation test checks this agreement (see
#' `vignette("RASbound-diagnostics")`). That is why this wrapper exists. It
#' is the one external result that should match this package's own engine
#' on a problem without bounds, so it is the most direct check against an
#' independent implementation.
#'
#' With bounds (what `survey` calls bounded raking), the bounds are
#' enforced as real constraints inside the same solve. As with
#' [bras_calibrate_linear()], they are optional and `Inf` is accepted, so
#' no replacement is needed, unlike [bras_calibrate_logit()]. See
#' `citation("survey")`.
#'
#' Rows and columns with a target of exactly 0 are removed before calling
#' `survey` and put back as zeros afterwards.
#'
#' @inheritParams bras_calibrate_logit
#' @return An object of class `bras_fit`; see `?bras_fit`. `converged` is
#'   computed from the margin error of the returned matrix. `feasible` is
#'   always `NA`, because the bounds are optional here and a successful
#'   return does not prove that they could be met. Check `converged`
#'   instead.
#' @examples
#' if (requireNamespace("survey", quietly = TRUE) &&
#'     requireNamespace("MASS", quietly = TRUE)) {
#'   x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
#'   fit <- bras_calibrate_raking(x0, c(20, 20), c(15, 25))
#'   fit$matrix
#'   # expected to match our own engine on the same problem:
#'   ras(x0, c(20, 20), c(15, 25))$matrix
#' }
#' @export
bras_calibrate_raking <- function(x0, row_targets, col_targets, lower = 0, upper = Inf,
                                   tol = 1e-6, ...) {
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package 'survey' is required for bras_calibrate_raking() but is not installed.",
         call. = FALSE)
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required by survey::calibrate() for bounded calibration ",
         "but is not installed. Install it with install.packages(\"MASS\").",
         call. = FALSE)
  }
  .validate_inputs(x0, row_targets, col_targets, lower, upper)
  m <- nrow(x0); n <- ncol(x0)
  lower_mat <- .as_bound_matrix(lower, m, n)
  upper_mat <- .as_bound_matrix(upper, m, n)
  if (any(lower_mat > upper_mat)) stop("`lower` must be <= `upper` in every cell.", call. = FALSE)
  x0 <- .clamp_seed_to_bounds(x0, lower_mat, upper_mat)
  m0 <- m; n0 <- n
  row_targets0 <- row_targets; col_targets0 <- col_targets
  drop <- .drop_zero_margins(x0, row_targets, col_targets, lower_mat, upper_mat)
  x0 <- drop$x0; row_targets <- drop$row_targets; col_targets <- drop$col_targets
  lower_mat <- drop$lower; upper_mat <- drop$upper
  m <- drop$m; n <- drop$n
  if (m <= 1L || n <= 1L) {
    return(.calibrate_trivial(drop, m0, n0, row_targets0, col_targets0,
                              lower_mat, upper_mat, tol, match.call()))
  }

  df <- .matrix_to_long(x0)
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)

  row_lv <- levels(df$row_id); col_lv <- levels(df$col_id)
  pop <- c(
    `(Intercept)` = sum(row_targets),
    stats::setNames(row_targets[-1], paste0("row_id", row_lv[-1])),
    stats::setNames(col_targets[-1], paste0("col_id", col_lv[-1]))
  )

  cal <- tryCatch(
    survey::calibrate(
      des, formula = ~row_id + col_id, population = pop,
      calfun = "raking",
      bounds = list(lower = as.vector(lower_mat), upper = as.vector(upper_mat)),
      bounds.const = TRUE,
      ...
    ),
    error = function(e) {
      stop("survey::calibrate(calfun = \"raking\") failed: ", conditionMessage(e), ". ",
           "See ?bras_calibrate_raking.", call. = FALSE)
    }
  )

  x_hat <- .reinsert_zero_margins(matrix(stats::weights(cal), nrow = m, ncol = n),
                                   drop, m0, n0)
  row_err <- abs(rowSums(x_hat) - row_targets0)
  col_err <- abs(colSums(x_hat) - col_targets0)

  .new_bras_fit(
    matrix = x_hat, converged = max(row_err, col_err) < tol, feasible = NA,
    n_outer_iter = 1L,
    fixed_cells = data.frame(row = integer(), col = integer(),
                              bound = character(), value = numeric(),
                              outer_pass = integer()),
    row_margin_error = row_err, col_margin_error = col_err,
    infeasible_at = NULL, method = "calibrate_raking (survey, wrapper)",
    call = match.call()
  )
}
