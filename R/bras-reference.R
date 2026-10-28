#' Classic RAS (no bounds)
#'
#' Plain biproportional RAS with no cell-level bounds; provided as the
#' baseline/control condition for comparison against the bounded variants.
#'
#' @inheritParams bras_freeze_sim
#' @return An object of class `bras_fit` (with an empty `fixed_cells`
#'   table, since there are no bounds to violate); see `?bras_fit`.
#' @examples
#' x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
#' fit <- ras(x0, row_targets = c(20, 20), col_targets = c(15, 25))
#' fit$matrix
#' rowSums(fit$matrix)  # matches row_targets exactly
#' @export
ras <- function(x0, row_targets, col_targets, tol = 1e-8, max_iter = 100) {
  .validate_inputs(x0, row_targets, col_targets, 0, Inf)
  core_result <- .ras_core(x0, row_targets, col_targets, tol = tol, max_iter = max_iter)
  .new_bras_fit(
    matrix = core_result$matrix,
    converged = core_result$converged,
    feasible = TRUE,
    n_outer_iter = 1L,
    fixed_cells = data.frame(row = integer(), col = integer(),
                              bound = character(), value = numeric(),
                              outer_pass = integer()),
    row_margin_error = core_result$row_margin_error,
    col_margin_error = core_result$col_margin_error,
    infeasible_at = NULL,
    method = "ras_classic",
    call = match.call()
  )
}

#' Bounded RAS: naive recompute-every-pass variant (baseline, not recommended)
#'
#' Recomputes the full RAS on **every** cell each pass, including cells
#' that already sit on a bound. After every row and column scaling step it
#' sets each cell outside `[lower, upper]` to the nearest bound, and never
#' removes a cell from further adjustment. The next row or column
#' correction can push such a cell outside its bound again. So this method
#' can oscillate and fail to reach the target margins, even when a
#' solution that respects the bounds exists.
#'
#' It is included as a reference baseline for comparing against
#' [bras_freeze_sim()], [bras_freeze_seq()] and [bras_freeze_hybrid()]. For
#' actual data adjustment, prefer one of the `bras_freeze_*()` functions,
#' which guarantee exact margins and respected bounds whenever a solution
#' exists.
#'
#' @inheritParams bras_freeze_sim
#' @return An object of class `bras_fit`; see `?bras_fit`. `feasible` is
#'   always `NA` here, since this variant does not attempt to diagnose
#'   infeasibility.
#' @examples
#' x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
#' fit <- bras_naive(x0, c(20, 20), c(15, 25),
#'                    upper = matrix(c(Inf, Inf, 12, Inf), nrow = 2))
#' fit$converged  # may be FALSE even though a solution exists
#' @export
bras_naive <- function(x0, row_targets, col_targets, lower = 0, upper = Inf,
                        tol = 1e-8, max_iter = 200, verbose = FALSE) {
  .validate_inputs(x0, row_targets, col_targets, lower, upper)
  m <- nrow(x0); n <- ncol(x0)
  lower <- .as_bound_matrix(lower, m, n)
  upper <- .as_bound_matrix(upper, m, n)
  x0 <- .clamp_seed_to_bounds(x0, lower, upper)

  x <- x0
  converged <- FALSE
  fixed_cells <- data.frame(row = integer(), col = integer(),
                             bound = character(), value = numeric(),
                             outer_pass = integer())

  for (iter in seq_len(max_iter)) {
    r <- ifelse(rowSums(x) > 0, row_targets / rowSums(x), 1)
    x <- x * r
    s <- ifelse(colSums(x) > 0, col_targets / colSums(x), 1)
    x <- sweep(x, 2, s, `*`)
    x <- pmin(pmax(x, lower), upper)  # clamp every cell, every pass, no freezing

    row_err <- abs(rowSums(x) - row_targets)
    col_err <- abs(colSums(x) - col_targets)
    if (max(row_err, col_err) < tol) { converged <- TRUE; break }
  }

  bound_hit <- (abs(x - lower) < tol) | (abs(x - upper) < tol & is.finite(upper))
  if (any(bound_hit)) {
    idx <- which(bound_hit, arr.ind = TRUE)
    fixed_cells <- data.frame(
      row = idx[, 1], col = idx[, 2],
      bound = ifelse(abs(x[idx] - lower[idx]) < tol, "lower", "upper"),
      value = x[idx], outer_pass = NA_integer_
    )
  }

  .new_bras_fit(
    matrix = x, converged = converged, feasible = NA,
    n_outer_iter = iter, fixed_cells = fixed_cells,
    row_margin_error = row_err, col_margin_error = col_err,
    infeasible_at = NULL, method = "naive_reclamp",
    call = match.call()
  )
}

#' Bounded RAS: soft proportional redistribution variant
#'
#' Sets free cells outside their bounds to the nearest bound, and shares
#' the excess or deficit among the other free cells of the same row and
#' column, in proportion to their remaining headroom. It is similar in
#' spirit to the bounded trimming of `survey::calibrate(calfun =
#' "linear")`. No cell is fixed for good. It is a comparator for
#' [bras_freeze_sim()] that needs no extra packages; see
#' `vignette("RASbound-intro")` for the trade-offs.
#'
#' @inheritParams bras_freeze_sim
#' @return An object of class `bras_fit`; see `?bras_fit`. `feasible` is
#'   always `NA` here, since this variant does not attempt to diagnose
#'   infeasibility.
#' @examples
#' x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
#' fit <- bras_soft(x0, c(20, 20), c(15, 25),
#'                   upper = matrix(c(Inf, Inf, 12, Inf), nrow = 2))
#' fit$matrix
#' @export
bras_soft <- function(x0, row_targets, col_targets, lower = 0, upper = Inf,
                       tol = 1e-8, max_iter = 200, verbose = FALSE) {
  .validate_inputs(x0, row_targets, col_targets, lower, upper)
  m <- nrow(x0); n <- ncol(x0)
  lower <- .as_bound_matrix(lower, m, n)
  upper <- .as_bound_matrix(upper, m, n)
  x0 <- .clamp_seed_to_bounds(x0, lower, upper)

  x <- x0
  converged <- FALSE

  # `room` (how much a cell could still absorb before hitting `upper`) is used
  # only to WEIGHT the proportional redistribution below, never to clamp the
  # final values. When upper == Inf, room is Inf, and an Inf weight times a
  # zero share (0 * Inf) is NaN in R, not 0 -- so for this weighting step only
  # we cap Inf entries at a large-but-finite sentinel, scaled to the problem,
  # to avoid that. The actual clamping a few lines below still uses the real
  # `upper` (with true Inf), so unbounded cells remain truly unbounded.
  scale_hint <- max(1, sum(abs(row_targets)), sum(abs(col_targets)))
  finite_cap <- 1e9 * scale_hint
  upper_room <- upper
  upper_room[!is.finite(upper_room)] <- finite_cap

  for (iter in seq_len(max_iter)) {
    # row step + soft redistribution
    r <- ifelse(rowSums(x) > 0, row_targets / rowSums(x), 1)
    x <- x * r
    over  <- pmax(x - upper, 0)
    under <- pmax(lower - x, 0)
    x <- pmin(pmax(x, lower), upper)
    excess <- rowSums(over) - rowSums(under)  # net mass to push back into the row
    room <- (upper_room - x)  # how much each free-ish cell could still absorb
    room[room < 0] <- 0
    room_row_sum <- rowSums(room)
    add <- ifelse(room_row_sum > 0, excess / room_row_sum, 0) * room
    x <- x + add
    x <- pmin(pmax(x, lower), upper)

    # column step + soft redistribution (mirrors the row step)
    s <- ifelse(colSums(x) > 0, col_targets / colSums(x), 1)
    x <- sweep(x, 2, s, `*`)
    over  <- pmax(x - upper, 0)
    under <- pmax(lower - x, 0)
    x <- pmin(pmax(x, lower), upper)
    excess <- colSums(over) - colSums(under)
    room <- (upper_room - x); room[room < 0] <- 0
    room_col_sum <- colSums(room)
    add <- sweep(room, 2, ifelse(room_col_sum > 0, excess / room_col_sum, 0), `*`)
    x <- x + add
    x <- pmin(pmax(x, lower), upper)

    row_err <- abs(rowSums(x) - row_targets)
    col_err <- abs(colSums(x) - col_targets)
    if (max(row_err, col_err) < tol) { converged <- TRUE; break }
  }

  .new_bras_fit(
    matrix = x, converged = converged, feasible = NA,
    n_outer_iter = iter,
    fixed_cells = data.frame(row = integer(), col = integer(),
                              bound = character(), value = numeric(),
                              outer_pass = integer()),
    row_margin_error = row_err, col_margin_error = col_err,
    infeasible_at = NULL, method = "soft_redistribute",
    call = match.call()
  )
}
