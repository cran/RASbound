# Shared outer-loop driver for the three hard-freeze variants: solve the
# free subproblem to convergence, THEN check bound violations among free
# cells, fix the offending ones, and repeat.

#' @keywords internal
#' @noRd
.bras_freeze_driver <- function(x0, row_targets, col_targets, lower, upper,
                                 tol, max_iter, max_outer, fixation,
                                 on_infeasible, verbose, call) {

  .validate_inputs(x0, row_targets, col_targets, lower, upper)
  m <- nrow(x0); n <- ncol(x0)
  lower <- .as_bound_matrix(lower, m, n)
  upper <- .as_bound_matrix(upper, m, n)
  if (any(lower > upper)) stop("`lower` must be <= `upper` in every cell.", call. = FALSE)
  x0 <- .clamp_seed_to_bounds(x0, lower, upper)

  free <- matrix(TRUE, m, n)
  fixed_values <- matrix(NA_real_, m, n)
  fixed_cells <- data.frame(row = integer(), col = integer(),
                             bound = character(), value = numeric(),
                             outer_pass = integer())

  core_result <- NULL
  feasible <- TRUE
  infeasible_at <- NULL
  outer_iter <- 0L
  exhausted <- FALSE
  stalled <- FALSE

  repeat {
    outer_iter <- outer_iter + 1L

    fc <- .check_feasibility(x0, row_targets, col_targets, lower, upper, free, fixed_values)
    if (!fc$feasible) {
      feasible <- FALSE
      infeasible_at <- list(rows = fc$offending_rows, cols = fc$offending_cols)
      break
    }

    # Recompute with the most recent free/fixed state BEFORE deciding whether
    # there is still budget to fix more cells, so the matrix and margin
    # errors we might return below always reflect the latest fixation.
    core_result <- .ras_core(x0, row_targets, col_targets, free, fixed_values,
                              tol = tol, max_iter = max_iter)
    x <- core_result$matrix

    viol_low  <- free & (x < lower - tol)
    viol_high <- free & (x > upper + tol)
    viol <- viol_low | viol_high

    if (!any(viol)) {
      if (core_result$converged) {
        if (verbose) message("Converged in ", outer_iter, " outer pass(es).")
      } else if (isTRUE(core_result$stalled)) {
        stalled <- TRUE
        if (verbose) {
          message("Inner solve stalled at a genuine fixed point after ", core_result$n_iter,
                   " iterations (margin error stopped changing) at outer pass ", outer_iter,
                   ". No bound violations remain, but the margins are not matched -- ",
                   "increasing max_iter will not help. This means the free cells left ",
                   "after fixation cannot jointly satisfy all margins, even though each ",
                   "margin's own aggregate capacity looked sufficient in isolation (the ",
                   "feasibility check is necessary but not sufficient). See ",
                   "vignette(\"RASbound-diagnostics\").")
        }
      } else if (verbose) {
        message("Inner solve did not converge within max_iter (", max_iter,
                 ") at outer pass ", outer_iter, "; consider increasing max_iter.")
      }
      break
    }

    if (outer_iter >= max_outer) {
      exhausted <- TRUE
      if (verbose) message("Reached max_outer (", max_outer,
                            ") with residual bound violations still unresolved.")
      break
    }

    idx <- which(viol, arr.ind = TRUE)

    to_fix <- if (fixation == "sequential") {
      dist <- ifelse(viol_low[idx], lower[idx] - x[idx], x[idx] - upper[idx])
      idx[which.max(dist), , drop = FALSE]
    } else {
      idx  # "simultaneous" and first attempt of "hybrid"
    }

    applied <- .apply_fixation(to_fix, viol_low, lower, upper,
                                fixed_values, free, outer_iter, fixed_cells)
    free <- applied$free
    fixed_values <- applied$fixed_values
    fixed_cells <- applied$fixed_cells

    if (fixation == "hybrid" && nrow(to_fix) > 1) {
      fc2 <- .check_feasibility(x0, row_targets, col_targets, lower, upper, free, fixed_values)
      if (!fc2$feasible) {
        if (verbose) message("Simultaneous fixation infeasible at pass ", outer_iter,
                              "; retrying that pass sequentially.")
        # undo this pass: go back to the free/fixed state before it
        this_pass <- fixed_cells$outer_pass == outer_iter
        undone <- fixed_cells[this_pass, ]
        for (r in seq_len(nrow(undone))) {
          free[undone$row[r], undone$col[r]] <- TRUE
          fixed_values[undone$row[r], undone$col[r]] <- NA_real_
        }
        fixed_cells <- fixed_cells[!this_pass, ]

        dist <- ifelse(viol_low[idx], lower[idx] - x[idx], x[idx] - upper[idx])
        one <- idx[which.max(dist), , drop = FALSE]
        applied <- .apply_fixation(one, viol_low, lower, upper,
                                    fixed_values, free, outer_iter, fixed_cells)
        free <- applied$free
        fixed_values <- applied$fixed_values
        fixed_cells <- applied$fixed_cells
      }
    }
  }

  .new_bras_fit(
    matrix = if (is.null(core_result)) x0 else core_result$matrix,
    converged = if (is.null(core_result)) FALSE else (core_result$converged && feasible && !exhausted),
    feasible = feasible,
    n_outer_iter = outer_iter,
    fixed_cells = fixed_cells,
    row_margin_error = if (is.null(core_result)) rep(NA_real_, m) else core_result$row_margin_error,
    col_margin_error = if (is.null(core_result)) rep(NA_real_, n) else core_result$col_margin_error,
    infeasible_at = infeasible_at,
    stalled = stalled,
    method = paste0("freeze_", fixation),
    call = call
  )
}

#' @keywords internal
#' @noRd
.apply_fixation <- function(to_fix, viol_low, lower, upper, fixed_values, free,
                             outer_iter, fixed_cells) {
  for (k in seq_len(nrow(to_fix))) {
    i <- to_fix[k, 1]; j <- to_fix[k, 2]
    is_low <- viol_low[i, j]
    bound_val <- if (is_low) lower[i, j] else upper[i, j]
    fixed_values[i, j] <- bound_val
    free[i, j] <- FALSE
    fixed_cells <- rbind(fixed_cells, data.frame(
      row = i, col = j,
      bound = if (is_low) "lower" else "upper",
      value = bound_val, outer_pass = outer_iter
    ))
  }
  list(free = free, fixed_values = fixed_values, fixed_cells = fixed_cells)
}

#' Bounded RAS: simultaneous hard-freeze variant
#'
#' Solves the free subproblem to convergence, fixes **all** free cells that
#' violate a bound in that pass, subtracts them from the residual margins,
#' and repeats. See `vignette("RASbound-intro")` for details and trade-offs
#' relative to [bras_freeze_seq()] and [bras_freeze_hybrid()].
#'
#' @param x0 Numeric matrix, the base/seed matrix (non-negative).
#' @param row_targets,col_targets Numeric vectors of target row/column totals.
#' @param lower,upper Scalar or matrix of per-cell lower/upper bounds.
#'   Defaults to 0 and `Inf`.
#' @param tol Convergence tolerance on margins for the inner RAS solve.
#' @param max_iter Maximum inner RAS iterations per outer pass.
#' @param max_outer Maximum number of outer fixation passes.
#' @param on_infeasible One of `"error"`, `"warning"`, `"force"`: what to do
#'   if no bound-respecting solution exists.
#' @param verbose Print progress messages.
#' @return An object of class `bras_fit`; see `?bras_fit` for the fields.
#' @examples
#' x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
#' row_targets <- c(20, 20)
#' col_targets <- c(15, 25)
#' upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)  # cap cell [1, 2] at 12
#'
#' fit <- bras_freeze_sim(x0, row_targets, col_targets, upper = upper)
#' fit$matrix
#' fit$fixed_cells
#' summary(fit)
#'
#' # An infeasible configuration is reported, not silently misapplied:
#' too_tight <- matrix(1, nrow = 2, ncol = 2)
#' tryCatch(
#'   bras_freeze_sim(matrix(1, 2, 2), c(100, 2), c(51, 51), upper = too_tight),
#'   error = function(e) conditionMessage(e)
#' )
#' @export
bras_freeze_sim <- function(x0, row_targets, col_targets, lower = 0, upper = Inf,
                             tol = 1e-8, max_iter = 100, max_outer = 50,
                             on_infeasible = c("error", "warning", "force"),
                             verbose = FALSE) {
  on_infeasible <- match.arg(on_infeasible)
  fit <- .bras_freeze_driver(x0, row_targets, col_targets, lower, upper,
                              tol, max_iter, max_outer, fixation = "simultaneous",
                              on_infeasible = on_infeasible, verbose = verbose,
                              call = match.call())
  .handle_infeasible(fit, on_infeasible)
}

#' Bounded RAS: sequential hard-freeze variant (most-extreme cell first)
#'
#' Same two-stage idea as [bras_freeze_sim()] (solve the free subproblem,
#' then fix violators, then repeat), but fixes only the single most
#' extreme violating cell per outer pass instead of all of them at once.
#' Slower, but never over-fixes: prefer this over [bras_freeze_sim()] when
#' you suspect simultaneous fixation is reporting a false infeasibility, or
#' when you want the most conservative, auditable fixation trail.
#'
#' @inheritParams bras_freeze_sim
#' @return An object of class `bras_fit`; see `?bras_fit` for the fields.
#' @examples
#' x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
#' fit <- bras_freeze_seq(x0, c(20, 20), c(15, 25),
#'                         upper = matrix(c(Inf, Inf, 12, Inf), nrow = 2))
#' fit$fixed_cells
#' @export
bras_freeze_seq <- function(x0, row_targets, col_targets, lower = 0, upper = Inf,
                             tol = 1e-8, max_iter = 100, max_outer = 50,
                             on_infeasible = c("error", "warning", "force"),
                             verbose = FALSE) {
  on_infeasible <- match.arg(on_infeasible)
  fit <- .bras_freeze_driver(x0, row_targets, col_targets, lower, upper,
                              tol, max_iter, max_outer, fixation = "sequential",
                              on_infeasible = on_infeasible, verbose = verbose,
                              call = match.call())
  .handle_infeasible(fit, on_infeasible)
}

#' Bounded RAS: hybrid hard-freeze variant (simultaneous with sequential fallback)
#'
#' Tries simultaneous fixation each pass; if that makes the residual
#' problem locally infeasible *within that one pass*, that pass is
#' retried fixing only the most extreme violating cell before continuing.
#' That per-pass safety net is not always enough: fixing slightly more than
#' necessary in several *different* passes can add up to an infeasible or
#' non-convergent run, even though no single pass looked infeasible (see
#' `vignette("RASbound-diagnostics")`). So if the per-pass approach ends
#' without a solution, the **entire problem is retried once, from
#' scratch, with [bras_freeze_seq()]'s one-cell-at-a-time fixation**
#' before giving up. This costs nothing in the common case, because it
#' only happens after the faster approach has failed. The result is as
#' robust as [bras_freeze_seq()] and as fast as [bras_freeze_sim()]
#' whenever that is safe. This is the recommended default when you are
#' unsure whether to pick [bras_freeze_sim()] or [bras_freeze_seq()].
#'
#' @inheritParams bras_freeze_sim
#' @return An object of class `bras_fit`; see `?bras_fit` for the fields.
#'   `method` is `"freeze_hybrid"` normally, or
#'   `"freeze_hybrid_retry_sequential"` if the global retry was the one
#'   that found the solution.
#' @examples
#' x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
#' fit <- bras_freeze_hybrid(x0, c(20, 20), c(15, 25),
#'                            upper = matrix(c(Inf, Inf, 12, Inf), nrow = 2))
#' fit$matrix
#' @export
bras_freeze_hybrid <- function(x0, row_targets, col_targets, lower = 0, upper = Inf,
                                tol = 1e-8, max_iter = 100, max_outer = 50,
                                on_infeasible = c("error", "warning", "force"),
                                verbose = FALSE) {
  on_infeasible <- match.arg(on_infeasible)
  fit <- .bras_freeze_driver(x0, row_targets, col_targets, lower, upper,
                              tol, max_iter, max_outer, fixation = "hybrid",
                              on_infeasible = on_infeasible, verbose = verbose,
                              call = match.call())

  if (!(isTRUE(fit$feasible) && isTRUE(fit$converged))) {
    if (verbose) {
      message("Per-pass hybrid fixation did not find a solution; retrying the ",
              "whole problem with pure sequential fixation before giving up.")
    }
    fit_seq <- .bras_freeze_driver(x0, row_targets, col_targets, lower, upper,
                                    tol, max_iter, max_outer, fixation = "sequential",
                                    on_infeasible = on_infeasible, verbose = verbose,
                                    call = match.call())
    if (isTRUE(fit_seq$feasible) && isTRUE(fit_seq$converged)) {
      fit_seq$method <- "freeze_hybrid_retry_sequential"
      fit <- fit_seq
    }
    # If the sequential retry ALSO fails, keep the original hybrid attempt's
    # diagnostics (they were the first to establish there is a real problem,
    # for whatever it is worth to the person reading them).
  }

  .handle_infeasible(fit, on_infeasible)
}

#' @keywords internal
#' @noRd
.handle_infeasible <- function(fit, on_infeasible) {
  if (!fit$feasible) {
    msg <- paste0("No bound-respecting solution found (infeasible at rows: ",
                   paste(fit$diagnostics$infeasible_at$rows, collapse = ", "),
                   "; cols: ", paste(fit$diagnostics$infeasible_at$cols, collapse = ", "), ").")
    if (on_infeasible == "error") stop(msg, call. = FALSE)
    if (on_infeasible == "warning") warning(msg, call. = FALSE)
  } else if (isTRUE(fit$diagnostics$stalled)) {
    # Distinct from both infeasibility and max_outer exhaustion: no bound
    # violations remain, but the inner solve reached a genuine fixed point
    # that does not satisfy the margins, confirmed by the margin error no
    # longer changing across iterations. Increasing max_iter/max_outer will
    # NOT help here -- say so explicitly, unlike the exhaustion message
    # below, so this isn't mistaken for "just needs more budget."
    warning("Stopped after ", fit$n_outer_iter, " outer pass(es): the remaining free ",
            "cells reached a fixed point that does not satisfy the margins, with no ",
            "bound violations left to fix. Increasing `max_iter` or `max_outer` will ",
            "NOT resolve this -- it means the free cells left after fixation cannot ",
            "jointly satisfy all margins, even though each margin's own aggregate ",
            "capacity looked sufficient in isolation. See ",
            "vignette(\"RASbound-diagnostics\") for how to investigate which rows/",
            "columns are involved.",
            call. = FALSE)
  } else if (!fit$converged) {
    # Distinct from infeasibility: a solution may still exist, but
    # `max_outer` ran out before it was found. Always surfaced, regardless
    # of `on_infeasible`, since it is not the same failure mode.
    warning("Stopped after ", fit$n_outer_iter, " outer pass(es) without fully ",
            "resolving bound violations. The returned matrix may still contain ",
            "cells outside their bounds; consider increasing `max_outer`.",
            call. = FALSE)
  }
  fit
}
