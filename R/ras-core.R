# Internal RAS engine, shared by every public function in the package.
# NOT exported. All public functions (ras(), bras_freeze_*(), bras_naive(),
# bras_soft()) call this to perform the alternating row/column scaling.
#

#' @keywords internal
#' @noRd
.ras_core <- function(x0,
                       row_targets,
                       col_targets,
                       free = matrix(TRUE, nrow(x0), ncol(x0)),
                       fixed_values = matrix(NA_real_, nrow(x0), ncol(x0)),
                       tol = 1e-8,
                       max_iter = 100,
                       stall_patience = 5L) {

  m <- nrow(x0); n <- ncol(x0)
  x <- x0
  x[!free] <- fixed_values[!free]

  converged <- FALSE
  stalled <- FALSE
  iter <- 0L
  prev_err <- Inf
  stall_count <- 0L

  for (iter in seq_len(max_iter)) {

    # --- Row step: scale free cells so each row sums to its target -------
    row_sum <- rowSums(x)
    free_row_mass <- rowSums(x * free)
    row_target_res <- row_targets - (row_sum - free_row_mass)  # target minus fixed mass
    r <- ifelse(free_row_mass > 0, row_target_res / free_row_mass, 1)
    scaled <- x * r
    x[free] <- scaled[free]

    # --- Column step: scale free cells so each column sums to its target -
    col_sum <- colSums(x)
    free_col_mass <- colSums(x * free)
    col_target_res <- col_targets - (col_sum - free_col_mass)
    s <- ifelse(free_col_mass > 0, col_target_res / free_col_mass, 1)
    scaled <- sweep(x, 2, s, `*`)
    x[free] <- scaled[free]

    # --- Convergence check on full margins (fixed + free) ----------------
    row_err <- abs(rowSums(x) - row_targets)
    col_err <- abs(colSums(x) - col_targets)
    cur_err <- max(row_err, col_err)

    if (cur_err < tol) {
      converged <- TRUE
      break
    }

    # --- Stall detection: a genuine fixed point of the alternating row/
    # column scaling that is NOT a solution. This happens when the free
    # cells left after fixation form a structurally degenerate subsystem
    # (row-wise pressure to scale up is exactly cancelled by column-wise
    # pressure to scale down, or vice versa) -- a real, known limitation of
    # alternating-projection RAS/IPF: the aggregate "sum(upper) >= residual
    # target" feasibility check is necessary but not sufficient (see
    # .check_feasibility()), and this is exactly the gap it cannot catch.
    # Unlike ordinary slow convergence, more iterations provably will not
    # help here -- the margin error stops changing entirely, iteration
    # after iteration. Detected as: negligible change in the worst-case
    # margin error for `stall_patience` consecutive iterations while still
    # above `tol`.
    if (abs(cur_err - prev_err) < 1e-9 * max(1, cur_err)) {
      stall_count <- stall_count + 1L
      if (stall_count >= stall_patience) {
        stalled <- TRUE
        break
      }
    } else {
      stall_count <- 0L
    }
    prev_err <- cur_err
  }

  list(
    matrix = x,
    converged = converged,
    stalled = stalled,
    n_iter = iter,
    row_margin_error = row_err,
    col_margin_error = col_err
  )
}

#' @keywords internal
#' @noRd
.check_feasibility <- function(x0, row_targets, col_targets, lower, upper, free, fixed_values,
                                eps = 1e-8) {

  fixed_row_sum <- rowSums(ifelse(!free, fixed_values, 0))
  fixed_col_sum <- colSums(ifelse(!free, fixed_values, 0))
  row_res <- row_targets - fixed_row_sum
  col_res <- col_targets - fixed_col_sum

  # A free cell seeded at exactly zero can never move away from zero under
  # multiplicative RAS scaling (x_ij <- x_ij * r_i * s_j; 0 times anything
  # is still 0), regardless of what `lower`/`upper` nominally allow. Its
  # true achievable range is {0}, not [lower, upper] -- reflect that here so
  # feasibility is not overstated for rows/columns that lean on such cells.
  achievable_lower <- ifelse(x0 == 0, 0, lower)
  achievable_upper <- ifelse(x0 == 0, 0, upper)

  row_lo <- rowSums(ifelse(free, achievable_lower, 0))
  row_hi <- rowSums(ifelse(free, achievable_upper, 0))
  col_lo <- colSums(ifelse(free, achievable_lower, 0))
  col_hi <- colSums(ifelse(free, achievable_upper, 0))

  bad_rows <- which(row_res < row_lo - eps | row_res > row_hi + eps)
  bad_cols <- which(col_res < col_lo - eps | col_res > col_hi + eps)

  list(
    feasible = length(bad_rows) == 0 && length(bad_cols) == 0,
    offending_rows = bad_rows,
    offending_cols = bad_cols
  )
}

#' @keywords internal
#' @noRd
.as_bound_matrix <- function(bound, m, n) {
  if (is.matrix(bound)) {
    if (!all(dim(bound) == c(m, n))) {
      stop("Bound matrix must have the same dimensions as `x0` (", m, " x ", n, ").",
           call. = FALSE)
    }
    return(bound)
  }
  if (length(bound) == 1L) {
    return(matrix(bound, nrow = m, ncol = n))
  }
  stop("`lower`/`upper` must be a single scalar or a matrix matching `x0`.", call. = FALSE)
}

#' @keywords internal
#' @noRd
.validate_inputs <- function(x0, row_targets, col_targets, lower, upper) {
  if (!is.matrix(x0) || !is.numeric(x0)) stop("`x0` must be a numeric matrix.", call. = FALSE)
  if (any(x0 < 0, na.rm = TRUE)) stop("`x0` must be non-negative.", call. = FALSE)
  m <- nrow(x0); n <- ncol(x0)
  if (length(row_targets) != m) stop("`row_targets` must have length nrow(x0).", call. = FALSE)
  if (length(col_targets) != n) stop("`col_targets` must have length ncol(x0).", call. = FALSE)
  if (abs(sum(row_targets) - sum(col_targets)) > 1e-6 * max(1, abs(sum(row_targets)))) {
    stop("`row_targets` and `col_targets` must have the same total.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.clamp_seed_to_bounds <- function(x0, lower, upper) {
  below <- x0 < lower
  above <- x0 > upper
  n_below <- sum(below); n_above <- sum(above)
  if (n_below == 0 && n_above == 0) return(x0)

  # x0 < lower (this includes x0 == 0 with lower > 0, otherwise unfixable:
  # 0 can never become positive under multiplicative RAS scaling, so
  # raising the seed up to `lower` is what makes the cell solvable at all,
  # not merely a speed optimization).
  x0[below] <- lower[below]
  # x0 > upper: clamping down does not change the final bound-respecting
  # result (the cell would be fixed there anyway once detected), just
  # saves the outer pass that would otherwise discover the violation.
  x0[above] <- upper[above]

  idx <- which(below | above, arr.ind = TRUE)
  preview <- utils::head(idx, 5)
  preview_str <- paste0("(", preview[, 1], ",", preview[, 2], ")", collapse = ", ")
  more <- if (nrow(idx) > 5) paste0(", and ", nrow(idx) - 5, " more") else ""
  warning(
    "`x0` had ", n_below + n_above, " cell(s) outside their own [lower, upper] ",
    "range before any adjustment (e.g. ", preview_str, more, "); imputed to the ",
    "nearest bound before solving, as if that value had been supplied as the seed. ",
    if (n_below > 0) paste0(
      "For cell(s) below `lower` (including any seeded at exactly 0 with a positive ",
      "`lower`), this substitution is what makes the cell solvable at all, not just ",
      "a speed optimization. "
    ) else "",
    "If this was unintentional -- a mismatch between `x0` and the bounds you meant to ",
    "set -- double-check your inputs; the final result may otherwise reflect a seed ",
    "value you did not intend.",
    call. = FALSE
  )
  x0
}
