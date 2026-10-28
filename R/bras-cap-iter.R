#' Bounded RAS: iterative capping with renormalization
#'
#' After each row and column scaling step, cells outside `[lower, upper]`
#' are set to the nearest bound ("capped"). Then only the cells of that
#' row or column that are not capped are rescaled, all by the same factor,
#' to make up for the mass lost or gained. This differs from [bras_naive()],
#' which leaves the margin off, and from [bras_soft()], which shares the
#' difference in proportion to the remaining headroom. No cell is removed
#' for good: a cell capped in one pass can move again in the next.
#'
#' This is RASbound's own implementation, without extra dependencies, of
#' the general "cap and iterate" idea used by some external weight-raking
#' tools. It is not a wrapper of, nor a copy of, any specific package. See
#' `vignette("RASbound-intro")` for how it compares to [bras_naive()] and
#' [bras_soft()].
#'
#' @inheritParams bras_freeze_sim
#' @return An object of class `bras_fit`; see `?bras_fit`. `feasible` is
#'   always `NA` here, since this variant does not attempt to diagnose
#'   infeasibility.
#' @examples
#' x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
#' fit <- bras_cap_iter(x0, c(20, 20), c(15, 25),
#'                       upper = matrix(c(Inf, Inf, 12, Inf), nrow = 2))
#' fit$matrix
#' @export
bras_cap_iter <- function(x0, row_targets, col_targets, lower = 0, upper = Inf,
                           tol = 1e-8, max_iter = 200, verbose = FALSE) {
  .validate_inputs(x0, row_targets, col_targets, lower, upper)
  m <- nrow(x0); n <- ncol(x0)
  lower <- .as_bound_matrix(lower, m, n)
  upper <- .as_bound_matrix(upper, m, n)
  x0 <- .clamp_seed_to_bounds(x0, lower, upper)

  x <- x0
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    # --- row step -----------------------------------------------------
    r <- ifelse(rowSums(x) > 0, row_targets / rowSums(x), 1)
    x <- x * r
    capped <- (x < lower) | (x > upper)
    x <- pmin(pmax(x, lower), upper)

    free_mass <- rowSums(x * !capped)
    deficit <- row_targets - rowSums(x)  # still missing/excess after capping
    corr <- ifelse(free_mass > 0, pmax(free_mass + deficit, 0) / free_mass, 1)
    scaled <- x * corr
    x[!capped] <- scaled[!capped]

    # --- column step (mirrors the row step) ----------------------------
    s <- ifelse(colSums(x) > 0, col_targets / colSums(x), 1)
    x <- sweep(x, 2, s, `*`)
    capped <- (x < lower) | (x > upper)
    x <- pmin(pmax(x, lower), upper)

    free_mass <- colSums(x * !capped)
    deficit <- col_targets - colSums(x)
    corr <- ifelse(free_mass > 0, pmax(free_mass + deficit, 0) / free_mass, 1)
    scaled <- sweep(x, 2, corr, `*`)
    x[!capped] <- scaled[!capped]

    row_err <- abs(rowSums(x) - row_targets)
    col_err <- abs(colSums(x) - col_targets)
    if (max(row_err, col_err) < tol) { converged <- TRUE; break }
  }

  bound_hit <- (abs(x - lower) < tol) | (abs(x - upper) < tol & is.finite(upper))
  fixed_cells <- data.frame(row = integer(), col = integer(),
                             bound = character(), value = numeric(),
                             outer_pass = integer())
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
    infeasible_at = NULL, method = "cap_iter_renorm",
    call = match.call()
  )
}
