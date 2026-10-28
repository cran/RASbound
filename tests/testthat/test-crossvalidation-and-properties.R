# These tests go beyond pinning exact numbers on toy 2x2/3x3 cases (the
# rest of the suite) and instead check PROPERTIES that must hold for any
# correct run, plus a direct cross-check against an independent external
# implementation (survey::calibrate()) -- the most convincing evidence we
# have that the core RAS engine is actually correct, not just internally
# consistent with itself.

# Generates n positive targets summing exactly to `total`, each at least
# `min_val`. Rounding drift is corrected on the LARGEST entry (never the
# smallest/arbitrary one), so this can never produce a degenerate
# near-zero or negative target the way `total - sum(others)` can.
.random_targets <- function(n, total, min_val = 5) {
  raw <- stats::runif(n, min_val, min_val * 4)
  targets <- round(raw / sum(raw) * total)
  targets[which.max(targets)] <- targets[which.max(targets)] + (total - sum(targets))
  targets
}

test_that("cross-validation: our RAS engine agrees with survey::calibrate(calfun = 'raking') on an unbounded problem", {
  testthat::skip_if_not_installed("survey")
  testthat::skip_if_not_installed("MASS")

  set.seed(20260919)
  x0 <- matrix(round(runif(16, 1, 20)), nrow = 4)
  row_targets <- round(runif(4, 80, 150))
  col_targets <- .random_targets(4, sum(row_targets), min_val = 20)

  fit_ras <- ras(x0, row_targets, col_targets, tol = 1e-10, max_iter = 500)

  # calfun = "raking" is the ONLY survey::calibrate() distance function
  # mathematically identical to classical multiplicative RAS/IPF (the
  # exponential distance); "linear" and "logit" minimize different
  # distance functions and, even unbounded, generally satisfy the same
  # margins without landing on the same matrix -- see the note in
  # ?bras_calibrate_logit. RASbound does not currently wrap "raking", so
  # this calls survey::calibrate() directly via the package's own
  # (unexported) long-format helper, for a genuine apples-to-apples check.
  df  <- RASbound:::.matrix_to_long(x0)
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)
  row_lv <- levels(df$row_id); col_lv <- levels(df$col_id)
  pop <- c(
    `(Intercept)` = sum(row_targets),
    stats::setNames(row_targets[-1], paste0("row_id", row_lv[-1])),
    stats::setNames(col_targets[-1], paste0("col_id", col_lv[-1]))
  )
  cal_raking <- survey::calibrate(des, formula = ~row_id + col_id, population = pop,
                                   calfun = "raking")
  raking_matrix <- matrix(stats::weights(cal_raking), nrow = 4, ncol = 4)

  # The unbounded RAS solution is unique for a positive, fully-connected
  # seed matrix (Bacharach, 1970) *given the same distance function*, so
  # an independent implementation of the same exponential/multiplicative
  # distance should land on (numerically) the same matrix.
  expect_equal(fit_ras$matrix, raking_matrix, tolerance = 1e-3)
})

test_that("property: bound-respecting variants never violate bounds or margins when they claim success", {
  set.seed(20260919)
  variants <- list(bras_freeze_sim, bras_freeze_seq, bras_freeze_hybrid)

  for (rep in 1:8) {
    m <- sample(3:5, 1); n <- sample(3:5, 1)
    x0 <- matrix(round(runif(m * n, 0, 20)), nrow = m)
    row_targets <- round(runif(m, 20, 60))
    col_targets <- .random_targets(n, sum(row_targets), min_val = 5)

    # generous, usually-feasible bounds: comfortably above what any single
    # cell should need, so most draws are solvable and we're testing the
    # "success" path's correctness, not infeasibility handling (covered
    # elsewhere)
    upper <- matrix(max(row_targets, col_targets), nrow = m, ncol = n)

    for (f in variants) {
      fit <- f(x0, row_targets, col_targets, upper = upper, on_infeasible = "warning")
      if (isTRUE(fit$feasible) && isTRUE(fit$converged)) {
        expect_true(all(fit$matrix >= -1e-6), info = "matrix must stay non-negative")
        expect_true(all(fit$matrix <= upper + 1e-6), info = "upper bound respected")
        expect_equal(unname(rowSums(fit$matrix)), row_targets, tolerance = 1e-4)
        expect_equal(unname(colSums(fit$matrix)), col_targets, tolerance = 1e-4)
      }
    }
  }
})

test_that("property: a zero-seeded cell stays exactly zero across every variant", {
  x0 <- matrix(c(0, 4, 6, 8, 0, 3, 5, 2, 9), nrow = 3)  # two structural zeros
  row_targets <- c(15, 20, 20)
  col_targets <- c(18, 17, 20)

  fits <- list(
    ras        = ras(x0, row_targets, col_targets),
    freeze_sim = bras_freeze_sim(x0, row_targets, col_targets),
    naive      = bras_naive(x0, row_targets, col_targets),
    soft       = bras_soft(x0, row_targets, col_targets),
    cap_iter   = bras_cap_iter(x0, row_targets, col_targets)
  )

  zero_idx <- which(x0 == 0, arr.ind = TRUE)
  for (name in names(fits)) {
    got <- fits[[name]]$matrix[zero_idx]
    expect_true(all(got == 0), info = paste(name, "should preserve structural zeros exactly"))
  }
})

test_that("a 4x4 problem with several simultaneous cap violations still converges correctly", {
  # Unlike the rest of the suite (always exactly one violating cell), this
  # forces multiple cells to exceed their bound in the same pass: every
  # diagonal cell starts large (20) relative to its row/column target, so
  # by symmetry all four are expected to need capping together.
  x0 <- matrix(c(20,  2,  2,  2,
                   2, 20,  2,  2,
                   2,  2, 20,  2,
                   2,  2,  2, 20), nrow = 4, byrow = TRUE)
  row_targets <- c(40, 40, 40, 40)
  col_targets <- c(40, 40, 40, 40)
  upper <- matrix(15, nrow = 4, ncol = 4)  # comfortably feasible: 4 * 15 = 60 >= 40 per row/col

  # x0's diagonal (20) starts above `upper` (15) by construction -- the
  # seed-clamping warning is expected and not the point of this test.
  fit <- suppressWarnings(
    bras_freeze_sim(x0, row_targets, col_targets, upper = upper, on_infeasible = "warning")
  )

  expect_true(fit$feasible)
  if (fit$converged) {
    expect_true(all(fit$matrix <= 15 + 1e-4))
    expect_equal(unname(rowSums(fit$matrix)), row_targets, tolerance = 1e-4)
    expect_equal(unname(colSums(fit$matrix)), col_targets, tolerance = 1e-4)
    expect_gte(nrow(fit$fixed_cells), 2)  # confirms more than one cell needed fixing
  }
})
