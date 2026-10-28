test_that("bras_calibrate_logit accepts Inf bounds via internal finite substitution", {
  testthat::skip_if_not_installed("survey")
  testthat::skip_if_not_installed("MASS")

  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  row_targets <- c(20, 20)
  col_targets <- c(15, 25)
  upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)

  fit <- bras_calibrate_logit(x0, row_targets, col_targets, upper = upper)

  expect_s3_class(fit, "bras_fit")
  expect_equal(unname(rowSums(fit$matrix)), row_targets, tolerance = 1e-4)
  expect_equal(unname(colSums(fit$matrix)), col_targets, tolerance = 1e-4)
  expect_true(fit$matrix[1, 2] <= 12 + 1e-4)
})

test_that("bras_calibrate_linear accepts Inf bounds natively", {
  testthat::skip_if_not_installed("survey")
  testthat::skip_if_not_installed("MASS")

  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  row_targets <- c(20, 20)
  col_targets <- c(15, 25)
  upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)

  fit <- bras_calibrate_linear(x0, row_targets, col_targets, upper = upper)

  expect_s3_class(fit, "bras_fit")
})

test_that("bras_calibrate_raking accepts Inf bounds natively and matches ras() when unbounded", {
  testthat::skip_if_not_installed("survey")
  testthat::skip_if_not_installed("MASS")

  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  row_targets <- c(20, 20)
  col_targets <- c(15, 25)

  fit_raking <- bras_calibrate_raking(x0, row_targets, col_targets)
  fit_ras    <- ras(x0, row_targets, col_targets, tol = 1e-10)

  expect_s3_class(fit_raking, "bras_fit")
  expect_true(fit_raking$converged)
  # the whole point of this wrapper: it should match our own engine
  # exactly (up to numerical tolerance) on an unbounded problem
  expect_equal(fit_raking$matrix, fit_ras$matrix, tolerance = 1e-3)

  # and it accepts Inf bounds natively too, like linear
  upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)
  fit_bounded <- bras_calibrate_raking(x0, row_targets, col_targets, upper = upper)
  expect_s3_class(fit_bounded, "bras_fit")
})

test_that("bras_calibrate_* reject lower > upper before calling survey", {
  testthat::skip_if_not_installed("survey")
  testthat::skip_if_not_installed("MASS")

  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  bad <- matrix(c(10, 0, 0, 0), nrow = 2)  # lower[1,1]=10 > upper[1,1]=0

  expect_error(
    bras_calibrate_logit(x0, c(20, 20), c(15, 25), lower = bad, upper = matrix(0, 2, 2)),
    "lower.*upper"
  )
  expect_error(
    bras_calibrate_linear(x0, c(20, 20), c(15, 25), lower = bad, upper = matrix(0, 2, 2)),
    "lower.*upper"
  )
  expect_error(
    bras_calibrate_raking(x0, c(20, 20), c(15, 25), lower = bad, upper = matrix(0, 2, 2)),
    "lower.*upper"
  )
})

test_that("bras_calibrate_logit's tryCatch wrapping is exercised by a real divergence", {
  testthat::skip_if_not_installed("survey")
  testthat::skip_if_not_installed("MASS")

  # Found while building the diagnostics vignette (not synthesized for the
  # test): a cell already clamped to a real, tight bound (upper[1,3] = 6,
  # with x0[1,3] = 9 starting above it) sitting next to cells substituted
  # to a much larger bound (unbounded_mult * max(targets) = 500) makes
  # survey's Newton-Raphson solve diverge outright, not just fail to
  # converge in time. This fills a previously-documented coverage gap:
  # the error-handling branch of bras_calibrate_logit()'s internal
  # tryCatch() had never been exercised by any test before this.
  x0 <- matrix(c(2, 4, 9,
                  6, 8, 1,
                  3, 5, 7), nrow = 3, byrow = TRUE)
  row_targets <- c(25, 20, 25)
  col_targets <- c(20, 25, 25)
  upper <- matrix(c(Inf, Inf, 6,
                     Inf, Inf, Inf,
                     Inf, Inf, Inf), nrow = 3, byrow = TRUE)

  expect_error(
    bras_calibrate_logit(x0, row_targets, col_targets, upper = upper),
    'survey::calibrate\\(calfun = "logit"\\) failed'
  )
})

test_that("bras_calibrate_logit/linear/raking compute converged from actual margin error", {
  testthat::skip_if_not_installed("survey")
  testthat::skip_if_not_installed("MASS")

  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  row_targets <- c(20, 20)
  col_targets <- c(15, 25)

  fit_logit  <- bras_calibrate_logit(x0, row_targets, col_targets)
  fit_linear <- bras_calibrate_linear(x0, row_targets, col_targets)
  fit_raking <- bras_calibrate_raking(x0, row_targets, col_targets)

  # A well-posed, unbounded problem should converge to the true margins.
  expect_true(fit_logit$converged)
  expect_true(fit_linear$converged)
  expect_true(fit_raking$converged)
  expect_equal(unname(rowSums(fit_logit$matrix)), row_targets, tolerance = 1e-4)
  expect_equal(unname(rowSums(fit_linear$matrix)), row_targets, tolerance = 1e-4)
  expect_equal(unname(rowSums(fit_raking$matrix)), row_targets, tolerance = 1e-4)

  # converged is not hard-coded: an unreasonably strict tol should flip it,
  # since floating-point margin error is never exactly 0.
  fit_strict <- bras_calibrate_logit(x0, row_targets, col_targets, tol = 0)
  expect_false(fit_strict$converged)
})

test_that("bras_calibrate_linear/raking's tryCatch wrapping is exercised by genuine infeasibility", {
  # Confirmed: survey::calibrate() errors by default (force = FALSE)
  # rather than returning an unconverged result when it cannot reach the
  # requested margins at all, for calfun = "linear" and "raking" too --
  # the same behaviour already confirmed for calfun = "logit" above.
  # Reuses the "too_tight" infeasibility pattern used to test
  # bras_freeze_sim()'s feasibility detection.
  testthat::skip_if_not_installed("survey")
  testthat::skip_if_not_installed("MASS")

  x0 <- matrix(1, nrow = 2, ncol = 2)
  too_tight <- matrix(1, nrow = 2, ncol = 2)  # each row can reach at most 2, but needs 100

  expect_error(
    bras_calibrate_linear(x0, c(100, 2), c(51, 51), upper = too_tight),
    'survey::calibrate\\(calfun = "linear"\\) failed', fixed = FALSE
  )
  expect_error(
    bras_calibrate_raking(x0, c(100, 2), c(51, 51), upper = too_tight),
    'survey::calibrate\\(calfun = "raking"\\) failed', fixed = FALSE
  )
})

test_that("zero-target rows/columns are dropped internally and reinserted as zeros", {
  testthat::skip_if_not_installed("survey")
  testthat::skip_if_not_installed("MASS")

  x0 <- matrix(c(2, 4, 6, 3, 5, 7, 8, 1, 9), nrow = 3)
  row_targets <- c(20, 0, 25)
  col_targets <- c(15, 0, 30)

  for (f in list(bras_calibrate_raking, bras_calibrate_linear)) {
    fit <- f(x0, row_targets, col_targets)
    expect_equal(dim(fit$matrix), c(3L, 3L))
    expect_true(all(fit$matrix[2, ] == 0))
    expect_true(all(fit$matrix[, 2] == 0))
    expect_equal(unname(rowSums(fit$matrix)), row_targets, tolerance = 1e-4)
    expect_equal(unname(colSums(fit$matrix)), col_targets, tolerance = 1e-4)
    expect_equal(length(fit$diagnostics$row_margin_error), 3L)
  }
})

test_that("zero-target row with a positive lower bound errors clearly", {
  testthat::skip_if_not_installed("survey")
  testthat::skip_if_not_installed("MASS")

  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  lower <- matrix(c(1, 0, 0, 0), nrow = 2)
  expect_error(
    bras_calibrate_raking(x0, c(0, 20), c(5, 15), lower = lower),
    "target 0 contains a cell with a positive `lower` bound"
  )
})

test_that(".drop_zero_margins / .reinsert_zero_margins round-trip", {
  x0 <- matrix(1:9, 3)
  d <- .drop_zero_margins(x0, c(5, 0, 7), c(0, 6, 6), matrix(0, 3, 3), matrix(Inf, 3, 3))
  expect_equal(c(d$m, d$n), c(2L, 2L))
  back <- .reinsert_zero_margins(matrix(1, 2, 2), d, 3, 3)
  expect_equal(back[2, ], c(0, 0, 0))
  expect_equal(back[, 2], c(1, 0, 1))
  expect_equal(back[, 1], c(0, 0, 0))
})

test_that("a single non-zero row/column is solved in closed form", {
  testthat::skip_if_not_installed("survey")
  testthat::skip_if_not_installed("MASS")

  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  fit <- bras_calibrate_raking(x0, c(0, 20), c(5, 15))
  expect_equal(unname(fit$matrix[1, ]), c(0, 0))
  expect_equal(unname(fit$matrix[2, ]), c(5, 15))
  expect_true(fit$converged)
  expect_error(
    bras_calibrate_raking(x0, c(0, 20), c(5, 15), upper = matrix(c(Inf, 4, Inf, Inf), 2)),
    "infeasible"
  )
})
