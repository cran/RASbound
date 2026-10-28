test_that("classic ras() matches target margins on a simple 2x2 case", {
  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  row_targets <- c(20, 20)
  col_targets <- c(15, 25)

  fit <- ras(x0, row_targets, col_targets, tol = 1e-10)

  expect_true(fit$converged)
  expect_equal(unname(rowSums(fit$matrix)), row_targets, tolerance = 1e-6)
  expect_equal(unname(colSums(fit$matrix)), col_targets, tolerance = 1e-6)
  expect_s3_class(fit, "bras_fit")
})

test_that("ras() rejects mismatched row/column totals", {
  x0 <- matrix(1:4, nrow = 2)
  expect_error(ras(x0, c(1, 2), c(1, 1)), "same total")
})

test_that("ras() rejects negative seed values", {
  x0 <- matrix(c(-1, 2, 3, 4), nrow = 2)
  expect_error(ras(x0, c(1, 8), c(2, 7)), "non-negative")
})

test_that("lower > 0 at a zero-seeded cell is imputed with a warning, not rejected", {
  x0 <- matrix(c(0, 4, 6, 8), nrow = 2)  # x0[1,1] == 0
  bad_lower <- matrix(c(5, 0, 0, 0), nrow = 2)  # lower[1,1] = 5, unsolvable as seeded

  expect_warning(
    fit <- bras_freeze_sim(x0, c(20, 20), c(15, 25), lower = bad_lower),
    "solvable at all"
  )
  expect_true(fit$matrix[1, 1] >= 5 - 1e-6)
  expect_equal(unname(rowSums(fit$matrix)), c(20, 20), tolerance = 1e-5)
})

test_that("x0 above upper is imputed with a warning and does not change the final result", {
  x0 <- matrix(c(2, 4, 20, 8), nrow = 2)  # x0[1,2] = 20, well above upper below
  upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)

  expect_warning(
    fit_seeded_high <- bras_freeze_sim(x0, c(20, 20), c(15, 25), upper = upper),
    "outside their own"
  )
  fit_seeded_ok <- bras_freeze_sim(matrix(c(2, 4, 6, 8), nrow = 2), c(20, 20), c(15, 25),
                                    upper = upper)
  expect_equal(fit_seeded_high$matrix, fit_seeded_ok$matrix, tolerance = 1e-6)
})

test_that("a row seeded entirely at zero is caught as infeasible, not overstated as feasible", {
  # row 1 is entirely zero-seeded: under multiplicative RAS it can only ever
  # sum to 0, but its target is 20 -- infeasible regardless of how loose the
  # nominal bounds are (default lower = 0, upper = Inf).
  x0 <- matrix(c(0, 4, 0, 8), nrow = 2)  # row 1: (0, 0); row 2: (4, 8)

  expect_error(
    bras_freeze_sim(x0, row_targets = c(20, 20), col_targets = c(15, 25)),
    "infeasible"
  )
})

test_that("bras_fit objects have the fields documented in ?bras_fit", {
  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  fit <- ras(x0, c(20, 20), c(15, 25))

  expect_named(fit, c("matrix", "converged", "feasible", "n_outer_iter",
                       "fixed_cells", "diagnostics", "method", "call"))
  expect_named(fit$diagnostics, c("row_margin_error", "col_margin_error", "infeasible_at", "stalled"))
  expect_named(fit$fixed_cells, c("row", "col", "bound", "value", "outer_pass"))
})
