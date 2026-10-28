test_that("bras_freeze_sim respects bounds and hits margins when feasible", {
  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  row_targets <- c(20, 20)
  col_targets <- c(15, 25)
  upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)  # force cell [1,2] to cap at 12

  fit <- bras_freeze_sim(x0, row_targets, col_targets, upper = upper, tol = 1e-8)

  expect_true(fit$feasible)
  expect_true(all(fit$matrix <= upper + 1e-6))
  expect_equal(unname(rowSums(fit$matrix)), row_targets, tolerance = 1e-5)
  expect_equal(unname(colSums(fit$matrix)), col_targets, tolerance = 1e-5)
})

test_that("bras_freeze_sim and bras_freeze_seq agree when only one cell violates", {
  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  row_targets <- c(20, 20)
  col_targets <- c(15, 25)
  upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)

  fit_sim <- bras_freeze_sim(x0, row_targets, col_targets, upper = upper)
  fit_seq <- bras_freeze_seq(x0, row_targets, col_targets, upper = upper)

  expect_equal(fit_sim$matrix, fit_seq$matrix, tolerance = 1e-5)
})

test_that("an infeasible bound configuration is detected and errors by default", {
  x0 <- matrix(c(1, 1, 1, 1), nrow = 2)
  row_targets <- c(100, 2)
  col_targets <- c(51, 51)
  upper <- matrix(1, nrow = 2, ncol = 2)  # row 1 can sum to at most 2, but needs 100

  expect_error(
    bras_freeze_sim(x0, row_targets, col_targets, upper = upper, on_infeasible = "error"),
    "infeasible"
  )

  expect_warning(
    fit <- bras_freeze_sim(x0, row_targets, col_targets, upper = upper, on_infeasible = "warning"),
    "infeasible"
  )
  expect_false(fit$feasible)
})

test_that("bras_naive and bras_soft run and return a bras_fit object", {
  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  row_targets <- c(20, 20)
  col_targets <- c(15, 25)
  upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)

  fit_naive <- bras_naive(x0, row_targets, col_targets, upper = upper)
  fit_soft  <- bras_soft(x0, row_targets, col_targets, upper = upper)

  expect_s3_class(fit_naive, "bras_fit")
  expect_s3_class(fit_soft, "bras_fit")
})

test_that("bras_cap_iter runs, respects bounds, and returns a bras_fit object", {
  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  row_targets <- c(20, 20)
  col_targets <- c(15, 25)
  upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)

  fit <- bras_cap_iter(x0, row_targets, col_targets, upper = upper)

  expect_s3_class(fit, "bras_fit")
  expect_true(all(fit$matrix <= upper + 1e-6))
})

test_that("exhausting max_outer is reported as not converged, not as infeasible", {
  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  row_targets <- c(20, 20)
  col_targets <- c(15, 25)
  upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)  # needs at least one fixation pass

  expect_warning(
    fit <- bras_freeze_sim(x0, row_targets, col_targets, upper = upper, max_outer = 1),
    "without fully"
  )
  expect_true(fit$feasible)       # not proven infeasible, just unresolved
  expect_false(fit$converged)     # honestly reports it did not finish
  expect_equal(fit$n_outer_iter, 1L)
})
