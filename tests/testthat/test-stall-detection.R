# Found via randomized search (trial 5, seed 1; see the project's ad hoc
# search script, not part of the installed package) after the real-world
# 55x55 case that first surfaced this failure mode -- a genuine fixed
# point of the alternating row/column scaling, distinct from both
# `max_outer` exhaustion and proven infeasibility: no bound violations
# remain, but the margins are not matched, and the margin error stops
# changing entirely across iterations (confirmed by tracing the real
# case: identical to 6 significant figures for 2000 consecutive
# iterations). See ?bras_fit and vignette("RASbound-diagnostics") for the
# explanation (the feasibility check is necessary but not sufficient).

test_that("a genuine stall (fixed point, not slow convergence) is detected and reported distinctly", {
  x0 <- matrix(c(4, 20,  5, 49,
                 10, 41, 44, 35,
                 29,  9, 24, 18,
                 12, 30, 16, 46,
                 46, 14,  6, 10), nrow = 5, byrow = TRUE)
  row_targets <- c(95, 120, 63, 67, 134)
  col_targets <- c(146, 215, 82, 36)
  upper <- matrix(c(33, 31, 13, 39,
                     35, 41, 13, 46,
                     17, 34, 48,  6,
                     23, 53, 51, 35,
                     53, 57, 50, 52), nrow = 5, byrow = TRUE)

  expect_warning(
    fit <- bras_freeze_sim(x0, row_targets, col_targets, upper = upper,
                            on_infeasible = "warning", max_outer = 30),
    "NOT resolve this"
  )

  expect_true(fit$feasible)          # not proven infeasible
  expect_false(fit$converged)        # but margins genuinely not matched
  expect_true(fit$diagnostics$stalled)
  # no bound violations remain, even though margins don't match:
  expect_true(all(fit$matrix <= upper + 1e-6))
})

test_that("stalled defaults to FALSE and does not appear for a well-posed problem", {
  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  fit <- bras_freeze_sim(x0, c(20, 20), c(15, 25))
  expect_false(fit$diagnostics$stalled)
})
