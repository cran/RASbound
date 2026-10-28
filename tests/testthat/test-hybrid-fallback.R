# The case below was found by randomized search (see find_hybrid_fallback.R
# in the project notes, not part of the installed package) after the first
# search attempt (looser bounds) found nothing in 2000 tries. This one --
# tight, unevenly distributed row-level bounds on a 5x3 matrix -- found a
# hit at iteration 3682 out of 5000, seed 1. Fixed here as a deterministic
# regression test so the sequential-fallback branch inside
# bras_freeze_hybrid() (simultaneous fixation found locally infeasible
# mid-pass, retried fixing only the most extreme cell) is actually
# exercised by the test suite, not just reasoned about.

test_that("bras_freeze_hybrid's sequential-fallback branch is actually exercised", {
  x0 <- matrix(c(29, 16, 22,
                 14,  4, 20,
                 27, 20, 23,
                 29, 23, 13,
                 26, 20, 12), nrow = 5, byrow = TRUE)
  row_targets <- c(91, 45, 67, 71, 70)
  col_targets <- c(133, 89, 122)
  upper <- matrix(c(40, 47, 27,
                     17, 19, 13,
                     44, 16, 17,
                     29, 19, 41,
                     25, 16, 33), nrow = 5, byrow = TRUE)

  # The key assertion: this specific input actually reaches the fallback
  # branch (survey/RASbound would otherwise have no way to know without
  # reading the source and reasoning about it, as we had to do before this
  # case was found).
  expect_message(
    suppressWarnings(
      fit <- bras_freeze_hybrid(x0, row_targets, col_targets, upper = upper,
                                 on_infeasible = "warning", verbose = TRUE,
                                 max_outer = 30)
    ),
    "infeasible at pass"
  )

  # Confirmed empirically (see vignette("RASbound-diagnostics")): on this
  # exact input, per-pass fallback alone is NOT enough -- bras_freeze_sim()
  # and a hybrid run WITHOUT the global retry both end up spuriously
  # infeasible here, while bras_freeze_seq() succeeds. This is exactly why
  # bras_freeze_hybrid() now retries the whole problem with pure sequential
  # fixation when the per-pass approach fails -- so it must succeed here,
  # via that retry.
  expect_true(fit$feasible)
  expect_true(fit$converged)
  expect_equal(fit$method, "freeze_hybrid_retry_sequential")
  expect_true(all(fit$matrix <= upper + 1e-4))
  expect_equal(unname(rowSums(fit$matrix)), row_targets, tolerance = 1e-3)
  expect_equal(unname(colSums(fit$matrix)), col_targets, tolerance = 1e-3)
})

test_that("bras_freeze_sim (no fallback) genuinely fails where hybrid's global retry succeeds", {
  # Same input as above: documents, with a fixed regression test, that
  # bras_freeze_sim() alone is NOT equivalent to bras_freeze_hybrid() here
  # -- over-fixation across several passes (not just within one) makes it
  # spuriously infeasible, confirmed empirically.
  x0 <- matrix(c(29, 16, 22,
                 14,  4, 20,
                 27, 20, 23,
                 29, 23, 13,
                 26, 20, 12), nrow = 5, byrow = TRUE)
  row_targets <- c(91, 45, 67, 71, 70)
  col_targets <- c(133, 89, 122)
  upper <- matrix(c(40, 47, 27,
                     17, 19, 13,
                     44, 16, 17,
                     29, 19, 41,
                     25, 16, 33), nrow = 5, byrow = TRUE)

  fit_sim <- suppressWarnings(
    bras_freeze_sim(x0, row_targets, col_targets, upper = upper,
                     on_infeasible = "warning", max_outer = 30)
  )
  expect_s3_class(fit_sim, "bras_fit")
  expect_false(fit_sim$feasible)

  fit_seq <- suppressWarnings(
    bras_freeze_seq(x0, row_targets, col_targets, upper = upper,
                     on_infeasible = "warning", max_outer = 30)
  )
  expect_true(fit_seq$feasible)
  expect_true(fit_seq$converged)
})
