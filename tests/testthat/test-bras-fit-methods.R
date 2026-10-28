test_that("print.bras_fit shows the expected summary information", {
  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  fit <- ras(x0, c(20, 20), c(15, 25))

  expect_output(print(fit), "bras_fit")
  expect_output(print(fit), "method")
})

test_that("print.bras_fit reports fixed cells when present", {
  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)
  fit <- bras_freeze_sim(x0, c(20, 20), c(15, 25), upper = upper)

  expect_output(print(fit), "cells fixed at a bound")
})

test_that("print.bras_fit reports infeasibility diagnostics when infeasible", {
  x0 <- matrix(1, nrow = 2, ncol = 2)
  too_tight <- matrix(1, nrow = 2, ncol = 2)
  fit <- suppressWarnings(
    bras_freeze_sim(x0, c(100, 2), c(51, 51), upper = too_tight, on_infeasible = "warning")
  )

  expect_output(print(fit), "infeasible at rows")
})

test_that("summary.bras_fit returns a one-row data frame with the documented columns", {
  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)
  fit <- bras_freeze_sim(x0, c(20, 20), c(15, 25), upper = upper)

  s <- summary(fit)
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 1L)
  expect_true(all(c("method", "converged", "feasible", "n_outer_iter",
                     "n_fixed_cells", "max_row_margin_error",
                     "max_col_margin_error") %in% names(s)))
  expect_equal(s$n_fixed_cells, nrow(fit$fixed_cells))
})

test_that("print.bras_fit does not error when feasible is NA", {
  # Regression test: print.bras_fit() used `if (!x$feasible)`, which errors
  # on NA (`!NA` is NA, and `if (NA)` is a hard error in R). Every variant
  # that sets feasible = NA (bras_naive, bras_soft, bras_cap_iter, and the
  # two survey wrappers with optional bounds) was completely unprintable
  # before this fix -- found on real data, not synthesized.
  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)

  fit_naive <- bras_naive(x0, c(20, 20), c(15, 25), upper = upper)
  fit_soft  <- bras_soft(x0, c(20, 20), c(15, 25), upper = upper)
  fit_cap   <- bras_cap_iter(x0, c(20, 20), c(15, 25), upper = upper)

  expect_true(is.na(fit_naive$feasible))
  expect_true(is.na(fit_soft$feasible))
  expect_true(is.na(fit_cap$feasible))

  expect_output(print(fit_naive), "bras_fit")
  expect_output(print(fit_soft), "bras_fit")
  expect_output(print(fit_cap), "bras_fit")
})

test_that("plot.bras_fit runs without error and returns its input invisibly", {
  x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
  upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)
  fit <- bras_freeze_sim(x0, c(20, 20), c(15, 25), upper = upper)

  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)

  result <- withVisible(plot(fit))
  expect_false(result$visible)
  expect_identical(result$value, fit)
})

test_that("plot.bras_fit draws labels from dimnames and works without them", {
  x0 <- matrix(c(2, 4, 6, 8), nrow = 2,
               dimnames = list(c("a", "b"), c("u", "v")))
  upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)
  fit <- bras_freeze_sim(x0, c(20, 20), c(15, 25), upper = upper)
  expect_false(is.null(rownames(fit$matrix)))

  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)

  expect_no_error(plot(fit))
  expect_no_error(plot(fit, labels = FALSE))

  fit_plain <- bras_freeze_sim(unname(x0), c(20, 20), c(15, 25), upper = upper)
  expect_no_error(plot(fit_plain))
})
