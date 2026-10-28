## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(RASbound)

## ----data---------------------------------------------------------------------
read_mat <- function(f) {
  as.matrix(utils::read.csv(system.file("extdata", f, package = "RASbound"),
                            header = FALSE))
}
x0 <- unname(read_mat("example_x0.csv"))
lower <- unname(read_mat("example_lower.csv"))
upper <- unname(read_mat("example_upper.csv"))
row_targets <- read_mat("example_row_targets.csv")[, 1]
col_targets <- read_mat("example_col_targets.csv")[, 1]

dim(x0)
c(sum(row_targets), sum(col_targets))
sum(x0 == 0)                  # exact zeros in the seed
sum(lower > 0)                # cells with a positive lower bound
tail(row_targets, 3)          # the last two rows (and columns) have target 0

## ----fit----------------------------------------------------------------------
args <- list(x0 = x0, row_targets = row_targets, col_targets = col_targets,
             lower = lower, upper = upper,
             tol = 1e-6, max_iter = 1000, max_outer = 500)
fit_sim <- suppressWarnings(do.call(bras_freeze_sim, args))
fit_seq <- suppressWarnings(do.call(bras_freeze_seq, args))
fit_hyb <- suppressWarnings(do.call(bras_freeze_hybrid, args))

sapply(list(sim = fit_sim, seq = fit_seq, hybrid = fit_hyb),
       function(f) c(converged = f$converged,
                     feasible = f$feasible,
                     n_fixed = nrow(f$fixed_cells)))

## ----check--------------------------------------------------------------------
ok <- function(f) c(
  within_bounds = all(f$matrix >= lower - 1e-6 & f$matrix <= upper + 1e-6),
  max_margin_error = max(f$diagnostics$row_margin_error,
                         f$diagnostics$col_margin_error)
)
sapply(list(sim = fit_sim, seq = fit_seq, hybrid = fit_hyb), ok)

## ----diff---------------------------------------------------------------------
max(abs(fit_sim$matrix - fit_seq$matrix))
max(abs(fit_sim$matrix - fit_hyb$matrix))

## ----plot, fig.width = 5, fig.height = 5--------------------------------------
plot(fit_seq)

## ----calibrate, eval = requireNamespace("survey", quietly = TRUE) && requireNamespace("MASS", quietly = TRUE)----
try_wrapper <- function(f) {
  tryCatch({
    fit <- suppressWarnings(f(x0, row_targets, col_targets,
                              lower = lower, upper = upper))
    paste("returned; converged =", fit$converged)
  }, error = function(e) paste("failed:", substr(conditionMessage(e), 1, 70)))
}
c(raking = try_wrapper(bras_calibrate_raking),
  linear = try_wrapper(bras_calibrate_linear))

