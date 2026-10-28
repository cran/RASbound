## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(RASbound)

## -----------------------------------------------------------------------------
x0 <- matrix(c(2, 4, 9,
               6, 8, 1,
               3, 5, 7), nrow = 3, byrow = TRUE)
row_targets <- c(25, 20, 25)
col_targets <- c(20, 25, 25)
upper <- matrix(c(Inf, Inf, 6,
                  Inf, Inf, Inf,
                  Inf, Inf, Inf), nrow = 3, byrow = TRUE)  # cap cell [1, 3] at 6

fit <- bras_freeze_hybrid(x0, row_targets, col_targets, upper = upper)
fit

## -----------------------------------------------------------------------------
fit$fixed_cells

## -----------------------------------------------------------------------------
too_tight <- matrix(1, nrow = 3, ncol = 3)
fit_bad <- suppressWarnings(  # x0 is also above `upper` here; that is not the point
  bras_freeze_sim(x0, row_targets, col_targets, upper = too_tight,
                  on_infeasible = "warning")
)
fit_bad$feasible
fit_bad$diagnostics$infeasible_at

## -----------------------------------------------------------------------------
rowSums(too_tight)  # every row can reach at most 3, but the targets are 25, 20, 25

## -----------------------------------------------------------------------------
fit_slow <- bras_freeze_sim(x0, row_targets, col_targets, upper = upper,
                            max_outer = 1, on_infeasible = "warning")
fit_slow$feasible    # TRUE: not proven infeasible
fit_slow$converged   # FALSE: it just did not finish

## -----------------------------------------------------------------------------
x0_with_zero <- x0
x0_with_zero[1, 1] <- 0
fit_zero <- bras_freeze_sim(x0_with_zero, row_targets, col_targets)
fit_zero$matrix[1, 1]  # exactly 0

## -----------------------------------------------------------------------------
bad_lower <- matrix(0, 3, 3); bad_lower[1, 1] <- 5
fit_imputed <- bras_freeze_sim(x0_with_zero, row_targets, col_targets,
                               lower = bad_lower)
fit_imputed$matrix[1, 1]  # at least 5, not 0

## -----------------------------------------------------------------------------
fits <- list(
  freeze_sim  = bras_freeze_sim(x0, row_targets, col_targets, upper = upper),
  freeze_seq  = bras_freeze_seq(x0, row_targets, col_targets, upper = upper),
  naive       = bras_naive(x0, row_targets, col_targets, upper = upper),
  soft        = bras_soft(x0, row_targets, col_targets, upper = upper),
  cap_iter    = bras_cap_iter(x0, row_targets, col_targets, upper = upper)
)

## ----eval = requireNamespace("survey", quietly = TRUE) && requireNamespace("MASS", quietly = TRUE)----
fits$calibrate_raking <- bras_calibrate_raking(x0, row_targets, col_targets, upper = upper)
fits$calibrate_linear <- bras_calibrate_linear(x0, row_targets, col_targets, upper = upper)
# bras_calibrate_logit() is left out on purpose. The section
# "When bras_calibrate_logit() fails" below uses this same problem.

## -----------------------------------------------------------------------------
data.frame(
  method    = names(fits),
  converged = vapply(fits, function(f) f$converged, logical(1)),
  max_row_err = vapply(fits, function(f) max(f$diagnostics$row_margin_error), numeric(1)),
  max_col_err = vapply(fits, function(f) max(f$diagnostics$col_margin_error), numeric(1)),
  dist_from_seed = vapply(fits, function(f) sum(abs(f$matrix - x0)), numeric(1))
)

## ----eval = requireNamespace("survey", quietly = TRUE) && requireNamespace("MASS", quietly = TRUE)----
set.seed(20260919)
x0_cmp <- matrix(round(runif(16, 1, 20)), nrow = 4)
row_targets_cmp <- round(runif(4, 80, 150))
props <- runif(4); props <- props / sum(props)
col_targets_cmp <- round(sum(row_targets_cmp) * props)
col_targets_cmp[4] <- sum(row_targets_cmp) - sum(col_targets_cmp[1:3])

fit_ras_cmp    <- ras(x0_cmp, row_targets_cmp, col_targets_cmp, tol = 1e-10, max_iter = 500)
fit_raking_cmp <- bras_calibrate_raking(x0_cmp, row_targets_cmp, col_targets_cmp)
fit_linear_cmp <- bras_calibrate_linear(x0_cmp, row_targets_cmp, col_targets_cmp)

# raking is very close to our own engine (same distance function):
max(abs(fit_ras_cmp$matrix - fit_raking_cmp$matrix))

# linear is not, although it meets the same margins to the same tolerance:
max(abs(fit_ras_cmp$matrix - fit_linear_cmp$matrix))
rbind(row_ras = rowSums(fit_ras_cmp$matrix), row_linear = rowSums(fit_linear_cmp$matrix))

## ----error = TRUE, eval = requireNamespace("survey", quietly = TRUE) && requireNamespace("MASS", quietly = TRUE)----
try({
# Self-contained: the running example again, so this chunk gives the same
# result in whatever order the document is run.
x0 <- matrix(c(2, 4, 9,
               6, 8, 1,
               3, 5, 7), nrow = 3, byrow = TRUE)
row_targets <- c(25, 20, 25)
col_targets <- c(20, 25, 25)
upper <- matrix(c(Inf, Inf, 6,
                  Inf, Inf, Inf,
                  Inf, Inf, Inf), nrow = 3, byrow = TRUE)
bras_calibrate_logit(x0, row_targets, col_targets, upper = upper)
})

## -----------------------------------------------------------------------------
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
  bras_freeze_sim(x0, row_targets, col_targets, upper = upper, on_infeasible = "warning")
)
fit_seq <- suppressWarnings(
  bras_freeze_seq(x0, row_targets, col_targets, upper = upper, on_infeasible = "warning")
)
fit_hybrid <- bras_freeze_hybrid(x0, row_targets, col_targets, upper = upper)

data.frame(
  method    = c("sim", "seq", "hybrid"),
  feasible  = c(fit_sim$feasible, fit_seq$feasible, fit_hybrid$feasible),
  converged = c(fit_sim$converged, fit_seq$converged, fit_hybrid$converged),
  n_fixed   = c(nrow(fit_sim$fixed_cells), nrow(fit_seq$fixed_cells), nrow(fit_hybrid$fixed_cells))
)

## -----------------------------------------------------------------------------
fit_hybrid$method

## -----------------------------------------------------------------------------
x0_stall <- matrix(c(4, 20,  5, 49,
                     10, 41, 44, 35,
                     29,  9, 24, 18,
                     12, 30, 16, 46,
                     46, 14,  6, 10), nrow = 5, byrow = TRUE)
row_targets_stall <- c(95, 120, 63, 67, 134)
col_targets_stall <- c(146, 215, 82, 36)
upper_stall <- matrix(c(33, 31, 13, 39,
                        35, 41, 13, 46,
                        17, 34, 48,  6,
                        23, 53, 51, 35,
                        53, 57, 50, 52), nrow = 5, byrow = TRUE)

fit_stall <- suppressWarnings(
  bras_freeze_sim(x0_stall, row_targets_stall, col_targets_stall,
                  upper = upper_stall, on_infeasible = "warning")
)

fit_stall$feasible             # TRUE: not proven infeasible
fit_stall$converged            # FALSE: the margins are not matched
fit_stall$diagnostics$stalled  # TRUE: a fixed point, not slow progress

## ----inf-data-----------------------------------------------------------------
read_mat <- function(f) {
  as.matrix(utils::read.csv(system.file("extdata", f, package = "RASbound"),
                            header = FALSE))
}
x0i <- unname(read_mat("infeasible_x0.csv"))
upperi <- unname(read_mat("infeasible_upper.csv"))
rti <- read_mat("infeasible_row_targets.csv")[, 1]
cti <- read_mat("infeasible_col_targets.csv")[, 1]

fit_i <- suppressWarnings(
  bras_freeze_seq(x0i, rti, cti, upper = upperi,
                  tol = 1e-6, max_iter = 1000, max_outer = 500,
                  on_infeasible = "warning")
)
fit_i$converged
fit_i$diagnostics$stalled

## ----inf-simple---------------------------------------------------------------
cap <- ifelse(x0i > 0, upperi, 0)   # effective capacity under RAS
any(rti > rowSums(cap))
any(cti > colSums(cap))

## ----inf-cert-----------------------------------------------------------------
i <- 13
rti[i]                         # what the row must receive
sum(cap[i, ])                  # its raw capacity: it looks sufficient
sum(pmin(cti, cap[i, ]))       # what it can really absorb

## ----inf-relaxed--------------------------------------------------------------
fit_r <- suppressWarnings(
  bras_freeze_seq(x0i, rti, cti, upper = round(upperi * 1.25),
                  tol = 1e-6, max_iter = 1000, max_outer = 500,
                  on_infeasible = "warning")
)
fit_r$converged

