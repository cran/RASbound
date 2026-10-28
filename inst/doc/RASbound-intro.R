## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----example------------------------------------------------------------------
library(RASbound)

x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
row_targets <- c(20, 20)
col_targets <- c(15, 25)
upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)  # cap cell [1, 2] at 12

fit <- bras_freeze_sim(x0, row_targets, col_targets, upper = upper)
fit
fit$matrix
summary(fit)

