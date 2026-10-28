# RASbound

Bounded RAS (biproportional matrix balancing, also called iterative
proportional fitting or raking) with lower and upper bounds on individual
cells.

Classic RAS adjusts a seed matrix so that its row and column sums match
given targets. RASbound adds cell-level bounds, for example when updating
or regionalizing an input-output table and some cells must stay within a
plausible range. If no matrix satisfies both the targets and the bounds,
the functions say so, and tell you which rows and columns are involved.

## Main functions

| Function | What it does |
|:---|:---|
| `ras()` | Classic RAS, no bounds |
| `bras_freeze_sim()`, `bras_freeze_seq()`, `bras_freeze_hybrid()` | Bounded RAS that solves, fixes violating cells at their bound, and repeats. Margins are matched exactly and bounds are respected. `bras_freeze_hybrid()` is the recommended default. |
| `bras_naive()`, `bras_soft()`, `bras_cap_iter()` | Simpler reference methods, for comparison |
| `bras_calibrate_logit()`, `bras_calibrate_linear()`, `bras_calibrate_raking()` | Wrappers around `survey::calibrate()` (need the `survey` and `MASS` packages) |

All functions take the same arguments and return an object of class
`bras_fit`.

## Example

```r
library(RASbound)

x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
upper <- matrix(c(Inf, Inf, 12, Inf), nrow = 2)  # cap cell [1, 2] at 12

fit <- bras_freeze_hybrid(x0, row_targets = c(20, 20),
                          col_targets = c(15, 25), upper = upper)
fit
fit$matrix
```

## Documentation

* `vignette("RASbound-intro")`: introduction and choice of function.
* `vignette("RASbound-example")`: a worked example with a real
  input-output-style table.
* `vignette("RASbound-diagnostics")`: reading the results and diagnosing
  failures.

## License

GPL-3. See `LICENSE.md`.
