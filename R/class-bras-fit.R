#' The `bras_fit` result object
#'
#' Every function in RASbound returns an object of class `bras_fit`. This
#' includes [ras()], the `bras_freeze_*()` functions, the reference
#' variants and the `bras_calibrate_*()` wrappers. It is a plain list with
#' the fields below, so code that uses one method's result also works with
#' the result of any other method.
#'
#' @section Fields:
#' \describe{
#'   \item{`matrix`}{The adjusted matrix. For the hard-freeze family this
#'     exactly matches `row_targets`/`col_targets` and respects
#'     `lower`/`upper` whenever `converged` is `TRUE`; for the reference
#'     variants and the wrappers, check `diagnostics` for how close it
#'     actually got.}
#'   \item{`converged`}{`TRUE` if the method finished having satisfied its
#'     own stopping criterion (exact margins for the hard-freeze family and
#'     `ras()`; margin tolerance `tol` for the reference variants).}
#'   \item{`feasible`}{For the hard-freeze family: `TRUE` unless a bound
#'     configuration was proven to admit no solution (see
#'     `diagnostics$infeasible_at`); can also be `TRUE` while `converged`
#'     is `FALSE` if `max_outer` ran out before a solution was found or
#'     ruled out. These are different situations; see Details. For
#'     methods that do not attempt to diagnose infeasibility (the
#'     reference variants, the wrappers), this is `NA`.}
#'   \item{`n_outer_iter`}{Number of outer fixation passes (hard-freeze
#'     family) or scaling iterations (reference variants) used.}
#'   \item{`fixed_cells`}{A data frame with one row per cell permanently
#'     fixed at a bound: `row`, `col`, `bound` (`"lower"`/`"upper"`),
#'     `value`, and `outer_pass` (which pass fixed it). Empty for methods
#'     that do not freeze cells.}
#'   \item{`diagnostics`}{A list with `row_margin_error` and
#'     `col_margin_error` (absolute deviation from `row_targets`/
#'     `col_targets` in the returned `matrix`), `infeasible_at`
#'     (`NULL` unless infeasibility was proven, in which case a list with
#'     `rows`/`cols` naming the offending margins), and `stalled` (hard-freeze
#'     family only; see Details).}
#'   \item{`method`}{A short string identifying which function produced
#'     this fit (e.g. `"freeze_simultaneous"`, `"ras_classic"`).}
#'   \item{`call`}{The original function call, as from [match.call()].}
#' }
#'
#' @details
#' `feasible = TRUE` together with `converged = FALSE` specifically means
#' the freeze algorithm ran out of `max_outer` passes before it found a
#' solution or proved that none exists. Increasing `max_outer` may solve
#' it. `feasible = FALSE` means that the bounds were shown to admit no
#' solution (see `diagnostics$infeasible_at`). More iterations cannot fix
#' that. You have to change the bounds or the targets.
#'
#' A third, distinct situation (hard-freeze family only):
#' `diagnostics$stalled = TRUE` means no bound violations remained, but the
#' inner alternating row/column solve reached a genuine fixed point that
#' does not satisfy the margins, and stayed there. The margin error stopped
#' changing across iterations; it was not just slow to converge. This happens when the free cells left after fixation form a
#' structurally degenerate subsystem: row-wise pressure to scale up is
#' exactly cancelled by column-wise pressure to scale down (or vice versa).
#' This is a known limitation of alternating-projection RAS/IPF methods in
#' general, not specific to this package. The feasibility check
#' (`sum(upper) >= residual target` for each row and column) is necessary
#' but not sufficient, and a stall is a case it cannot detect. A stall
#' often means that no matrix satisfies the targets and the bounds
#' together. Unlike `max_outer` exhaustion, increasing `max_iter` or
#' `max_outer` will not help; see `vignette("RASbound-diagnostics")`.
#'
#' @name bras_fit
#' @aliases bras_fit-class
NULL

#' @keywords internal
#' @noRd
.new_bras_fit <- function(matrix, converged, feasible, n_outer_iter,
                           fixed_cells, row_margin_error, col_margin_error,
                           infeasible_at, method, call, stalled = FALSE) {
  structure(
    list(
      matrix = matrix,
      converged = converged,
      feasible = feasible,
      n_outer_iter = n_outer_iter,
      fixed_cells = fixed_cells,
      diagnostics = list(
        row_margin_error = row_margin_error,
        col_margin_error = col_margin_error,
        infeasible_at = infeasible_at,
        stalled = stalled
      ),
      method = method,
      call = call
    ),
    class = "bras_fit"
  )
}

#' Print a fitted bounded-RAS object
#'
#' @param x A `bras_fit` object.
#' @param ... Unused, present for S3 consistency.
#' @return `x`, invisibly. The function is called for its side effect: it
#'   prints the method, whether the fit converged and was feasible, the
#'   number of outer passes, the number of cells fixed at a bound, and, if
#'   the problem was found infeasible, the rows and columns involved.
#' @export
print.bras_fit <- function(x, ...) {
  cat("<bras_fit> method:", x$method, "\n")
  cat("  converged:", x$converged, " | feasible:", x$feasible,
      " | outer passes:", x$n_outer_iter, "\n")
  if (nrow(x$fixed_cells) > 0) {
    cat("  cells fixed at a bound:", nrow(x$fixed_cells), "\n")
  }
  if (isFALSE(x$feasible)) {
    cat("  infeasible at rows:", paste(x$diagnostics$infeasible_at$rows, collapse = ", "),
        " cols:", paste(x$diagnostics$infeasible_at$cols, collapse = ", "), "\n")
  }
  invisible(x)
}

#' Summarize a fitted bounded-RAS object
#'
#' @param object A `bras_fit` object.
#' @param ... Unused, present for S3 consistency.
#' @return A `data.frame` with one row and the columns `method` (the
#'   method used), `converged` (logical), `feasible` (logical, `NA` for
#'   methods that do not check feasibility), `n_outer_iter` (number of outer
#'   passes), `n_fixed_cells` (number of cells fixed at a bound) and
#'   `max_row_margin_error` and `max_col_margin_error` (the largest absolute
#'   error in the row and column totals of the fitted matrix).
#' @export
summary.bras_fit <- function(object, ...) {
  data.frame(
    method = object$method,
    converged = object$converged,
    feasible = object$feasible,
    n_outer_iter = object$n_outer_iter,
    n_fixed_cells = nrow(object$fixed_cells),
    max_row_margin_error = suppressWarnings(max(object$diagnostics$row_margin_error)),
    max_col_margin_error = suppressWarnings(max(object$diagnostics$col_margin_error))
  )
}

#' Plot which cells were fixed at a bound
#'
#' Draws the fitted matrix as an image and marks the cells that were fixed
#' at a bound.
#'
#' @section How to read the plot:
#' * Each square is one cell of the matrix. Row 1 is at the top and column
#'   1 is on the left. If the fitted matrix has row or column names
#'   (`rownames()` and `colnames()`), they are shown as axis labels. Without
#'   names there are no labels.
#' * The colour shows the value of the cell in the fitted matrix. Lower
#'   values are pale yellow and higher values are dark red. A few very
#'   large cells can make all the others look similar, since the colours
#'   follow the raw values.
#' * An upward triangle marks a cell fixed at its **upper** bound. A
#'   downward triangle marks a cell fixed at its **lower** bound. Cells
#'   without a marker were never fixed.
#' * Methods that do not fix cells (for example [ras()] or
#'   [bras_naive()]) show no markers.
#'
#' The plot gives a quick view of where the bounds are binding. The exact
#' cells, bounds and values are in `fit$fixed_cells`.
#'
#' @param x A `bras_fit` object.
#' @param labels Logical. If `TRUE` (the default), show the row and column
#'   names of the matrix as axis labels when it has them. Use `FALSE` to
#'   hide them.
#' @param cex.axis Size of the axis labels. Reduce it for large matrices.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @examples
#' x0 <- matrix(c(2, 4, 6, 8), nrow = 2)
#' fit <- bras_freeze_sim(x0, c(20, 20), c(15, 25),
#'                        upper = matrix(c(Inf, Inf, 12, Inf), nrow = 2))
#' plot(fit)
#' @export
plot.bras_fit <- function(x, labels = TRUE, cex.axis = 0.6, ...) {
  m <- x$matrix
  nr <- nrow(m); nc <- ncol(m)
  xs <- if (nc > 1) (seq_len(nc) - 1) / (nc - 1) else 0.5
  ys <- if (nr > 1) 1 - (seq_len(nr) - 1) / (nr - 1) else 0.5

  rn <- rownames(m); cn <- colnames(m)
  show_labels <- isTRUE(labels) && (!is.null(rn) || !is.null(cn))
  if (show_labels) {
    op <- graphics::par(mar = c(6, 6, 4, 1))
    on.exit(graphics::par(op), add = TRUE)
  }

  graphics::image(t(m[nr:1, , drop = FALSE]), axes = FALSE,
                   col = grDevices::hcl.colors(12, "YlOrRd", rev = TRUE),
                   main = paste0("RASbound fit (", x$method, ")"))
  if (show_labels) {
    if (!is.null(cn)) {
      graphics::axis(1, at = xs, labels = cn, las = 2, tick = FALSE,
                     cex.axis = cex.axis)
    }
    if (!is.null(rn)) {
      graphics::axis(2, at = ys, labels = rn, las = 2, tick = FALSE,
                     cex.axis = cex.axis)
    }
  }
  if (nrow(x$fixed_cells) > 0) {
    graphics::points(
      xs[x$fixed_cells$col], ys[x$fixed_cells$row],
      pch = ifelse(x$fixed_cells$bound == "upper", 24, 25),
      bg = "white"
    )
  }
  invisible(x)
}
