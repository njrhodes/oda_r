###############################################################################
# R/utils.R
# Shared micro-utilities for odacore.
# Nothing here is exported except via NAMESPACE; this file is loaded first
# because R sources files in the R/ directory in alphabetical order.
###############################################################################

# Suppress R CMD check visible-binding notes for ggplot2 aes() column names.
# .eid   — ellipse grouping column in graphics_v3.R geom_polygon calls
# status_lbl — annotation column in plot_cta_balance_effects()
utils::globalVariables(c(".eid", "status_lbl"))

# Null-coalescing operator (defined once here; not re-defined elsewhere)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Tight integer quantisation used for floating-point tie comparisons.
# Must NEVER return NA_integer_; uses a safe sentinel for non-finite input.
tick <- function(z) {
  if (is.na(z) || !is.finite(z)) return(-2147483647L)
  v    <- round(as.numeric(z) * 1e9)
  vmax <- .Machine$integer.max
  vmin <- -vmax                         # NOT -vmax-1 (overflow hazard)
  if (v > vmax) v <- vmax
  if (v < vmin) v <- vmin
  as.integer(v)
}

# Simple formatting helpers used in harness output
fmt2 <- function(x) sprintf("%.2f", as.numeric(x))
fmt6 <- function(x) sprintf("%.6f", as.numeric(x))

identical_fmt2 <- function(x, y) identical(fmt2(x), fmt2(y))
identical_fmt6 <- function(x, y) identical(fmt6(x), fmt6(y))

# Internal case-weight validator  -  call at every public fit entrypoint.
# NULL w is always accepted (interpreted as unit weights by the caller).
# Non-NULL w must be numeric/integer, length n, finite, non-missing, > 0.
.validate_case_weights <- function(w, n, arg = "w") {
  if (is.null(w)) return(invisible(NULL))
  if (!is.numeric(w) && !is.integer(w))
    stop(sprintf("'%s' must be numeric or integer.", arg), call. = FALSE)
  if (length(w) != n)
    stop(sprintf("'%s' must have length %d (got %d).", arg, n, length(w)),
         call. = FALSE)
  if (anyNA(w))
    stop(sprintf("'%s' must not contain NA or NaN.", arg), call. = FALSE)
  if (any(!is.finite(w)))
    stop(sprintf("'%s' must be finite (no Inf or -Inf).", arg), call. = FALSE)
  if (any(w <= 0))
    stop(sprintf("'%s' must be strictly positive (no zeros or negatives).", arg),
         call. = FALSE)
  invisible(NULL)
}

# p-value bucket (for deterministic MC significance category checks)
p_bucket <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<.001")
  if (p < 0.01)  return("<.01")
  if (p < 0.05)  return("<.05")
  if (p < 0.10)  return("<.10")
  ">=.10"
}
