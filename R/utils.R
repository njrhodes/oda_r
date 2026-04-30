###############################################################################
# R/utils.R
# Shared micro-utilities for odacore.
# Nothing here is exported except via NAMESPACE; this file is loaded first
# because R sources files in the R/ directory in alphabetical order.
###############################################################################

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

# p-value bucket (for deterministic MC significance category checks)
p_bucket <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<.001")
  if (p < 0.01)  return("<.01")
  if (p < 0.05)  return("<.05")
  if (p < 0.10)  return("<.10")
  ">=.10"
}
