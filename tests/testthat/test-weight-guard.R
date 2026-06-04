###############################################################################
# tests/testthat/test-weight-guard.R
# Fast tests for .validate_case_weights() wired into public fit entrypoints.
# All tests use small synthetic data so no canon-fixture machinery is needed.
###############################################################################

# Minimal binary dataset for oda_fit tests
x2 <- c(1, 2, 3, 4, 5, 6)
y2 <- c(0L, 0L, 0L, 1L, 1L, 1L)
n2 <- length(x2)

# Minimal dataset for cta fits
X_cta <- data.frame(a = c(1, 2, 3, 4, 5, 6, 7, 8))
y_cta <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
n_cta <- nrow(X_cta)

# ------------------------------------------------------------------
# .validate_case_weights - direct unit tests (internal helper)
# ------------------------------------------------------------------

test_that("valid positive weights pass silently", {
  expect_silent(.validate_case_weights(c(1, 2, 1, 2, 1, 2), 6L))
})

test_that("NULL weights pass silently (unit-weight sentinel)", {
  expect_silent(.validate_case_weights(NULL, 6L))
})

test_that("NA weights are rejected", {
  expect_error(
    .validate_case_weights(c(1, NA, 1, 1, 1, 1), 6L),
    regexp = "NA or NaN"
  )
})

test_that("NaN weights are rejected", {
  expect_error(
    .validate_case_weights(c(1, NaN, 1, 1, 1, 1), 6L),
    regexp = "NA or NaN"
  )
})

test_that("Inf weights are rejected", {
  expect_error(
    .validate_case_weights(c(1, Inf, 1, 1, 1, 1), 6L),
    regexp = "finite"
  )
})

test_that("negative-Inf weights are rejected", {
  expect_error(
    .validate_case_weights(c(1, -Inf, 1, 1, 1, 1), 6L),
    regexp = "finite"
  )
})

test_that("zero weights are rejected", {
  expect_error(
    .validate_case_weights(c(1, 0, 1, 1, 1, 1), 6L),
    regexp = "strictly positive"
  )
})

test_that("negative weights are rejected", {
  expect_error(
    .validate_case_weights(c(1, -0.5, 1, 1, 1, 1), 6L),
    regexp = "strictly positive"
  )
})

test_that("wrong-length weights are rejected", {
  expect_error(
    .validate_case_weights(c(1, 1, 1), 6L),
    regexp = "length 6"
  )
})

test_that("non-numeric weights are rejected", {
  expect_error(
    .validate_case_weights(as.character(1:6), 6L),
    regexp = "numeric or integer"
  )
})

# ------------------------------------------------------------------
# oda_fit() wiring
# ------------------------------------------------------------------

test_that("oda_fit rejects NA weights", {
  w_bad <- c(1, NA, 1, 1, 1, 1)
  expect_error(
    oda_fit(x2, y2, w = w_bad, mcarlo = FALSE),
    regexp = "NA or NaN"
  )
})

test_that("oda_fit rejects zero weights", {
  w_bad <- c(1, 0, 1, 1, 1, 1)
  expect_error(
    oda_fit(x2, y2, w = w_bad, mcarlo = FALSE),
    regexp = "strictly positive"
  )
})

test_that("oda_fit accepts valid weights", {
  w_ok <- rep(2, n2)
  expect_no_error(oda_fit(x2, y2, w = w_ok, mcarlo = FALSE))
})

# ------------------------------------------------------------------
# oda_cta_fit() wiring
# ------------------------------------------------------------------

test_that("oda_cta_fit rejects NA weights", {
  w_bad <- c(1, NA, rep(1, n_cta - 2L))
  expect_error(
    oda_cta_fit(X_cta, y_cta, w = w_bad, mindenom = 2L, mc_iter = 100L),
    regexp = "NA or NaN"
  )
})

test_that("oda_cta_fit rejects zero weights", {
  w_bad <- c(0, rep(1, n_cta - 1L))
  expect_error(
    oda_cta_fit(X_cta, y_cta, w = w_bad, mindenom = 2L, mc_iter = 100L),
    regexp = "strictly positive"
  )
})

test_that("oda_cta_fit accepts valid weights", {
  w_ok <- rep(1.5, n_cta)
  expect_no_error(
    oda_cta_fit(X_cta, y_cta, w = w_ok, mindenom = 2L, mc_iter = 100L,
                mc_seed = 1L)
  )
})

# ------------------------------------------------------------------
# cta_descendant_family() wiring
# ------------------------------------------------------------------

test_that("cta_descendant_family rejects negative weights", {
  w_bad <- c(-1, rep(1, n_cta - 1L))
  expect_error(
    cta_descendant_family(X_cta, y_cta, w = w_bad,
                          start_mindenom = 2L, mc_iter = 100L),
    regexp = "strictly positive"
  )
})
