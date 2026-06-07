###############################################################################
# test-miss-codes.R
#
# Contract tests for miss_codes behavior across the four public entry points:
#   M-1  oda_fit:              miss-coded obs predicted as NA_integer_
#   M-2  oda_cta_fit:          miss_codes stored; miss-coded obs excluded from fitting
#   M-3  predict.cta_tree:     path-local miss (missing_action="na") -> NA prediction
#   M-4  cta_assign_endpoints: path-local miss -> NA endpoint_id and endpoint_node_id
#
# All tests are CRAN-safe (no heavy MC, no slow fixtures).
###############################################################################

# Shared small dataset used by M-2 through M-4.
# x1: sentinel -9 at positions 4 and 8 (out of 12 rows).
# v2 is a clean second attribute that lets the tree split on v1 if available.
.mc_X <- data.frame(
  v1 = c(1, 2, 3, -9, 5, 6, 7, -9, 9, 10, 4, 8),
  v2 = c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 6, 4)
)
.mc_y <- c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 0L, 1L)

# =============================================================================
# M-1  oda_fit: miss-coded obs return NA from predict()
# =============================================================================

test_that("M-1: oda_fit miss_codes are stored and predict returns NA for sentinel values", {
  x <- c(1, 2, 3, -9, 5, 6, 7, -9, 9, 10)
  y <- c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L)
  fit <- oda_fit(x, y, miss_codes = -9, mcarlo = FALSE)

  expect_equal(fit$miss_codes, -9)
  preds <- predict(fit, x)
  expect_true(is.na(preds[4]))   # position of -9
  expect_true(is.na(preds[8]))   # position of -9
  # Non-missing positions must be non-NA
  expect_false(any(is.na(preds[-c(4, 8)])))
})

# =============================================================================
# M-2  oda_cta_fit: miss_codes stored in fitted tree
# =============================================================================

test_that("M-2: oda_cta_fit stores miss_codes in the returned cta_tree", {
  tree <- oda_cta_fit(.mc_X, .mc_y, miss_codes = -9, mc_iter = 50L, mindenom = 3L)
  expect_equal(tree$miss_codes, -9)
})

# =============================================================================
# M-3  predict.cta_tree: path-local NA for miss-coded split attribute
# =============================================================================

test_that("M-3: predict.cta_tree returns NA for miss-coded split attribute (missing_action='na')", {
  tree <- oda_cta_fit(.mc_X, .mc_y, miss_codes = -9, mc_iter = 50L, mindenom = 3L)
  # Three new rows: clean low, miss-coded, clean high.
  newX <- data.frame(v1 = c(2, -9, 8), v2 = c(9, 5, 3))
  preds <- predict(tree, newX, missing_action = "na")
  expect_length(preds, 3L)
  expect_true(is.na(preds[2]))          # -9 sentinel -> NA
  expect_false(is.na(preds[1]))         # clean obs -> non-NA
  expect_false(is.na(preds[3]))         # clean obs -> non-NA
})

# =============================================================================
# M-4  cta_assign_endpoints: path-local NA endpoint for miss-coded obs
# =============================================================================

test_that("M-4: cta_assign_endpoints returns NA endpoint for miss-coded obs (missing_action='na')", {
  tree <- oda_cta_fit(.mc_X, .mc_y, miss_codes = -9, mc_iter = 50L, mindenom = 3L)
  newX <- data.frame(v1 = c(2, -9, 8), v2 = c(9, 5, 3))
  ep   <- cta_assign_endpoints(tree, newX, missing_action = "na")
  expect_true(is.data.frame(ep))
  expect_true(is.na(ep$endpoint_node_id[2]))
  expect_true(is.na(ep$endpoint_id[2]))
  expect_false(is.na(ep$endpoint_node_id[1]))
  expect_false(is.na(ep$endpoint_node_id[3]))
})
