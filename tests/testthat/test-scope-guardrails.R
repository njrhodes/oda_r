###############################################################################
# test-scope-guardrails.R
#
# Structural guardrail tests documenting and protecting:
#   1. Lean-fit invariant: cta_tree / cta_ort objects must not accumulate
#      large fit-time data stores (training rows, caches, permutations).
#   2. SORT/GORT export absence: reserved taxonomy names must not be exported.
#
# Note: LORT taxonomy metadata (method == "lort", global_optimization == FALSE,
# sda_anchored == FALSE) is already asserted in test-cta-ort.R under the
# "LORT taxonomy" block.  These tests add the lean-fit and export dimensions.
#
# Tier: CRAN-safe.  All fits use small mc_iter and loo = "off".
###############################################################################

# ---- Minimal fixtures -------------------------------------------------------
# B perfectly separates: B > 0.5 -> y = 1; B <= 0.5 -> y = 0.
.sg_X <- data.frame(
  A = c(rep(0L, 15L), rep(1L, 15L), rep(1L, 15L)),
  B = c(rep(0L, 15L), rep(0L, 15L), rep(1L, 15L))
)
.sg_y <- c(rep(0L, 15L), rep(0L, 15L), rep(1L, 15L))

.sg_tree <- cta_fit(.sg_X, .sg_y,
                    mindenom = 1L, mc_iter = 50L, mc_seed = 42L, loo = "off")

.sg_ort  <- cta_fit(.sg_X, .sg_y,
                    recursive = TRUE,
                    mc_iter = 50L, mc_seed = 42L, loo = "off", min_n = 5L)

# Known-forbidden field names: large fit-time data stores that must never be
# added to cta_tree or cta_ort objects.
.lean_forbidden <- c(
  "training_data", "X", "y", "frame_data",
  "node_rows", "row_indices_by_node",
  "full_fit_cache", "bootstrap_cache", "permutations"
)

# =============================================================================
# 1.  Lean-fit invariant
# =============================================================================

test_that("lean-fit: cta_tree top-level has no forbidden storage fields", {
  bad <- intersect(names(.sg_tree), .lean_forbidden)
  expect_equal(bad, character(0L),
               info = paste("forbidden fields found:", paste(bad, collapse = ", ")))
})

test_that("lean-fit: cta_ort top-level has no forbidden storage fields", {
  bad <- intersect(names(.sg_ort), .lean_forbidden)
  expect_equal(bad, character(0L),
               info = paste("forbidden fields found:", paste(bad, collapse = ", ")))
})

test_that("lean-fit: cta_ort ort_nodes have no forbidden storage fields", {
  for (nd in .sg_ort$ort_nodes) {
    if (is.null(nd)) next
    bad <- intersect(names(nd), .lean_forbidden)
    expect_equal(bad, character(0L),
                 info = paste("forbidden fields in ort_node:", paste(bad, collapse = ", ")))
  }
})

# =============================================================================
# 2.  SORT/GORT export absence
# =============================================================================

test_that("SORT/GORT reserved names are not exported from the oda namespace", {
  exports  <- getNamespaceExports("oda")
  reserved <- c("sort_fit", "gort_fit",
                "cta_sort_fit", "cta_gort_fit",
                "lort_sort_fit", "oda_sort_fit", "oda_gort_fit")
  found <- intersect(exports, reserved)
  expect_equal(found, character(0L),
               info = paste("reserved names unexpectedly exported:",
                            paste(found, collapse = ", ")))
})
