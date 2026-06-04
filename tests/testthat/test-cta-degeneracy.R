###############################################################################
# test-cta-degeneracy.R
#
# Regression: CTA must never accept a degenerate all-same-class endpoint tree.
#
# CLAUDE.md canon: "CTA never produces degenerate trees: a candidate split that
# predicts only one class is ineligible at the node gate. A CTA or recursive
# CTA (LORT) result in which all terminal endpoints predict the same class is
# not a valid tree -- it is a model failure and must surface as no_tree."
#
# Two layers of protection are tested:
#
# B-1 / B-2  -- public no_tree path: no admissible root split is constructible
#   (mindenom guard or constant predictor).  Result: no_tree = TRUE.
#
# B-3 / B-4  -- internal degeneracy guard: .cta_predictions_degenerate().
#   This guard fires post-pruning inside the ENUMERATE expanded candidate loop
#   (cta_core.R line ~1469) and in the stump phase (line ~1541).  It rejects
#   any candidate tree whose terminal prediction vector collapses to a single
#   class after pruning.
#
#   A natural public-entrypoint repro of the post-pruning path would require
#   an expanded candidate whose sub-split MC p falls inside a narrow Sidak
#   window [alpha_k, alpha_split), making it simultaneously MC-significant and
#   Sidak-prunable.  That is too stochastic for a CRAN-safe regression.
#   The guard is therefore tested deterministically via the internal function.
#
# All paths: CRAN-safe, no skip gate, no MC, fast tier.
###############################################################################

# =============================================================================
# B-1  n < 2*mindenom: no cut can satisfy MINDENOM on both sides
# =============================================================================
# n=3, mindenom=2: cut at x<=1 leaves left with 1 obs < 2; cut at x<=2 leaves
# right with 1 obs < 2.  Both candidates rejected by the mindenom guard in
# oda_univariate_core (line 1066).  oda_fit returns ok=FALSE for attribute x.
# fast_screen produces no candidates -> root_cands is empty -> no_tree = TRUE.

test_that("B-1: n < 2*mindenom returns no_tree (degenerate tree rejected)", {
  X   <- data.frame(x = 1:3)
  y   <- c(0L, 1L, 0L)
  fit <- suppressMessages(
    oda_cta_fit(X, y, mindenom = 2L, mc_iter = 200L, mc_seed = 1L, loo = "off")
  )

  expect_true(isTRUE(fit$no_tree))
  # Single root leaf: attribute is NA_character_ (leaf sentinel), not a split.
  expect_length(fit$nodes[!vapply(fit$nodes, is.null, logical(1L))], 1L)
  expect_true(isTRUE(fit$nodes[[1L]]$leaf))
  expect_true(is.na(fit$nodes[[1L]]$attribute))
})

# =============================================================================
# B-2  Constant predictor
# =============================================================================
# y has both classes; x is constant (all 5).
# oda_make_blocks_ordered returns m = 1 block (one unique x value) ->
# oda_univariate_core returns ok=FALSE ("no_blocks") for this attribute.
# fast_screen produces no candidates -> root_cands is empty -> no_tree = TRUE.

test_that("B-2: constant x returns no_tree (no admissible cut, degenerate tree rejected)", {
  X   <- data.frame(x = rep(5, 8))
  y   <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  fit <- suppressMessages(
    oda_cta_fit(X, y, mindenom = 1L, mc_iter = 200L, mc_seed = 1L, loo = "off")
  )

  expect_true(isTRUE(fit$no_tree))
  expect_length(fit$nodes[!vapply(fit$nodes, is.null, logical(1L))], 1L)
  expect_true(isTRUE(fit$nodes[[1L]]$leaf))
  expect_true(is.na(fit$nodes[[1L]]$attribute))
})

# =============================================================================
# B-3  .cta_predictions_degenerate(): TRUE when all non-NA preds are one class
# B-4  .cta_predictions_degenerate(): FALSE when preds contain two classes
# =============================================================================
# This guard (cta_core.R) is called in two places:
#   1. Expanded ENUMERATE post-prune: rejects candidate if pruning collapses
#      all terminal endpoints to one class (e.g. sl_A artifact for direction
#      "1->0" causes reversed leaf pred_class assignments; pruning a sub-split
#      back to its local majority can yield all-same-class terminals).
#   2. Stump phase: defensive guard before stump is recorded.
#
# A public MC-based repro of path 1 requires the sub-split MC p to land in
# [alpha_k, alpha_split) -- a Sidak window too narrow for a CRAN-safe fixture.
# The guard is tested directly via the internal function.

test_that("B-3: .cta_predictions_degenerate() returns TRUE for all-same-class predictions", {
  degen <- odacore:::.cta_predictions_degenerate

  # All class 1
  expect_true(degen(c(1L, 1L, 1L)))
  # All class 0
  expect_true(degen(c(0L, 0L, 0L)))
  # Single observation
  expect_true(degen(c(1L)))
  # NAs excluded; remaining obs all same class
  expect_true(degen(c(1L, NA_integer_, 1L)))
  # All NA: no non-NA preds -> 0 distinct classes < 2 -> degenerate
  expect_true(degen(c(NA_integer_, NA_integer_)))
})

test_that("B-4: .cta_predictions_degenerate() returns FALSE when two classes present", {
  degen <- odacore:::.cta_predictions_degenerate

  # Two classes
  expect_false(degen(c(0L, 1L, 0L)))
  # Two classes with NAs interspersed
  expect_false(degen(c(0L, NA_integer_, 1L)))
  # Longer vector
  expect_false(degen(c(1L, 1L, 0L, 1L, 0L)))
})
