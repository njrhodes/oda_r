###############################################################################
# test-cta-degeneracy.R
#
# Regression: CTA must never accept a degenerate all-same-class endpoint tree.
# When no non-degenerate tree is constructible, the result must be no_tree=TRUE
# (single leaf node, no split).
#
# CLAUDE.md canon: "CTA never produces degenerate trees: a candidate split that
# predicts only one class is ineligible at the node gate. A CTA or recursive
# CTA (LORT) result in which all terminal endpoints predict the same class is
# not a valid tree -- it is a model failure and must surface as no_tree."
#
# Two deterministic paths to no_tree are exercised:
#   B-1  n < 2*mindenom  -> no admissible ordered cut (mindenom guard)
#   B-2  constant x      -> no admissible ordered cut (no_blocks guard)
# Both paths are CRAN-safe, no skip gate, no MC, fast tier.
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
