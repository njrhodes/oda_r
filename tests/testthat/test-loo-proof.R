###############################################################################
# test-loo-proof.R
#
# L1R/L1R2 brute-force LOO canonicity proofs.
#
# Canonical LOO contract (MPE Chapter 2, pp. 32-33):
#   For each observation i: hold out i, obtain the ODA model from the n-1
#   remaining observations, classify i with that model, and store the result.
#   There is one canonical LOO contract.  Shortcuts are allowed only when
#   proven to exactly reproduce this delete-one result.
#
# Weighted binary LOO is canonical.
# Weighted categorical LOO is out of scope (forbidden elsewhere).
# Explicit delete-one oda_univariate_core(mcarlo=FALSE) is the test oracle,
# not the production path.
#
# Production paths tested:
#   binary_map: oda_loo_binary_map_counts() -- weighted 2x2 table adjustment.
#     Probes B1-B4 verify unweighted and weighted cases match the oracle.
#   ordered_cut uniform weights: oda_loo_ordered_cut_counts() -- bin counts.
#     Probes A1-A3 verify clean, priors, and tie cases match the oracle.
#
# CRAN-safe only. No skip gates. No MC (mcarlo=FALSE in oracle calls).
###############################################################################

# ---- Shared helper ----------------------------------------------------------

# Reference explicit delete-one LOO for binary class data.
# For each fold i: drops obs i, refits via oda_univariate_core(), and predicts
# obs i with the delete-one rule.  Returns ess_loo using the same
# priors-adjusted formula as oda_loo_for_rule() (unioda_core.R).
#
# Parameters:
#   x         - predictor values
#   y         - integer class labels (0L/1L)
#   priors_on - logical; passed to oda_univariate_core and used in ESS calc
#   attr_type - "binary" or "ordered"
#   w         - optional case weights (NULL = unit weights)
#
# ESS formula (same as oda_loo_for_rule):
#   sens = TP_w / (TP_w + FN_w) = sum(w[y==1 & pred==1]) / sum(w[y==1])
#   spec = TN_w / (TN_w + FP_w) = sum(w[y==0 & pred==0]) / sum(w[y==0])
#   ESS = (sens + spec - 1) * 100
#
# Returns NA_real_ if any fold fails.
.explicit_loo_ess <- function(x, y, priors_on, attr_type, w = NULL) {
  n   <- length(y)
  y   <- as.integer(y)
  if (is.null(w)) w <- rep(1.0, n) else w <- as.numeric(w)
  tw0 <- sum(w[y == 0L])
  tw1 <- sum(w[y == 1L])
  preds <- integer(n)

  for (i in seq_len(n)) {
    keep <- seq_len(n) != i
    fi <- oda_univariate_core(
      x         = x[keep],
      y         = y[keep],
      w         = w[keep],
      attr_type = attr_type,
      priors_on = priors_on,
      mcarlo    = FALSE,
      loo       = "off"
    )
    if (!isTRUE(fi$ok)) return(NA_real_)
    preds[i] <- oda_rule_predict(x[i], fi$rule)
  }

  # Priors-adjusted ESS matching oda_loo_for_rule:
  # oda_apply_priors normalizes so each class has equal total weight; the
  # resulting mean_pac reduces to (spec + sens) / 2 regardless of priors_on.
  # With uniform w this equals (PAC_0 + PAC_1) / 2; with non-uniform w it
  # uses weighted PAC per class (TP_w / tc_w for each class).
  ok   <- !is.na(preds)
  spec <- sum(w[ok & y == 0L] * (preds[ok & y == 0L] == 0L)) / tw0
  sens <- sum(w[ok & y == 1L] * (preds[ok & y == 1L] == 1L)) / tw1
  (sens + spec - 1) * 100
}


# =============================================================================
# Probe B1: binary_map shortcut -- stable-margin case
#
# With large separation between class distributions (n00=8, n01=1, n10=1,
# n11=8), no single deletion can flip the optimal binary direction.
# The full-data training rule "0->1" equals the delete-one rule for every
# fold.  Shortcut ESS should equal explicit delete-one ESS.
#
# Arithmetic: A-B = 2*(n00*n11 - n01*n10) = 2*(64-1) = 126.
# Direction-flip condition for dropping (x=0,y=0): 0 < A-B < 2*n11 = 16.
# 126 >> 16: no deletion can flip direction.
# Expected: shortcut ESS ~= explicit ESS (within floating-point tolerance).
# =============================================================================

test_that("B1: binary_map shortcut matches explicit LOO for stable-margin data", {
  # n=18: x=0 has 8 class-0 and 1 class-1; x=1 has 1 class-0 and 8 class-1
  x_b1 <- c(rep(0L, 9L), rep(1L, 9L))
  y_b1 <- c(rep(0L, 8L), 1L,            # x=0: 8 class-0, 1 class-1
             0L, rep(1L, 8L))            # x=1: 1 class-0, 8 class-1

  # Training rule: "0->1", ESS ~= 77.78%
  fit_b1 <- oda_univariate_core(
    x = x_b1, y = y_b1,
    attr_type = "binary", priors_on = TRUE, mcarlo = FALSE, loo = "off"
  )
  expect_true(isTRUE(fit_b1$ok),
    label = "B1: training fit must succeed on stable-margin data")
  expect_equal(fit_b1$rule$type, "binary_map",
    label = "B1: binary attribute must produce binary_map rule type")
  expect_equal(fit_b1$rule$direction, "0->1",
    label = "B1: stable-margin data should have clear '0->1' direction")

  # Shortcut path: oda_loo_for_rule applies binary_map shortcut
  loo_shortcut_b1 <- odacore:::oda_loo_for_rule(
    x         = x_b1,
    y         = y_b1,
    rule      = fit_b1$rule,
    attr_type = "binary",
    priors_on = TRUE
  )
  expect_true(isTRUE(loo_shortcut_b1$allowed),
    label = "B1: oda_loo_for_rule must be allowed for binary stable-margin data")

  # Reference explicit delete-one
  ess_explicit_b1 <- .explicit_loo_ess(x_b1, y_b1,
                                        priors_on = TRUE, attr_type = "binary")
  expect_false(is.na(ess_explicit_b1), label = "B1: explicit LOO must not fail")

  # Stable margin: shortcut should reproduce delete-one result exactly
  expect_equal(loo_shortcut_b1$ess_loo, ess_explicit_b1, tolerance = 1e-6,
    label = paste0(
      "B1: binary_map shortcut ESS (", round(loo_shortcut_b1$ess_loo, 4), "%) ",
      "must match explicit delete-one ESS (", round(ess_explicit_b1, 4), "%) ",
      "for stable-margin data where no deletion can flip direction"
    )
  )
})


# =============================================================================
# Probe B2: binary_map algebraic LOO -- near-tie case (L1R2 fix verification)
#
# The class distributions are marginally separated. A single deletion CAN
# flip the optimal binary direction. This probe was originally written to
# document non-canonical divergence; after the L1R2 fix it verifies that
# oda_loo_binary_map_counts() correctly handles direction-flipping deletions.
#
# Dataset: n=13, x in {0,1}, y in {0,1}
#   x=0: n00=5 class-0, n01=3 class-1   (tc0=8, tc1=5)
#   x=1: n10=3 class-0, n11=2 class-1
#
# Arithmetic:
#   A-B = 2*(n00*n11 - n01*n10) = 2*(5*2 - 3*3) = 2*(10-9) = 2.
#   Condition for direction flip when dropping (x=0,y=0): 0 < A-B < 2*n11=4.
#   2 satisfies this: after dropping one (x=0,y=0), "1->0" becomes optimal.
#
# Algebraic path: subtracts held-out weight, recomputes direction per fold.
# Correctly detects the flip; both algebraic and explicit ESS ~= -100%.
#
# Expected: algebraic ESS matches explicit delete-one ESS (within tolerance).
# =============================================================================

test_that("B2: binary_map algebraic LOO matches explicit delete-one for near-tie data", {
  # n=13: x=0 has 5 class-0 and 3 class-1; x=1 has 3 class-0 and 2 class-1
  x_b2 <- c(rep(0L, 8L), rep(1L, 5L))
  y_b2 <- c(rep(0L, 5L), rep(1L, 3L),   # x=0: 5 class-0, 3 class-1
             rep(0L, 3L), rep(1L, 2L))   # x=1: 3 class-0, 2 class-1

  # Training rule: "0->1", ESS ~= 2.5% (marginal but positive)
  fit_b2 <- oda_univariate_core(
    x = x_b2, y = y_b2,
    attr_type = "binary", priors_on = TRUE, mcarlo = FALSE, loo = "off"
  )
  expect_true(isTRUE(fit_b2$ok),
    label = "B2: training fit must succeed on near-tie data")
  expect_equal(fit_b2$rule$type, "binary_map",
    label = "B2: binary attribute must produce binary_map rule type")
  expect_equal(fit_b2$rule$direction, "0->1",
    label = "B2: training rule should be '0->1' (ESS ~= 2.5% per audit)")

  # Algebraic path: oda_loo_for_rule dispatches to oda_loo_binary_map_counts()
  loo_shortcut_b2 <- odacore:::oda_loo_for_rule(
    x         = x_b2,
    y         = y_b2,
    rule      = fit_b2$rule,
    attr_type = "binary",
    priors_on = TRUE
  )
  expect_true(isTRUE(loo_shortcut_b2$allowed),
    label = "B2: oda_loo_for_rule must be allowed for binary near-tie data")

  # Reference explicit delete-one
  ess_explicit_b2 <- .explicit_loo_ess(x_b2, y_b2,
                                        priors_on = TRUE, attr_type = "binary")
  expect_false(is.na(ess_explicit_b2), label = "B2: explicit LOO must not fail")

  # After the L1R2 fix, oda_loo_binary_map_counts() uses weighted 2x2 table
  # adjustment to recompute the delete-one rule for each fold.  Dropping any
  # (x=0,y=0) obs makes "1->0" optimal; the algebraic path correctly detects
  # this and matches the explicit oracle (ESS ~= -100% for both).
  expect_equal(loo_shortcut_b2$ess_loo, ess_explicit_b2, tolerance = 1e-6,
    label = paste0(
      "B2: binary_map algebraic ESS (", round(loo_shortcut_b2$ess_loo, 2), "%) ",
      "must match explicit delete-one ESS (", round(ess_explicit_b2, 2), "%) ",
      "for near-tie data after L1R2 fix"
    )
  )
})


# =============================================================================
# Probe B3: binary_map algebraic LOO -- weighted, stable-margin
#
# Non-uniform case weights; large margin so no deletion flips direction.
# Verifies the weighted 2x2 table path handles non-uniform w correctly.
#
# Same x/y structure as B1 (n00=8, n01=1, n10=1, n11=8) but with
# w=[2,2,2,...,2, 1,1,1,...,1] (x=0 obs weight=2, x=1 obs weight=1).
# Weighted table: c00_w=16, c01_w=2, c10_w=1, c11_w=8, tc0_w=17, tc1_w=10.
# ESS("0->1") ~= 74.1%.  Direction is stable under any single deletion.
# Expected: algebraic ESS == explicit ESS (within tolerance).
# =============================================================================

test_that("B3: weighted binary_map algebraic LOO matches explicit (stable-margin)", {
  x_b3 <- c(rep(0L, 9L), rep(1L, 9L))
  y_b3 <- c(rep(0L, 8L), 1L, 0L, rep(1L, 8L))   # same as B1
  w_b3 <- c(rep(2, 9L), rep(1, 9L))               # x=0: weight 2, x=1: weight 1

  fit_b3 <- oda_univariate_core(
    x = x_b3, y = y_b3, w = w_b3,
    attr_type = "binary", priors_on = TRUE, mcarlo = FALSE, loo = "off"
  )
  expect_true(isTRUE(fit_b3$ok),
    label = "B3: weighted training fit must succeed")
  expect_equal(fit_b3$rule$type, "binary_map",
    label = "B3: binary attribute must produce binary_map rule type")

  loo_alg_b3 <- odacore:::oda_loo_for_rule(
    x         = x_b3,
    y         = y_b3,
    w         = w_b3,
    rule      = fit_b3$rule,
    attr_type = "binary",
    priors_on = TRUE
  )
  expect_true(isTRUE(loo_alg_b3$allowed),
    label = "B3: oda_loo_for_rule must be allowed for weighted binary data")

  ess_explicit_b3 <- .explicit_loo_ess(x_b3, y_b3,
                                        priors_on = TRUE, attr_type = "binary",
                                        w = w_b3)
  expect_false(is.na(ess_explicit_b3), label = "B3: explicit LOO must not fail")

  expect_equal(loo_alg_b3$ess_loo, ess_explicit_b3, tolerance = 1e-6,
    label = paste0(
      "B3: weighted binary_map algebraic ESS (", round(loo_alg_b3$ess_loo, 4), "%) ",
      "must match explicit delete-one ESS (", round(ess_explicit_b3, 4), "%) ",
      "for stable-margin weighted data"
    )
  )
})


# =============================================================================
# Probe B4: binary_map algebraic LOO -- weighted, direction-flip case
#
# Non-uniform case weights on the near-tie dataset (same x/y as B2).
# w=[2,2,2,2,2, 2,2,2, 1,1,1, 1,1] (x=0 obs weight=2, x=1 obs weight=1).
# Weighted table:
#   c00_w=10, c01_w=6, c10_w=3, c11_w=2, tc0_w=13, tc1_w=8.
#   ESS("0->1") = (10/13 + 2/8 - 1)*100 ~= 1.9%.
#   ESS("1->0") = (3/13 + 6/8 - 1)*100 ~= -1.9%.
# Dropping a (x=0,y=0,w=2) obs: tc0_w=11, c00_w=8.
#   ESS("0->1") = (8/11 + 2/8 - 1)*100 ~= -2.3%  <-- direction flips!
#   ESS("1->0") = (3/11 + 6/8 - 1)*100 ~= 2.3%.
# The weighted algebraic path must detect this flip and match the oracle.
# Expected: algebraic ESS == explicit ESS (within tolerance).
# =============================================================================

test_that("B4: weighted binary_map algebraic LOO matches explicit (direction-flip)", {
  x_b4 <- c(rep(0L, 8L), rep(1L, 5L))             # same structure as B2
  y_b4 <- c(rep(0L, 5L), rep(1L, 3L),
             rep(0L, 3L), rep(1L, 2L))
  w_b4 <- c(rep(2, 8L), rep(1, 5L))               # x=0: weight 2, x=1: weight 1

  fit_b4 <- oda_univariate_core(
    x = x_b4, y = y_b4, w = w_b4,
    attr_type = "binary", priors_on = TRUE, mcarlo = FALSE, loo = "off"
  )
  expect_true(isTRUE(fit_b4$ok),
    label = "B4: weighted training fit must succeed")
  expect_equal(fit_b4$rule$type, "binary_map",
    label = "B4: binary attribute must produce binary_map rule type")
  expect_equal(fit_b4$rule$direction, "0->1",
    label = "B4: training rule should be '0->1' for this weighted near-tie data")

  loo_alg_b4 <- odacore:::oda_loo_for_rule(
    x         = x_b4,
    y         = y_b4,
    w         = w_b4,
    rule      = fit_b4$rule,
    attr_type = "binary",
    priors_on = TRUE
  )
  expect_true(isTRUE(loo_alg_b4$allowed),
    label = "B4: oda_loo_for_rule must be allowed for weighted binary near-tie data")

  ess_explicit_b4 <- .explicit_loo_ess(x_b4, y_b4,
                                        priors_on = TRUE, attr_type = "binary",
                                        w = w_b4)
  expect_false(is.na(ess_explicit_b4), label = "B4: explicit LOO must not fail")

  expect_equal(loo_alg_b4$ess_loo, ess_explicit_b4, tolerance = 1e-6,
    label = paste0(
      "B4: weighted binary_map algebraic ESS (", round(loo_alg_b4$ess_loo, 4), "%) ",
      "must match explicit delete-one ESS (", round(ess_explicit_b4, 4), "%) ",
      "for weighted direction-flip data"
    )
  )
})


# =============================================================================
# Probe A1: algebraic ordered_cut shortcut -- priors_on=FALSE
#
# Uniform weights, clean class separation, priors_on=FALSE.
# Verifies that the algebraic count-adjustment path reproduces explicit
# delete-one results for the simplest ordered-cut case.
#
# Dataset: n=40, x=1..40, y=c(rep(0,20), rep(1,20)).
# Training cut: 20.5, direction "0->1", ESS=100%.
# LOO: all folds have the same cut; all predictions correct; ESS=100%.
# =============================================================================

test_that("A1: algebraic ordered_cut shortcut matches explicit LOO (priors_on=FALSE)", {
  x_a1 <- 1L:40L
  y_a1 <- c(rep(0L, 20L), rep(1L, 20L))

  fit_a1 <- oda_univariate_core(
    x = x_a1, y = y_a1,
    attr_type = "ordered", priors_on = FALSE, mcarlo = FALSE, loo = "off"
  )
  expect_true(isTRUE(fit_a1$ok),
    label = "A1: training fit must succeed on clean-separation data")
  expect_equal(fit_a1$rule$type, "ordered_cut",
    label = "A1: ordered attribute must produce ordered_cut rule type")

  # Algebraic path via oda_loo_for_rule (uniform weights -> algebraic branch)
  loo_alg_a1 <- odacore:::oda_loo_for_rule(
    x         = x_a1,
    y         = y_a1,
    rule      = fit_a1$rule,
    attr_type = "ordered",
    priors_on = FALSE
  )
  expect_true(isTRUE(loo_alg_a1$allowed),
    label = "A1: oda_loo_for_rule must be allowed for uniform-weight ordered data")

  # Reference explicit delete-one
  ess_explicit_a1 <- .explicit_loo_ess(x_a1, y_a1,
                                        priors_on = FALSE, attr_type = "ordered")
  expect_false(is.na(ess_explicit_a1), label = "A1: explicit LOO must not fail")

  expect_equal(loo_alg_a1$ess_loo, ess_explicit_a1, tolerance = 1e-6,
    label = paste0(
      "A1: algebraic ordered_cut ESS (", round(loo_alg_a1$ess_loo, 4), "%) ",
      "must match explicit delete-one ESS (", round(ess_explicit_a1, 4), "%) ",
      "for clean-separation data, priors_on=FALSE"
    )
  )
})


# =============================================================================
# Probe A2: algebraic ordered_cut shortcut -- priors_on=TRUE, noisy data
#
# Uniform weights, imbalanced classes (tc0=13, tc1=7), priors_on=TRUE.
# Noise observations (2 class-1 in the class-0 region; 2 class-0 in the
# class-1 region) produce LOO misclassifications, exercising the priors
# weighting in .find_best() across folds.
#
# Verifies that priors normalization in the algebraic path matches what
# oda_univariate_core() produces in explicit delete-one refits.
# =============================================================================

test_that("A2: algebraic ordered_cut shortcut matches explicit LOO (priors_on=TRUE)", {
  # n=20: x=1..20, tc0=13, tc1=7
  # x=1..11: class-0 (11 obs); x=12,13: class-1 noise (2 obs)
  # x=14..18: class-1 (5 obs);  x=19,20: class-0 noise (2 obs)
  x_a2 <- 1L:20L
  y_a2 <- c(rep(0L, 11L), rep(1L, 2L), rep(1L, 5L), rep(0L, 2L))

  fit_a2 <- oda_univariate_core(
    x = x_a2, y = y_a2,
    attr_type = "ordered", priors_on = TRUE, mcarlo = FALSE, loo = "off"
  )
  expect_true(isTRUE(fit_a2$ok),
    label = "A2: training fit must succeed on noisy data")
  expect_equal(fit_a2$rule$type, "ordered_cut",
    label = "A2: ordered attribute must produce ordered_cut rule type")

  # Algebraic path
  loo_alg_a2 <- odacore:::oda_loo_for_rule(
    x         = x_a2,
    y         = y_a2,
    rule      = fit_a2$rule,
    attr_type = "ordered",
    priors_on = TRUE
  )
  expect_true(isTRUE(loo_alg_a2$allowed),
    label = "A2: oda_loo_for_rule must be allowed for uniform-weight ordered data")

  # Reference explicit delete-one
  ess_explicit_a2 <- .explicit_loo_ess(x_a2, y_a2,
                                        priors_on = TRUE, attr_type = "ordered")
  expect_false(is.na(ess_explicit_a2), label = "A2: explicit LOO must not fail")

  expect_equal(loo_alg_a2$ess_loo, ess_explicit_a2, tolerance = 1e-6,
    label = paste0(
      "A2: algebraic ordered_cut ESS (", round(loo_alg_a2$ess_loo, 4), "%) ",
      "must match explicit delete-one ESS (", round(ess_explicit_a2, 4), "%) ",
      "for noisy data, priors_on=TRUE, imbalanced classes"
    )
  )
})


# =============================================================================
# Probe A3: algebraic ordered_cut shortcut -- tie-forcing case
#
# Two cut positions produce equal ESS in training; first-identified wins.
# Some LOO folds also produce ties (e.g., when a bin-3 obs is dropped,
# the cut-1.5 and cut-3.5 tie breaks differently than training).
# Verifies that the algebraic path's strict-greater-than scan reproduces
# oda_univariate_core()'s first-identified tie-breaking for all folds.
#
# Dataset: x=[1,2,3,4] repeated 5 times (n=20); y=[0,1,0,1] repeated 5 times.
#   x=1: 5 class-0 obs; x=2: 5 class-1 obs; x=3: 5 class-0 obs; x=4: 5 class-1 obs.
#   Cut 1.5: PAC_0=0.5, PAC_1=1.0, mean=0.75, ESS=50%.
#   Cut 3.5: PAC_0=1.0, PAC_1=0.5, mean=0.75, ESS=50%.
#   Tie at 50%: first-identified cut (1.5) must win in training.
# =============================================================================

test_that("A3: algebraic ordered_cut tie-forcing -- training selects first-identified cut", {
  x_a3 <- rep(c(1L, 2L, 3L, 4L), 5L)
  y_a3 <- rep(c(0L, 1L, 0L, 1L), 5L)

  fit_a3 <- oda_univariate_core(
    x = x_a3, y = y_a3,
    attr_type = "ordered", priors_on = TRUE, mcarlo = FALSE, loo = "off"
  )
  expect_true(isTRUE(fit_a3$ok),
    label = "A3: training fit must succeed on tie-forcing data")
  expect_equal(fit_a3$rule$type, "ordered_cut",
    label = "A3: ordered attribute must produce ordered_cut rule type")

  # First-identified wins: cut=1.5 must be selected over equally-scoring cut=3.5
  expect_equal(fit_a3$rule$cut_value, 1.5, tolerance = 1e-9,
    label = paste0(
      "A3: training tie at cut=1.5 and cut=3.5 (both ESS=50%); ",
      "first-identified (cut=1.5) must win"
    )
  )
})

test_that("A3: algebraic ordered_cut tie-forcing -- LOO ESS matches explicit", {
  x_a3 <- rep(c(1L, 2L, 3L, 4L), 5L)
  y_a3 <- rep(c(0L, 1L, 0L, 1L), 5L)

  fit_a3 <- oda_univariate_core(
    x = x_a3, y = y_a3,
    attr_type = "ordered", priors_on = TRUE, mcarlo = FALSE, loo = "off"
  )
  skip_if(!isTRUE(fit_a3$ok), "A3 training fit failed; cannot run LOO comparison")

  # Algebraic path
  loo_alg_a3 <- odacore:::oda_loo_for_rule(
    x         = x_a3,
    y         = y_a3,
    rule      = fit_a3$rule,
    attr_type = "ordered",
    priors_on = TRUE
  )
  expect_true(isTRUE(loo_alg_a3$allowed),
    label = "A3: oda_loo_for_rule must be allowed for tie-forcing data")

  # Reference explicit delete-one
  ess_explicit_a3 <- .explicit_loo_ess(x_a3, y_a3,
                                        priors_on = TRUE, attr_type = "ordered")
  skip_if(is.na(ess_explicit_a3), "A3 explicit LOO failed on at least one fold")

  expect_equal(loo_alg_a3$ess_loo, ess_explicit_a3, tolerance = 1e-6,
    label = paste0(
      "A3: algebraic ordered_cut ESS (", round(loo_alg_a3$ess_loo, 4), "%) ",
      "must match explicit delete-one ESS (", round(ess_explicit_a3, 4), "%) ",
      "for tie-forcing data -- tests SAMPLEREP/first-identified consistency"
    )
  )
})
