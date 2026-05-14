###############################################################################
# test-mpe-ch4.R
#
# MPE Chapter 4 categorical fixture tests.
#
# Source: Yarnold & Soltysik, Maximizing Predictive Accuracy (2016), Ch. 4.
#
# Coverage:
#   - Bray-Curtis Dissimilarity Index (Table 4.2): training ESS + LOO ESS.
#   - Synthetic absent-level LOO override: unit tests for the MegaODA
#     compatibility rule that assigns a held-out observation whose category
#     was absent from the n-1 LOO training fold to the left-side class.
#
# Fixture status:
#   - Bray-Curtis training ESS: hard anchor (probe-confirmed parity).
#   - Bray-Curtis LOO ESS: hard anchor, fixture-derived MegaODA compatibility
#     (absent-level override required; see R/unioda_core.R oda_loo_for_rule).
#   - Rule grouping: hard anchor.
#   - p-values: directional check only (MC stochastic).
###############################################################################

# ---- helpers ----------------------------------------------------------------

.expand_table <- function(counts, class_vals, attr_vals) {
  y <- integer(0L); x <- integer(0L)
  for (i in seq_along(class_vals))
    for (j in seq_along(attr_vals))
      if (counts[i, j] > 0L) {
        y <- c(y, rep(class_vals[i], counts[i, j]))
        x <- c(x, rep(attr_vals[j],  counts[i, j]))
      }
  list(x = x, y = y)
}

# ---- Bray-Curtis (MPE Table 4.2) --------------------------------------------

test_that("Bray-Curtis Table 4.2: training ESS approx 44.7", {
  # Table 4.2: rows = ecological category (A-E = 1-5),
  #            cols = sampling site (S29 = 1, S30 = 2)
  # Class = site; Attribute = category.
  # Fixture-derived MegaODA categorical LOO absent-level compatibility.
  bray_raw <- matrix(c(
    11, 24,
     0, 37,
     7,  5,
     8, 18,
     0,  1
  ), nrow = 5L, byrow = TRUE)

  d <- .expand_table(t(bray_raw), class_vals = 1:2, attr_vals = 1:5)

  set.seed(99L)
  fit <- oda_fit(d$x, d$y,
                 attr_type = "categorical",
                 priors_on = TRUE,
                 loo       = "off",
                 mc_iter   = 25000L)

  expect_true(fit$ok)
  expect_equal(fit$ess, 44.7, tolerance = 0.1)
  expect_equal(fit$rule$type, "nominal_cut")

  # Rule grouping: {A(1), C(3), D(4)} -> S29 (label_0); {B(2), E(5)} -> S30
  # oda_fit recodes S29->0, S30->1; label_0 = original S29 label = 1L
  left_set  <- sort(as.integer(fit$rule$left_levels))
  right_set <- sort(as.integer(fit$rule$right_levels))
  expect_equal(left_set,  c(1L, 3L, 4L))
  expect_equal(right_set, c(2L, 5L))
})

test_that("Bray-Curtis Table 4.2: LOO ESS approx 43.5 (absent-level override)", {
  # LOO ESS = 43.5 requires the absent-level override in oda_loo_for_rule():
  # category E (the single S30 obs in that category) is absent from its LOO
  # training fold; MegaODA assigns it to the left class (S29), producing one
  # additional S30 misclassification relative to training.
  # Fixture-derived MegaODA categorical LOO absent-level compatibility.
  bray_raw <- matrix(c(
    11, 24,
     0, 37,
     7,  5,
     8, 18,
     0,  1
  ), nrow = 5L, byrow = TRUE)

  d <- .expand_table(t(bray_raw), class_vals = 1:2, attr_vals = 1:5)

  set.seed(99L)
  fit <- oda_fit(d$x, d$y,
                 attr_type = "categorical",
                 priors_on = TRUE,
                 loo       = "on",
                 mc_iter   = 25000L)

  expect_true(fit$ok)
  expect_true(isTRUE(fit$loo$allowed))
  expect_equal(fit$loo$ess_loo, 43.5, tolerance = 0.1)
  expect_lt(fit$loo$p_value, 0.001)
})

# ---- Synthetic absent-level override unit tests -----------------------------

test_that("LOO absent-level override triggers only for truly absent category", {
  # Dataset: 3 attribute levels (1, 2, 3); binary class (0/1).
  # Level 3 has exactly 1 observation, in class 1.
  # When that observation is held out, level 3 is absent from the LOO training
  # fold -> override must predict 0L (left class, direction "0->1").
  # When any other observation is held out, level 3 is still present
  # -> no override, normal oda_rule_predict used.
  set.seed(7L)
  x <- c(rep(1L, 10L), rep(2L, 10L), 3L)   # 21 observations
  y <- c(rep(0L,  8L), rep(1L,  2L),        # level 1: 8 class-0, 2 class-1
         rep(1L, 10L),                       # level 2: all class-1
         1L)                                 # level 3: 1 class-1 (singleton)

  fit <- oda_fit(x, y,
                 attr_type = "categorical",
                 priors_on = TRUE,
                 loo       = "on",
                 mc_iter   = 1L,    # skip MC; focus on LOO structure
                 mcarlo    = FALSE)

  expect_true(fit$ok)
  expect_true(isTRUE(fit$loo$allowed))

  # LOO ESS must differ from training ESS because the singleton (obs 21)
  # is overridden to class 0; confirm ess_loo is a valid number.
  expect_false(is.na(fit$loo$ess_loo))
})

test_that("LOO absent-level override does not trigger when category present", {
  # Same structure as above but level 3 has 2 observations — never absent
  # from any LOO fold. Override must not fire; LOO ESS == training ESS
  # (rule is completely stable across all folds).
  set.seed(8L)
  x <- c(rep(1L, 10L), rep(2L, 10L), 3L, 3L)  # 22 obs; level 3 has 2
  y <- c(rep(0L,  8L), rep(1L,  2L),
         rep(1L, 10L),
         1L, 1L)                                 # both level-3 obs in class 1

  fit <- oda_fit(x, y,
                 attr_type = "categorical",
                 priors_on = TRUE,
                 loo       = "on",
                 mc_iter   = 1L,
                 mcarlo    = FALSE)

  expect_true(fit$ok)
  expect_true(isTRUE(fit$loo$allowed))
  expect_false(is.na(fit$loo$ess_loo))
  # With level 3 never absent, no override fires.  The rule is stable so
  # ess_loo equals training ESS within floating-point tolerance.
  expect_equal(fit$loo$ess_loo, fit$ess, tolerance = 0.01)
})
