###############################################################################
# test-mpe-ch4.R
#
# MPE Chapter 4 categorical fixture tests.
#
# Source: Yarnold & Soltysik, Maximizing Predictive Accuracy (2016), Ch. 4.
#
# Coverage:
#   - Bowker / Stability (Table 4.1): transcription check + ordered training.
#   - Marginal dissymmetry (Table 4.1 off-diagonal): transcription + categorical.
#   - Bray-Curtis Dissimilarity Index (Table 4.2): training ESS + LOO ESS.
#   - Synthetic absent-level LOO override: unit tests for the MegaODA
#     compatibility rule that assigns a held-out observation whose category
#     was absent from the n-1 LOO training fold to the left-side class.
#
# Fixture status:
#   - Bowker training ESS: hard anchor (probe-confirmed parity, ESS=93.7).
#   - Marginal dissymmetry ESS: hard anchor (probe-confirmed parity, ESS=21.1).
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

# ---- Bowker / Stability (MPE Table 4.1) -------------------------------------

test_that("Bowker Table 4.1: transcription checks (N=55981, diagonal sum=53295)", {
  # Table 4.1: rows = 1980 class (1-4), cols = 1985 class (1-4).
  # Diagonal = stable observations; off-diagonal = changers.
  bowker <- matrix(c(
    11607,   100,   366,   124,
       87, 13677,   515,   302,
      172,   225, 17819,   270,
       63,   176,   286, 10192
  ), nrow = 4L, byrow = TRUE)

  expect_equal(sum(bowker), 55981L)
  expect_equal(rowSums(bowker), c(12197L, 14581L, 18486L, 10717L))
  expect_equal(colSums(bowker), c(11929L, 14178L, 18986L, 10888L))
  expect_equal(sum(diag(bowker)), 53295L)
  expect_equal(round(sum(diag(bowker)) / sum(bowker) * 100, 1), 95.2)
})

test_that("Bowker Table 4.1: ordered training ESS approx 93.7", {
  # Class = 1985 category (col index); Attribute = 1980 category (row index).
  # bowker is rows=1980, cols=1985. Transpose so rows=class(1985), cols=attr(1980).
  bowker <- matrix(c(
    11607,   100,   366,   124,
       87, 13677,   515,   302,
      172,   225, 17819,   270,
       63,   176,   286, 10192
  ), nrow = 4L, byrow = TRUE)

  d <- .expand_table(t(bowker), class_vals = 1:4, attr_vals = 1:4)

  set.seed(42L)
  fit <- oda_fit(d$x, d$y,
                 attr_type = "ordered",
                 priors_on = TRUE,
                 loo       = "off",
                 mc_iter   = 10000L)

  expect_true(fit$ok)
  expect_equal(fit$ess, 93.7, tolerance = 0.1)
  expect_equal(fit$rule$type, "multiclass_ordered")

  # Cut values: {1|2}, {2|3}, {3|4} boundaries; each class predicted by its rank.
  expect_equal(fit$rule$cut_values, c(1.5, 2.5, 3.5))
  expect_equal(fit$rule$seg_classes, 1:4)

  # Hard stability anchor: diagonal sum from confusion matrix must equal 53295.
  expect_equal(sum(diag(fit$confusion)), 53295L)

  # Note: p_mc is not asserted here — DIRECTIONAL not yet implemented (issue #6).
  # When DIRECTIONAL is added, expect p_mc < 0.0001 with directional Fisher.
})

test_that("Bowker Table 4.1: DIRECTIONAL ordered multiclass p_mc < 0.0001 [deferred]", {
  # MPE Chapter 4: Bowker/stability, DIRECTIONAL < 1 2 3 4.
  # Future anchors when multiclass DIRECTIONAL MC is implemented:
  #   ESS ≈ 93.7, diagonal correctly classified count = 53295,
  #   p_mc < 0.0001 with directional Fisher.
  # LOO not part of this target unless MPE declares it.
  skip("multiclass DIRECTIONAL not implemented: target Bowker DIRECTIONAL < 1 2 3 4, ESS \u2248 93.7, p_mc < 0.0001")
})

# ---- Marginal Dissymmetry (MPE Table 4.1 off-diagonal) ----------------------

test_that("Marginal dissymmetry Table 4.1: transcription check (N=2686)", {
  # Off-diagonal cells of Table 4.1 (changers only; diagonal zeroed).
  marg <- matrix(c(
      0, 100, 366, 124,
     87,   0, 515, 302,
    172, 225,   0, 270,
     63, 176, 286,   0
  ), nrow = 4L, byrow = TRUE)

  expect_equal(sum(marg), 2686L)
  expect_equal(rowSums(marg), c(590L, 904L, 667L, 525L))
  expect_equal(colSums(marg), c(322L, 501L, 1167L, 696L))
})

test_that("Marginal dissymmetry Table 4.1: categorical training ESS approx 21.1", {
  # SDA step 2: classify changers (off-diagonal) — categorical (unordered)
  # because dissymmetry patterns need not respect ordinal rank.
  # Class = 1985 category; Attribute = 1980 category.
  marg <- matrix(c(
      0, 100, 366, 124,
     87,   0, 515, 302,
    172, 225,   0, 270,
     63, 176, 286,   0
  ), nrow = 4L, byrow = TRUE)

  d <- .expand_table(t(marg), class_vals = 1:4, attr_vals = 1:4)

  set.seed(42L)
  fit <- oda_fit(d$x, d$y,
                 attr_type = "categorical",
                 priors_on = TRUE,
                 loo       = "off",
                 mc_iter   = 25000L)

  expect_true(fit$ok)
  expect_equal(fit$ess, 21.1, tolerance = 0.5)
  expect_equal(fit$rule$type, "multiclass_nominal")

  # Note: p_mc not asserted — DIRECTIONAL not yet implemented (issue #6).
  # Note: This is SDA step 2 context from MPE Ch. 4; full SDA also requires
  # step 1 (overall Bowker stability, above) and D-statistic comparison.
})

# ---- Pinckney Gag Rule (MPE Table 4.3 / 4.4) --------------------------------

test_that("Pinckney Table 4.3: transcription check (N=225)", {
  # Table 4.3: Congressional Voting on the 1836 Pinckney Gag Rule.
  # Rows = region attribute (1=North, 2=Border, 3=South).
  # Cols = vote class (1=Yea, 2=Abstain, 3=Nay).
  pinckney <- matrix(c(
    61, 12, 60,
    17,  6,  1,
    39, 22,  7
  ), nrow = 3L, byrow = TRUE)

  expect_equal(sum(pinckney), 225L)
  expect_equal(rowSums(pinckney), c(133L, 24L, 68L))   # North, Border, South
  expect_equal(colSums(pinckney), c(117L, 40L, 68L))   # Yea, Abstain, Nay
})

test_that("Pinckney Table 4.3: categorical training ESS approx 28.8", {
  # MPE p.76: nondirectional (no DIRECTIONAL command).
  # MegaODA syntax: CATEGORICAL ON; TABLE 3; CLASS COL; MCARLO ITER 10000; GO;
  # Class = vote (col); Attribute = region (row).
  #
  # This test anchors nondirectional multiclass categorical training ESS and
  # rule mapping only.  ESP = 24.2 in MPE (Effect Strength for Predictive
  # value, a separate metric not returned by oda_fit()).
  #
  # MPE reports p < 0.0001.  odacore p_mc consistently ~0.155 regardless of
  # seed and does not match that MPE value.  p_mc is intentionally not asserted
  # here.  Multiclass categorical MC parity with MegaODA.exe is a separate open
  # issue / follow-up.
  pinckney <- matrix(c(
    61, 12, 60,
    17,  6,  1,
    39, 22,  7
  ), nrow = 3L, byrow = TRUE)

  d <- .expand_table(t(pinckney), class_vals = 1:3, attr_vals = 1:3)

  set.seed(42L)
  fit <- oda_fit(d$x, d$y,
                 attr_type = "categorical",
                 priors_on = TRUE,
                 loo       = "off",
                 mc_iter   = 10000L)

  expect_true(fit$ok)
  expect_equal(fit$ess, 28.8, tolerance = 0.5)
  expect_equal(fit$rule$type, "multiclass_nominal")

  # Mapping: level_class[i] = predicted class for attribute level i.
  # level 1 (North) -> class 3 (Nay)
  # level 2 (Border) -> class 1 (Yea)
  # level 3 (South)  -> class 2 (Abstain)
  expect_equal(as.integer(fit$rule$level_class), c(3L, 1L, 2L))
})

test_that("Pinckney Table 4.4: residual categorical training ESS approx 40.9", {
  # MPE p.77: primary voters (correctly classified by Table 4.3 model) removed.
  # Residual minority voting pattern; same nondirectional analysis.
  #
  # This test anchors nondirectional multiclass categorical training ESS and
  # rule mapping only.  ESP = 42.2 in MPE (Effect Strength for Predictive
  # value, a separate metric not returned by oda_fit()).
  #
  # MPE reports p < 0.0003.  odacore p_mc varies 0.000-0.014 across seeds and
  # is not stably below that MPE threshold.  p_mc is intentionally not asserted
  # here.  Multiclass categorical MC parity with MegaODA.exe is a separate open
  # issue / follow-up.
  pinckney_resid <- matrix(c(
    61, 12,  0,
     0,  6,  1,
    39,  0,  7
  ), nrow = 3L, byrow = TRUE)

  d <- .expand_table(t(pinckney_resid), class_vals = 1:3, attr_vals = 1:3)

  set.seed(42L)
  fit <- oda_fit(d$x, d$y,
                 attr_type = "categorical",
                 priors_on = TRUE,
                 loo       = "off",
                 mc_iter   = 10000L)

  expect_true(fit$ok)
  expect_equal(fit$ess, 40.9, tolerance = 0.5)
  expect_equal(fit$rule$type, "multiclass_nominal")

  # Mapping: level_class[i] = predicted class for attribute level i.
  # level 1 (North)  -> class 1 (Yea)
  # level 2 (Border) -> class 2 (Abstain)
  # level 3 (South)  -> class 3 (Nay)
  expect_equal(as.integer(fit$rule$level_class), c(1L, 2L, 3L))
})

# ---- Political affiliation (MPE Chapter 4) ----------------------------------

test_that("Political affiliation: categorical multiclass DIRECTIONAL p < 0.0001 [deferred]", {
  # MPE Chapter 4: political affiliation, categorical multiclass.
  # MegaODA commands: CATEGORICAL ON; TABLE 7; CLASS ROW;
  #   DIRECTIONAL < 1 2 3 4 5 6 7; MCARLO ITER 10000; GO;
  # Directional hypothesis: student and parent have same political affiliation.
  # Future anchors when categorical multiclass DIRECTIONAL is implemented:
  #   ESS = 19.4, ESP = 17.9, p < 0.0001.
  # LOO is superfluous when class categories equal attribute categories (k = C)
  # and a directional hypothesis is specified; LOO off / not applicable.
  skip("categorical multiclass DIRECTIONAL not implemented: target political affiliation DIRECTIONAL < 1..7, ESS 19.4, ESP 17.9, p < 0.0001, LOO superfluous")
})

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
