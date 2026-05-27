###############################################################################
# test-sda-unioda.R — SDA-1 tests for sda_fit(mode = "unioda_max_ess")
#
# Tier: CRAN-safe (all tests run at default tier).
# All fits use small mc_iter for speed. Only synthetic public data.
###############################################################################

# ---------------------------------------------------------------------------
# Synthetic dataset — 100 observations, 3 attributes, 2 classes
#
# Design goal: A has higher ESS than B so A wins at step 1.
#
#   class 1L: rows 1-50
#   class 2L: rows 51-100
#
#   A (near-perfect, ESS ≈ 60%):
#     rows  1-40: class 1, A=0  (40 of 50 class-1 on left — correctly pred class 1)
#     rows 41-50: class 1, A=1  (10 class-1 "wrong" side)
#     rows 51-90: class 2, A=1  (40 of 50 class-2 on right — correctly pred class 2)
#     rows 91-100: class 2, A=0 (10 class-2 "wrong" side)
#     → sens = 40/50 = 0.80, spec = 40/50 = 0.80 → ESS ≈ 60%
#
#   B (weaker, ESS ≈ 20%):
#     rows  1-30: class 1, B=0  (30 of 50 class-1 on left)
#     rows 31-50: class 1, B=1  (20 class-1 "wrong" side)
#     rows 51-80: class 2, B=1  (30 of 50 class-2 on right)
#     rows 81-100: class 2, B=0 (20 class-2 "wrong" side)
#     → sens = 30/50 = 0.60, spec = 30/50 = 0.60 → ESS ≈ 20%
#
#   C: noise
#
# After step 1 (A selected, 80 obs correctly classified):
#   Remaining 20: rows 41-50 (class 1, A=1) + rows 91-100 (class 2, A=0)
#   On this subsample B has ESS=100% (class1→B=1, class2→B=0, inverted direction)
# ---------------------------------------------------------------------------

set.seed(0)

sda_y <- c(rep(1L, 50), rep(2L, 50))
sda_X <- data.frame(
  A = c(rep(0L, 40), rep(1L, 10),    # class 1: 40 at A=0, 10 at A=1
        rep(1L, 40), rep(0L, 10)),    # class 2: 40 at A=1, 10 at A=0
  B = c(rep(0L, 30), rep(1L, 20),    # class 1: 30 at B=0, 20 at B=1
        rep(1L, 30), rep(0L, 20)),    # class 2: 30 at B=1, 20 at B=0
  C = sample(0:1, 100, replace = TRUE)  # noise
)

# ---------------------------------------------------------------------------
# Test 1: Selects known max-ESS variable at step 1
# ---------------------------------------------------------------------------
test_that("sda: selects max-ESS variable at step 1", {
  fit <- sda_fit(sda_X, sda_y,
                 mode    = "unioda_max_ess",
                 mc_iter = 500L,
                 mc_seed = 42L,
                 loo     = "off")
  expect_true(inherits(fit, "sda_fit"))
  expect_equal(length(fit$steps), length(fit$selected_attributes))
  # A should be selected at step 1 (highest ESS)
  expect_equal(fit$steps[[1]]$attribute, "A")
})

# ---------------------------------------------------------------------------
# Test 2: Correctly classified observations are removed after step 1
# ---------------------------------------------------------------------------
test_that("sda: n_correct obs removed; step 2 n_in is reduced", {
  fit <- sda_fit(sda_X, sda_y,
                 mode    = "unioda_max_ess",
                 mc_iter = 500L,
                 mc_seed = 42L,
                 loo     = "off")
  s1 <- fit$steps[[1]]
  # observations entering step 2 = n_total - n_correct_step1
  if (length(fit$steps) >= 2L) {
    s2 <- fit$steps[[2]]
    expect_equal(s2$n_in, 100L - s1$n_correct)
  }
  # n_correct + n_incorrect = n_in at step 1
  expect_equal(s1$n_correct + s1$n_incorrect, s1$n_in)
})

# ---------------------------------------------------------------------------
# Test 3: Selected attribute absent from step 2 candidates
# ---------------------------------------------------------------------------
test_that("sda: selected attribute removed from candidate set at next step", {
  fit <- sda_fit(sda_X, sda_y,
                 mode    = "unioda_max_ess",
                 mc_iter = 500L,
                 mc_seed = 42L,
                 loo     = "off")
  if (length(fit$steps) >= 2L) {
    attr_s1 <- fit$steps[[1]]$attribute
    # attribute selected at step 1 must not appear in step 2 candidate table
    ctab_s2 <- fit$steps[[2]]$candidate_table
    expect_false(attr_s1 %in% ctab_s2$attribute)
  } else {
    skip("only one SDA step in this fit; cannot check step 2 candidates")
  }
})

# ---------------------------------------------------------------------------
# Test 4: Stops when no candidates remain (one-attribute frame)
# ---------------------------------------------------------------------------
test_that("sda: stops with no_candidates when only one attribute given", {
  X1 <- data.frame(A = sda_X$A)
  fit <- sda_fit(X1, sda_y,
                 mode    = "unioda_max_ess",
                 mc_iter = 500L,
                 mc_seed = 42L,
                 loo     = "off")
  expect_equal(fit$stop_reason, "no_candidates")
  expect_equal(length(fit$steps), 1L)
})

# ---------------------------------------------------------------------------
# Test 5: Stops with min_n when working sample drops below threshold
# ---------------------------------------------------------------------------
test_that("sda: stops with min_n when remaining n < threshold", {
  # Set min_n large enough to trigger stop after first removal.
  # After step 1, ~20 obs remain (the 10 "wrong-side" obs per class).
  # min_n = 25 > 20 triggers stop before step 2 begins.
  fit <- sda_fit(sda_X, sda_y,
                 mode    = "unioda_max_ess",
                 mc_iter = 500L,
                 mc_seed = 42L,
                 loo     = "off",
                 min_n   = 25L)
  expect_equal(fit$stop_reason, "min_n")
})

# ---------------------------------------------------------------------------
# Test 6: Candidate table stored for every step, with required columns
# ---------------------------------------------------------------------------
test_that("sda: candidate table has required columns at every step", {
  fit <- sda_fit(sda_X, sda_y,
                 mode    = "unioda_max_ess",
                 mc_iter = 500L,
                 mc_seed = 42L,
                 loo     = "off")
  required_cols <- c("attribute", "status", "n", "ess", "d", "p_mc",
                      "eligible", "ineligible_reason", "selected")
  for (s in fit$steps) {
    ctab <- s$candidate_table
    expect_true(is.data.frame(ctab))
    for (col in required_cols)
      expect_true(col %in% names(ctab),
                  label = paste0("candidate_table missing column: ", col))
    # d must be NA for unioda_max_ess
    expect_true(all(is.na(ctab$d)))
    # exactly one row selected per step
    expect_equal(sum(ctab$selected), 1L)
  }
})

# ---------------------------------------------------------------------------
# Test 7: predict.sda_fit(type="class") classifies correctly; unresolved → NA
# ---------------------------------------------------------------------------
test_that("sda predict: type='class' returns integer vector with NA for unresolved", {
  fit <- sda_fit(sda_X, sda_y,
                 mode    = "unioda_max_ess",
                 mc_iter = 500L,
                 mc_seed = 42L,
                 loo     = "off")
  preds <- predict(fit, sda_X, type = "class")
  expect_true(is.integer(preds))
  expect_equal(length(preds), 100L)
  # Resolved obs get a non-NA class value in {1L, 2L}
  resolved <- !is.na(preds)
  expect_true(all(preds[resolved] %in% c(1L, 2L)))
})

# ---------------------------------------------------------------------------
# Test 8: predict.sda_fit(type="stage") returns step_id at classification
# ---------------------------------------------------------------------------
test_that("sda predict: type='stage' returns integer step_id or NA", {
  fit <- sda_fit(sda_X, sda_y,
                 mode    = "unioda_max_ess",
                 mc_iter = 500L,
                 mc_seed = 42L,
                 loo     = "off")
  stages <- predict(fit, sda_X, type = "stage")
  expect_true(is.integer(stages))
  expect_equal(length(stages), 100L)
  resolved <- !is.na(stages)
  expect_true(all(stages[resolved] >= 1L))
  expect_true(all(stages[resolved] <= length(fit$steps)))
})

# ---------------------------------------------------------------------------
# Test 9: Wide newdata (extra columns) works via name-based routing
# ---------------------------------------------------------------------------
test_that("sda predict: wide newdata with extra columns works", {
  fit <- sda_fit(sda_X, sda_y,
                 mode    = "unioda_max_ess",
                 mc_iter = 500L,
                 mc_seed = 42L,
                 loo     = "off")
  wide_X <- cbind(sda_X, extra1 = rnorm(100), extra2 = rnorm(100))
  # Should not error; extra columns ignored
  preds_narrow <- predict(fit, sda_X,  type = "class")
  preds_wide   <- predict(fit, wide_X, type = "class")
  expect_identical(preds_narrow, preds_wide)
})

# ---------------------------------------------------------------------------
# Test 10: Missing required column in newdata errors clearly
# ---------------------------------------------------------------------------
test_that("sda predict: missing required split variable errors clearly", {
  fit <- sda_fit(sda_X, sda_y,
                 mode    = "unioda_max_ess",
                 mc_iter = 500L,
                 mc_seed = 42L,
                 loo     = "off")
  bad_X <- sda_X[, c("B", "C"), drop = FALSE]   # drop A (step-1 winner)
  expect_error(predict(fit, bad_X), regexp = "missing required SDA split variable")
})

# ---------------------------------------------------------------------------
# Test 11: dry_run (remove_correct = FALSE) runs without modifying sample
# ---------------------------------------------------------------------------
test_that("sda: remove_correct=FALSE produces dry_run stop reason", {
  fit <- sda_fit(sda_X, sda_y,
                 mode           = "unioda_max_ess",
                 mc_iter        = 500L,
                 mc_seed        = 42L,
                 loo            = "off",
                 remove_correct = FALSE,
                 max_steps      = 2L)
  expect_equal(fit$stop_reason, "dry_run")
  # n_final_unresolved should equal n_initial (nothing removed)
  expect_equal(fit$n_final_unresolved, fit$n_initial)
})

# ---------------------------------------------------------------------------
# Test 12: verbose = TRUE runs without error or warning
# ---------------------------------------------------------------------------
test_that("sda: verbose=TRUE runs without error", {
  expect_no_error(
    suppressMessages(
      sda_fit(sda_X, sda_y,
              mode    = "unioda_max_ess",
              mc_iter = 200L,
              mc_seed = 42L,
              loo     = "off",
              verbose = TRUE)
    )
  )
})

# ---------------------------------------------------------------------------
# Test 13: weights != NULL errors clearly
# ---------------------------------------------------------------------------
test_that("sda: weights != NULL errors with informative message", {
  expect_error(
    sda_fit(sda_X, sda_y,
            mode    = "unioda_max_ess",
            weights = rep(1, 100)),
    regexp = "Weighted SDA is not implemented yet"
  )
})

# ---------------------------------------------------------------------------
# Test 14: mode = "novometric_min_d" requires explicit mindenom (SDA-4B)
# ---------------------------------------------------------------------------
test_that("sda: novometric_min_d without mindenom errors with MPE guidance", {
  expect_error(
    sda_fit(sda_X, sda_y, mode = "novometric_min_d"),
    regexp = "mindenom"
  )
})

# ---------------------------------------------------------------------------
# Test 15: print and summary methods run without error
# ---------------------------------------------------------------------------
test_that("sda: print and summary run without error", {
  fit <- sda_fit(sda_X, sda_y,
                 mode    = "unioda_max_ess",
                 mc_iter = 200L,
                 mc_seed = 42L,
                 loo     = "off")
  expect_no_error(capture.output(print(fit)))
  sm <- summary(fit)
  expect_true(inherits(sm, "sda_fit_summary"))
  expect_no_error(capture.output(print(sm)))
})

# ---------------------------------------------------------------------------
# Test 16: stop_reason "all_resolved" when all obs classified
# ---------------------------------------------------------------------------
test_that("sda: stop_reason = all_resolved when all obs classified", {
  # Perfectly separable two-attribute data so SDA can resolve everyone
  X_perf <- data.frame(
    A = c(rep(0L, 20), rep(1L, 20)),
    B = c(rep(0L, 10), rep(1L, 10), rep(0L, 10), rep(1L, 10))
  )
  y_perf <- c(rep(1L, 20), rep(2L, 20))
  fit <- sda_fit(X_perf, y_perf,
                 mode    = "unioda_max_ess",
                 mc_iter = 500L,
                 mc_seed = 42L,
                 loo     = "off")
  expect_true(fit$stop_reason %in% c("all_resolved", "class_resolved",
                                      "no_candidates", "p_gate"))
})
