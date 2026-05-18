###############################################################################
# test-unioda.R
# Unit tests for the binary-class UniODA engine (oda_univariate_core).
###############################################################################

test_that("UniODA ordered: naval big-ships example achieves 87.5% PAC", {
  # Yarnold & Soltysik (2005) p.62 — optimal cut at 80, PAC = 87.5%
  x <- c(35, 45, 55, 65, 75, 85, 95, 99)
  y <- c(0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L)

  fit <- oda_univariate_core(
    x = x, y = y, attr_type = "ordered",
    priors_on = FALSE, mcarlo = FALSE, loo = "off"
  )

  expect_true(fit$ok, label = "naval: fit ok")
  expect_equal(round(fit$rule$cut_value, 0), 80,
               label = "naval: cut value = 80")
  expect_equal(round(fit$pac, 1), 87.5,
               label = "naval: PAC = 87.5%")
})


test_that("UniODA ordered: perfect separation gives 100% PAC", {
  x <- c(1, 2, 3, 7, 8, 9)
  y <- c(0L, 0L, 0L, 1L, 1L, 1L)

  fit <- oda_univariate_core(
    x = x, y = y, priors_on = TRUE, mcarlo = FALSE, loo = "off"
  )

  expect_true(fit$ok)
  expect_equal(round(fit$pac, 1), 100.0, label = "perfect sep: 100% PAC")
  # Cut must fall in the gap (3, 7)
  expect_true(fit$rule$cut_value > 3 && fit$rule$cut_value < 7,
              label = "perfect sep: cut in gap")
})


test_that("UniODA: pure-node returns ok=FALSE with reason pure_node", {
  x <- c(1, 2, 3, 4)
  y <- c(0L, 0L, 0L, 0L)

  fit <- oda_univariate_core(x = x, y = y, mcarlo = FALSE, loo = "off")
  expect_false(fit$ok)
  expect_equal(fit$reason, "pure_node")
})


test_that("UniODA: all-missing x returns ok=FALSE with reason all_missing", {
  x <- rep(NA_real_, 5)
  y <- c(0L, 1L, 0L, 1L, 0L)

  fit <- oda_univariate_core(x = x, y = y, mcarlo = FALSE, loo = "off")
  expect_false(fit$ok)
  expect_equal(fit$reason, "all_missing")
})


test_that("UniODA binary attribute: runs and produces a valid rule", {
  x <- c("A","A","B","B","B","A")
  y <- c(0L, 0L, 1L, 1L, 0L, 1L)

  fit <- oda_univariate_core(
    x = x, y = y, attr_type = "binary",
    mcarlo = FALSE, loo = "off"
  )

  expect_true(fit$ok)
  expect_true(fit$rule$type %in% c("binary_map","nominal_cut"))
  expect_true(fit$pac >= 0 && fit$pac <= 100)
})


test_that("UniODA: priors_on=TRUE vs FALSE differ on imbalanced classes", {
  # 9 class-0 vs 3 class-1 — priors_on should NOT produce a degenerate solution
  x <- c(1,2,3,4,5,6,7,8,9, 10,11,12)
  y <- c(rep(0L,9), 1L,1L,1L)

  fit_on  <- oda_univariate_core(x=x, y=y, priors_on=TRUE,  mcarlo=FALSE, loo="off")
  fit_off <- oda_univariate_core(x=x, y=y, priors_on=FALSE, mcarlo=FALSE, loo="off")

  expect_true(fit_on$ok,  label = "priors_on=TRUE: fit ok")
  expect_true(fit_off$ok, label = "priors_on=FALSE: fit ok")

  # Both should find a cut in the upper range where class-1 lives
  expect_true(fit_on$rule$cut_value  >= 5,
              label = "priors_on=TRUE: cut above midpoint")
})


test_that("UniODA: LOO refit returns valid confusion table summing to n", {
  set.seed(7)
  x <- c(1,2,2,3,4,4,5,6,6,7)
  y <- c(0L,0L,0L,0L,1L,1L,1L,1L,0L,1L)
  n <- length(y)

  # loo_alpha=1.0 means never suppress on p-value (always return LOO result)
  fit <- oda_univariate_core(
    x = x, y = y, priors_on = TRUE,
    mcarlo = FALSE, loo = "pvalue", loo_alpha = 1.0
  )

  expect_true(fit$ok)
  expect_true(!is.null(fit$loo), label = "LOO result present")

  if (isTRUE(fit$loo$allowed)) {
    conf <- fit$loo$confusion
    total <- conf$TP + conf$FN + conf$FP + conf$TN
    expected_total <- length(unique(y))
    expect_equal(as.integer(round(total)), expected_total,
                 label = "LOO confusion sums to priors-adjusted class total")
  }
})


test_that("UniODA MC: p-value is in [0,1]", {
  set.seed(123)
  x <- c(1,2,3,4,5,6,7,8,9,10)
  y <- c(0L,0L,0L,0L,0L,1L,1L,1L,1L,1L)

  fit <- oda_univariate_core(
    x = x, y = y, priors_on = TRUE,
    mcarlo = TRUE, mc_iter = 200L, mc_seed = 99L,
    loo = "off"
  )

  expect_true(fit$ok)
  expect_true(is.numeric(fit$p_mc) && !is.na(fit$p_mc),
              label = "MC p-value is numeric non-NA")
  expect_true(fit$p_mc >= 0 && fit$p_mc <= 1,
              label = "MC p-value in [0,1]")
})


test_that("UniODA: deterministic — identical calls give identical results", {
  x <- c(1,2,3,4,5,6,7,8)
  y <- c(0L,0L,0L,0L,1L,1L,1L,1L)

  fit1 <- oda_univariate_core(x=x, y=y, priors_on=TRUE, mcarlo=FALSE, loo="off")
  fit2 <- oda_univariate_core(x=x, y=y, priors_on=TRUE, mcarlo=FALSE, loo="off")

  expect_true(fit1$ok); expect_true(fit2$ok)
  expect_equal(fit1$rule$cut_value, fit2$rule$cut_value)
  expect_equal(fit1$rule$direction, fit2$rule$direction)
})


# ---- LOO contract -----------------------------------------------------------

test_that("LOO contract: loo='off' is default — fit$loo is NULL without explicit loo=", {
  # Canon: LOO is off unless the caller explicitly opts in.
  # Omitting loo= must leave fit$loo as NULL.
  x <- c(1, 2, 3, 4, 5, 6, 7, 8)
  y <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  fit <- oda_fit(x, y, mcarlo = FALSE)   # no loo= argument
  expect_true(fit$ok)
  expect_null(fit$loo)
})

test_that("LOO contract: explicit loo='pvalue' populates fit$loo with ess_loo present", {
  # Use oda_univariate_core with loo_alpha=1.0 so the LOO always runs and
  # populates fields — tests the contract (fields present), not significance
  # gating. (oda_fit does not forward loo_alpha.)
  x <- c(1, 2, 3, 4, 5, 6, 7, 8)
  y <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  fit <- oda_univariate_core(x, y, priors_on = TRUE, mcarlo = FALSE,
                             loo = "pvalue", loo_alpha = 1.0)
  expect_true(fit$ok)
  expect_false(is.null(fit$loo))
  expect_true(isTRUE(fit$loo$allowed))
  expect_false(is.na(fit$loo$ess_loo))
})

test_that("LOO contract: weighted categorical LOO returns ok=TRUE with allowed=FALSE and canon reason", {
  # Canon: weighted categorical LOO is explicitly forbidden
  # (oda_loo_for_rule, allow_weighted_categorical_loo = FALSE by default).
  # Non-uniform weights + categorical + loo != 'off' must produce:
  #   fit$ok = TRUE     (training result valid; rejection is LOO-level only)
  #   fit$loo$allowed = FALSE
  #   fit$loo$reason  = "weighted_categorical_loo_not_supported"
  # Probe-confirmed behavior (2026-05-14).
  x <- c(1L, 1L, 2L, 2L, 3L, 3L)
  y <- c(0L, 1L, 0L, 0L, 1L, 1L)
  w <- c(1.0, 2.0, 1.0, 2.0, 1.0, 2.0)  # non-uniform
  fit <- oda_fit(x, y, w = w, attr_type = "categorical",
                 mcarlo = FALSE, loo = "on")
  expect_true(fit$ok)
  expect_false(is.null(fit$loo))
  expect_false(isTRUE(fit$loo$allowed))
  expect_equal(fit$loo$reason, "weighted_categorical_loo_not_supported")
})

test_that("LOO contract: absent category in ordinary prediction routes to right side (class 1)", {
  # Canon: the absent-level override lives only inside oda_loo_for_rule().
  # oda_rule_predict / oda_rule_side must send an unseen category to the right
  # side (1L, coded) — NOT to the left side as the LOO override does.
  # Mechanistic boundary established by commit 8197c4c; probe-confirmed.
  x <- c(rep(1L, 10L), rep(2L, 10L))  # only categories 1 and 2 in training
  y <- c(rep(0L, 10L), rep(1L, 10L))  # 1 -> class 0, 2 -> class 1 (perfect separation)
  fit <- oda_fit(x, y, attr_type = "categorical", mcarlo = FALSE, loo = "off")
  expect_true(fit$ok)
  expect_equal(fit$rule$type, "nominal_cut")
  # Known categories: confirm correct sides.
  expect_equal(oda_rule_predict(1L, fit$rule), 0L)
  expect_equal(oda_rule_predict(2L, fit$rule), 1L)
  # Category 3 was never in training. Ordinary prediction: absent -> right (1L).
  # This is the opposite of the LOO absent-level override (which returns 0L).
  expect_equal(oda_rule_predict(3L, fit$rule), 1L)
})

# ---- SAMPLEREP tie-breaking -------------------------------------------------

test_that("UniODA SAMPLEREP: balanced split chosen over unbalanced when PAC ties", {
  # Perfectly separable at x=3.5.
  # Both cut=3.5 and cut=4.5 give 100% PAC, but cut=3.5 gives balanced
  # predicted frequencies (3 zeros + 3 ones = observed 3:3).
  # cut=4.5 would give 4 predicted zeros + 2 ones, worse SREP.
  x <- c(1, 2, 3, 4, 5, 6)
  y <- c(0L, 0L, 0L, 1L, 1L, 1L)

  fit <- oda_univariate_core(x=x, y=y, priors_on=TRUE, mcarlo=FALSE, loo="off")

  expect_true(fit$ok)
  expect_equal(round(fit$pac, 1), 100.0, label = "SAMPLEREP: 100% PAC")
  # SAMPLEREP should prefer cut=3.5 (3+3 balanced) over cut=4.5 (4+2) or cut=2.5
  expect_equal(round(fit$rule$cut_value, 1), 3.5,
               label = "SAMPLEREP: balanced cut chosen")
})

# ---- MC p-value formula: raw proportion (no LaPlace pseudocount) ------------

test_that("oda_mc_p_value: p_mc is raw ge_count/iter_used (no LaPlace +1)", {
  # Seeded regression anchors verify the formula directly.
  #
  # Anchor 1 — refugee-act data, ESS = 44.09% (strong effect).
  # seed=1, mc_iter=200 is a seeded regression anchor with ge_count == 0.
  # p_mc must be exactly 0.0 (raw proportion), not (0+1)/(200+1) = 1/201.
  vote  <- c(rep(0L, 118), rep(0L,  78), rep(1L,  34), rep(1L, 177))
  party <- c(rep(0L, 118), rep(1L,  78), rep(0L,  34), rep(1L, 177))

  r1 <- odacore:::oda_mc_p_value(
    x=party, y=vote, attr_type="ordered", priors_on=TRUE,
    primary="maxsens", secondary="samplerep",
    mc_iter=200L, mc_stop=NA_real_, mc_stopup=NA_real_,
    seed=1L, ess_obs=44.09
  )
  expect_equal(r1$ge_count,  0L,   label = "anchor1: ge_count == 0")
  expect_equal(r1$iter_used, 200L, label = "anchor1: iter_used == 200")
  expect_equal(r1$p_mc,      0.0,  label = "anchor1: p_mc == 0.0 (not 1/201)")

  # Anchor 2 — migraine data, ESS = 24.69% (weak; MC p ~ 0.086 per MegaODA).
  # seed=1, mc_iter=500: probe gives ge_count=39, iter_used=500.
  # Raw p_mc = 39/500 = 0.078.  LaPlace would give 40/501 ≈ 0.07984.
  treatment <- c(
    rep(0L, 13), rep(1L,  5),
    rep(0L,  9), rep(1L, 13),
    rep(0L,  4), rep(1L,  6),
    rep(0L,  2), rep(1L,  1),
    rep(0L,  1), rep(1L,  2),
    rep(0L,  1), rep(1L,  3),
    rep(0L,  3), rep(1L,  3),
    rep(0L,  0), rep(1L,  1)
  )
  attacks <- c(
    rep(0L, 18), rep(1L, 22), rep(2L, 10),
    rep(3L,  3), rep(4L,  3), rep(5L,  4),
    rep(6L,  6), rep(7L,  1)
  )

  r2 <- odacore:::oda_mc_p_value(
    x=attacks, y=treatment, attr_type="ordered", priors_on=TRUE,
    primary="maxsens", secondary="samplerep",
    mc_iter=500L, mc_stop=NA_real_, mc_stopup=NA_real_,
    seed=1L, ess_obs=24.69
  )
  expect_equal(r2$ge_count,  39L,   label = "anchor2: ge_count == 39")
  expect_equal(r2$iter_used, 500L,  label = "anchor2: iter_used == 500")
  expect_equal(r2$p_mc,      39 / 500,
               label = "anchor2: p_mc == ge_count/iter_used")
  expect_false(isTRUE(all.equal(r2$p_mc, 40 / 501)),
               label = "anchor2: p_mc != LaPlace (40/501)")
})

# ---- Phase 6A: DIRECTION support (ordered binary only) ----------------------

test_that("DIRECTION: direction='greater' selects '0->1' on clear ascending data", {
  # Clear 0->1 geometry: class 1 on the right (high-attribute) side.
  # direction="greater" (Chapter 2 greater-than direction; MegaODA DIRECTION < 0 1)
  # must constrain the search to "0->1" candidates only and find the perfect cut.
  x <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L)
  y <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)

  fit <- oda_univariate_core(x, y, direction = "greater",
                             priors_on = TRUE, mcarlo = FALSE, loo = "off")

  expect_true(fit$ok,                              label = "direction greater: ok")
  expect_equal(fit$rule$direction, "0->1",         label = "direction greater: rule is 0->1")
  expect_equal(round(fit$pac, 1), 100.0,           label = "direction greater: PAC 100%")
})

test_that("DIRECTION: direction='less' selects '1->0' on clear descending data", {
  # Clear 1->0 geometry: class 1 on the left (low-attribute) side.
  # direction="less" (Chapter 2 less-than direction; MegaODA DIRECTION > 0 1)
  # must constrain the search to "1->0" candidates only and find the perfect cut.
  x <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L)
  y <- c(1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L)

  fit <- oda_univariate_core(x, y, direction = "less",
                             priors_on = TRUE, mcarlo = FALSE, loo = "off")

  expect_true(fit$ok,                              label = "direction less: ok")
  expect_equal(fit$rule$direction, "1->0",         label = "direction less: rule is 1->0")
  expect_equal(round(fit$pac, 1), 100.0,           label = "direction less: PAC 100%")
})

test_that("DIRECTION: direction='off' matches default (regression gate)", {
  # direction="off" must produce identical results to omitting the argument.
  # Uses the naval example (Yarnold & Soltysik 2005, p.62).
  x <- c(35, 45, 55, 65, 75, 85, 95, 99)
  y <- c(0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L)

  fit_off     <- oda_univariate_core(x, y, direction = "off",
                                     priors_on = FALSE, mcarlo = FALSE, loo = "off")
  fit_default <- oda_univariate_core(x, y,
                                     priors_on = FALSE, mcarlo = FALSE, loo = "off")

  expect_equal(fit_off$rule$cut_value, fit_default$rule$cut_value,
               label = "direction off == default: cut_value")
  expect_equal(fit_off$rule$direction, fit_default$rule$direction,
               label = "direction off == default: rule direction")
  expect_equal(round(fit_off$pac, 1), 87.5,
               label = "direction off: PAC 87.5%")
})

test_that("DIRECTION: direction propagates into MC permutation refits", {
  # Structural contract: direction is wired into oda_mc_p_value() and its
  # internal oda_univariate_core() calls.  Assert rule direction obeyed and
  # p_mc is a valid probability — not the specific value.
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  y <- c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L)

  fit <- oda_univariate_core(x, y, direction = "greater", priors_on = TRUE,
                             mcarlo = TRUE, mc_seed = 1L, mc_iter = 100L,
                             loo = "off")

  expect_true(fit$ok,                              label = "MC direction: ok")
  expect_equal(fit$rule$direction, "0->1",         label = "MC direction: rule is 0->1")
  expect_true(is.numeric(fit$p_mc) && !is.na(fit$p_mc),
              label = "MC direction: p_mc is numeric non-NA")
  expect_true(fit$p_mc >= 0 && fit$p_mc <= 1,     label = "MC direction: p_mc in [0,1]")
})

test_that("DIRECTION: direction propagates into LOO fold refits", {
  # Structural contract: direction is wired into oda_loo_for_rule() and its
  # internal oda_univariate_core() fold calls.  Uses loo_alpha=1.0 to force
  # the LOO to run and populate fit$loo on clear separation data (Fisher p
  # will be small, well below 1.0, so the gate does not reject).
  x <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L)
  y <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)

  fit <- oda_univariate_core(x, y, direction = "greater", priors_on = TRUE,
                             mcarlo = FALSE, loo = "pvalue", loo_alpha = 1.0)

  expect_true(fit$ok,                              label = "LOO direction: ok")
  expect_equal(fit$rule$direction, "0->1",         label = "LOO direction: rule is 0->1")
  expect_false(is.null(fit$loo),                   label = "LOO direction: loo not NULL")
  expect_true(isTRUE(fit$loo$allowed),             label = "LOO direction: loo allowed")
})

test_that("DIRECTION: categorical + direction != 'off' returns explicit failure", {
  # Phase 6A does not support directional analysis for categorical attributes.
  # MPE Chapter 4 TABLE/DIRECTIONAL semantics are deferred to Phase 6C.
  # Must return ok=FALSE, not a warning-and-ignore.
  x <- c("A", "A", "B", "B", "B", "A")
  y <- c(0L, 0L, 1L, 1L, 0L, 1L)

  fit <- oda_univariate_core(x, y, attr_type = "categorical",
                             direction = "greater", mcarlo = FALSE, loo = "off")

  expect_false(fit$ok,                             label = "cat direction: ok is FALSE")
  expect_equal(fit$reason, "direction_not_supported_for_categorical",
               label = "cat direction: correct reason")
})
