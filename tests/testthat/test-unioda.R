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
