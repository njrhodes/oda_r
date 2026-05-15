###############################################################################
# test-tie-breaking.R
#
# Targeted tests for the three-layer tie-breaking spec:
#   PRIMARY   = MAXSENS (overall PAC in priors-weighted space)
#   SECONDARY = SAMPLEREP (L1 distance between predicted/observed class freqs)
#   TERTIARY  = FIRST IDENTIFIED (enumeration order)
###############################################################################

# ---- Layer 1: PRIMARY (MAXSENS) ------------------------------------------ #

test_that("layer1 PRIMARY: perfectly separated 3-class data gets 100% mean PAC", {
  x <- c(1,1,1, 5,5, 9,9,9)
  y <- c(1L,1L,1L, 2L,2L, 3L,3L,3L)

  fit <- oda_multiclass_unioda_core(
    x=x, y=y, attr_type="ordered",
    priors_on=TRUE, K_segments=3L, degen=FALSE,
    mcarlo=FALSE, loo="off"
  )

  expect_true(fit$ok)
  expect_equal(round(fit$mean_pac, 1), 100.0,
               label = "perfect separation: 100% mean PAC")
  expect_equal(length(unique(as.integer(fit$rule$seg_classes))), 3L,
               label = "all 3 classes predicted")
})


# ---- Layer 2: SAMPLEREP decisive ----------------------------------------- #

test_that("layer2 SAMPLEREP: seg=(2,1,3) is the only valid nondegenerate at cut (1.5,2.5)", {
  # Data: nc1=4,nc2=3,nc3=4, n=11
  # At cut (1.5,2.5):
  #   seg1 argmax=c2 (unique: c2_pw=1/3 > c1_pw=1/4)
  #   seg2 argmax=c1 (unique: c1_pw=0.5 > c3_pw=0.25)
  #   seg3 argmax=c3 (unique: c3_pw=0.75 > c2_pw=0.667 > c1_pw=0.25)
  # → only (2,1,3) is nondegenerate. First-identified-cut also picks this cut.
  # At cut (2.5,3.5): assignment (1,3,2) also ties on primary.
  # Two-level: cut (1.5,2.5) enumerated first → wins → seg=(2,1,3).
  x <- c(1L, 1L,  2L, 2L, 2L,  3L, 3L, 3L, 3L, 3L, 3L)
  y <- c(1L, 2L,  1L, 1L, 3L,  1L, 2L, 2L, 3L, 3L, 3L)

  fit <- oda_multiclass_unioda_core(
    x=x, y=y, attr_type="ordered",
    priors_on=TRUE, K_segments=3L, degen=FALSE,
    mcarlo=FALSE, loo="off"
  )

  expect_true(fit$ok)
  expect_equal(round(fit$rule$cut_values, 2), c(1.5, 2.5),
               label = "first enumerated cut (1.5,2.5) wins")
  expect_equal(as.integer(fit$rule$seg_classes), c(2L, 1L, 3L),
               label = "only valid nondegenerate assignment at (1.5,2.5) is (2,1,3)")
})


test_that("layer2 SAMPLEREP: balanced predicted distribution preferred", {
  # 6 obs, 3 classes, 2 each — any valid cut achieves perfect PAC,
  # but SAMPLEREP should ensure predicted counts match observed (2,2,2)
  x <- c(1, 1, 3, 3, 5, 5)
  y <- c(1L,1L, 2L,2L, 3L,3L)

  fit <- oda_multiclass_unioda_core(
    x=x, y=y, attr_type="ordered",
    priors_on=TRUE, K_segments=3L, degen=FALSE,
    mcarlo=FALSE, loo="off"
  )

  expect_true(fit$ok)
  yhat <- oda_rule_predict_multiclass(x, fit$rule)
  pred_counts <- sort(as.integer(table(yhat)))
  expect_equal(pred_counts, c(2L,2L,2L),
               label = "SAMPLEREP: predicted (2,2,2) = observed")
})


# ---- Layer 3: FIRST IDENTIFIED (deterministic) --------------------------- #

test_that("layer3 FIRST IDENTIFIED: same call → identical rule every time", {
  x <- c(1,2,3,4,5,6,7,8,9)
  y <- c(1L,1L,1L,2L,2L,2L,3L,3L,3L)

  fit1 <- oda_multiclass_unioda_core(x=x, y=y, attr_type="ordered",
    priors_on=TRUE, K_segments=3L, degen=FALSE, mcarlo=FALSE, loo="off")
  fit2 <- oda_multiclass_unioda_core(x=x, y=y, attr_type="ordered",
    priors_on=TRUE, K_segments=3L, degen=FALSE, mcarlo=FALSE, loo="off")

  expect_true(fit1$ok); expect_true(fit2$ok)
  expect_equal(fit1$rule$cut_values,  fit2$rule$cut_values,
               label = "deterministic: same cuts")
  expect_equal(fit1$rule$seg_classes, fit2$rule$seg_classes,
               label = "deterministic: same seg classes")
})


# ---- PRIORS ON vs OFF on imbalanced data --------------------------------- #

test_that("priors_on=TRUE avoids degenerate solution on heavily imbalanced data", {
  # 8 class-1, 2 class-2, 2 class-3
  x <- c(1:8, 9,9, 10,10)
  y <- c(rep(1L,8), 2L,2L, 3L,3L)

  fit_on <- oda_multiclass_unioda_core(
    x=x, y=y, priors_on=TRUE, K_segments=3L, degen=FALSE,
    mcarlo=FALSE, loo="off"
  )

  expect_true(fit_on$ok, label = "priors_on: fit ok")
  expect_equal(length(unique(as.integer(fit_on$rule$seg_classes))), 3L,
               label = "priors_on: 3-class non-degenerate solution")
})


test_that("degen=FALSE and only 2 classes in data returns not_multiclass", {
  x <- c(1,2,3,4,5)
  y <- c(1L,1L,2L,2L,2L)   # only 2 distinct classes → multiclass core rejects

  fit <- oda_multiclass_unioda_core(
    x=x, y=y, attr_type="ordered",
    priors_on=TRUE, K_segments=3L, degen=FALSE,
    mcarlo=FALSE, loo="off"
  )

  expect_false(fit$ok)
  expect_equal(fit$reason, "not_multiclass",
               label = "2-class data rejected by multiclass core")
})


# ---- LOO internal consistency -------------------------------------------- #

test_that("LOO: core confusion_raw equals harness refit (internal consistency)", {
  x <- c(1,2,3,4,5,6,7,8,9)
  y <- c(1L,1L,1L,2L,2L,2L,3L,3L,3L)

  fit <- oda_multiclass_unioda_core(
    x=x, y=y, attr_type="ordered",
    priors_on=TRUE, K_segments=3L, degen=FALSE,
    mcarlo=FALSE, loo="on", boundary_mode="right_closed",
    loo_opts=list(grid_mode="refit")
  )

  expect_true(fit$ok)
  expect_true(isTRUE(fit$loo$allowed), label = "LOO: allowed")

  h <- harness_loo_refit_ordered_raw(
    x=x, y=y, priors_on=TRUE, degen=FALSE, K_segments=3L,
    boundary_mode="right_closed"
  )
  expect_true(h$ok, label = "harness LOO: ok")

  core_conf <- unname(fit$loo$confusion_raw$confusion)
  storage.mode(core_conf) <- "integer"
  expect_equal(core_conf, h$confusion,
               label = "core LOO == harness refit (internal consistency)")
})


# ---- ESS sanity ---------------------------------------------------------- #

test_that("ESS is in (-100, 100] for finite-PAC results", {
  x <- c(1,2,3,4,5,6)
  y <- c(1L,1L,2L,2L,3L,3L)

  fit <- oda_multiclass_unioda_core(
    x=x, y=y, attr_type="ordered",
    priors_on=TRUE, K_segments=3L, degen=FALSE,
    mcarlo=FALSE, loo="off"
  )

  expect_true(fit$ok)
  expect_true(is.finite(fit$ess), label = "ESS is finite")
  expect_true(fit$ess > -100 && fit$ess <= 100,
               label = "ESS in (-100,100]")
})
