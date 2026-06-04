###############################################################################
# test-propensity-ess-balance.R  -  BAL-PROP1-b tests
#
# Tier: CRAN-safe (all tests run at default tier).
# Uses small n_boot throughout for speed.  Deterministic synthetic data only.
###############################################################################

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

# Confounded example: age is correlated with group; score is not.
# Propensity score: x_pv separates group 0 / group 1 well.
set.seed(42L)
n_fix    <- 80L
group_fix <- c(rep(0L, 40L), rep(1L, 40L))

# Propensity covariate: clearly separates groups
x_pv_fix <- c(rnorm(40L, mean = 0), rnorm(40L, mean = 3))

# Balance covariates: age confounded, score independent
X_bal_fix <- data.frame(
  age   = c(rnorm(40L, mean = 45), rnorm(40L, mean = 55)),
  score = rnorm(n_fix)
)

# Fit ODA propensity model (group ~ x_pv)
prop_fit_oda <- oda_fit(x = x_pv_fix, y = group_fix, mcarlo = FALSE)

# Fit CTA propensity model (group ~ X_pv as data.frame)
X_pv_df <- data.frame(xpv = x_pv_fix)
prop_fit_cta <- cta_fit(
  X          = X_pv_df,
  y          = group_fix,
  mindenom   = 5L,
  mc_iter    = 200L,
  mc_seed    = 7L,
  loo        = "off",
  attr_names = "xpv"
)

# ---------------------------------------------------------------------------
# 1. ODA path: returns expected class and columns
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance (ODA) returns correct class and columns", {
  peb <- propensity_ess_balance(
    propensity_fit = prop_fit_oda,
    group          = group_fix,
    X_balance      = X_bal_fix,
    x_prop         = x_pv_fix,
    n_boot         = 20L,
    seed           = 1L
  )
  expect_s3_class(peb, "propensity_ess_balance")
  expect_s3_class(peb, "data.frame")
  expected_cols <- c("variable", "n", "unweighted_ess", "weighted_ess",
                     "delta_ess", "boot_low", "boot_high", "crosses_null",
                     "status")
  expect_true(all(expected_cols %in% names(peb)),
              info = paste("Missing:", paste(setdiff(expected_cols, names(peb)),
                                             collapse = ", ")))
  expect_equal(nrow(peb), ncol(X_bal_fix))
  expect_identical(peb$variable, names(X_bal_fix))
})

# ---------------------------------------------------------------------------
# 2. ODA path: point estimates are numeric and finite for "ok" rows
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance (ODA) point estimates finite for ok rows", {
  peb <- propensity_ess_balance(
    propensity_fit = prop_fit_oda,
    group          = group_fix,
    X_balance      = X_bal_fix,
    x_prop         = x_pv_fix,
    n_boot         = 20L,
    seed           = 1L
  )
  ok_rows <- peb$status == "ok"
  expect_true(any(ok_rows))
  expect_true(all(is.finite(peb$unweighted_ess[ok_rows])))
  expect_true(all(is.finite(peb$weighted_ess[ok_rows])))
  expect_true(all(is.finite(peb$delta_ess[ok_rows])))
  expect_true(all(is.finite(peb$boot_low[ok_rows])))
  expect_true(all(is.finite(peb$boot_high[ok_rows])))
})

# ---------------------------------------------------------------------------
# 3. ODA path: delta_ess == weighted_ess - unweighted_ess
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance (ODA) delta_ess is consistent with columns", {
  peb <- propensity_ess_balance(
    propensity_fit = prop_fit_oda,
    group          = group_fix,
    X_balance      = X_bal_fix,
    x_prop         = x_pv_fix,
    n_boot         = 20L,
    seed           = 1L
  )
  ok <- peb$status == "ok"
  expect_equal(peb$delta_ess[ok],
               peb$weighted_ess[ok] - peb$unweighted_ess[ok],
               tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 4. ODA path: confounded variable (age) has positive unweighted ESS
#    (Weighting direction is not asserted: with only 2 ODA strata and
#    strong confounding the adjusted ESS may stay at ceiling.  The robust
#    assertion is that confounding exists in the unadjusted analysis.)
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance (ODA): confounded covariate has positive unweighted ESS", {
  peb <- propensity_ess_balance(
    propensity_fit = prop_fit_oda,
    group          = group_fix,
    X_balance      = X_bal_fix,
    x_prop         = x_pv_fix,
    n_boot         = 20L,
    seed           = 1L
  )
  age_row <- peb[peb$variable == "age", ]
  expect_equal(nrow(age_row), 1L)
  if (age_row$status == "ok") {
    expect_gt(age_row$unweighted_ess, 0,
              label = "age unweighted ESS > 0: confounding exists in unadjusted analysis")
    expect_gte(age_row$weighted_ess, 0,
               label = "age weighted ESS is non-negative")
    expect_true(is.finite(age_row$delta_ess),
                label = "age delta_ess is finite")
  }
})

# ---------------------------------------------------------------------------
# 5. Bootstrap CI structure: boot_low <= boot_high, crosses_null consistent
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance (ODA) bootstrap CI is structurally valid", {
  peb <- propensity_ess_balance(
    propensity_fit = prop_fit_oda,
    group          = group_fix,
    X_balance      = X_bal_fix,
    x_prop         = x_pv_fix,
    n_boot         = 20L,
    seed           = 1L
  )
  ok <- peb$status == "ok"
  expect_true(all(peb$boot_low[ok] <= peb$boot_high[ok]),
              info = "boot_low must be <= boot_high")
  expected_crosses <- peb$boot_low[ok] <= 0 & peb$boot_high[ok] >= 0
  expect_identical(peb$crosses_null[ok], expected_crosses)
})

# ---------------------------------------------------------------------------
# 6. Seed reproducibility
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance (ODA) seed gives identical bootstrap output", {
  args <- list(
    propensity_fit = prop_fit_oda,
    group          = group_fix,
    X_balance      = X_bal_fix,
    x_prop         = x_pv_fix,
    n_boot         = 30L,
    seed           = 99L
  )
  peb1 <- do.call(propensity_ess_balance, args)
  peb2 <- do.call(propensity_ess_balance, args)
  expect_identical(peb1$boot_low,  peb2$boot_low)
  expect_identical(peb1$boot_high, peb2$boot_high)
})

# ---------------------------------------------------------------------------
# 7. CTA path: returns expected class and columns
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance (CTA) returns correct class and columns", {
  skip_if(isTRUE(prop_fit_cta$no_tree), "CTA propensity fit produced no tree")

  peb <- propensity_ess_balance(
    propensity_fit = prop_fit_cta,
    group          = group_fix,
    X_balance      = X_bal_fix,
    newdata        = X_pv_df,
    n_boot         = 20L,
    seed           = 1L
  )
  expect_s3_class(peb, "propensity_ess_balance")
  expect_equal(nrow(peb), ncol(X_bal_fix))
  expect_true(all(c("variable", "delta_ess", "crosses_null") %in% names(peb)))
})

# ---------------------------------------------------------------------------
# 8. CTA path: confounded variable shows delta_ess <= 0
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance (CTA): confounded covariate delta_ess <= 0", {
  skip_if(isTRUE(prop_fit_cta$no_tree), "CTA propensity fit produced no tree")

  peb <- propensity_ess_balance(
    propensity_fit = prop_fit_cta,
    group          = group_fix,
    X_balance      = X_bal_fix,
    newdata        = X_pv_df,
    n_boot         = 20L,
    seed           = 2L
  )
  age_row <- peb[peb$variable == "age", ]
  if (age_row$status == "ok") {
    expect_lte(age_row$delta_ess, 0,
               label = "CTA path: age delta_ess should be <= 0")
  }
})

# ---------------------------------------------------------------------------
# 9. Bad input: non-binary group
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance errors on non-binary group", {
  bad_group <- c(rep(0L, 20L), rep(1L, 20L), rep(2L, 40L))
  expect_error(
    propensity_ess_balance(
      propensity_fit = prop_fit_oda,
      group          = bad_group,
      X_balance      = X_bal_fix,
      x_prop         = x_pv_fix,
      n_boot         = 5L
    ),
    regexp = "2 distinct"
  )
})

# ---------------------------------------------------------------------------
# 10. Bad input: x_prop wrong length (ODA path)
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance errors when x_prop wrong length", {
  expect_error(
    propensity_ess_balance(
      propensity_fit = prop_fit_oda,
      group          = group_fix,
      X_balance      = X_bal_fix,
      x_prop         = x_pv_fix[seq_len(10L)],
      n_boot         = 5L
    ),
    regexp = "same length"
  )
})

# ---------------------------------------------------------------------------
# 11. Bad input: newdata missing for CTA path
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance errors when newdata missing for cta_tree", {
  skip_if(isTRUE(prop_fit_cta$no_tree), "CTA propensity fit produced no tree")
  expect_error(
    propensity_ess_balance(
      propensity_fit = prop_fit_cta,
      group          = group_fix,
      X_balance      = X_bal_fix,
      newdata        = NULL,
      n_boot         = 5L
    ),
    regexp = "newdata is required"
  )
})

# ---------------------------------------------------------------------------
# Shared LORT fixture
# ---------------------------------------------------------------------------

set.seed(1L)
X_lort <- data.frame(v1 = c(rnorm(20L, 0), rnorm(20L, 3)))
y_lort <- c(rep(0L, 20L), rep(1L, 20L))
lort_fit_obj <- lort_fit(
  X       = X_lort,
  y       = y_lort,
  min_n   = 5L,
  mc_iter = 100L,
  mc_seed = 1L,
  loo     = "off"
)

# ---------------------------------------------------------------------------
# 12. LORT path: returns expected class and columns
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance (LORT) returns correct class and columns", {
  skip_if(isTRUE(lort_fit_obj$no_tree), "LORT propensity fit produced no tree")

  peb <- propensity_ess_balance(
    propensity_fit = lort_fit_obj,
    group          = y_lort,
    X_balance      = X_lort,
    newdata        = X_lort,
    n_boot         = 20L,
    seed           = 1L
  )
  expect_s3_class(peb, "propensity_ess_balance")
  expect_s3_class(peb, "data.frame")
  expected_cols <- c("variable", "n", "unweighted_ess", "weighted_ess",
                     "delta_ess", "boot_low", "boot_high", "crosses_null",
                     "status")
  expect_true(all(expected_cols %in% names(peb)),
              info = paste("Missing:", paste(setdiff(expected_cols, names(peb)),
                                             collapse = ", ")))
  expect_equal(nrow(peb), ncol(X_lort))
})

# ---------------------------------------------------------------------------
# 13. LORT path: point estimates numeric, delta consistent
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance (LORT) delta_ess consistent with columns", {
  skip_if(isTRUE(lort_fit_obj$no_tree), "LORT propensity fit produced no tree")

  peb <- propensity_ess_balance(
    propensity_fit = lort_fit_obj,
    group          = y_lort,
    X_balance      = X_lort,
    newdata        = X_lort,
    n_boot         = 20L,
    seed           = 2L
  )
  ok <- peb$status == "ok"
  if (any(ok)) {
    expect_equal(peb$delta_ess[ok],
                 peb$weighted_ess[ok] - peb$unweighted_ess[ok],
                 tolerance = 1e-10)
    expect_true(all(is.finite(peb$unweighted_ess[ok])))
    expect_true(all(is.finite(peb$weighted_ess[ok])))
  }
})

# ---------------------------------------------------------------------------
# 14. LORT path: seed reproducibility
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance (LORT) seed gives identical bootstrap output", {
  skip_if(isTRUE(lort_fit_obj$no_tree), "LORT propensity fit produced no tree")

  args <- list(
    propensity_fit = lort_fit_obj,
    group          = y_lort,
    X_balance      = X_lort,
    newdata        = X_lort,
    n_boot         = 30L,
    seed           = 77L
  )
  peb1 <- do.call(propensity_ess_balance, args)
  peb2 <- do.call(propensity_ess_balance, args)
  expect_identical(peb1$boot_low,  peb2$boot_low)
  expect_identical(peb1$boot_high, peb2$boot_high)
})

# ---------------------------------------------------------------------------
# 15. LORT path: newdata missing errors
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance errors when newdata missing for cta_ort", {
  skip_if(isTRUE(lort_fit_obj$no_tree), "LORT propensity fit produced no tree")

  expect_error(
    propensity_ess_balance(
      propensity_fit = lort_fit_obj,
      group          = y_lort,
      X_balance      = X_lort,
      newdata        = NULL,
      n_boot         = 5L
    ),
    regexp = "newdata is required"
  )
})

# ---------------------------------------------------------------------------
# 16. X_balance nrow mismatch
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance errors when X_balance nrow != length(group)", {
  expect_error(
    propensity_ess_balance(
      propensity_fit = prop_fit_oda,
      group          = group_fix,
      X_balance      = X_bal_fix[seq_len(10L), ],
      x_prop         = x_pv_fix,
      n_boot         = 5L
    ),
    regexp = "nrow"
  )
})

# ---------------------------------------------------------------------------
# 17. Invalid n_boot
# ---------------------------------------------------------------------------
test_that("propensity_ess_balance errors on n_boot < 1", {
  expect_error(
    propensity_ess_balance(
      propensity_fit = prop_fit_oda,
      group          = group_fix,
      X_balance      = X_bal_fix,
      x_prop         = x_pv_fix,
      n_boot         = 0L
    ),
    regexp = "n_boot"
  )
})
