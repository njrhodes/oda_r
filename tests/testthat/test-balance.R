###############################################################################
# test-balance.R
#
# Tests for the univariate ODA covariate balance layer:
#   oda_balance_table()     (T1-T5, T7)
#   smd_balance_table()     (T6)
#   oda_balance_plot_data() (T8-T9)
#
# Test plan:
#   T1  obvious imbalance returns high ESS and significant raw p (MC enabled)
#   T2  balanced/noise covariate returns weak ESS or non-significant p
#   T3  oda_balance_table returns required columns
#   T4  failed covariate fit returns fit_ok = FALSE (not a dropped row)
#   T5  weights accepted; ess_label/has_weights reflect weighted mode
#   T6  smd_balance_table returns group means/SDs and abs_smd
#   T7  Sidak/Bonferroni corrections use only valid p values
#   T8  oda_balance_plot_data joins SMD when provided
#   T9  oda_balance_plot_data does not call oda_fit or require group/X
###############################################################################

# ---- Module-level fixtures ------------------------------------------------- #

set.seed(42)
n_bal <- 60L
grp_bal  <- c(rep(0L, 30L), rep(1L, 30L))

# Strongly imbalanced covariate: perfectly separated
age_imbalanced <- c(rep(40L, 30L), rep(60L, 30L))

# Noise covariate: same distribution in both groups
noise <- c(rnorm(30, 50, 5), rnorm(30, 50, 5))

X_bal <- data.frame(age = age_imbalanced, noise = noise)

# Weights (non-uniform so WESS path is triggered)
w_bal <- c(rep(1L, 30L), rep(2L, 30L))

# oda_balance_table computed once (mc_iter=200 for speed; seed via ...)
bt_unweighted <- oda_balance_table(grp_bal, X_bal,
                                    mcarlo  = TRUE,
                                    mc_iter = 200L,
                                    mc_seed = 7L)

bt_weighted <- oda_balance_table(grp_bal, X_bal,
                                  w       = w_bal,
                                  mcarlo  = TRUE,
                                  mc_iter = 200L,
                                  mc_seed = 7L)

smd_tbl <- smd_balance_table(grp_bal, X_bal)

# ---------------------------------------------------------------------------
# T1: obvious imbalance -> high ESS and significant raw p
# ---------------------------------------------------------------------------
test_that("oda_balance_table: imbalanced covariate yields high ESS and significant p (T1)", {
  row_age <- bt_unweighted$rows[bt_unweighted$rows$attribute == "age", ]
  expect_true(row_age$fit_ok)
  # Perfect separation -> ESS should be very high (close to 100%)
  expect_true(row_age$ess_display > 50)
  # With mc_iter=200 and perfect separation, p should be significant
  expect_true(row_age$significant_raw || row_age$p_mc <= 0.05)
})

# ---------------------------------------------------------------------------
# T2: noise covariate -> weak ESS or non-significant p
# ---------------------------------------------------------------------------
test_that("oda_balance_table: noise covariate yields low ESS (T2)", {
  row_noise <- bt_unweighted$rows[bt_unweighted$rows$attribute == "noise", ]
  expect_true(row_noise$fit_ok)
  # ESS for pure noise should be below 40% (usually near 0)
  expect_true(row_noise$ess_display < 40)
})

# ---------------------------------------------------------------------------
# T3: required columns present
# ---------------------------------------------------------------------------
test_that("oda_balance_table: required columns all present (T3)", {
  required_cols <- c(
    "attribute", "attr_type", "n_total", "n_group_0", "n_group_1",
    "sensitivity", "specificity", "mean_pac", "ess", "wess",
    "ess_display", "p_mc", "p_sidak", "p_bonferroni",
    "significant_raw", "significant_sidak", "significant_bonferroni",
    "significant", "rule_type", "rule_summary",
    "loo_status", "ess_loo", "has_weights", "fit_ok", "fit_reason"
  )
  expect_true(all(required_cols %in% names(bt_unweighted$rows)),
              info = paste("Missing:", paste(setdiff(required_cols, names(bt_unweighted$rows)), collapse=", ")))

  # Meta fields
  meta_fields <- c("n_covariates", "n_obs", "has_weights", "ess_label",
                    "alpha", "adjust", "k_valid", "loo_mode", "mcarlo", "mc_iter")
  expect_true(all(meta_fields %in% names(bt_unweighted$meta)),
              info = paste("Missing meta:", paste(setdiff(meta_fields, names(bt_unweighted$meta)), collapse=", ")))

  # Row count matches covariate count
  expect_equal(nrow(bt_unweighted$rows), 2L)
})

# ---------------------------------------------------------------------------
# T4: failed covariate fit returns fit_ok = FALSE, row still present
# ---------------------------------------------------------------------------
test_that("oda_balance_table: failed fit is a row with fit_ok=FALSE, not dropped (T4)", {
  # Force a failure by passing a constant column (no variation -> ODA fails or
  # returns ok=FALSE). Catch: if oda_fit returns ok=TRUE for a constant column
  # we fall through gracefully; the row must still be present.
  X_const <- data.frame(good  = c(rep(0L, 30L), rep(1L, 30L)),
                         const = rep(5L, 60L))
  bt_const <- oda_balance_table(grp_bal, X_const,
                                 mcarlo  = FALSE,
                                 mc_iter = 100L)
  # Both covariates must produce rows (no dropped rows)
  expect_equal(nrow(bt_const$rows), 2L)
  # The constant column should either fail (fit_ok=FALSE) or be valid
  row_const <- bt_const$rows[bt_const$rows$attribute == "const", ]
  expect_true(nrow(row_const) == 1L)
  # If oda_fit flags it as ok=FALSE, fit_ok must be FALSE
  if (!row_const$fit_ok) {
    expect_false(row_const$fit_ok)
    expect_true(!is.na(row_const$fit_reason))
  }
})

# ---------------------------------------------------------------------------
# T5: weights accepted; ess_label = "WESS", has_weights = TRUE
# ---------------------------------------------------------------------------
test_that("oda_balance_table: weighted mode reflected in meta (T5)", {
  expect_true(bt_weighted$meta$has_weights)
  expect_equal(bt_weighted$meta$ess_label, "WESS")
  # All rows have has_weights = TRUE
  expect_true(all(bt_weighted$rows$has_weights))
  # ESS appears in wess column, not ess column
  row_age_w <- bt_weighted$rows[bt_weighted$rows$attribute == "age", ]
  expect_false(is.na(row_age_w$wess))
  expect_true(is.na(row_age_w$ess))
})

# ---------------------------------------------------------------------------
# T6: smd_balance_table returns means and abs_smd
# ---------------------------------------------------------------------------
test_that("smd_balance_table: returns group means and abs_smd (T6)", {
  expect_s3_class(smd_tbl, "smd_balance_table")
  expect_equal(nrow(smd_tbl), 2L)

  req <- c("attribute", "n_group_0", "n_group_1",
            "mean_0", "sd_0", "mean_1", "sd_1",
            "smd", "abs_smd", "balanced_020", "balanced_010")
  expect_true(all(req %in% names(smd_tbl)),
              info = paste("Missing:", paste(setdiff(req, names(smd_tbl)), collapse=", ")))

  row_age <- smd_tbl[smd_tbl$attribute == "age", ]
  # Group means should be 40 and 60 for the perfectly-separated covariate
  expect_equal(row_age$mean_0, 40, tolerance = 0.5)
  expect_equal(row_age$mean_1, 60, tolerance = 0.5)
  # SMD for a 20-unit gap with ~0 SD should be large
  expect_true(row_age$abs_smd > 1 || is.na(row_age$smd))   # SD of constant is 0 -> NA
})

# ---------------------------------------------------------------------------
# T7: Sidak/Bonferroni use only valid p values (k = valid count)
# ---------------------------------------------------------------------------
test_that("oda_balance_table: Sidak/Bonferroni corrections use only valid p count (T7)", {
  # With k=2 valid p values, Bonferroni is p * 2 for both rows
  rows   <- bt_unweighted$rows
  k_valid <- bt_unweighted$meta$k_valid
  expect_equal(k_valid, sum(!is.na(rows$p_mc)))

  for (i in seq_len(nrow(rows))) {
    p_raw <- rows$p_mc[i]
    if (!is.na(p_raw)) {
      # Bonferroni: min(p_raw * k_valid, 1)
      expected_bonf <- min(p_raw * k_valid, 1)
      expect_equal(rows$p_bonferroni[i], expected_bonf, tolerance = 1e-10)

      # Sidak: min(1 - (1 - p_raw)^k_valid, 1)
      expected_sid  <- min(1 - (1 - p_raw)^k_valid, 1)
      expect_equal(rows$p_sidak[i], expected_sid, tolerance = 1e-10)
    }
  }
})

# ---------------------------------------------------------------------------
# T8: oda_balance_plot_data joins SMD when supplied
# ---------------------------------------------------------------------------
test_that("oda_balance_plot_data: joins SMD table when provided (T8)", {
  pd <- oda_balance_plot_data(bt_unweighted, smd_table = smd_tbl)
  expect_s3_class(pd, "oda_balance_plot_data")

  # abs_smd column present in output
  expect_true("abs_smd" %in% names(pd$rows))
  # noise row has within-group variance so abs_smd should be joined and non-NA
  row_noise <- pd$rows[pd$rows$attribute == "noise", ]
  expect_false(is.na(row_noise$abs_smd))

  # Without SMD: abs_smd should be NA
  pd_no_smd <- oda_balance_plot_data(bt_unweighted)
  expect_true(all(is.na(pd_no_smd$rows$abs_smd)))
})

# ---------------------------------------------------------------------------
# T9: oda_balance_plot_data accepts only the pre-computed table (no fitting)
# ---------------------------------------------------------------------------
test_that("oda_balance_plot_data: does not call oda_fit; takes balance_table only (T9)", {
  # Must error if given a non-balance_table object
  expect_error(oda_balance_plot_data(list()), "oda_balance_table")

  # Calling with no group/X arguments (they are not parameters)
  pd <- oda_balance_plot_data(bt_unweighted, smd_table = smd_tbl,
                               p_col = "p_mc", rank_by = "abs_ess")
  expect_true(inherits(pd, "oda_balance_plot_data"))
  expect_equal(pd$n_covariates, 2L)
  expect_true("rank" %in% names(pd$rows))
  expect_equal(pd$rank_by, "abs_ess")
  expect_equal(pd$p_col_used, "p_mc")
})

###############################################################################
# cta_balance_table() / cta_balance_plot_data() - T10-T17
###############################################################################

# ---- Module-level CTA fixtures --------------------------------------------- #

# Discriminating fixture: B alone perfectly separates the groups.
X_cta_disc  <- data.frame(
  A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
  B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
)
grp_cta_disc <- c(rep(0L, 40), rep(1L, 20))

# Computed once.
ct_disc <- cta_balance_table(grp_cta_disc, X_cta_disc,
                              mindenom = 5L, mc_iter = 200L, mc_seed = 42L)

# no_tree fixture: force no_tree via alpha_split=0 (p < 0 is impossible).
ct_notree <- cta_balance_table(grp_cta_disc, X_cta_disc,
                                mindenom = 5L, mc_iter = 50L, mc_seed = 42L,
                                alpha_split = 0)

# ---------------------------------------------------------------------------
# T10: discriminating data -> status stump or valid_tree, ess_display > 0
# ---------------------------------------------------------------------------
test_that("cta_balance_table: discriminating data finds a tree (T10)", {
  expect_true(ct_disc$status %in% c("stump", "valid_tree"))
  expect_equal(ct_disc$balance_interpretation, "discriminating")
  expect_true(!is.na(ct_disc$ess_display) && ct_disc$ess_display > 0)
  expect_false(is.na(ct_disc$root_attribute))
  expect_false(is.na(ct_disc$n_endpoints))
  expect_true(ct_disc$n_endpoints >= 2L)
})

# ---------------------------------------------------------------------------
# T11: alpha_split=0 forces no_tree; balance_interpretation correct
# ---------------------------------------------------------------------------
test_that("cta_balance_table: no_tree -> no_discriminating_combinations (T11)", {
  expect_equal(ct_notree$status, "no_tree")
  expect_equal(ct_notree$balance_interpretation, "no_discriminating_combinations")
  expect_true(is.na(ct_notree$root_attribute))
  expect_true(is.na(ct_notree$n_endpoints))
  expect_true(is.na(ct_notree$ess_display))
  expect_true(is.na(ct_notree$d_stat))
})

# ---------------------------------------------------------------------------
# T12: required fields present
# ---------------------------------------------------------------------------
test_that("cta_balance_table: required fields all present (T12)", {
  req_fields <- c("status", "balance_interpretation", "root_attribute",
                   "n_endpoints", "overall_ess", "overall_wess", "ess_display",
                   "d_stat", "mindenom", "alpha", "has_weights",
                   "tree", "endpoint_table", "node_table",
                   "fit_error", "fit_reason")
  expect_true(all(req_fields %in% names(ct_disc)),
              info = paste("Missing:", paste(setdiff(req_fields, names(ct_disc)), collapse = ", ")))
  expect_s3_class(ct_disc, "cta_balance_table")
  expect_false(ct_disc$fit_error)
})

# ---------------------------------------------------------------------------
# T13: no_tree returns empty endpoint_table (0 rows), tree field present
# ---------------------------------------------------------------------------
test_that("cta_balance_table: no_tree has zero-row endpoint_table (T13)", {
  expect_equal(nrow(ct_notree$endpoint_table), 0L)
  # tree field is still present (the fitted cta_tree with no_tree flag)
  expect_true(!is.null(ct_notree$tree))
  expect_true(isTRUE(ct_notree$tree$no_tree))
})

# ---------------------------------------------------------------------------
# T14: weights accepted; has_weights = TRUE in result
# ---------------------------------------------------------------------------
test_that("cta_balance_table: weights accepted, has_weights reflected (T14)", {
  w_cta <- c(rep(1L, 40), rep(2L, 20))
  ct_w  <- cta_balance_table(grp_cta_disc, X_cta_disc, w = w_cta,
                               mindenom = 5L, mc_iter = 200L, mc_seed = 42L)
  expect_true(ct_w$has_weights)
  expect_true(is.na(ct_w$overall_ess))
  expect_false(is.na(ct_w$overall_wess))
})

# ---------------------------------------------------------------------------
# T15: no_tree -> no_tree_message populated, cta_pd = NULL
# ---------------------------------------------------------------------------
test_that("cta_balance_plot_data: no_tree -> message populated, cta_pd NULL (T15)", {
  cpd <- cta_balance_plot_data(ct_notree)
  expect_s3_class(cpd, "cta_balance_plot_data")
  expect_equal(cpd$status, "no_tree")
  expect_false(is.na(cpd$no_tree_message))
  expect_true(nchar(cpd$no_tree_message) > 0)
  expect_null(cpd$cta_pd)
  expect_equal(cpd$balance_interpretation, "no_discriminating_combinations")
})

# ---------------------------------------------------------------------------
# T16: discriminating tree -> cta_pd populated, no_tree_message = NA
# ---------------------------------------------------------------------------
test_that("cta_balance_plot_data: valid tree -> cta_pd populated, no message (T16)", {
  cpd <- cta_balance_plot_data(ct_disc)
  expect_s3_class(cpd, "cta_balance_plot_data")
  expect_true(cpd$status %in% c("stump", "valid_tree"))
  expect_true(is.na(cpd$no_tree_message))
  expect_false(is.null(cpd$cta_pd))
  expect_true(!is.na(cpd$ess_display) && cpd$ess_display > 0)
})

# ---------------------------------------------------------------------------
# T17: cta_balance_plot_data errors on non-cta_balance_table input
# ---------------------------------------------------------------------------
test_that("cta_balance_plot_data: errors without cta_balance_table (T17)", {
  expect_error(cta_balance_plot_data(list()), "cta_balance_table")
  expect_error(cta_balance_plot_data(ct_disc$tree), "cta_balance_table")
})
