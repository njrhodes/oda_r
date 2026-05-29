# tests/testthat/test-live-fire.R
# Slice T — Live-fire public API release rehearsal
#
# End-to-end operator workflows using tiny synthetic data.
# Each test_that block is a complete workflow, not an isolated unit test.
# Assert object shapes, field presence, and non-NA evidence — not exact values.
#
# All tests skip at cran tier; gated at fast tier and above.

LIVEFIRE_SKIP <- Sys.getenv("ODACORE_TEST_TIER") %in% c("", "cran")
SKIP_MSG      <- "live-fire workflow: skip at cran tier"

# ---------------------------------------------------------------------------
# Synthetic data generator (shared across workflows)
# ---------------------------------------------------------------------------

local({
  set.seed(42L)
  n <- 60L

  # Binary outcome
  lf_y <<- c(rep(0L, 30L), rep(1L, 30L))

  # Three ordered/binary attributes (all have signal)
  lf_X <<- data.frame(
    V1 = c(rnorm(30L, mean = -1), rnorm(30L, mean = 1)),
    V2 = c(rnorm(30L, mean =  0), rnorm(30L, mean = 2)),
    V3 = sample(c(0L, 1L), n, replace = TRUE, prob = c(0.4, 0.6))
  )

  # Case weights (mild variation)
  lf_w <<- c(rep(1.0, 20L), rep(1.5, 20L), rep(2.0, 20L))
})

# ---------------------------------------------------------------------------
# T1: ODA workflow
# ---------------------------------------------------------------------------

test_that("T1: ODA operator workflow completes end-to-end", {
  skip_if(LIVEFIRE_SKIP, SKIP_MSG)

  # --- Step 1: Readiness check ---
  rc <- oda_readiness_check(lf_X, lf_y, w = lf_w)
  expect_type(rc, "list")
  expect_true("ok" %in% names(rc))

  # --- Step 2: Fit binary ODA ---
  fit <- oda_fit(lf_X[["V1"]], lf_y, attr_type = "ordered",
                 mcarlo = FALSE, loo = "off")
  expect_s3_class(fit, "oda_fit")
  expect_true(fit$ok)
  expect_s3_class(fit, "oda_fit_binary")

  # --- Step 3: Summary / print ---
  sm <- summary(fit)
  expect_s3_class(sm, "oda_fit_summary")
  out <- capture.output(print(fit))
  expect_true(length(out) > 0L)

  # --- Step 4: Evidence accessors ---
  conf <- oda_confusion(fit)
  expect_true(all(c("TP", "TN", "FP", "FN") %in% names(conf)))
  expect_true(sum(c(conf$TP, conf$TN, conf$FP, conf$FN)) == 60L)

  met <- oda_metrics(fit)
  expect_true("ess" %in% names(met))
  expect_false(is.na(met$ess))

  d_val <- oda_d_stat(fit)
  expect_true(is.numeric(d_val))

  # --- Step 5: Predict ---
  preds <- predict(fit, lf_X[["V1"]])
  expect_length(preds, 60L)
  expect_true(all(preds %in% c(0L, 1L, NA_integer_)))

  # --- Step 6: NOVOboot ---
  ci <- novo_boot_ci(fit, nboot = 100L, seed = 1L)
  expect_s3_class(ci, "novo_boot_ci")
  expect_equal(ci$n, 60L)
  expect_true(is.logical(ci$significant))
})

# ---------------------------------------------------------------------------
# T2: CTA workflow
# ---------------------------------------------------------------------------

test_that("T2: CTA operator workflow completes end-to-end", {
  skip_if(LIVEFIRE_SKIP, SKIP_MSG)

  # --- Step 1: Readiness check ---
  rc <- oda_readiness_check(lf_X, lf_y)
  expect_true(rc$ok)

  # --- Step 2: Fit CTA ---
  tree <- cta_fit(lf_X, lf_y, mc_iter = 200L, mc_seed = 1L,
                  loo = "off", mindenom = 5L)
  expect_s3_class(tree, "cta_tree")

  # --- Step 3: Summary / print ---
  sm <- summary(tree)
  out <- capture.output(print(tree))
  expect_true(length(out) > 0L)

  # --- Step 4: Translation functions ---
  nt <- cta_node_table(tree)
  expect_s3_class(nt, "data.frame")

  st <- cta_strata(tree)
  # st is an integer (strata count) or NA for no_tree
  expect_true(is.integer(st) || is.na(st))

  d_val <- cta_d_stat(tree)
  expect_true(is.numeric(d_val))

  ct <- cta_confusion_table(tree)
  expect_s3_class(ct, "data.frame")

  # --- Step 5: predict ---
  preds <- predict(tree, lf_X, missing_action = "na")
  expect_length(preds, 60L)

  # --- Step 6: Plot-data evidence fields ---
  pd <- cta_plot_data(tree)
  expect_true(all(c("overall_ess", "ess_label", "d", "model_label", "training_n")
                  %in% names(pd)))
  expect_equal(pd$model_label, "CTA")
  if (!isTRUE(pd$no_tree)) {
    expect_false(is.na(pd$overall_ess))
    expect_true(pd$ess_label %in% c("ESS", "WESS"))
  }

  # --- Step 7: NOVOboot (skip for no_tree) ---
  if (!isTRUE(tree$no_tree)) {
    ci <- novo_boot_ci(tree, nboot = 100L, seed = 1L)
    expect_s3_class(ci, "novo_boot_ci")
    expect_equal(ci$confusion, tree$training_confusion)
  }

  # --- Step 8: ggplot2 renderer with show_metrics ---
  skip_if(!requireNamespace("ggplot2", quietly = TRUE), "ggplot2 not available")
  if (!isTRUE(tree$no_tree)) {
    p <- plot_cta_tree(tree, show_metrics = TRUE, color_by = "none")
    expect_s3_class(p, "ggplot")
    # Subtitle should contain ESS label
    if (!is.na(pd$overall_ess)) {
      expect_true(grepl("ESS|WESS", p$labels$subtitle %||% ""))
    }
  }
})

# ---------------------------------------------------------------------------
# T3: LORT workflow
# ---------------------------------------------------------------------------

test_that("T3: LORT operator workflow completes end-to-end", {
  skip_if(LIVEFIRE_SKIP, SKIP_MSG)

  # --- Step 1: Fit LORT ---
  lort <- lort_fit(lf_X, lf_y, mc_iter = 200L, mc_seed = 1L,
                   loo = "off", min_n = 5L)
  expect_s3_class(lort, "cta_ort")
  expect_s3_class(lort, "cta_tree")

  # Confirm LORT metadata invariants
  expect_equal(lort$ort_settings$method, "lort")
  expect_false(lort$ort_settings$global_optimization)
  expect_false(lort$ort_settings$sda_anchored)

  # --- Step 2: Summary / print ---
  sm <- summary(lort)
  out <- capture.output(print(lort))
  expect_true(length(out) > 0L)

  # --- Step 3: Plot-data evidence fields ---
  pd <- ort_plot_data(lort)
  expect_true(all(c("overall_ess", "ess_label", "d", "model_label", "training_n")
                  %in% names(pd)))
  expect_equal(pd$model_label, "LORT")
  expect_true(is.numeric(pd$overall_ess))

  # --- Step 4: NOVOboot from strata ---
  st <- lort$strata
  if (!is.null(st) && nrow(st) > 0L) {
    ci <- novo_boot_ci(lort, nboot = 100L, seed = 1L)
    expect_s3_class(ci, "novo_boot_ci")
    total_from_strata <- sum(vapply(st$class_counts, sum, integer(1L)))
    expect_equal(ci$n, total_from_strata)
  }

  # --- Step 5: ggplot2 renderer with show_metrics ---
  skip_if(!requireNamespace("ggplot2", quietly = TRUE), "ggplot2 not available")
  p <- plot_lort_tree(lort, show_metrics = TRUE, color_by = "none")
  expect_s3_class(p, "ggplot")
  if (!is.na(pd$overall_ess)) {
    expect_true(grepl("ESS|WESS", p$labels$subtitle %||% ""))
  }
})

# ---------------------------------------------------------------------------
# T4: Balance + propensity workflow
# ---------------------------------------------------------------------------

test_that("T4: Balance and propensity workflow completes end-to-end", {
  skip_if(LIVEFIRE_SKIP, SKIP_MSG)

  # --- Step 1: ODA balance table (group vs covariates) ---
  # group = lf_y (binary), X = lf_X (3 covariates)
  bal_oda <- oda_balance_table(lf_y, lf_X, w = lf_w,
                                mc_iter = 100L, mc_seed = 1L)
  expect_s3_class(bal_oda, "oda_balance_table")
  expect_s3_class(bal_oda, "data.frame")
  expect_true(nrow(bal_oda) == 3L)  # one row per attribute
  expect_true(all(c("attribute", "ess", "p_mc") %in% names(bal_oda)))

  # --- Step 2: SMD balance table ---
  smd_tbl <- smd_balance_table(lf_y, lf_X, w = lf_w)
  expect_s3_class(smd_tbl, "smd_balance_table")
  expect_true(nrow(smd_tbl) == 3L)
  expect_true("smd" %in% names(smd_tbl))

  # --- Step 3: ODA balance plot data ---
  bpd <- oda_balance_plot_data(bal_oda)
  expect_s3_class(bpd, "oda_balance_plot_data")
  expect_true(is.data.frame(bpd$rows))

  # --- Step 4: ODA propensity weights ---
  fit_oda <- oda_fit(lf_X[["V1"]], lf_y, attr_type = "ordered",
                     mcarlo = FALSE, loo = "off")
  expect_true(fit_oda$ok)
  prop_oda <- oda_propensity_weights(fit_oda)
  expect_s3_class(prop_oda, "data.frame")
  expect_equal(nrow(prop_oda), 4L)  # 2 strata x 2 classes
  expect_equal(prop_oda$model_family[1L], "oda")
  expect_true("weight" %in% names(prop_oda))
  expect_true(all(is.finite(prop_oda$weight[prop_oda$weight != Inf])))

  # --- Step 5: LORT propensity weights ---
  lort <- lort_fit(lf_X, lf_y, mc_iter = 200L, mc_seed = 1L,
                   loo = "off", min_n = 5L)
  prop_lort <- lort_propensity_weights(lort, target_class = 1L)
  expect_s3_class(prop_lort, "data.frame")
  expect_true(nrow(prop_lort) > 0L)
  expect_equal(prop_lort$model_family[1L], "lort")
  expect_false(prop_lort$global_optimization[1L])
  expect_false(prop_lort$sda_anchored[1L])
  expect_true("weight" %in% names(prop_lort))

  # --- Step 6: CTA multivariate balance ---
  bal_cta <- cta_balance_table(lf_y, lf_X, w = lf_w,
                                mc_iter = 100L, mc_seed = 1L, mindenom = 5L)
  expect_s3_class(bal_cta, "cta_balance_table")
  expect_true("status" %in% names(bal_cta))
  expect_true(bal_cta$status %in% c("valid_tree", "stump", "no_tree", "fit_error"))
})

# ---------------------------------------------------------------------------
# T5: SDA anchor workflow
# ---------------------------------------------------------------------------

test_that("T5: SDA anchor operator workflow completes end-to-end", {
  skip_if(LIVEFIRE_SKIP, SKIP_MSG)

  # --- Step 1: Fit SDA ---
  sda <- sda_fit(lf_X, lf_y, mc_iter = 100L, mc_seed = 1L,
                 mode = "novometric_min_d")
  expect_s3_class(sda, "sda_fit")

  # --- Step 2: SDA accessors ---
  sel <- sda_selected_attributes(sda)
  expect_type(sel, "character")
  expect_true(length(sel) >= 0L)  # may be 0 if none significant

  step_tbl <- sda_step_table(sda)
  expect_s3_class(step_tbl, "data.frame")

  # --- Step 3: Convert to anchor ---
  anchor <- as_sda_anchor(sda)
  expect_s3_class(anchor, "sda_anchor")

  # --- Step 4: Validate anchor ---
  vr <- validate_sda_anchor(anchor)
  expect_type(vr, "list")
  expect_true("ok" %in% names(vr))

  # --- Step 5: Anchor invariants ---
  expect_true("prohibited_downstream" %in% names(anchor))
  expect_true("propensity_weighting" %in% anchor$prohibited_downstream)
  expect_true("fraud_demo"           %in% anchor$prohibited_downstream)

  # --- Step 6: print / summary ---
  out_print   <- capture.output(print(anchor))
  out_summary <- capture.output(summary(anchor))
  expect_true(length(out_print) > 0L)
  expect_true(length(out_summary) > 0L)

  # --- Step 7: Explicit anchor path (manual construction) ---
  if (length(sel) > 0L) {
    manual_anchor <- sda_anchor(
      selected_attributes = sel,
      candidate_universe  = names(lf_X),
      group_levels        = c(0L, 1L)
    )
    expect_s3_class(manual_anchor, "sda_anchor")
    expect_true("propensity_weighting" %in% manual_anchor$prohibited_downstream)
  }
})
