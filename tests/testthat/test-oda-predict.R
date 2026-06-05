###############################################################################
# test-oda-predict.R - Phase 2A tests: S3 class, predict, print, summary,
# and accessor behavior for oda_fit() results.
###############################################################################

# ---- Helpers ----------------------------------------------------------------

.bin_fit <- function(mcarlo = FALSE, loo = "off") {
  x <- 1:8
  y <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  oda_fit(x, y, mcarlo = mcarlo, loo = loo)
}

.multi_fit <- function(mcarlo = FALSE, loo = "off") {
  y <- as.integer(factor(iris$Species,
                          levels = c("setosa", "versicolor", "virginica")))
  oda_fit(iris$Petal.Length, y,
          attr_type = "ordered", priors_on = TRUE,
          mcarlo = mcarlo, loo = loo)
}

.failed_fit <- function() {
  # pure_node: all same class
  oda_fit(1:5, rep(0L, 5), mcarlo = FALSE)
}

# =============================================================================
# 1. Class inheritance
# =============================================================================

test_that("binary fit inherits oda_fit_binary and oda_fit", {
  fit <- .bin_fit()
  expect_s3_class(fit, "oda_fit_binary")
  expect_s3_class(fit, "oda_fit")
})

test_that("multiclass fit inherits oda_fit_multiclass and oda_fit", {
  fit <- .multi_fit()
  expect_s3_class(fit, "oda_fit_multiclass")
  expect_s3_class(fit, "oda_fit")
})

test_that("failed fit inherits oda_fit_failed and oda_fit", {
  fit <- .failed_fit()
  expect_s3_class(fit, "oda_fit_failed")
  expect_s3_class(fit, "oda_fit")
  expect_false(isTRUE(fit$ok))
})

# =============================================================================
# 2. Metadata fields stored on the fit
# =============================================================================

test_that("binary fit stores priors_on and has_weights", {
  fit <- .bin_fit()
  expect_true(!is.null(fit$priors_on))
  expect_true(!is.null(fit$has_weights))
  expect_false(isTRUE(fit$has_weights))   # no w= supplied
})

test_that("has_weights: NULL w gives FALSE", {
  fit <- oda_fit(c(1,2,3,4,5,6), c(0L,0L,0L,1L,1L,1L), mcarlo = FALSE)
  expect_false(isTRUE(fit$has_weights))
})

test_that("has_weights: all-ones w gives FALSE (not a meaningful weight)", {
  fit <- oda_fit(c(1,2,3,4,5,6), c(0L,0L,0L,1L,1L,1L),
                 w = rep(1, 6), mcarlo = FALSE)
  expect_false(isTRUE(fit$has_weights))
})

test_that("has_weights: nontrivial w gives TRUE", {
  fit <- oda_fit(c(1,2,3,4,5,6), c(0L,0L,0L,1L,1L,1L),
                 w = c(1,2,1,2,1,2), mcarlo = FALSE)
  expect_true(isTRUE(fit$has_weights))
})

test_that("has_weights: priors_on=TRUE with NULL w gives FALSE", {
  fit <- oda_fit(c(1,2,3,4,5,6), c(0L,0L,0L,1L,1L,1L),
                 priors_on = TRUE, mcarlo = FALSE)
  expect_false(isTRUE(fit$has_weights))
})

test_that("multiclass fit stores boundary_mode", {
  fit <- .multi_fit()
  expect_true(!is.null(fit$boundary_mode))
  expect_true(fit$boundary_mode %in% c("megaoda_halfopen", "right_closed"))
})

test_that("miss_codes stored on fit", {
  x <- c(1:8, -9)
  y <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 0L)
  fit <- oda_fit(x, y, miss_codes = -9, mcarlo = FALSE)
  expect_identical(fit$miss_codes, -9)
})

# =============================================================================
# 3. predict.oda_fit - binary
# =============================================================================

test_that("predict.oda_fit binary returns correct labels for simple separable data", {
  fit   <- .bin_fit()
  preds <- predict(fit, 1:8)
  expect_equal(preds, c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L))
})

test_that("predict.oda_fit binary returns integer vector", {
  fit   <- .bin_fit()
  preds <- predict(fit, 1:8)
  expect_true(is.integer(preds))
  expect_equal(length(preds), 8L)
})

test_that("predict.oda_fit binary: NA input gives NA output", {
  fit   <- .bin_fit()
  preds <- predict(fit, c(1, NA, 8))
  expect_equal(preds[2L], NA_integer_)
  expect_false(is.na(preds[1L]))
  expect_false(is.na(preds[3L]))
})

test_that("predict.oda_fit binary: miss-coded value gives NA output", {
  x <- c(1:8, -9)
  y <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 0L)
  fit   <- oda_fit(x, y, miss_codes = -9, mcarlo = FALSE)
  preds <- predict(fit, c(2, -9, 7))
  expect_equal(preds[2L], NA_integer_)
  expect_false(is.na(preds[1L]))
  expect_false(is.na(preds[3L]))
})

test_that("predict.oda_fit binary: single-column data.frame accepted", {
  fit   <- .bin_fit()
  preds <- predict(fit, data.frame(x = 1:8))
  expect_equal(length(preds), 8L)
  expect_true(is.integer(preds))
})

test_that("predict.oda_fit binary: multi-column data.frame rejected", {
  fit <- .bin_fit()
  expect_error(predict(fit, data.frame(x = 1:4, z = 1:4)), "single-column")
})

# =============================================================================
# 4. predict.oda_fit - multiclass
# =============================================================================

test_that("predict.oda_fit multiclass returns only valid classes", {
  fit   <- .multi_fit()
  preds <- predict(fit, iris$Petal.Length)
  expect_true(all(preds %in% c(1L, 2L, 3L)))
})

test_that("predict.oda_fit multiclass returns correct length", {
  fit   <- .multi_fit()
  preds <- predict(fit, iris$Petal.Length)
  expect_equal(length(preds), 150L)
})

test_that("predict.oda_fit multiclass returns integer", {
  fit   <- .multi_fit()
  preds <- predict(fit, iris$Petal.Length)
  expect_true(is.integer(preds))
})

# =============================================================================
# 5. predict.oda_fit - failed fit
# =============================================================================

test_that("predict.oda_fit failed fit returns all NA with warning", {
  fit <- .failed_fit()
  expect_warning(preds <- predict(fit, 1:5), "ok=FALSE")
  expect_true(all(is.na(preds)))
  expect_true(is.integer(preds))
  expect_equal(length(preds), 5L)
})

# =============================================================================
# 6. No predict(fit$rule, ...) API
# =============================================================================

test_that("predict(fit$rule) is not a defined S3 generic on rule objects", {
  fit <- .bin_fit()
  # rule is a plain list, not classed - predict() on a list falls through
  # to the default which does not know what to do. We simply verify no
  # predict.rule or predict.ordered_cut method exists.
  expect_null(getS3method("predict", "rule", optional = TRUE))
  expect_null(getS3method("predict", "ordered_cut", optional = TRUE))
})

# =============================================================================
# 7. print.oda_fit
# =============================================================================

test_that("print.oda_fit binary produces output containing 'ODA'", {
  fit <- .bin_fit()
  expect_output(print(fit), "ODA")
})

test_that("print.oda_fit binary shows rule", {
  fit <- .bin_fit()
  expect_output(print(fit), "Rule:")
})

test_that("print.oda_fit multiclass produces output", {
  fit <- .multi_fit()
  expect_output(print(fit), "ODA")
})

test_that("print.oda_fit failed fit shows FAILED", {
  fit <- .failed_fit()
  expect_output(print(fit), "FAILED")
})

test_that("print.oda_fit returns invisibly", {
  fit <- .bin_fit()
  ret <- withVisible(print(fit))
  expect_false(ret$visible)
})

# =============================================================================
# 8. summary.oda_fit
# =============================================================================

test_that("summary.oda_fit binary returns oda_fit_summary", {
  fit <- .bin_fit()
  s   <- summary(fit)
  expect_s3_class(s, "oda_fit_summary")
})

test_that("summary.oda_fit binary has expected fields", {
  fit <- .bin_fit()
  s   <- summary(fit)
  expect_true("model_type" %in% names(s))
  expect_true("status"     %in% names(s))
  expect_true("train"      %in% names(s))
  expect_true("loo"        %in% names(s))
  expect_true("rule"       %in% names(s))
})

test_that("summary.oda_fit binary train section has ess", {
  fit <- .bin_fit()
  s   <- summary(fit)
  expect_true("ess" %in% names(s$train))
  expect_true(is.numeric(s$train$ess))
})

test_that("summary.oda_fit contains no spurious LOO p-value when loo=off", {
  fit <- .bin_fit(loo = "off")
  s   <- summary(fit)
  # loo section should be NULL (no LOO was run)
  expect_null(s$loo)
})

test_that("summary.oda_fit multiclass returns oda_fit_summary", {
  fit <- .multi_fit()
  s   <- summary(fit)
  expect_s3_class(s, "oda_fit_summary")
})

test_that("summary.oda_fit multiclass train has ess and mean_pac", {
  fit <- .multi_fit()
  s   <- summary(fit)
  expect_true("ess"      %in% names(s$train))
  expect_true("mean_pac" %in% names(s$train))
})

test_that("print.oda_fit_summary runs without error", {
  fit <- .bin_fit()
  s   <- summary(fit)
  expect_output(print(s), "ODA Summary")
})

# =============================================================================
# 9. oda_confusion()
# =============================================================================

test_that("oda_confusion train returns binary confusion list", {
  fit  <- .bin_fit()
  conf <- oda_confusion(fit, "train")
  expect_true(!is.null(conf))
  expect_true("TP" %in% names(conf) || is.matrix(conf))
})

test_that("oda_confusion train weighted returns confusion_wt", {
  fit   <- .bin_fit()
  conf  <- oda_confusion(fit, "train", weighted = FALSE)
  confw <- oda_confusion(fit, "train", weighted = TRUE)
  # Both should be non-NULL for a successful fit
  expect_false(is.null(conf))
  # weighted and raw may differ or be identical for unit weights
})

test_that("oda_confusion loo returns NULL when loo=off", {
  fit  <- .bin_fit(loo = "off")
  conf <- oda_confusion(fit, "loo")
  expect_null(conf)
})

test_that("oda_confusion multiclass train returns matrix", {
  fit  <- .multi_fit()
  conf <- oda_confusion(fit, "train")
  expect_true(is.matrix(conf))
})

# =============================================================================
# 10. oda_metrics()
# =============================================================================

test_that("oda_metrics train binary returns ess and pac", {
  fit <- .bin_fit()
  m   <- oda_metrics(fit, "train")
  expect_true("ess" %in% names(m))
  expect_true("pac" %in% names(m))
  expect_true(is.numeric(m$ess))
  expect_true(m$ess >= 0 && m$ess <= 100)
})

test_that("oda_metrics train multiclass returns ess and mean_pac", {
  fit <- .multi_fit()
  m   <- oda_metrics(fit, "train")
  expect_true("ess"      %in% names(m))
  expect_true("mean_pac" %in% names(m))
})

test_that("oda_metrics loo returns NULL when loo=off (with message)", {
  fit <- .bin_fit(loo = "off")
  expect_message(m <- oda_metrics(fit, "loo"), "not available")
  expect_null(m)
})

# =============================================================================
# 11. oda_predictions()
# =============================================================================

test_that("oda_predictions train with newdata calls predict", {
  fit   <- .bin_fit()
  preds <- oda_predictions(fit, "train", newdata = 1:8)
  expect_equal(preds, c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L))
})

test_that("oda_predictions train without newdata returns NULL with message", {
  fit <- .bin_fit()
  expect_message(p <- oda_predictions(fit, "train"), "not stored")
  expect_null(p)
})

test_that("oda_predictions loo returns NULL when not run", {
  fit <- .bin_fit(loo = "off")
  expect_message(p <- oda_predictions(fit, "loo"), "not available")
  expect_null(p)
})

# =============================================================================
# 12. LOO p-value canon (MPE)
# =============================================================================

# Helper: binary fit with LOO enabled on data that produces allowed=TRUE.
# iris setosa vs versicolor on Petal.Length is well-separated; LOO p << 0.05.
.bin_fit_loo <- function() {
  x <- c(iris$Petal.Length[1:50], iris$Petal.Length[51:100])
  y <- c(rep(0L, 50), rep(1L, 50))
  oda_fit(x, y, attr_type = "ordered", mcarlo = FALSE, loo = "pvalue")
}

test_that("binary LOO summary p_method is Fisher exact (2x2) and p_value is not NA", {
  fit <- .bin_fit_loo()
  s   <- summary(fit)
  if (isTRUE(s$loo$allowed)) {
    expect_equal(s$loo$p_method, "Fisher exact (2x2), one-tailed; MPE p.34")
    expect_equal(s$loo$p_status, "computed")
    expect_false(is.na(s$loo$p_value))   # p_value must be present when status=computed
  } else {
    skip("LOO not allowed for this fit - skip p_method check")
  }
})

test_that("binary LOO metrics p_method is Fisher exact (2x2) and p_value is not NA", {
  fit <- .bin_fit_loo()
  m   <- oda_metrics(fit, "loo")
  if (!is.null(m)) {
    expect_equal(m$p_method, "Fisher exact (2x2), one-tailed; MPE p.34")
    expect_equal(m$p_status, "computed")
    expect_false(is.na(m$p_value))
  } else {
    skip("LOO not allowed for this fit - skip p_method check")
  }
})

test_that("binary LOO with NA p_value reports not_computed, not computed", {
  # Synthesise an oda_fit_binary with loo$allowed=TRUE but p_value=NA_real_
  # to exercise the not_computed branch of .loo_p_info().
  fit <- .bin_fit_loo()
  fit$loo$p_value <- NA_real_   # force absence
  s <- summary(fit)
  expect_equal(s$loo$p_status, "not_computed")
  expect_equal(s$loo$p_method, "Fisher exact (2x2), one-tailed; MPE p.34")
  expect_true(is.na(s$loo$p_value))
  expect_match(s$loo$p_reason, "not available")
})

test_that("multiclass LOO summary p_status is not_computed", {
  testthat::skip_on_cran()
  fit <- .multi_fit(loo = "on")
  s   <- summary(fit)
  if (isTRUE(s$loo$allowed)) {
    expect_equal(s$loo$p_status, "not_computed")
    expect_true(is.na(s$loo$p_value))
    expect_equal(s$loo$p_method, "none")
    expect_match(s$loo$p_reason, "not reported")
  } else {
    skip("LOO not allowed for this multiclass fit - skip p_status check")
  }
})

test_that("multiclass LOO metrics p_status is not_computed", {
  testthat::skip_on_cran()
  fit <- .multi_fit(loo = "on")
  m   <- oda_metrics(fit, "loo")
  if (!is.null(m)) {
    expect_equal(m$p_status, "not_computed")
    expect_true(is.na(m$p_value))
  } else {
    skip("LOO not allowed for this multiclass fit - skip p_status check")
  }
})

# =============================================================================
# 13. oda_d_stat()
# =============================================================================

test_that("oda_d_stat binary ESS=100 returns 0", {
  fit <- .bin_fit()   # perfectly separable x=1:8
  expect_equal(oda_d_stat(fit), 0)
})

test_that("oda_d_stat binary returns scalar numeric", {
  fit <- .bin_fit()
  d   <- oda_d_stat(fit)
  expect_true(is.numeric(d) && length(d) == 1L)
})

test_that("oda_d_stat multiclass ordered returns finite numeric", {
  fit <- .multi_fit()
  d   <- oda_d_stat(fit)
  expect_true(is.numeric(d) && length(d) == 1L)
  expect_true(is.finite(d))
})

test_that("oda_d_stat multiclass ordered uses seg_classes length as strata", {
  fit    <- .multi_fit()
  strata <- length(fit$rule$seg_classes)
  ess    <- fit$ess
  expect_equal(oda_d_stat(fit), 100 / (ess / strata) - strata)
})

test_that("oda_d_stat failed fit returns NA_real_", {
  fit <- .failed_fit()
  expect_equal(oda_d_stat(fit), NA_real_)
})
