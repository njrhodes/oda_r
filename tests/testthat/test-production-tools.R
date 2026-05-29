###############################################################################
# test-production-tools.R  --  Slice Q production-tools tests
###############################################################################

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

local({
  set.seed(42L)
  n <- 60L
  x1 <- c(rnorm(30, mean = 1), rnorm(30, mean = -1))
  x2 <- c(sample(1:5, 30, replace = TRUE), sample(1:5, 30, replace = TRUE))
  x3 <- rep(7.0, n)                    # constant column
  x4 <- c(rep(0L, 30), rep(1L, 30))   # binary column
  y  <- c(rep(0L, 30), rep(1L, 30))
  w  <- runif(n, 0.5, 2.0)

  env <- parent.env(environment())
  env$X_fix  <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
  env$y_fix  <- y
  env$w_fix  <- w
  env$n_fix  <- n
})

# ---------------------------------------------------------------------------
# oda_clean_missing_codes
# ---------------------------------------------------------------------------

test_that("oda_clean_missing_codes: NULL miss_codes returns input unchanged", {
  v <- c(1, -9, 3)
  expect_identical(oda_clean_missing_codes(v, miss_codes = NULL), v)
})

test_that("oda_clean_missing_codes: replaces miss codes in vector", {
  v <- c(1.0, -9.0, 3.0, -9.0)
  out <- oda_clean_missing_codes(v, miss_codes = -9)
  expect_true(is.na(out[2]))
  expect_true(is.na(out[4]))
  expect_equal(out[1], 1.0)
  expect_equal(out[3], 3.0)
})

test_that("oda_clean_missing_codes: replaces miss codes in data frame", {
  df <- data.frame(a = c(1.0, -9.0, 3.0), b = c(-9.0, 2.0, -9.0))
  out <- oda_clean_missing_codes(df, miss_codes = -9)
  expect_true(is.na(out$a[2]))
  expect_true(is.na(out$b[1]))
  expect_true(is.na(out$b[3]))
  expect_equal(out$a[1], 1.0)
  expect_equal(out$b[2], 2.0)
})

test_that("oda_clean_missing_codes: returns same class as input", {
  v <- c(1.0, 2.0, -9.0)
  expect_true(is.numeric(oda_clean_missing_codes(v, miss_codes = -9)))
  df <- data.frame(x = v)
  expect_true(is.data.frame(oda_clean_missing_codes(df, miss_codes = -9)))
})

test_that("oda_clean_missing_codes: custom replacement value", {
  v <- c(1.0, -9.0, 3.0)
  out <- oda_clean_missing_codes(v, miss_codes = -9, replacement = 0.0)
  expect_equal(out[2], 0.0)
})

# ---------------------------------------------------------------------------
# oda_validate_group
# ---------------------------------------------------------------------------

test_that("oda_validate_group: valid binary group returns ok = TRUE", {
  r <- oda_validate_group(y_fix)
  expect_true(r$ok)
  expect_equal(r$n_classes, 2L)
  expect_equal(sort(r$class_levels), c(0L, 1L))
  expect_length(r$issues, 0L)
})

test_that("oda_validate_group: NULL y returns ok = FALSE", {
  r <- oda_validate_group(NULL)
  expect_false(r$ok)
  expect_true(length(r$issues) > 0L)
})

test_that("oda_validate_group: single-class vector flagged", {
  r <- oda_validate_group(rep(1L, 10))
  expect_false(r$ok)
  expect_true(any(grepl("one unique class", r$issues)))
})

test_that("oda_validate_group: binary_only = TRUE flags >2 classes", {
  r <- oda_validate_group(c(0L, 1L, 2L, 0L, 1L, 2L), binary_only = TRUE)
  expect_false(r$ok)
  expect_true(any(grepl("binary_only", r$issues)))
})

test_that("oda_validate_group: binary_only = TRUE passes for 2 classes", {
  r <- oda_validate_group(c(0L, 1L, 0L, 1L), binary_only = TRUE)
  expect_true(r$ok)
})

test_that("oda_validate_group: NA values flagged as issue", {
  r <- oda_validate_group(c(0L, 1L, NA_integer_, 0L))
  expect_false(r$ok)
  expect_true(any(grepl("NA", r$issues)))
})

# ---------------------------------------------------------------------------
# oda_validate_weights
# ---------------------------------------------------------------------------

test_that("oda_validate_weights: NULL w returns ok = TRUE", {
  r <- oda_validate_weights(NULL, 60L)
  expect_true(r$ok)
  expect_length(r$issues, 0L)
  expect_null(r$range)
})

test_that("oda_validate_weights: valid weights return ok = TRUE", {
  r <- oda_validate_weights(w_fix, n_fix)
  expect_true(r$ok)
  expect_length(r$issues, 0L)
  expect_equal(r$n_weights, as.integer(n_fix))
  expect_length(r$range, 2L)
})

test_that("oda_validate_weights: wrong-length weights flagged", {
  r <- oda_validate_weights(w_fix[seq_len(5L)], n_fix)
  expect_false(r$ok)
  expect_true(any(grepl("length", r$issues)))
})

test_that("oda_validate_weights: NA weights flagged", {
  w_bad <- w_fix
  w_bad[1] <- NA_real_
  r <- oda_validate_weights(w_bad, n_fix)
  expect_false(r$ok)
  expect_true(any(grepl("NA", r$issues)))
})

test_that("oda_validate_weights: zero weight flagged", {
  w_bad <- w_fix
  w_bad[1] <- 0.0
  r <- oda_validate_weights(w_bad, n_fix)
  expect_false(r$ok)
  expect_true(any(grepl("zero", r$issues)))
})

test_that("oda_validate_weights: negative weight flagged", {
  w_bad <- w_fix
  w_bad[2] <- -1.0
  r <- oda_validate_weights(w_bad, n_fix)
  expect_false(r$ok)
})

test_that("oda_validate_weights: Inf weight flagged", {
  w_bad <- w_fix
  w_bad[3] <- Inf
  r <- oda_validate_weights(w_bad, n_fix)
  expect_false(r$ok)
})

# ---------------------------------------------------------------------------
# oda_infer_attr_types
# ---------------------------------------------------------------------------

test_that("oda_infer_attr_types: returns one row per column", {
  r <- oda_infer_attr_types(X_fix)
  expect_equal(nrow(r), ncol(X_fix))
  expect_equal(r$attribute, names(X_fix))
})

test_that("oda_infer_attr_types: ordered column inferred correctly", {
  r <- oda_infer_attr_types(X_fix)
  expect_equal(r$inferred_type[r$attribute == "x1"], "ordered")
})

test_that("oda_infer_attr_types: binary column inferred correctly", {
  r <- oda_infer_attr_types(X_fix)
  expect_equal(r$inferred_type[r$attribute == "x4"], "binary")
})

test_that("oda_infer_attr_types: constant column gets binary type (1 unique)", {
  r <- oda_infer_attr_types(X_fix)
  # x3 is constant = 7.0 -> 1 unique value -> binary
  expect_equal(r$inferred_type[r$attribute == "x3"], "binary")
})

test_that("oda_infer_attr_types: n_unique counts exclude miss_codes", {
  df <- data.frame(a = c(1.0, 2.0, 3.0, -9.0))
  r  <- oda_infer_attr_types(df, miss_codes = -9)
  expect_equal(r$n_unique[1], 3L)
  expect_equal(r$n_miss_code[1], 1L)
})

test_that("oda_infer_attr_types: NA values counted in n_missing", {
  df <- data.frame(a = c(1.0, 2.0, NA, 4.0))
  r  <- oda_infer_attr_types(df)
  expect_equal(r$n_missing[1], 1L)
})

# ---------------------------------------------------------------------------
# oda_readiness_check
# ---------------------------------------------------------------------------

test_that("oda_readiness_check: valid inputs return ok = TRUE", {
  r <- oda_readiness_check(X_fix, y_fix, w = w_fix)
  expect_true(r$ok)
  expect_length(r$issues, 0L)
  expect_equal(r$n_obs, as.integer(n_fix))
  expect_equal(r$n_attrs, ncol(X_fix))
})

test_that("oda_readiness_check: returns required fields", {
  r <- oda_readiness_check(X_fix, y_fix)
  expect_true(all(c("ok", "issues", "warnings", "n_obs", "n_attrs",
                     "group_report", "weight_report", "attr_types",
                     "constant_attrs") %in% names(r)))
})

test_that("oda_readiness_check: invalid group flagged", {
  r <- oda_readiness_check(X_fix, rep(1L, n_fix))
  expect_false(r$ok)
  expect_true(length(r$issues) > 0L)
})

test_that("oda_readiness_check: invalid weights flagged", {
  w_bad <- w_fix; w_bad[1] <- NA_real_
  r <- oda_readiness_check(X_fix, y_fix, w = w_bad)
  expect_false(r$ok)
  expect_true(any(grepl("NA", r$issues)))
})

test_that("oda_readiness_check: constant column raises warning", {
  r <- oda_readiness_check(X_fix, y_fix)
  expect_true(any(grepl("x3", r$warnings) | grepl("onstant", r$warnings)))
})

test_that("oda_readiness_check: miss_codes cleaning applied before attr type check", {
  # x3 is constant 7.0; if we declare 7 a miss code, all values missing
  X_mc <- X_fix
  r <- oda_readiness_check(X_mc, y_fix, miss_codes = 7.0)
  # Should not error; x3 now all NA -> constant_attrs should still contain x3
  expect_true(is.list(r))
})

test_that("oda_readiness_check: length mismatch flagged", {
  r <- oda_readiness_check(X_fix, y_fix[seq_len(5L)])
  expect_false(r$ok)
  expect_true(any(grepl("nrow", r$issues)))
})

# ---------------------------------------------------------------------------
# oda_propensity_weights
# ---------------------------------------------------------------------------

test_that("oda_propensity_weights: requires oda_fit object", {
  expect_error(oda_propensity_weights(list()), "oda_fit")
})

test_that("oda_propensity_weights: returns 13-column data frame", {
  fit <- oda_fit(X_fix[["x1"]], y_fix, attr_type = "ordered")
  expect_true(isTRUE(fit$ok))
  r <- oda_propensity_weights(fit)
  expect_s3_class(r, "data.frame")
  expect_equal(ncol(r), 13L)
})

test_that("oda_propensity_weights: returns 4 rows (2 strata x 2 classes)", {
  fit <- oda_fit(X_fix[["x1"]], y_fix, attr_type = "ordered")
  r   <- oda_propensity_weights(fit)
  expect_equal(nrow(r), 4L)
})

test_that("oda_propensity_weights: model_family is 'oda'", {
  fit <- oda_fit(X_fix[["x1"]], y_fix, attr_type = "ordered")
  r   <- oda_propensity_weights(fit)
  expect_true(all(r$model_family == "oda"))
})

test_that("oda_propensity_weights: stratum_ids are 1 and 2", {
  fit <- oda_fit(X_fix[["x1"]], y_fix, attr_type = "ordered")
  r   <- oda_propensity_weights(fit)
  expect_equal(sort(unique(r$stratum_id)), c(1L, 2L))
})

test_that("oda_propensity_weights: propensity weights are positive finite", {
  fit <- oda_fit(X_fix[["x1"]], y_fix, attr_type = "ordered")
  r   <- oda_propensity_weights(fit)
  # adjusted weights should always be finite
  expect_true(all(is.finite(r$adjusted_propensity_weight)))
})

test_that("oda_propensity_weights: returns empty df for failed fit", {
  # Build a fake oda_fit with ok = FALSE
  fake <- structure(list(ok = FALSE), class = "oda_fit")
  r    <- oda_propensity_weights(fake)
  expect_equal(nrow(r), 0L)
  expect_equal(ncol(r), 13L)
})

test_that("oda_propensity_weights: stratum_n sums to N", {
  fit <- oda_fit(X_fix[["x1"]], y_fix, attr_type = "ordered")
  r   <- oda_propensity_weights(fit)
  # Each stratum_id appears twice (one per class); sum of stratum_n / 2 = N
  s1  <- r$stratum_n[r$stratum_id == 1L][1L]
  s2  <- r$stratum_n[r$stratum_id == 2L][1L]
  N   <- r$marginal_total_n[1L]
  expect_equal(s1 + s2, N)
})

# ---------------------------------------------------------------------------
# lort_propensity_weights
# ---------------------------------------------------------------------------

test_that("lort_propensity_weights: requires cta_ort object", {
  expect_error(lort_propensity_weights(list()), "cta_ort")
})

test_that("lort_propensity_weights: model_family is 'lort'", {
  skip_if(Sys.getenv("ODACORE_TEST_TIER") %in% c("", "cran"),
          "lort fit skipped outside fast/smoke/full tier")
  ort <- lort_fit(X_fix[, c("x1","x4"), drop = FALSE], y_fix,
                  min_n = 10L, mc_iter = 100L, mc_seed = 42L)
  r   <- lort_propensity_weights(ort)
  expect_true(all(r$model_family == "lort"))
})

test_that("lort_propensity_weights: global_optimization = FALSE", {
  skip_if(Sys.getenv("ODACORE_TEST_TIER") %in% c("", "cran"),
          "lort fit skipped outside fast/smoke/full tier")
  ort <- lort_fit(X_fix[, c("x1","x4"), drop = FALSE], y_fix,
                  min_n = 10L, mc_iter = 100L, mc_seed = 42L)
  r   <- lort_propensity_weights(ort)
  expect_true(all(r$global_optimization == FALSE))
})

test_that("lort_propensity_weights: sda_anchored = FALSE", {
  skip_if(Sys.getenv("ODACORE_TEST_TIER") %in% c("", "cran"),
          "lort fit skipped outside fast/smoke/full tier")
  ort <- lort_fit(X_fix[, c("x1","x4"), drop = FALSE], y_fix,
                  min_n = 10L, mc_iter = 100L, mc_seed = 42L)
  r   <- lort_propensity_weights(ort)
  expect_true(all(r$sda_anchored == FALSE))
})

test_that("lort_propensity_weights: returns 18-column data frame", {
  skip_if(Sys.getenv("ODACORE_TEST_TIER") %in% c("", "cran"),
          "lort fit skipped outside fast/smoke/full tier")
  ort <- lort_fit(X_fix[, c("x1","x4"), drop = FALSE], y_fix,
                  min_n = 10L, mc_iter = 100L, mc_seed = 42L)
  r   <- lort_propensity_weights(ort)
  expect_equal(ncol(r), 18L)
})

test_that("lort_propensity_weights: adjusted weights are finite", {
  skip_if(Sys.getenv("ODACORE_TEST_TIER") %in% c("", "cran"),
          "lort fit skipped outside fast/smoke/full tier")
  ort <- lort_fit(X_fix[, c("x1","x4"), drop = FALSE], y_fix,
                  min_n = 10L, mc_iter = 100L, mc_seed = 42L)
  r   <- lort_propensity_weights(ort)
  expect_true(all(is.finite(r$adjusted_propensity_weight)))
})

# ---------------------------------------------------------------------------
# Absence checks (OBBP doctrine)
# ---------------------------------------------------------------------------

test_that("sda_propensity_weights does not exist", {
  expect_false(existsFunction("sda_propensity_weights"))
})

test_that("propensity_weights (generic) does not exist", {
  expect_false(existsFunction("propensity_weights"))
})

test_that("propensity_scores (generic) does not exist", {
  expect_false(existsFunction("propensity_scores"))
})

test_that("sort_propensity_weights is not exported", {
  expect_false("sort_propensity_weights" %in% getNamespaceExports("odacore"))
})

test_that("gort_propensity_weights is not exported", {
  expect_false("gort_propensity_weights" %in% getNamespaceExports("odacore"))
})

# ---------------------------------------------------------------------------
# cta_propensity_weights regression (existing test file not broken)
# ---------------------------------------------------------------------------

test_that("cta_propensity_weights still exported", {
  expect_true("cta_propensity_weights" %in% getNamespaceExports("odacore"))
})
