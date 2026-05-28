###############################################################################
# test-lort-fit.R
#
# Tests for lort_fit() — the preferred explicit entry point for the LORT
# workflow layer.
#
# Verifies:
#   T1  lort_fit() returns cta_ort / cta_tree classes
#   T2  method metadata is "lort"; global_optimization and sda_anchored FALSE
#   T3  lort_fit() output equals cta_fit(..., recursive=TRUE) on same args
#   T4  cta_fit(..., recursive=TRUE) still works after lort_fit added
#   T5  passing mindenom to lort_fit() errors
#   T6  passing recursive to lort_fit() errors
###############################################################################

# Minimal deterministic synthetic fixture shared across tests.
lf_X <- data.frame(
  A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
  B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
)
lf_y <- c(rep(0L, 20), rep(0L, 20), rep(1L, 20))

lf_args <- list(X = lf_X, y = lf_y,
                mc_iter = 100L, mc_seed = 42L,
                loo = "off", min_n = 5L)

# Module-level fixture computed once.
lf_fit <- do.call(lort_fit, lf_args)

# ---------------------------------------------------------------------------
# T1: class contract
# ---------------------------------------------------------------------------
test_that("lort_fit: returns cta_ort and cta_tree (T1)", {
  expect_s3_class(lf_fit, "cta_ort")
  expect_s3_class(lf_fit, "cta_tree")
  expect_true(isTRUE(lf_fit$recursive))
})

# ---------------------------------------------------------------------------
# T2: LORT method metadata
# ---------------------------------------------------------------------------
test_that("lort_fit: method metadata is 'lort' (T2)", {
  expect_equal(lf_fit$ort_settings$method, "lort")
  expect_equal(lf_fit$ort_settings$method_label, "Locally Optimal Recursive Tree")
  expect_false(isTRUE(lf_fit$ort_settings$global_optimization))
  expect_false(isTRUE(lf_fit$ort_settings$sda_anchored))
})

# ---------------------------------------------------------------------------
# T3: equivalence to cta_fit(..., recursive = TRUE)
# ---------------------------------------------------------------------------
test_that("lort_fit: output equals cta_fit(recursive=TRUE) on same args (T3)", {
  cta_equivalent <- cta_fit(lf_X, lf_y, recursive = TRUE,
                             mc_iter = 100L, mc_seed = 42L,
                             loo = "off", min_n = 5L)

  expect_equal(lf_fit$strata$n,            cta_equivalent$strata$n)
  expect_equal(lf_fit$strata$path,          cta_equivalent$strata$path)
  expect_equal(lf_fit$n_strata,             cta_equivalent$n_strata)
  expect_equal(lf_fit$ort_settings$method,  cta_equivalent$ort_settings$method)

  pred_lort <- predict(lf_fit,         lf_X, type = "class")
  pred_cta  <- predict(cta_equivalent, lf_X, type = "class")
  expect_equal(pred_lort, pred_cta)
})

# ---------------------------------------------------------------------------
# T4: cta_fit(recursive=TRUE) still works (legacy-compat)
# ---------------------------------------------------------------------------
test_that("cta_fit(recursive=TRUE) still works after lort_fit added (T4)", {
  fit <- cta_fit(lf_X, lf_y, recursive = TRUE,
                 mc_iter = 100L, mc_seed = 42L,
                 loo = "off", min_n = 5L)
  expect_s3_class(fit, "cta_ort")
  expect_true(isTRUE(fit$recursive))
  expect_equal(fit$ort_settings$method, "lort")
})

# ---------------------------------------------------------------------------
# T5: mindenom is not a lort_fit parameter — errors
# ---------------------------------------------------------------------------
test_that("lort_fit: passing mindenom errors (T5)", {
  expect_error(
    lort_fit(lf_X, lf_y, mindenom = 1L,
             mc_iter = 100L, mc_seed = 42L, loo = "off", min_n = 5L)
  )
})

# ---------------------------------------------------------------------------
# T6: recursive is not a lort_fit parameter — errors
# ---------------------------------------------------------------------------
test_that("lort_fit: passing recursive errors (T6)", {
  expect_error(
    lort_fit(lf_X, lf_y, recursive = FALSE,
             mc_iter = 100L, mc_seed = 42L, loo = "off", min_n = 5L)
  )
})
