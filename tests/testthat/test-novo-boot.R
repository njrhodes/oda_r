# tests/testthat/test-novo-boot.R
# Tests for novo_boot_ci() — novometric bootstrap CI from a fixed 2x2 confusion.

# ---- helpers -----------------------------------------------------------------

good_conf <- function() matrix(c(146, 40, 36, 33), nrow = 2L, byrow = TRUE)

# ---- input validation --------------------------------------------------------

test_that("non-matrix input is rejected", {
  expect_error(novo_boot_ci(c(10, 5, 3, 8)), "2x2 numeric matrix")
})

test_that("non-2x2 matrix is rejected", {
  expect_error(novo_boot_ci(matrix(1:9, 3, 3)), "2x2 numeric matrix")
})

test_that("non-numeric matrix is rejected", {
  m <- matrix(c("a","b","c","d"), 2, 2)
  expect_error(novo_boot_ci(m), "2x2 numeric matrix")
})

test_that("NA in confusion is rejected", {
  m <- matrix(c(10L, NA_integer_, 5L, 8L), 2, 2)
  expect_error(novo_boot_ci(m), "finite")
})

test_that("Inf in confusion is rejected", {
  m <- matrix(c(10, Inf, 5, 8), 2, 2)
  expect_error(novo_boot_ci(m), "finite")
})

test_that("negative counts are rejected", {
  m <- matrix(c(10, -1, 5, 8), 2, 2)
  expect_error(novo_boot_ci(m), "non-negative")
})

test_that("zero-total confusion is rejected", {
  m <- matrix(0L, 2, 2)
  expect_error(novo_boot_ci(m), "n > 0")
})

test_that("nboot < 1 is rejected", {
  expect_error(novo_boot_ci(good_conf(), nboot = 0L), "positive integer")
})

test_that("nboot = NA is rejected", {
  expect_error(novo_boot_ci(good_conf(), nboot = NA_integer_), "positive integer")
})

test_that("sample_frac too small is rejected", {
  expect_error(novo_boot_ci(good_conf(), sample_frac = 1e-6), "k =")
})

test_that("non-integer counts warn", {
  m <- matrix(c(10.7, 5.1, 3.2, 8.9), 2, 2)
  expect_warning(novo_boot_ci(m, nboot = 50L, seed = 1L), "non-integer")
})

# ---- output class and required fields ----------------------------------------

test_that("output class is c('novo_boot_ci', 'list')", {
  ci <- novo_boot_ci(good_conf(), nboot = 100L, seed = 1L)
  expect_s3_class(ci, "novo_boot_ci")
  expect_true(inherits(ci, "list"))
})

test_that("output has all required fields", {
  ci <- novo_boot_ci(good_conf(), nboot = 100L, seed = 1L)
  expected_names <- c("call", "confusion", "n", "k", "nboot", "sample_frac",
                      "probs", "alternative", "has_zero_cells",
                      "observed", "model", "chance", "quantiles", "ci",
                      "significant")
  expect_true(all(expected_names %in% names(ci)))
})

test_that("n and k are correct", {
  conf <- good_conf()
  ci <- novo_boot_ci(conf, nboot = 100L, seed = 1L, sample_frac = 0.5)
  expect_equal(ci$n, sum(conf))
  expect_equal(ci$k, as.integer(round(0.5 * sum(conf))))
})

test_that("confusion field is an integer matrix matching input", {
  conf <- good_conf()
  ci <- novo_boot_ci(conf, nboot = 50L, seed = 1L)
  expect_true(is.matrix(ci$confusion))
  expect_equal(storage.mode(ci$confusion), "integer")
  expect_equal(as.integer(ci$confusion), as.integer(conf))
  expect_identical(dim(ci$confusion), c(2L, 2L))
})

test_that("has_zero_cells is FALSE for good_conf", {
  ci <- novo_boot_ci(good_conf(), nboot = 50L, seed = 1L)
  expect_false(ci$has_zero_cells)
})

test_that("has_zero_cells is TRUE when a cell is zero", {
  m <- matrix(c(10L, 0L, 5L, 8L), 2, 2)
  ci <- novo_boot_ci(m, nboot = 50L, seed = 1L)
  expect_true(ci$has_zero_cells)
})

test_that("significant is a logical scalar", {
  ci <- novo_boot_ci(good_conf(), nboot = 100L, seed = 1L)
  expect_type(ci$significant, "logical")
  expect_length(ci$significant, 1L)
})

# ---- observed data frame -----------------------------------------------------

test_that("observed has metric and value columns", {
  ci <- novo_boot_ci(good_conf(), nboot = 50L, seed = 1L)
  expect_true(is.data.frame(ci$observed))
  expect_true(all(c("metric", "value") %in% names(ci$observed)))
})

test_that("observed has a row for each expected metric", {
  ci <- novo_boot_ci(good_conf(), nboot = 50L, seed = 1L)
  expected_metrics <- c("sensitivity", "specificity", "mean_pac", "ess",
                        "odds_ratio", "risk_ratio")
  expect_true(all(expected_metrics %in% ci$observed$metric))
})

test_that("observed ESS is consistent with confusion arithmetic", {
  conf <- good_conf()
  ci   <- novo_boot_ci(conf, nboot = 50L, seed = 1L)
  TN <- conf[1,1]; FP <- conf[1,2]; FN <- conf[2,1]; TP <- conf[2,2]
  sens_obs <- TP / (TP + FN) * 100
  spec_obs <- TN / (TN + FP) * 100
  mpac_obs <- (sens_obs + spec_obs) / 2
  ess_obs  <- 100 * (mpac_obs / 100 - 0.5) / 0.5
  expect_equal(ci$observed$value[ci$observed$metric == "ess"], ess_obs,
               tolerance = 1e-10)
})

# ---- model and chance data frames --------------------------------------------

test_that("model and chance have nboot rows", {
  ci <- novo_boot_ci(good_conf(), nboot = 200L, seed = 1L)
  expect_equal(nrow(ci$model),  200L)
  expect_equal(nrow(ci$chance), 200L)
})

test_that("model and chance have required columns", {
  ci <- novo_boot_ci(good_conf(), nboot = 50L, seed = 1L)
  expected_cols <- c("sensitivity", "specificity", "mean_pac", "ess",
                     "odds_ratio", "risk_ratio", "p_value")
  expect_true(all(expected_cols %in% names(ci$model)))
  expect_true(all(expected_cols %in% names(ci$chance)))
})

test_that("p_value column is present and numeric in model and chance", {
  ci <- novo_boot_ci(good_conf(), nboot = 100L, seed = 1L)
  expect_type(ci$model$p_value,  "double")
  expect_type(ci$chance$p_value, "double")
})

test_that("sensitivity and specificity are in [0, 100]", {
  ci <- novo_boot_ci(good_conf(), nboot = 200L, seed = 1L)
  expect_true(all(ci$model$sensitivity >= 0 & ci$model$sensitivity <= 100,
                  na.rm = TRUE))
  expect_true(all(ci$model$specificity >= 0 & ci$model$specificity <= 100,
                  na.rm = TRUE))
})

test_that("ESS is consistent with mean_pac in model distribution", {
  ci <- novo_boot_ci(good_conf(), nboot = 200L, seed = 1L)
  computed_ess <- 100 * (ci$model$mean_pac / 100 - 0.5) / 0.5
  expect_equal(ci$model$ess, computed_ess, tolerance = 1e-10)
})

test_that("no Inf in odds_ratio or risk_ratio (NA instead)", {
  m <- matrix(c(50L, 0L, 0L, 50L), 2, 2)
  ci <- novo_boot_ci(m, nboot = 50L, seed = 1L)
  expect_false(any(is.infinite(ci$model$odds_ratio), na.rm = FALSE))
  expect_false(any(is.infinite(ci$model$risk_ratio), na.rm = FALSE))
})

# ---- quantiles data frame ----------------------------------------------------

test_that("quantiles has prob column and correct nrow", {
  probs <- c(0, .025, .5, .975, 1)
  ci <- novo_boot_ci(good_conf(), nboot = 100L, seed = 1L, probs = probs)
  expect_equal(nrow(ci$quantiles), length(probs))
  expect_true("prob" %in% names(ci$quantiles))
  expect_equal(ci$quantiles$prob, probs)
})

test_that("quantiles contains p_value columns for model and chance", {
  ci <- novo_boot_ci(good_conf(), nboot = 100L, seed = 1L)
  expect_true("p_value_model"  %in% names(ci$quantiles))
  expect_true("p_value_chance" %in% names(ci$quantiles))
})

# ---- ci data frame -----------------------------------------------------------

test_that("ci has exact required column names", {
  ci <- novo_boot_ci(good_conf(), nboot = 100L, seed = 1L)
  expect_identical(names(ci$ci),
                   c("metric", "model_lower", "model_upper",
                     "chance_lower", "chance_upper", "overlap"))
})

test_that("ci has one row per metric including ess and p_value", {
  ci <- novo_boot_ci(good_conf(), nboot = 100L, seed = 1L)
  expect_true("ess"     %in% ci$ci$metric)
  expect_true("p_value" %in% ci$ci$metric)
})

test_that("overlap is logical", {
  ci <- novo_boot_ci(good_conf(), nboot = 100L, seed = 1L)
  expect_type(ci$ci$overlap, "logical")
})

test_that("model_lower <= model_upper and chance_lower <= chance_upper", {
  ci <- novo_boot_ci(good_conf(), nboot = 100L, seed = 1L)
  ok_m <- ci$ci$model_lower  <= ci$ci$model_upper  | is.na(ci$ci$model_lower)
  ok_c <- ci$ci$chance_lower <= ci$ci$chance_upper | is.na(ci$ci$chance_lower)
  expect_true(all(ok_m))
  expect_true(all(ok_c))
})

# ---- significance logic ------------------------------------------------------

test_that("significant = TRUE for perfect-separation confusion", {
  m <- matrix(c(50L, 0L, 0L, 50L), 2, 2)
  ci <- novo_boot_ci(m, nboot = 500L, seed = 42L)
  expect_true(ci$significant)
})

test_that("significant = FALSE for chance-level confusion", {
  m <- matrix(c(25L, 25L, 25L, 25L), 2, 2)
  ci <- novo_boot_ci(m, nboot = 500L, seed = 42L)
  expect_false(ci$significant)
})

test_that("significant is consistent with ci ESS rows", {
  ci <- novo_boot_ci(good_conf(), nboot = 200L, seed = 1L)
  ess_row <- ci$ci[ci$ci$metric == "ess", ]
  expected_sig <- isTRUE(ess_row$model_lower > ess_row$chance_upper)
  expect_identical(ci$significant, expected_sig)
})

# ---- seed reproducibility ----------------------------------------------------

test_that("same seed produces identical results", {
  ci1 <- novo_boot_ci(good_conf(), nboot = 200L, seed = 99L)
  ci2 <- novo_boot_ci(good_conf(), nboot = 200L, seed = 99L)
  expect_identical(ci1$model,  ci2$model)
  expect_identical(ci1$chance, ci2$chance)
  expect_identical(ci1$significant, ci2$significant)
})

test_that("different seeds produce different model distributions", {
  ci1 <- novo_boot_ci(good_conf(), nboot = 200L, seed = 1L)
  ci2 <- novo_boot_ci(good_conf(), nboot = 200L, seed = 2L)
  expect_false(identical(ci1$model$ess, ci2$model$ess))
})

# ---- zero-cell no crash ------------------------------------------------------

test_that("zero-cell confusion runs without error", {
  m <- matrix(c(80L, 0L, 20L, 0L), 2, 2)
  expect_no_error(novo_boot_ci(m, nboot = 100L, seed = 1L))
})

# ---- print method ------------------------------------------------------------

test_that("print returns invisibly and produces output", {
  ci <- novo_boot_ci(good_conf(), nboot = 100L, seed = 1L)
  expect_output(ret <- print(ci), "Novometric")
  expect_identical(ret, ci)
})

test_that("print output contains observed ESS and Mean PAC lines", {
  ci <- novo_boot_ci(good_conf(), nboot = 100L, seed = 1L)
  expect_output(print(ci), "Observed")
  expect_output(print(ci), "ESS")
})

# ---- myeloma anchor (slow) ---------------------------------------------------

test_that("myeloma MINDENOM=1 anchor: structure, observed, and significance", {
  conf <- matrix(c(146, 40,
                    36, 33), nrow = 2L, byrow = TRUE)
  ci <- novo_boot_ci(conf, nboot = 5000L, seed = 2024L)

  expect_s3_class(ci, "novo_boot_ci")
  expect_equal(ci$n, 255L)
  expect_false(ci$has_zero_cells)

  # observed fields present and non-NA
  expect_true(all(!is.na(ci$observed$value[ci$observed$metric %in%
                         c("ess", "mean_pac", "sensitivity", "specificity")])))

  # ci column names exact
  expect_identical(names(ci$ci),
                   c("metric", "model_lower", "model_upper",
                     "chance_lower", "chance_upper", "overlap"))

  # ESS model 95% CI lower should be positive (real signal)
  ess_row <- ci$ci[ci$ci$metric == "ess", ]
  expect_true(ess_row$model_lower > 0,
              info = "model ESS 2.5th pctile should be positive for myeloma")
})
