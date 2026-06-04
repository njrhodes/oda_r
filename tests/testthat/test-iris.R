###############################################################################
# test-iris.R
# MegaODA gold regression: iris dataset, all 4 attributes, K=3.
#
# Canon source:
#   tests/testthat/fixtures/iris/IRIS_CANON_ALL4.OUT
#   tests/testthat/fixtures/iris/IRIS_EXPECTED_GOLD.csv
#
# MegaODA run: IRIS_CANON_ALL4.PGM, MCARLO ITER 25000, LOO.
# Status: 0, elapsed: ~160 s.  Attributes V2-V5 (MegaODA naming);
# V1 = class (Species recoded: setosa=1, versicolor=2, virginica=3).
#
# Performance: fits are cached at file scope -- each attribute is fitted
# exactly once (4 fits total, one per attribute).  Each fit runs LOO
# refit (loo="on", grid_mode="refit") on n=150, K=3.  mcarlo=FALSE at
# fast tier: LOO coverage only; MC results are in the canon fixture.
###############################################################################

iris_y <- function() {
  as.integer(factor(iris$Species,
                    levels = c("setosa", "versicolor", "virginica")))
}

run_iris_attr <- function(col_name) {
  x <- iris[[col_name]]
  y <- iris_y()
  oda_multiclass_unioda_core(
    x             = x, y = y,
    attr_type     = "ordered",
    priors_on     = TRUE,
    K_segments    = 3L,
    degen         = FALSE,
    mcarlo        = FALSE,
    loo           = "on",
    boundary_mode = "right_closed",
    loo_opts      = list(grid_mode    = "refit",
                         use_samplerep = FALSE,
                         priors_mode  = "fold")
  )
}

# ---------------------------------------------------------------------------
# File-scope cache: each attribute is fitted exactly once.
# V2=Sepal.Length, V3=Sepal.Width, V4=Petal.Length, V5=Petal.Width
# (MegaODA variable naming; V1 = class).
# ---------------------------------------------------------------------------
iris_fits <- local({
  list(
    V2 = run_iris_attr("Sepal.Length"),
    V3 = run_iris_attr("Sepal.Width"),
    V4 = run_iris_attr("Petal.Length"),
    V5 = run_iris_attr("Petal.Width")
  )
})

# ---------------------------------------------------------------------------
# Load canonical expected values from MegaODA fixture.
# Columns: oda_attribute, source_name, expected_cut_1, expected_cut_2,
#          expected_segment_classes, expected_train_overall_PAC,
#          expected_LOO_overall_PAC
# ---------------------------------------------------------------------------
iris_gold <- local({
  p <- testthat::test_path("fixtures", "iris", "IRIS_EXPECTED_GOLD.csv")
  g <- read.csv(p, stringsAsFactors = FALSE)
  g$seg_classes_int <- lapply(
    strsplit(g$expected_segment_classes, ","),
    function(x) as.integer(trimws(x))
  )
  g
})

gold_row <- function(attr_id) {
  iris_gold[iris_gold$oda_attribute == attr_id, ]
}

# ---- V2: Sepal.Length -------------------------------------------------------

test_that("iris V2 (Sepal.Length): rule matches MegaODA gold", {
  testthat::skip_on_cran()
  fit <- iris_fits$V2
  g   <- gold_row("V2")
  expect_true(fit$ok)
  expect_equal(round(fit$rule$cut_values, 2),
               c(g$expected_cut_1, g$expected_cut_2),
               label = "V2: cuts")
  expect_equal(as.integer(fit$rule$seg_classes), g$seg_classes_int[[1L]],
               label = "V2: seg classes")
})

test_that("iris V2 (Sepal.Length): training confusion matches gold", {
  testthat::skip_on_cran()
  fit  <- iris_fits$V2
  gold <- matrix(c(45, 5, 0,  6, 28, 16,  1, 10, 39), nrow = 3L, byrow = TRUE)
  yhat <- oda_rule_predict_multiclass(iris$Sepal.Length, fit$rule,
                                      boundary = "right_closed")
  conf <- confusion_raw(iris_y(), yhat, C = 3L)
  expect_equal(conf, gold, label = "V2: train confusion")
  g <- gold_row("V2")
  expect_equal(round(pac_overall(conf), 2), g$expected_train_overall_PAC,
               label = "V2: train PAC")
})

test_that("iris V2 (Sepal.Length): LOO confusion matches gold", {
  testthat::skip_on_cran()
  fit  <- iris_fits$V2
  gold <- matrix(c(45, 5, 0,  6, 28, 16,  1, 12, 37), nrow = 3L, byrow = TRUE)
  g    <- gold_row("V2")
  expect_true(isTRUE(fit$loo$allowed), label = "V2: LOO allowed")
  conf_l <- unname(fit$loo$confusion_raw$confusion)
  storage.mode(conf_l) <- "integer"
  expect_equal(conf_l, gold, label = "V2: LOO confusion")
  expect_equal(round(pac_overall(conf_l), 2), g$expected_LOO_overall_PAC,
               label = "V2: LOO PAC")
})

# ---- V3: Sepal.Width --------------------------------------------------------

test_that("iris V3 (Sepal.Width): rule matches MegaODA gold", {
  testthat::skip_on_cran()
  fit <- iris_fits$V3
  g   <- gold_row("V3")
  expect_true(fit$ok)
  expect_equal(round(fit$rule$cut_values, 2),
               c(g$expected_cut_1, g$expected_cut_2),
               label = "V3: cuts")
  expect_equal(as.integer(fit$rule$seg_classes), g$seg_classes_int[[1L]],
               label = "V3: seg classes")
})

test_that("iris V3 (Sepal.Width): training PAC matches gold", {
  testthat::skip_on_cran()
  fit  <- iris_fits$V3
  yhat <- oda_rule_predict_multiclass(iris$Sepal.Width, fit$rule,
                                      boundary = "right_closed")
  conf <- confusion_raw(iris_y(), yhat, C = 3L)
  g    <- gold_row("V3")
  expect_equal(round(pac_overall(conf), 2), g$expected_train_overall_PAC,
               label = "V3: train PAC")
})

test_that("iris V3 (Sepal.Width): LOO PAC matches gold", {
  testthat::skip_on_cran()
  fit  <- iris_fits$V3
  g    <- gold_row("V3")
  expect_true(isTRUE(fit$loo$allowed), label = "V3: LOO allowed")
  conf_l <- unname(fit$loo$confusion_raw$confusion)
  storage.mode(conf_l) <- "integer"
  expect_equal(round(pac_overall(conf_l), 2), g$expected_LOO_overall_PAC,
               label = "V3: LOO PAC")
})

# ---- V4: Petal.Length -------------------------------------------------------

test_that("iris V4 (Petal.Length): rule matches MegaODA gold", {
  testthat::skip_on_cran()
  fit <- iris_fits$V4
  g   <- gold_row("V4")
  expect_true(fit$ok)
  expect_equal(round(fit$rule$cut_values, 2),
               c(g$expected_cut_1, g$expected_cut_2),
               label = "V4: cuts")
  expect_equal(as.integer(fit$rule$seg_classes), g$seg_classes_int[[1L]],
               label = "V4: seg classes")
})

test_that("iris V4 (Petal.Length): training confusion matches gold", {
  testthat::skip_on_cran()
  fit  <- iris_fits$V4
  gold <- matrix(c(50, 0, 0,  0, 44, 6,  0, 1, 49), nrow = 3L, byrow = TRUE)
  yhat <- oda_rule_predict_multiclass(iris$Petal.Length, fit$rule,
                                      boundary = "right_closed")
  conf <- confusion_raw(iris_y(), yhat, C = 3L)
  expect_equal(conf, gold, label = "V4: train confusion")
  g <- gold_row("V4")
  expect_equal(round(pac_overall(conf), 2), g$expected_train_overall_PAC,
               label = "V4: train PAC")
})

test_that("iris V4 (Petal.Length): LOO confusion matches gold", {
  testthat::skip_on_cran()
  fit  <- iris_fits$V4
  gold <- matrix(c(50, 0, 0,  0, 44, 6,  0, 3, 47), nrow = 3L, byrow = TRUE)
  g    <- gold_row("V4")
  conf_l <- unname(fit$loo$confusion_raw$confusion)
  storage.mode(conf_l) <- "integer"
  expect_equal(conf_l, gold, label = "V4: LOO confusion")
  expect_equal(round(pac_overall(conf_l), 2), g$expected_LOO_overall_PAC,
               label = "V4: LOO PAC")
})

# ---- V5: Petal.Width --------------------------------------------------------

test_that("iris V5 (Petal.Width): rule matches MegaODA gold", {
  testthat::skip_on_cran()
  fit <- iris_fits$V5
  g   <- gold_row("V5")
  expect_true(fit$ok)
  expect_equal(round(fit$rule$cut_values, 2),
               c(g$expected_cut_1, g$expected_cut_2),
               label = "V5: cuts")
  expect_equal(as.integer(fit$rule$seg_classes), g$seg_classes_int[[1L]],
               label = "V5: seg classes")
})

test_that("iris V5 (Petal.Width): training confusion matches gold", {
  testthat::skip_on_cran()
  fit  <- iris_fits$V5
  gold <- matrix(c(50, 0, 0,  0, 48, 2,  0, 4, 46), nrow = 3L, byrow = TRUE)
  yhat <- oda_rule_predict_multiclass(iris$Petal.Width, fit$rule,
                                      boundary = "right_closed")
  conf <- confusion_raw(iris_y(), yhat, C = 3L)
  expect_equal(conf, gold, label = "V5: train confusion")
  g <- gold_row("V5")
  expect_equal(round(pac_overall(conf), 2), g$expected_train_overall_PAC,
               label = "V5: train PAC")
})

test_that("iris V5 (Petal.Width): LOO confusion matches gold", {
  testthat::skip_on_cran()
  fit  <- iris_fits$V5
  gold <- matrix(c(50, 0, 0,  0, 48, 2,  0, 5, 45), nrow = 3L, byrow = TRUE)
  g    <- gold_row("V5")
  conf_l <- unname(fit$loo$confusion_raw$confusion)
  storage.mode(conf_l) <- "integer"
  expect_equal(conf_l, gold, label = "V5: LOO confusion")
  expect_equal(round(pac_overall(conf_l), 2), g$expected_LOO_overall_PAC,
               label = "V5: LOO PAC")
})
