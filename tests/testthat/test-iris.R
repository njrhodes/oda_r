###############################################################################
# test-iris.R
# MegaODA.exe gold regression - iris dataset, all 4 attributes, K=3.
#
# Gold values from MegaODA.exe output. Species recoded: setosa=1,
# versicolor=2, virginica=3.  Boundary = right_closed, LOO = refit.
###############################################################################

iris_y <- function() {
  as.integer(factor(iris$Species,
                    levels = c("setosa","versicolor","virginica")))
}

run_iris_attr <- function(col_name) {
  x <- iris[[col_name]]
  y <- iris_y()
  oda_multiclass_unioda_core(
    x = x, y = y,
    attr_type         = "ordered",
    priors_on         = TRUE,
    K_segments        = 3L,
    degen             = FALSE,
    mcarlo            = FALSE,
    loo               = "on",
    boundary_mode     = "right_closed",
    loo_opts          = list(grid_mode = "refit", use_samplerep = FALSE, priors_mode = "fold")
  )
}

# ---- V1: Sepal.Length -------------------------------------------------------
test_that("iris V1 (Sepal.Length): rule matches MegaODA gold", {
  testthat::skip_on_cran()
  fit <- run_iris_attr("Sepal.Length")
  expect_true(fit$ok)
  expect_equal(round(fit$rule$cut_values, 2), c(5.45, 6.15),
               label = "V1: cuts")
  expect_equal(as.integer(fit$rule$seg_classes), c(1L, 2L, 3L),
               label = "V1: seg classes")
})

test_that("iris V1 (Sepal.Length): training confusion matches gold", {
  testthat::skip_on_cran()
  fit <- run_iris_attr("Sepal.Length")
  gold <- matrix(c(45,5,0, 6,28,16, 1,10,39), nrow=3, byrow=TRUE)
  yhat <- oda_rule_predict_multiclass(iris$Sepal.Length, fit$rule,
                                      boundary="right_closed")
  conf <- confusion_raw(iris_y(), yhat, C=3L)
  expect_equal(conf, gold, label = "V1: train confusion")
  expect_equal(round(pac_overall(conf), 2), 74.67, label = "V1: train PAC")
})

test_that("iris V1 (Sepal.Length): LOO confusion matches gold", {
  testthat::skip_on_cran()
  fit <- run_iris_attr("Sepal.Length")
  gold <- matrix(c(45,5,0, 6,28,16, 1,12,37), nrow=3, byrow=TRUE)
  expect_true(isTRUE(fit$loo$allowed), label = "V1: LOO allowed")
  conf_l <- unname(fit$loo$confusion_raw$confusion)
  storage.mode(conf_l) <- "integer"
  expect_equal(conf_l, gold, label = "V1: LOO confusion")
  expect_equal(round(pac_overall(conf_l), 2), 73.33, label = "V1: LOO PAC")
})

# ---- V2: Sepal.Width --------------------------------------------------------
test_that("iris V2 (Sepal.Width): rule matches MegaODA gold", {
  testthat::skip_on_cran()
  fit <- run_iris_attr("Sepal.Width")
  expect_true(fit$ok)
  expect_equal(round(fit$rule$cut_values, 2), c(2.95, 3.35),
               label = "V2: cuts")
  expect_equal(as.integer(fit$rule$seg_classes), c(2L, 3L, 1L),
               label = "V2: seg classes")
})

test_that("iris V2 (Sepal.Width): training PAC matches gold", {
  testthat::skip_on_cran()
  fit <- run_iris_attr("Sepal.Width")
  yhat <- oda_rule_predict_multiclass(iris$Sepal.Width, fit$rule,
                                      boundary="right_closed")
  conf <- confusion_raw(iris_y(), yhat, C=3L)
  expect_equal(round(pac_overall(conf), 2), 59.33, label = "V2: train PAC")
})

test_that("iris V2 (Sepal.Width): LOO PAC matches gold", {
  testthat::skip_on_cran()
  fit <- run_iris_attr("Sepal.Width")
  expect_true(isTRUE(fit$loo$allowed), label = "V2: LOO allowed")
  conf_l <- unname(fit$loo$confusion_raw$confusion)
  storage.mode(conf_l) <- "integer"
  expect_equal(round(pac_overall(conf_l), 2), 59.33, label = "V2: LOO PAC")
})

# ---- V3: Petal.Length -------------------------------------------------------
test_that("iris V3 (Petal.Length): rule matches MegaODA gold", {
  testthat::skip_on_cran()
  fit <- run_iris_attr("Petal.Length")
  expect_true(fit$ok)
  expect_equal(round(fit$rule$cut_values, 2), c(2.45, 4.75),
               label = "V3: cuts")
  expect_equal(as.integer(fit$rule$seg_classes), c(1L, 2L, 3L),
               label = "V3: seg classes")
})

test_that("iris V3 (Petal.Length): training confusion matches gold", {
  testthat::skip_on_cran()
  fit <- run_iris_attr("Petal.Length")
  gold <- matrix(c(50,0,0, 0,44,6, 0,1,49), nrow=3, byrow=TRUE)
  yhat <- oda_rule_predict_multiclass(iris$Petal.Length, fit$rule,
                                      boundary="right_closed")
  conf <- confusion_raw(iris_y(), yhat, C=3L)
  expect_equal(conf, gold, label = "V3: train confusion")
  expect_equal(round(pac_overall(conf), 2), 95.33, label = "V3: train PAC")
})

test_that("iris V3 (Petal.Length): LOO confusion matches gold", {
  testthat::skip_on_cran()
  fit <- run_iris_attr("Petal.Length")
  gold <- matrix(c(50,0,0, 0,44,6, 0,3,47), nrow=3, byrow=TRUE)
  conf_l <- unname(fit$loo$confusion_raw$confusion)
  storage.mode(conf_l) <- "integer"
  expect_equal(conf_l, gold, label = "V3: LOO confusion")
  expect_equal(round(pac_overall(conf_l), 2), 94.00, label = "V3: LOO PAC")
})

# ---- V4: Petal.Width --------------------------------------------------------
test_that("iris V4 (Petal.Width): rule matches MegaODA gold", {
  testthat::skip_on_cran()
  fit <- run_iris_attr("Petal.Width")
  expect_true(fit$ok)
  expect_equal(round(fit$rule$cut_values, 2), c(0.80, 1.65),
               label = "V4: cuts")
  expect_equal(as.integer(fit$rule$seg_classes), c(1L, 2L, 3L),
               label = "V4: seg classes")
})

test_that("iris V4 (Petal.Width): training confusion matches gold", {
  testthat::skip_on_cran()
  fit <- run_iris_attr("Petal.Width")
  gold <- matrix(c(50,0,0, 0,48,2, 0,4,46), nrow=3, byrow=TRUE)
  yhat <- oda_rule_predict_multiclass(iris$Petal.Width, fit$rule,
                                      boundary="right_closed")
  conf <- confusion_raw(iris_y(), yhat, C=3L)
  expect_equal(conf, gold, label = "V4: train confusion")
  expect_equal(round(pac_overall(conf), 2), 96.00, label = "V4: train PAC")
})

test_that("iris V4 (Petal.Width): LOO confusion matches gold", {
  testthat::skip_on_cran()
  fit <- run_iris_attr("Petal.Width")
  gold <- matrix(c(50,0,0, 0,48,2, 0,5,45), nrow=3, byrow=TRUE)
  conf_l <- unname(fit$loo$confusion_raw$confusion)
  storage.mode(conf_l) <- "integer"
  expect_equal(conf_l, gold, label = "V4: LOO confusion")
  expect_equal(round(pac_overall(conf_l), 2), 95.33, label = "V4: LOO PAC")
})
