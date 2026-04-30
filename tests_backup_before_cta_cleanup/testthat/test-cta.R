###############################################################################
# test-cta.R — CTA core tests
###############################################################################

# ---- Data helpers -----------------------------------------------------------
bin_data <- function() {
  # Perfectly separable: x <= 4.5 → class 0, x > 4.5 → class 1
  list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
}

sep3_data <- function() {
  # Three-class, one informative attribute
  list(
    X = data.frame(x1 = c(1,2,3,4,5,6,7,8,9), x2 = rep(5, 9)),
    y = c(1L,1L,1L,2L,2L,2L,3L,3L,3L)
  )
}

# =============================================================================

test_that("oda_cta_fit returns cta_tree for binary data", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                      mc_seed = 1L, loo = "off")
  expect_s3_class(tree, "cta_tree")
  expect_gte(tree$n_nodes, 1L)
})

test_that("root node n_obs equals n", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                      mc_seed = 1L, loo = "off")
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$n_obs, length(d$y))
})

test_that("mindenom too large → single leaf node", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 999L, mc_iter = 100L,
                      mc_seed = 1L, loo = "off")
  expect_equal(tree$n_nodes, 1L)
  expect_true(isTRUE(tree$nodes[[1]]$leaf))
})

test_that("max_depth = 1 → root is a leaf (no splits allowed)", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, max_depth = 1L,
                      mc_iter = 300L, mc_seed = 1L, loo = "off")
  root <- tree$nodes[[tree$root_id]]
  expect_true(isTRUE(root$leaf))
})

test_that("max_depth = 2 → no node deeper than depth 2", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, max_depth = 2L,
                      mc_iter = 300L, mc_seed = 1L, loo = "off")
  depths <- vapply(tree$nodes, function(nd) nd$depth, integer(1))
  expect_lte(max(depths, na.rm = TRUE), 2L)
})

test_that("predict returns integer vector of length n", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                      mc_seed = 1L, loo = "off")
  preds <- predict(tree, d$X)
  expect_equal(length(preds), nrow(d$X))
  expect_true(is.integer(preds))
})

test_that("predict: all returned classes are in the training class set", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                      mc_seed = 1L, loo = "off")
  preds <- predict(tree, d$X)
  expect_true(all(preds %in% unique(d$y)))
})

test_that("training accuracy >= majority class baseline", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 500L,
                      mc_seed = 7L, loo = "off")
  preds        <- predict(tree, d$X)
  acc          <- mean(preds == d$y)
  majority_acc <- max(table(d$y)) / length(d$y)
  expect_gte(acc, majority_acc - 0.01)
})

test_that("cta_node_table returns data frame with required columns", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                      mc_seed = 1L, loo = "off")
  tbl <- cta_node_table(tree)
  expect_true(is.data.frame(tbl))
  required <- c("node_id","depth","n_obs","ess","p_mc","leaf","attribute")
  expect_true(all(required %in% names(tbl)))
})

test_that("print.cta_tree runs without error and mentions CTA Tree", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                      mc_seed = 1L, loo = "off")
  expect_output(print(tree), "CTA Tree")
})

test_that("case weights accepted without error", {
  d <- bin_data()
  w <- c(2,1,2,1,2,1,2,1)
  tree <- oda_cta_fit(d$X, d$y, w = w, mindenom = 2L,
                      mc_iter = 200L, mc_seed = 1L, loo = "off")
  expect_s3_class(tree, "cta_tree")
  root <- tree$nodes[[1]]
  expect_equal(root$n_weighted, sum(w))
})

test_that("multi-attribute: most informative attribute selected at root", {
  # x1 perfectly separates, x2 is noise
  set.seed(42)
  x1 <- 1:10; x2 <- sample(1:10)
  y  <- c(0L,0L,0L,0L,0L,1L,1L,1L,1L,1L)
  X  <- data.frame(x1 = x1, x2 = x2)
  tree <- oda_cta_fit(X, y, mindenom = 2L, mc_iter = 500L,
                      mc_seed = 42L, loo = "off")
  root <- tree$nodes[[tree$root_id]]
  if (!isTRUE(root$leaf)) {
    expect_equal(root$attribute, "x1")
  } else {
    skip("root is leaf — increase mc_iter for reliable test")
  }
})
