###############################################################################
# test-cta-endpoint-counts.R — per-endpoint class-count storage
#
# Verifies that .leaf_nd() stores class_counts_raw and class_counts_weighted
# on every terminal node at fit time.  These fields are the foundation for
# future staging-table and event-rate reporting.
#
# No public accessor yet — tests reach node fields directly via $nodes.
###############################################################################

# ---- Helpers -----------------------------------------------------------------

.counts_load_data <- function() {
  list(X = data.frame(x = 1:8), y = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L))
}

.counts_no_tree_fit <- function() {
  d <- .counts_load_data()
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 999L, mc_iter = 100L,
                mc_seed = 1L, loo = "off")
  )
}

.counts_stump_fit <- function() {
  d <- .counts_load_data()
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off")
  )
}

# Weighted stump: four obs class-0 weight 2, four obs class-1 weight 3.
.counts_weighted_fit <- function() {
  d <- .counts_load_data()
  w <- c(2, 2, 2, 2, 3, 3, 3, 3)
  suppressMessages(
    oda_cta_fit(d$X, d$y, w = w, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off")
  )
}

.leaf_nodes <- function(tree) {
  Filter(function(nd) !is.null(nd) && isTRUE(nd$leaf), tree$nodes)
}

# =============================================================================
# No-tree fit: single root leaf
# =============================================================================

test_that("counts: no-tree leaf has class_counts_raw", {
  tree   <- .counts_no_tree_fit()
  leaves <- .leaf_nodes(tree)
  expect_true(length(leaves) >= 1L)
  expect_true(!is.null(leaves[[1L]]$class_counts_raw))
})

test_that("counts: no-tree leaf has class_counts_weighted", {
  tree   <- .counts_no_tree_fit()
  leaves <- .leaf_nodes(tree)
  expect_true(!is.null(leaves[[1L]]$class_counts_weighted))
})

test_that("counts: no-tree class_counts_raw names are class labels", {
  tree   <- .counts_no_tree_fit()
  leaves <- .leaf_nodes(tree)
  expect_equal(names(leaves[[1L]]$class_counts_raw), c("0", "1"))
})

test_that("counts: no-tree class_counts_weighted names are class labels", {
  tree   <- .counts_no_tree_fit()
  leaves <- .leaf_nodes(tree)
  expect_equal(names(leaves[[1L]]$class_counts_weighted), c("0", "1"))
})

test_that("counts: no-tree class_counts_raw sums to n_obs", {
  tree   <- .counts_no_tree_fit()
  leaves <- .leaf_nodes(tree)
  nd     <- leaves[[1L]]
  expect_equal(sum(nd$class_counts_raw), nd$n_obs)
})

test_that("counts: no-tree class_counts_weighted sums to n_weighted", {
  tree   <- .counts_no_tree_fit()
  leaves <- .leaf_nodes(tree)
  nd     <- leaves[[1L]]
  expect_equal(sum(nd$class_counts_weighted), nd$n_weighted)
})

test_that("counts: no-tree class_counts_raw is integer", {
  tree   <- .counts_no_tree_fit()
  leaves <- .leaf_nodes(tree)
  expect_type(leaves[[1L]]$class_counts_raw, "integer")
})

test_that("counts: no-tree class_counts_weighted is double", {
  tree   <- .counts_no_tree_fit()
  leaves <- .leaf_nodes(tree)
  expect_type(leaves[[1L]]$class_counts_weighted, "double")
})

test_that("counts: no-tree unweighted fit has class_counts_raw == class_counts_weighted", {
  tree   <- .counts_no_tree_fit()
  leaves <- .leaf_nodes(tree)
  nd     <- leaves[[1L]]
  expect_equal(unname(as.numeric(nd$class_counts_raw)), unname(nd$class_counts_weighted))
})

# =============================================================================
# Stump: 2 leaves
# =============================================================================

test_that("counts: stump all leaves have class_counts_raw", {
  tree <- .counts_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  leaves <- .leaf_nodes(tree)
  has_field <- vapply(leaves, function(nd) !is.null(nd$class_counts_raw), logical(1L))
  expect_true(all(has_field))
})

test_that("counts: stump all leaves have class_counts_weighted", {
  tree <- .counts_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  leaves <- .leaf_nodes(tree)
  has_field <- vapply(leaves, function(nd) !is.null(nd$class_counts_weighted), logical(1L))
  expect_true(all(has_field))
})

test_that("counts: stump class_counts_raw names are c('0','1') on every leaf", {
  tree <- .counts_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  leaves <- .leaf_nodes(tree)
  ok <- vapply(leaves, function(nd) identical(names(nd$class_counts_raw), c("0","1")), logical(1L))
  expect_true(all(ok))
})

test_that("counts: stump per-leaf class_counts_raw sums to n_obs", {
  tree <- .counts_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  leaves <- .leaf_nodes(tree)
  ok <- vapply(leaves, function(nd) sum(nd$class_counts_raw) == nd$n_obs, logical(1L))
  expect_true(all(ok))
})

test_that("counts: stump per-leaf class_counts_weighted sums to n_weighted", {
  tree <- .counts_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  leaves <- .leaf_nodes(tree)
  ok <- vapply(leaves,
               function(nd) isTRUE(all.equal(sum(nd$class_counts_weighted), nd$n_weighted)),
               logical(1L))
  expect_true(all(ok))
})

test_that("counts: stump total class_counts_raw sums to 8 (n total)", {
  tree <- .counts_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  leaves <- .leaf_nodes(tree)
  total  <- sum(vapply(leaves, function(nd) sum(nd$class_counts_raw), integer(1L)))
  expect_equal(total, 8L)
})

test_that("counts: stump unweighted leaves have class_counts_raw == class_counts_weighted", {
  tree <- .counts_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  leaves <- .leaf_nodes(tree)
  ok <- vapply(leaves,
               function(nd) isTRUE(all.equal(unname(as.numeric(nd$class_counts_raw)),
                                             unname(nd$class_counts_weighted))),
               logical(1L))
  expect_true(all(ok))
})

test_that("counts: stump zero-count class preserved in names", {
  # y = c(0,0,0,0,1,1,1,1) with cut ~4.5 → left leaf all class-0,
  # right leaf all class-1.  Each leaf should have both "0" and "1" names
  # with one entry = 0.
  tree <- .counts_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  leaves <- .leaf_nodes(tree)
  all_names_ok <- vapply(leaves,
                         function(nd) setequal(names(nd$class_counts_raw), c("0","1")),
                         logical(1L))
  expect_true(all(all_names_ok))
  # At least one leaf must have a zero count for one class.
  any_zero <- any(vapply(leaves,
                         function(nd) any(nd$class_counts_raw == 0L),
                         logical(1L)))
  expect_true(any_zero)
})

# =============================================================================
# Weighted stump: raw vs weighted counts differ
# =============================================================================

test_that("counts: weighted fit class_counts_raw are integer raw obs counts", {
  tree <- .counts_weighted_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  leaves <- .leaf_nodes(tree)
  ok <- vapply(leaves, function(nd) is.integer(nd$class_counts_raw), logical(1L))
  expect_true(all(ok))
})

test_that("counts: weighted fit class_counts_weighted differ from class_counts_raw when mixed", {
  tree <- .counts_weighted_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  leaves <- .leaf_nodes(tree)
  # At least one leaf should have raw != weighted (because weights differ by class).
  differs <- vapply(leaves,
                    function(nd) !isTRUE(all.equal(as.numeric(nd$class_counts_raw),
                                                   nd$class_counts_weighted)),
                    logical(1L))
  expect_true(any(differs))
})

test_that("counts: weighted fit class_counts_raw sums to n_obs on each leaf", {
  tree <- .counts_weighted_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  leaves <- .leaf_nodes(tree)
  ok <- vapply(leaves, function(nd) sum(nd$class_counts_raw) == nd$n_obs, logical(1L))
  expect_true(all(ok))
})

test_that("counts: weighted fit class_counts_weighted sums to n_weighted on each leaf", {
  tree <- .counts_weighted_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  leaves <- .leaf_nodes(tree)
  ok <- vapply(leaves,
               function(nd) isTRUE(all.equal(sum(nd$class_counts_weighted), nd$n_weighted)),
               logical(1L))
  expect_true(all(ok))
})

# =============================================================================
# Reconciliation with cta_confusion_table()
# =============================================================================

test_that("counts: stump actual-class totals from leaf counts match confusion table", {
  # cta_confusion_table() returns a tidy data frame: columns actual, predicted, n.
  # Actual class totals: sum(n) where actual == class.
  # These must match sum of class_counts_raw[class] across all leaves.
  tree <- .counts_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  leaves  <- .leaf_nodes(tree)
  ct      <- cta_confusion_table(tree)
  # Actual class totals from confusion table.
  ct_total_0 <- sum(ct$n[ct$actual == 0L])
  ct_total_1 <- sum(ct$n[ct$actual == 1L])
  # Actual class totals from stored leaf counts.
  stored_0   <- sum(vapply(leaves, function(nd) nd$class_counts_raw[["0"]], integer(1L)))
  stored_1   <- sum(vapply(leaves, function(nd) nd$class_counts_raw[["1"]], integer(1L)))
  expect_equal(stored_0, ct_total_0)
  expect_equal(stored_1, ct_total_1)
})

# =============================================================================
# Slow fixture tests — myeloma canon
# =============================================================================

.counts_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.counts_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.counts_myeloma_fit <- function(mindenom) {
  df <- .counts_load_myeloma()
  suppressMessages(
    oda_cta_fit(
      X           = df[, .counts_myeloma_attrs],
      y           = as.integer(df$V1),
      w           = df$V2,
      priors_on   = TRUE,
      miss_codes  = -9,
      alpha_split = 0.05,
      mindenom    = mindenom,
      prune_alpha = 0.05,
      max_depth   = 20L,
      mc_iter     = 5000L,
      mc_target   = 0.05,
      mc_stop     = 99.9,
      mc_stopup   = 99.9,
      mc_seed     = NULL,
      loo         = "stable",
      attr_names  = .counts_myeloma_attrs
    )
  )
}

test_that("counts: myeloma MINDENOM=1 total class_counts_raw sums to 255", {
  skip_if_slow_tests_disabled("cta-endpoint-counts")
  tree   <- .counts_myeloma_fit(1L)
  leaves <- .leaf_nodes(tree)
  total  <- sum(vapply(leaves, function(nd) sum(nd$class_counts_raw), integer(1L)))
  expect_equal(total, 255L)
})

test_that("counts: myeloma MINDENOM=30 total class_counts_raw sums to 186", {
  skip_if_slow_tests_disabled("cta-endpoint-counts")
  tree   <- .counts_myeloma_fit(30L)
  leaves <- .leaf_nodes(tree)
  total  <- sum(vapply(leaves, function(nd) sum(nd$class_counts_raw), integer(1L)))
  expect_equal(total, 186L)
})

test_that("counts: myeloma MINDENOM=56 no-tree leaf has class_counts_raw safely", {
  skip_if_slow_tests_disabled("cta-endpoint-counts")
  tree   <- .counts_myeloma_fit(56L)
  expect_true(isTRUE(tree$no_tree))
  leaves <- .leaf_nodes(tree)
  expect_equal(length(leaves), 1L)
  expect_true(!is.null(leaves[[1L]]$class_counts_raw))
  expect_equal(names(leaves[[1L]]$class_counts_raw), c("0", "1"))
})

test_that("counts: myeloma MINDENOM=1 all leaf names are c('0','1')", {
  skip_if_slow_tests_disabled("cta-endpoint-counts")
  tree   <- .counts_myeloma_fit(1L)
  leaves <- .leaf_nodes(tree)
  ok <- vapply(leaves, function(nd) identical(names(nd$class_counts_raw), c("0","1")), logical(1L))
  expect_true(all(ok))
})

test_that("counts: myeloma MINDENOM=1 per-leaf class_counts_raw sums to n_obs", {
  skip_if_slow_tests_disabled("cta-endpoint-counts")
  tree   <- .counts_myeloma_fit(1L)
  leaves <- .leaf_nodes(tree)
  ok <- vapply(leaves, function(nd) sum(nd$class_counts_raw) == nd$n_obs, logical(1L))
  expect_true(all(ok))
})

test_that("counts: myeloma MINDENOM=1 class-0 total matches confusion table", {
  skip_if_slow_tests_disabled("cta-endpoint-counts")
  tree   <- .counts_myeloma_fit(1L)
  leaves <- .leaf_nodes(tree)
  ct     <- cta_confusion_table(tree)
  stored_0   <- sum(vapply(leaves, function(nd) nd$class_counts_raw[["0"]], integer(1L)))
  ct_total_0 <- sum(ct$n[ct$actual == 0L])
  expect_equal(stored_0, ct_total_0)
})

test_that("counts: myeloma MINDENOM=1 class-1 total matches confusion table", {
  skip_if_slow_tests_disabled("cta-endpoint-counts")
  tree   <- .counts_myeloma_fit(1L)
  leaves <- .leaf_nodes(tree)
  ct     <- cta_confusion_table(tree)
  stored_1   <- sum(vapply(leaves, function(nd) nd$class_counts_raw[["1"]], integer(1L)))
  ct_total_1 <- sum(ct$n[ct$actual == 1L])
  expect_equal(stored_1, ct_total_1)
})

test_that("counts: myeloma MINDENOM=30 per-leaf class_counts_weighted sums to n_weighted", {
  skip_if_slow_tests_disabled("cta-endpoint-counts")
  tree   <- .counts_myeloma_fit(30L)
  leaves <- .leaf_nodes(tree)
  ok <- vapply(leaves,
               function(nd) isTRUE(all.equal(sum(nd$class_counts_weighted), nd$n_weighted)),
               logical(1L))
  expect_true(all(ok))
})
