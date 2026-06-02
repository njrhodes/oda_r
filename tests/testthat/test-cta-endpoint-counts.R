###############################################################################
# test-cta-endpoint-counts.R — per-leaf class-count storage contract
#
# Verifies that .leaf_nd() stores class_counts_raw and class_counts_weighted
# on every terminal node at fit time.
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
# Contract tests
# =============================================================================

test_that("counts: no-tree leaf stores class_counts_raw and _weighted correctly", {
  tree   <- .counts_no_tree_fit()
  leaves <- .leaf_nodes(tree)
  nd     <- leaves[[1L]]
  # Both fields present with correct names and types
  expect_equal(names(nd$class_counts_raw),      c("0", "1"))
  expect_equal(names(nd$class_counts_weighted),  c("0", "1"))
  expect_type(nd$class_counts_raw,      "integer")
  expect_type(nd$class_counts_weighted, "double")
  # Sums reconcile with node counts
  expect_equal(sum(nd$class_counts_raw), nd$n_obs)
  expect_equal(sum(nd$class_counts_weighted), nd$n_weighted)
  # Unweighted fit: raw and weighted agree
  expect_equal(unname(as.numeric(nd$class_counts_raw)), unname(nd$class_counts_weighted))
})

test_that("counts: stump leaves store class_counts correctly; zero-count preserved; reconciles with confusion", {
  tree <- .counts_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  leaves <- .leaf_nodes(tree)
  # All leaves have correct names and per-leaf sums
  for (nd in leaves) {
    expect_equal(names(nd$class_counts_raw), c("0", "1"))
    expect_equal(sum(nd$class_counts_raw), nd$n_obs)
    expect_true(isTRUE(all.equal(sum(nd$class_counts_weighted), nd$n_weighted)))
  }
  # Zero-count class preserved: perfectly separated data -> one zero per leaf
  expect_true(any(vapply(leaves, function(nd) any(nd$class_counts_raw == 0L), logical(1L))))
  # Totals reconcile with cta_confusion_table()
  ct       <- cta_confusion_table(tree)
  stored_0 <- sum(vapply(leaves, function(nd) nd$class_counts_raw[["0"]], integer(1L)))
  stored_1 <- sum(vapply(leaves, function(nd) nd$class_counts_raw[["1"]], integer(1L)))
  expect_equal(stored_0, sum(ct$n[ct$actual == 0L]))
  expect_equal(stored_1, sum(ct$n[ct$actual == 1L]))
})

test_that("counts: weighted leaves — raw counts are integer raw obs; differ from weighted", {
  tree <- .counts_weighted_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  leaves <- .leaf_nodes(tree)
  for (nd in leaves) {
    expect_type(nd$class_counts_raw, "integer")
    expect_equal(sum(nd$class_counts_raw), nd$n_obs)
    expect_true(isTRUE(all.equal(sum(nd$class_counts_weighted), nd$n_weighted)))
  }
  # At least one leaf has raw != weighted (weights differ by class)
  differs <- vapply(leaves,
                    function(nd) !isTRUE(all.equal(as.numeric(nd$class_counts_raw),
                                                   nd$class_counts_weighted)),
                    logical(1L))
  expect_true(any(differs))
})

# =============================================================================
# Smoke: myeloma canon — internal storage reconciles with confusion table
# =============================================================================

.counts_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.counts_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.counts_myeloma_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      df <- .counts_load_myeloma()
      fit <<- suppressMessages(oda_cta_fit(
        X = df[, .counts_myeloma_attrs], y = as.integer(df$V1), w = df$V2,
        priors_on = TRUE, miss_codes = -9, alpha_split = 0.05, mindenom = 1L,
        prune_alpha = 0.05, max_depth = 20L, mc_iter = 5000L,
        mc_target = 0.05, mc_stop = 99.9, mc_stopup = 99.9, mc_seed = NULL,
        loo = "stable", attr_names = .counts_myeloma_attrs
      ))
    }
    fit
  }
})

test_that("counts: myeloma MINDENOM=1 — leaf storage names, sums, reconcile confusion (n=255)", {
  skip_if_slow_tests_disabled("cta-endpoint-counts")
  tree   <- .counts_myeloma_fit()
  leaves <- .leaf_nodes(tree)
  # All leaves have canonical names and sums reconcile
  for (nd in leaves) {
    expect_equal(names(nd$class_counts_raw), c("0", "1"))
    expect_equal(sum(nd$class_counts_raw), nd$n_obs)
  }
  # Total classified = 255
  expect_equal(sum(vapply(leaves, function(nd) sum(nd$class_counts_raw), integer(1L))), 255L)
  # Class totals match confusion table marginals
  ct       <- cta_confusion_table(tree)
  stored_0 <- sum(vapply(leaves, function(nd) nd$class_counts_raw[["0"]], integer(1L)))
  stored_1 <- sum(vapply(leaves, function(nd) nd$class_counts_raw[["1"]], integer(1L)))
  expect_equal(stored_0, sum(ct$n[ct$actual == 0L]))
  expect_equal(stored_1, sum(ct$n[ct$actual == 1L]))
})
