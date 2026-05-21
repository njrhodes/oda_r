###############################################################################
# test-cta-endpoint-table.R — cta_endpoint_table() (Phase 2B)
###############################################################################

# ---- Helpers ----------------------------------------------------------------

.eptest_no_tree_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 999L, mc_iter = 100L,
                mc_seed = 1L, loo = "off")
  )
}

.eptest_stump_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off")
  )
}

# Synthetic 3-leaf tree (2 splits) — mirrors .synthetic_valid_tree() in test-cta-summary.R
# but with parent linkage and child_ids for path reconstruction.
.synthetic_3leaf_tree <- function() {
  # Layout:
  #   node 1 (split): x<=4.5 → node 2 (leaf, class 0)
  #                   x>4.5  → node 3 (split): x<=6.5 → node 4 (leaf, class 0)
  #                                             x>6.5  → node 5 (leaf, class 1)
  node1 <- list(
    node_id = 1L, parent_id = 0L, depth = 1L, leaf = FALSE,
    attribute = "x", attr_col = 1L, attr_type = "ordered",
    n_obs = 10L, n_weighted = 10,
    rule = list(type = "ordered_cut", cut_value = 4.5, direction = "0->1"),
    ess = 60.0, ess_weighted = NA_real_, p_mc = 0.01,
    loo_status = "STABLE", loo_ess = 58.0, loo_p = NA_real_,
    confusion = NULL,
    split_labels = c(0L, 1L),
    child_ids    = c(2L, 3L)
  )
  node2 <- list(
    node_id = 2L, parent_id = 1L, depth = 2L, leaf = TRUE,
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    n_obs = 4L, n_weighted = 4, majority_class = 0L,
    rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  node3 <- list(
    node_id = 3L, parent_id = 1L, depth = 2L, leaf = FALSE,
    attribute = "x", attr_col = 1L, attr_type = "ordered",
    n_obs = 6L, n_weighted = 6,
    rule = list(type = "ordered_cut", cut_value = 6.5, direction = "0->1"),
    ess = 50.0, ess_weighted = NA_real_, p_mc = 0.02,
    loo_status = "STABLE", loo_ess = 48.0, loo_p = NA_real_,
    confusion = NULL,
    split_labels = c(0L, 1L),
    child_ids    = c(4L, 5L)
  )
  node4 <- list(
    node_id = 4L, parent_id = 3L, depth = 3L, leaf = TRUE,
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    n_obs = 3L, n_weighted = 3, majority_class = 0L,
    rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  node5 <- list(
    node_id = 5L, parent_id = 3L, depth = 3L, leaf = TRUE,
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    n_obs = 3L, n_weighted = 3, majority_class = 1L,
    rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  structure(
    list(nodes       = list(node1, node2, node3, node4, node5),
         no_tree     = FALSE,
         overall_ess = 75.0,
         n_nodes     = 5L,
         root_id     = 1L,
         has_weights = FALSE,
         mindenom    = 2L,
         alpha_split = 0.05,
         prune_alpha = 1.0,
         loo         = "off"),
    class = "cta_tree"
  )
}

# =============================================================================
# Return type and column structure
# =============================================================================

test_that("cta_endpoint_table: returns a data.frame", {
  df <- cta_endpoint_table(.eptest_no_tree_fit())
  expect_s3_class(df, "data.frame")
})

test_that("cta_endpoint_table: no-tree has zero rows", {
  df <- cta_endpoint_table(.eptest_no_tree_fit())
  expect_equal(nrow(df), 0L)
})

test_that("cta_endpoint_table: required columns present", {
  required <- c("node_id", "parent_id", "depth", "path",
                "majority_class", "n_obs", "n_weighted",
                "ess", "ess_weighted", "loo_status", "loo_ess")
  df <- cta_endpoint_table(.eptest_no_tree_fit())
  expect_true(all(required %in% names(df)),
              info = paste("missing:", paste(setdiff(required, names(df)),
                                             collapse = ", ")))
})

# =============================================================================
# Stump: 2 leaves
# =============================================================================

test_that("cta_endpoint_table: stump has 2 rows", {
  tree <- .eptest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed — stump tests need a split")
  df <- cta_endpoint_table(tree)
  expect_equal(nrow(df), 2L)
})

test_that("cta_endpoint_table: stump node_ids are integers", {
  tree <- .eptest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_table(tree)
  expect_type(df$node_id, "integer")
})

test_that("cta_endpoint_table: stump n_obs are positive integers", {
  tree <- .eptest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_table(tree)
  expect_true(all(df$n_obs > 0L))
  expect_type(df$n_obs, "integer")
})

test_that("cta_endpoint_table: stump n_obs sums to <= total obs", {
  tree <- .eptest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_table(tree)
  # Both leaves together should sum to 8 (perfectly separable 8-obs dataset)
  expect_equal(sum(df$n_obs), 8L)
})

test_that("cta_endpoint_table: stump paths are non-empty strings", {
  tree <- .eptest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_table(tree)
  expect_true(all(nzchar(df$path)))
})

test_that("cta_endpoint_table: stump paths contain attribute name", {
  tree <- .eptest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_table(tree)
  # The attribute is 'x'; both branch strings should contain it
  expect_true(all(grepl("x", df$path, fixed = TRUE)))
})

test_that("cta_endpoint_table: stump majority_class values are 0 or 1", {
  tree <- .eptest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_table(tree)
  expect_true(all(df$majority_class %in% c(0L, 1L)))
})

test_that("cta_endpoint_table: stump ess columns are NA (leaf nodes)", {
  tree <- .eptest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_table(tree)
  expect_true(all(is.na(df$ess)))
  expect_true(all(is.na(df$ess_weighted)))
})

test_that("cta_endpoint_table: stump rows sorted by node_id", {
  tree <- .eptest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_table(tree)
  expect_equal(df$node_id, sort(df$node_id))
})

# =============================================================================
# Synthetic 3-leaf tree — path reconstruction and structural checks
# =============================================================================

test_that("cta_endpoint_table: synthetic 3-leaf tree has 3 rows", {
  df <- cta_endpoint_table(.synthetic_3leaf_tree())
  expect_equal(nrow(df), 3L)
})

test_that("cta_endpoint_table: synthetic 3-leaf node_ids are 2, 4, 5", {
  df <- cta_endpoint_table(.synthetic_3leaf_tree())
  expect_equal(sort(df$node_id), c(2L, 4L, 5L))
})

test_that("cta_endpoint_table: synthetic 3-leaf node2 path = 'x<=4.5'", {
  df <- cta_endpoint_table(.synthetic_3leaf_tree())
  nd2 <- df[df$node_id == 2L, , drop = FALSE]
  expect_equal(nd2$path, "x<=4.5")
})

test_that("cta_endpoint_table: synthetic 3-leaf node4 path = 'x>4.5 AND x<=6.5'", {
  df <- cta_endpoint_table(.synthetic_3leaf_tree())
  nd4 <- df[df$node_id == 4L, , drop = FALSE]
  expect_equal(nd4$path, "x>4.5 AND x<=6.5")
})

test_that("cta_endpoint_table: synthetic 3-leaf node5 path = 'x>4.5 AND x>6.5'", {
  df <- cta_endpoint_table(.synthetic_3leaf_tree())
  nd5 <- df[df$node_id == 5L, , drop = FALSE]
  expect_equal(nd5$path, "x>4.5 AND x>6.5")
})

test_that("cta_endpoint_table: synthetic 3-leaf depths correct", {
  df <- cta_endpoint_table(.synthetic_3leaf_tree())
  expect_equal(df[df$node_id == 2L, "depth"], 2L)
  expect_equal(df[df$node_id == 4L, "depth"], 3L)
  expect_equal(df[df$node_id == 5L, "depth"], 3L)
})

test_that("cta_endpoint_table: synthetic 3-leaf majority classes", {
  df <- cta_endpoint_table(.synthetic_3leaf_tree())
  expect_equal(df[df$node_id == 2L, "majority_class"], 0L)
  expect_equal(df[df$node_id == 4L, "majority_class"], 0L)
  expect_equal(df[df$node_id == 5L, "majority_class"], 1L)
})

test_that("cta_endpoint_table: synthetic 3-leaf n_obs correct", {
  df <- cta_endpoint_table(.synthetic_3leaf_tree())
  expect_equal(df[df$node_id == 2L, "n_obs"], 4L)
  expect_equal(df[df$node_id == 4L, "n_obs"], 3L)
  expect_equal(df[df$node_id == 5L, "n_obs"], 3L)
})

# =============================================================================
# Slow fixture tests — gated (canonical myeloma settings)
# =============================================================================

.eptest_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.eptest_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.eptest_myeloma_fit <- function(mindenom) {
  df <- .eptest_load_myeloma()
  suppressMessages(
    oda_cta_fit(
      X           = df[, .eptest_myeloma_attrs],
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
      attr_names  = .eptest_myeloma_attrs
    )
  )
}

test_that("cta_endpoint_table: myeloma MINDENOM=1 has 3 rows, correct paths", {
  skip_if_slow_tests_disabled("cta-endpoint-fixture")
  tree <- .eptest_myeloma_fit(1L)
  df   <- cta_endpoint_table(tree)
  # LOO STABLE selects V14->V15 (strata=3).
  expect_equal(nrow(df), 3L)
  expect_equal(nrow(df), cta_strata(tree))
  # All paths must contain "V14" (root attribute)
  expect_true(all(grepl("V14", df$path, fixed = TRUE)))
  # At least one path must contain "V15" (child attribute)
  expect_true(any(grepl("V15", df$path, fixed = TRUE)))
})

test_that("cta_endpoint_table: myeloma MINDENOM=56 returns zero-row df", {
  skip_if_slow_tests_disabled("cta-endpoint-fixture")
  tree <- .eptest_myeloma_fit(56L)
  df   <- cta_endpoint_table(tree)
  expect_equal(nrow(df), 0L)
  expect_true("path" %in% names(df))
})

test_that("cta_endpoint_table: myeloma MINDENOM=30 stump has 2 rows, V17 in paths", {
  skip_if_slow_tests_disabled("cta-endpoint-fixture")
  tree <- .eptest_myeloma_fit(30L)
  df   <- cta_endpoint_table(tree)
  expect_equal(nrow(df), 2L)
  expect_true(all(grepl("V17", df$path, fixed = TRUE)))
})
