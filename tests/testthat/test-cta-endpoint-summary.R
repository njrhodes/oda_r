###############################################################################
# test-cta-endpoint-summary.R — cta_endpoint_summary()
#
# Conservative endpoint reporting accessor.
# Returns one row per terminal leaf with stored structural fields only.
# Does NOT include class counts, target-class proportions, or staging order.
###############################################################################

# ---- Helpers ----------------------------------------------------------------

.epsum_no_tree_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 999L, mc_iter = 100L,
                mc_seed = 1L, loo = "off")
  )
}

.epsum_stump_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off")
  )
}

# Synthetic 3-leaf tree — matches .synthetic_3leaf_tree() in
# test-cta-endpoint-table.R; defined locally to keep this file self-contained.
.epsum_3leaf_tree <- function() {
  node1 <- list(
    node_id = 1L, parent_id = 0L, depth = 1L, leaf = FALSE,
    attribute = "x", attr_col = 1L, attr_type = "ordered",
    n_obs = 10L, n_weighted = 10,
    rule = list(type = "ordered_cut", cut_value = 4.5, direction = "0->1"),
    ess = 60.0, ess_weighted = NA_real_, p_mc = 0.01,
    loo_status = "STABLE", loo_ess = 58.0, loo_p = NA_real_,
    confusion = NULL, split_labels = c(0L, 1L), child_ids = c(2L, 3L)
  )
  node2 <- list(
    node_id = 2L, parent_id = 1L, depth = 2L, leaf = TRUE,
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    n_obs = 4L, n_weighted = 4, majority_class = 0L, rule = NULL,
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
    confusion = NULL, split_labels = c(0L, 1L), child_ids = c(4L, 5L)
  )
  node4 <- list(
    node_id = 4L, parent_id = 3L, depth = 3L, leaf = TRUE,
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    n_obs = 3L, n_weighted = 3, majority_class = 0L, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  node5 <- list(
    node_id = 5L, parent_id = 3L, depth = 3L, leaf = TRUE,
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    n_obs = 3L, n_weighted = 3, majority_class = 1L, rule = NULL,
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

test_that("cta_endpoint_summary: returns a data.frame", {
  expect_s3_class(cta_endpoint_summary(.epsum_no_tree_fit()), "data.frame")
})

test_that("cta_endpoint_summary: no-tree returns zero rows", {
  expect_equal(nrow(cta_endpoint_summary(.epsum_no_tree_fit())), 0L)
})

test_that("cta_endpoint_summary: no-tree has exactly the required columns", {
  expected <- c("endpoint_id", "endpoint_node_id", "path", "depth",
                "terminal_prediction", "n_obs", "n_weighted", "denominator")
  df <- cta_endpoint_summary(.epsum_no_tree_fit())
  expect_equal(names(df), expected)
})

test_that("cta_endpoint_summary: no-tree column types are correct", {
  df <- cta_endpoint_summary(.epsum_no_tree_fit())
  expect_type(df$endpoint_id,         "integer")
  expect_type(df$endpoint_node_id,    "integer")
  expect_type(df$path,                "character")
  expect_type(df$depth,               "integer")
  expect_type(df$terminal_prediction, "integer")
  expect_type(df$n_obs,               "integer")
  expect_type(df$n_weighted,          "double")
  expect_type(df$denominator,         "integer")
})

# =============================================================================
# Stump: 2 leaves
# =============================================================================

test_that("cta_endpoint_summary: stump returns 2 rows", {
  tree <- .epsum_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed — stump tests need a split")
  df <- cta_endpoint_summary(tree)
  expect_equal(nrow(df), 2L)
})

test_that("cta_endpoint_summary: stump endpoint_id is 1..2", {
  tree <- .epsum_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_summary(tree)
  expect_equal(df$endpoint_id, c(1L, 2L))
})

test_that("cta_endpoint_summary: stump denominator equals n_obs", {
  tree <- .epsum_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_summary(tree)
  expect_equal(df$denominator, df$n_obs)
})

test_that("cta_endpoint_summary: stump terminal_prediction values are 0 or 1", {
  tree <- .epsum_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_summary(tree)
  expect_true(all(df$terminal_prediction %in% c(0L, 1L)))
})

test_that("cta_endpoint_summary: stump path is non-empty string", {
  tree <- .epsum_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_summary(tree)
  expect_true(all(nzchar(df$path)))
})

test_that("cta_endpoint_summary: stump n_obs are positive integers", {
  tree <- .epsum_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_summary(tree)
  expect_type(df$n_obs, "integer")
  expect_true(all(df$n_obs > 0L))
})

test_that("cta_endpoint_summary: stump rows sorted by endpoint_node_id", {
  tree <- .epsum_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_summary(tree)
  expect_equal(df$endpoint_node_id, sort(df$endpoint_node_id))
})

# =============================================================================
# Synthetic 3-leaf tree
# =============================================================================

test_that("cta_endpoint_summary: synthetic 3-leaf returns 3 rows", {
  df <- cta_endpoint_summary(.epsum_3leaf_tree())
  expect_equal(nrow(df), 3L)
})

test_that("cta_endpoint_summary: synthetic 3-leaf endpoint_id is 1..3", {
  df <- cta_endpoint_summary(.epsum_3leaf_tree())
  expect_equal(df$endpoint_id, c(1L, 2L, 3L))
})

test_that("cta_endpoint_summary: synthetic 3-leaf endpoint_node_ids are 2, 4, 5", {
  df <- cta_endpoint_summary(.epsum_3leaf_tree())
  expect_equal(df$endpoint_node_id, c(2L, 4L, 5L))
})

test_that("cta_endpoint_summary: synthetic 3-leaf denominator equals n_obs", {
  df <- cta_endpoint_summary(.epsum_3leaf_tree())
  expect_equal(df$denominator, df$n_obs)
})

test_that("cta_endpoint_summary: synthetic 3-leaf terminal_prediction matches stored majority_class", {
  df <- cta_endpoint_summary(.epsum_3leaf_tree())
  expect_equal(df$terminal_prediction[df$endpoint_node_id == 2L], 0L)
  expect_equal(df$terminal_prediction[df$endpoint_node_id == 4L], 0L)
  expect_equal(df$terminal_prediction[df$endpoint_node_id == 5L], 1L)
})

test_that("cta_endpoint_summary: synthetic 3-leaf n_obs correct", {
  df <- cta_endpoint_summary(.epsum_3leaf_tree())
  expect_equal(df$n_obs[df$endpoint_node_id == 2L], 4L)
  expect_equal(df$n_obs[df$endpoint_node_id == 4L], 3L)
  expect_equal(df$n_obs[df$endpoint_node_id == 5L], 3L)
})

test_that("cta_endpoint_summary: synthetic 3-leaf depth correct", {
  df <- cta_endpoint_summary(.epsum_3leaf_tree())
  expect_equal(df$depth[df$endpoint_node_id == 2L], 2L)
  expect_equal(df$depth[df$endpoint_node_id == 4L], 3L)
  expect_equal(df$depth[df$endpoint_node_id == 5L], 3L)
})

test_that("cta_endpoint_summary: synthetic 3-leaf node2 path is 'x<=4.5'", {
  df <- cta_endpoint_summary(.epsum_3leaf_tree())
  expect_equal(df$path[df$endpoint_node_id == 2L], "x<=4.5")
})

test_that("cta_endpoint_summary: synthetic 3-leaf paths are non-empty", {
  df <- cta_endpoint_summary(.epsum_3leaf_tree())
  expect_true(all(nzchar(df$path)))
})

test_that("cta_endpoint_summary: no target_class column present", {
  df <- cta_endpoint_summary(.epsum_3leaf_tree())
  expect_false("target_class" %in% names(df))
})

test_that("cta_endpoint_summary: no staging or event-rate columns present", {
  df <- cta_endpoint_summary(.epsum_3leaf_tree())
  forbidden <- c("target_class", "target_event_n", "target_n",
                 "target_proportion", "odds", "staging_tier", "status")
  expect_true(!any(forbidden %in% names(df)))
})

# =============================================================================
# Slow fixture tests — myeloma canon
# =============================================================================

.epsum_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.epsum_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.epsum_myeloma_fit <- function(mindenom) {
  df <- .epsum_load_myeloma()
  suppressMessages(
    oda_cta_fit(
      X           = df[, .epsum_myeloma_attrs],
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
      attr_names  = .epsum_myeloma_attrs
    )
  )
}

test_that("cta_endpoint_summary: myeloma MINDENOM=1 has 3 rows", {
  skip_if_slow_tests_disabled("cta-endpoint-summary")
  df <- cta_endpoint_summary(.epsum_myeloma_fit(1L))
  expect_equal(nrow(df), 3L)
})

test_that("cta_endpoint_summary: myeloma MINDENOM=30 has 2 rows", {
  skip_if_slow_tests_disabled("cta-endpoint-summary")
  df <- cta_endpoint_summary(.epsum_myeloma_fit(30L))
  expect_equal(nrow(df), 2L)
})

test_that("cta_endpoint_summary: myeloma MINDENOM=56 returns zero rows", {
  skip_if_slow_tests_disabled("cta-endpoint-summary")
  df <- cta_endpoint_summary(.epsum_myeloma_fit(56L))
  expect_equal(nrow(df), 0L)
  expected <- c("endpoint_id", "endpoint_node_id", "path", "depth",
                "terminal_prediction", "n_obs", "n_weighted", "denominator")
  expect_equal(names(df), expected)
})

test_that("cta_endpoint_summary: myeloma MINDENOM=1 denominators match cta_endpoint_denominators()", {
  skip_if_slow_tests_disabled("cta-endpoint-summary")
  tree <- .epsum_myeloma_fit(1L)
  df   <- cta_endpoint_summary(tree)
  canon_denoms <- sort(cta_endpoint_denominators(tree))
  expect_equal(sort(df$denominator), canon_denoms)
})

test_that("cta_endpoint_summary: myeloma MINDENOM=30 denominators match cta_endpoint_denominators()", {
  skip_if_slow_tests_disabled("cta-endpoint-summary")
  tree <- .epsum_myeloma_fit(30L)
  df   <- cta_endpoint_summary(tree)
  canon_denoms <- sort(cta_endpoint_denominators(tree))
  expect_equal(sort(df$denominator), canon_denoms)
})

test_that("cta_endpoint_summary: myeloma MINDENOM=1 paths contain V14 (root attribute)", {
  skip_if_slow_tests_disabled("cta-endpoint-summary")
  df <- cta_endpoint_summary(.epsum_myeloma_fit(1L))
  expect_true(all(grepl("V14", df$path, fixed = TRUE)))
})

test_that("cta_endpoint_summary: myeloma MINDENOM=30 paths all contain V17", {
  skip_if_slow_tests_disabled("cta-endpoint-summary")
  df <- cta_endpoint_summary(.epsum_myeloma_fit(30L))
  expect_true(all(grepl("V17", df$path, fixed = TRUE)))
})
