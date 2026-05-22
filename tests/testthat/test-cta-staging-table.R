###############################################################################
# test-cta-staging-table.R — cta_staging_table()
#
# Per-endpoint staging table ordered by ascending target-class propensity.
# Consumes cta_endpoint_summary() + cta_endpoint_counts() only.
# No refitting, no prediction, no recomputation of tree metrics.
###############################################################################

# ---- Helpers ----------------------------------------------------------------

.est_no_tree_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 999L, mc_iter = 100L,
                mc_seed = 1L, loo = "off")
  )
}

.est_stump_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off")
  )
}

.est_weighted_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  w <- c(2, 2, 2, 2, 3, 3, 3, 3)
  suppressMessages(
    oda_cta_fit(d$X, d$y, w = w, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off")
  )
}

# Synthetic 2-leaf tree with known class counts for exact arithmetic checks.
#   node2 (x <= 6.5): class0=6, class1=2 → target_n=2, denom=8, prop=0.25
#   node3 (x >  6.5): class0=1, class1=7 → target_n=7, denom=8, prop=0.875
.est_2leaf_known <- function() {
  node1 <- list(
    node_id = 1L, parent_id = 0L, depth = 1L, leaf = FALSE,
    attribute = "x", attr_col = 1L, attr_type = "ordered",
    n_obs = 16L, n_weighted = 16,
    rule = list(type = "ordered_cut", cut_value = 6.5, direction = "0->1"),
    ess = 60.0, ess_weighted = NA_real_, p_mc = 0.01,
    loo_status = "STABLE", loo_ess = 58.0, loo_p = NA_real_,
    confusion = NULL, split_labels = c(0L, 1L), child_ids = c(2L, 3L)
  )
  node2 <- list(
    node_id = 2L, parent_id = 1L, depth = 2L, leaf = TRUE,
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    n_obs = 8L, n_weighted = 8, majority_class = 0L, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    class_counts_raw      = c("0" = 6L, "1" = 2L),
    class_counts_weighted = c("0" = 6.0, "1" = 2.0),
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  node3 <- list(
    node_id = 3L, parent_id = 1L, depth = 2L, leaf = TRUE,
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    n_obs = 8L, n_weighted = 8, majority_class = 1L, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    class_counts_raw      = c("0" = 1L, "1" = 7L),
    class_counts_weighted = c("0" = 1.0, "1" = 7.0),
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  structure(
    list(nodes       = list(node1, node2, node3),
         no_tree     = FALSE,
         overall_ess = 65.0,
         n_nodes     = 3L,
         root_id     = 1L,
         has_weights = FALSE,
         mindenom    = 2L,
         alpha_split = 0.05,
         prune_alpha = 1.0,
         loo         = "off"),
    class = "cta_tree"
  )
}

# Synthetic 2-leaf tree: one normal endpoint, one perfectly predicted endpoint.
#   node2 (x <= 6.5): class0=6, class1=2 → not perfectly predicted
#   node3 (x >  6.5): class0=0, class1=4 → perfectly predicted (non_target_n=0)
.est_perfect_ep_tree <- function() {
  node1 <- list(
    node_id = 1L, parent_id = 0L, depth = 1L, leaf = FALSE,
    attribute = "x", attr_col = 1L, attr_type = "ordered",
    n_obs = 12L, n_weighted = 12,
    rule = list(type = "ordered_cut", cut_value = 6.5, direction = "0->1"),
    ess = 55.0, ess_weighted = NA_real_, p_mc = 0.01,
    loo_status = "STABLE", loo_ess = 53.0, loo_p = NA_real_,
    confusion = NULL, split_labels = c(0L, 1L), child_ids = c(2L, 3L)
  )
  node2 <- list(
    node_id = 2L, parent_id = 1L, depth = 2L, leaf = TRUE,
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    n_obs = 8L, n_weighted = 8, majority_class = 0L, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    class_counts_raw      = c("0" = 6L, "1" = 2L),
    class_counts_weighted = c("0" = 6.0, "1" = 2.0),
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  node3 <- list(
    node_id = 3L, parent_id = 1L, depth = 2L, leaf = TRUE,
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    n_obs = 4L, n_weighted = 4, majority_class = 1L, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    class_counts_raw      = c("0" = 0L, "1" = 4L),
    class_counts_weighted = c("0" = 0.0, "1" = 4.0),
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  structure(
    list(nodes       = list(node1, node2, node3),
         no_tree     = FALSE,
         overall_ess = 60.0,
         n_nodes     = 3L,
         root_id     = 1L,
         has_weights = FALSE,
         mindenom    = 2L,
         alpha_split = 0.05,
         prune_alpha = 1.0,
         loo         = "off"),
    class = "cta_tree"
  )
}

# Synthetic 2-leaf tree with weighted counts differing from raw counts.
#   node2: class0_raw=4, class1_raw=2; w0=8.0, w1=6.0
#   node3: class0_raw=1, class1_raw=3; w0=2.0, w1=9.0
.est_weighted_known <- function() {
  node1 <- list(
    node_id = 1L, parent_id = 0L, depth = 1L, leaf = FALSE,
    attribute = "x", attr_col = 1L, attr_type = "ordered",
    n_obs = 10L, n_weighted = 25,
    rule = list(type = "ordered_cut", cut_value = 6.5, direction = "0->1"),
    ess = 55.0, ess_weighted = NA_real_, p_mc = 0.01,
    loo_status = "STABLE", loo_ess = 53.0, loo_p = NA_real_,
    confusion = NULL, split_labels = c(0L, 1L), child_ids = c(2L, 3L)
  )
  node2 <- list(
    node_id = 2L, parent_id = 1L, depth = 2L, leaf = TRUE,
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    n_obs = 6L, n_weighted = 14, majority_class = 0L, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    class_counts_raw      = c("0" = 4L, "1" = 2L),
    class_counts_weighted = c("0" = 8.0, "1" = 6.0),
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  node3 <- list(
    node_id = 3L, parent_id = 1L, depth = 2L, leaf = TRUE,
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    n_obs = 4L, n_weighted = 11, majority_class = 1L, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    class_counts_raw      = c("0" = 1L, "1" = 3L),
    class_counts_weighted = c("0" = 2.0, "1" = 9.0),
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  structure(
    list(nodes       = list(node1, node2, node3),
         no_tree     = FALSE,
         overall_ess = 60.0,
         n_nodes     = 3L,
         root_id     = 1L,
         has_weights = TRUE,
         mindenom    = 2L,
         alpha_split = 0.05,
         prune_alpha = 1.0,
         loo         = "off"),
    class = "cta_tree"
  )
}

# Synthetic 3-class tree (for error-on-NULL target_class check).
.est_3class_tree <- function() {
  leaf <- list(
    node_id = 1L, parent_id = 0L, depth = 1L, leaf = TRUE,
    majority_class = 1L, n_obs = 9L, n_weighted = 9,
    class_counts_raw      = c("0" = 3L, "1" = 3L, "2" = 3L),
    class_counts_weighted = c("0" = 3.0, "1" = 3.0, "2" = 3.0),
    attribute = NA_character_, attr_col = NA_integer_,
    attr_type = NA_character_, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  structure(
    list(nodes = list(leaf), no_tree = FALSE, overall_ess = NA_real_,
         n_nodes = 1L, root_id = 1L, has_weights = FALSE,
         mindenom = 2L, alpha_split = 0.05, prune_alpha = 1.0, loo = "off"),
    class = "cta_tree"
  )
}

# =============================================================================
# Fast tests
# =============================================================================

# ---- No-tree ----------------------------------------------------------------

test_that("est: no-tree fit returns empty df with 21 columns", {
  df <- cta_staging_table(.est_no_tree_fit())
  expect_equal(nrow(df), 0L)
  expected_cols <- c(
    "stage", "endpoint_id", "endpoint_node_id", "path",
    "terminal_prediction", "target_class", "target_n", "denominator",
    "target_proportion", "non_target_n", "odds", "perfectly_predicted",
    "adjusted", "adjusted_target_n", "adjusted_denominator",
    "adjusted_target_proportion", "adjusted_non_target_n", "adjusted_odds",
    "weighted", "n_obs", "n_weighted"
  )
  expect_equal(names(df), expected_cols)
})

# ---- Column types -----------------------------------------------------------

test_that("est: column types are correct for 21-column schema", {
  df <- cta_staging_table(.est_2leaf_known())
  expect_true(is.integer(df$stage))
  expect_true(is.integer(df$endpoint_id))
  expect_true(is.integer(df$endpoint_node_id))
  expect_true(is.character(df$path))
  expect_true(is.integer(df$terminal_prediction))
  expect_true(is.integer(df$target_class))
  expect_true(is.numeric(df$target_n))
  expect_true(is.numeric(df$denominator))
  expect_true(is.numeric(df$target_proportion))
  expect_true(is.numeric(df$non_target_n))
  expect_true(is.numeric(df$odds))
  expect_true(is.logical(df$perfectly_predicted))
  expect_true(is.logical(df$adjusted))
  expect_true(is.numeric(df$adjusted_target_n))
  expect_true(is.numeric(df$adjusted_denominator))
  expect_true(is.numeric(df$adjusted_target_proportion))
  expect_true(is.numeric(df$adjusted_non_target_n))
  expect_true(is.numeric(df$adjusted_odds))
  expect_true(is.logical(df$weighted))
  expect_true(is.integer(df$n_obs))
  expect_true(is.numeric(df$n_weighted))
})

# ---- Row count and stage sequence ------------------------------------------

test_that("est: 2-leaf tree returns 2 rows", {
  df <- cta_staging_table(.est_2leaf_known())
  expect_equal(nrow(df), 2L)
})

test_that("est: stage is always 1..nrow", {
  df <- cta_staging_table(.est_2leaf_known())
  expect_equal(df$stage, seq_len(nrow(df)))
})

# ---- target_class resolution ------------------------------------------------

test_that("est: NULL target_class defaults to class 1 for binary tree", {
  df <- cta_staging_table(.est_2leaf_known())
  expect_equal(unique(df$target_class), 1L)
})

test_that("est: NULL target_class stops for 3-class tree", {
  expect_error(
    cta_staging_table(.est_3class_tree()),
    regexp = "target_class must be specified"
  )
})

test_that("est: unknown target_class stops with informative error", {
  expect_error(
    cta_staging_table(.est_2leaf_known(), target_class = 9L),
    regexp = "not found in tree classes"
  )
})

test_that("est: explicit target_class=0 uses class 0 as target", {
  df <- cta_staging_table(.est_2leaf_known(), target_class = 0L)
  expect_equal(unique(df$target_class), 0L)
  # class-0 counts: node2=6, node3=1
  expect_equal(sort(df$target_n), c(1.0, 6.0))
})

# ---- Arithmetic: known tree -------------------------------------------------

test_that("est: target_n and denominator correct for known tree (weighted=FALSE)", {
  df <- cta_staging_table(.est_2leaf_known())
  # stage 1 = lower proportion endpoint (node2: 2/8=0.25)
  # stage 2 = higher proportion endpoint (node3: 7/8=0.875)
  expect_equal(df$target_n[df$endpoint_node_id == 2L], 2.0)
  expect_equal(df$denominator[df$endpoint_node_id == 2L], 8.0)
  expect_equal(df$target_n[df$endpoint_node_id == 3L], 7.0)
  expect_equal(df$denominator[df$endpoint_node_id == 3L], 8.0)
})

test_that("est: target_proportion = target_n / denominator", {
  df <- cta_staging_table(.est_2leaf_known())
  expect_equal(df$target_proportion,
               df$target_n / df$denominator,
               tolerance = 1e-12)
})

test_that("est: non_target_n = denominator - target_n", {
  df <- cta_staging_table(.est_2leaf_known())
  expect_equal(df$non_target_n, df$denominator - df$target_n)
})

test_that("est: odds = target_n / non_target_n for non-perfect endpoints", {
  df <- cta_staging_table(.est_2leaf_known())
  expect_false(any(df$perfectly_predicted))
  expect_equal(df$odds, df$target_n / df$non_target_n, tolerance = 1e-12)
})

# ---- Stage ordering ---------------------------------------------------------

test_that("est: rows ordered by ascending adjusted_target_proportion (default)", {
  df <- cta_staging_table(.est_2leaf_known())
  expect_true(all(diff(df$adjusted_target_proportion) >= 0))
})

test_that("est: rows ordered by ascending target_proportion when adjust_perfect=FALSE", {
  df <- cta_staging_table(.est_2leaf_known(), adjust_perfect = FALSE)
  expect_true(all(diff(df$target_proportion) >= 0))
})

test_that("est: lower-proportion endpoint is stage 1", {
  df <- cta_staging_table(.est_2leaf_known())
  expect_equal(df$endpoint_node_id[1L], 2L)   # node2: prop 0.25
  expect_equal(df$endpoint_node_id[2L], 3L)   # node3: prop 0.875
})

# ---- Perfectly predicted detection -----------------------------------------

test_that("est: non_target_n==0 endpoint detected as perfectly_predicted", {
  df <- cta_staging_table(.est_perfect_ep_tree())
  # node3: class0=0, class1=4 → non_target_n=0 → perfectly_predicted
  expect_true(df$perfectly_predicted[df$endpoint_node_id == 3L])
})

test_that("est: endpoint with both classes present is not perfectly_predicted", {
  df <- cta_staging_table(.est_perfect_ep_tree())
  expect_false(df$perfectly_predicted[df$endpoint_node_id == 2L])
})

test_that("est: odds is NA for perfectly_predicted endpoint", {
  df <- cta_staging_table(.est_perfect_ep_tree())
  expect_true(is.na(df$odds[df$endpoint_node_id == 3L]))
  expect_false(is.na(df$odds[df$endpoint_node_id == 2L]))
})

# ---- Adjustment (adjust_perfect=TRUE default) --------------------------------

test_that("est: adjusted=TRUE for perfectly_predicted when adjust_perfect=TRUE", {
  df <- cta_staging_table(.est_perfect_ep_tree())
  expect_true(df$adjusted[df$endpoint_node_id == 3L])
  expect_false(df$adjusted[df$endpoint_node_id == 2L])
})

test_that("est: adjusted=FALSE for all endpoints when adjust_perfect=FALSE", {
  df <- cta_staging_table(.est_perfect_ep_tree(), adjust_perfect = FALSE)
  expect_true(all(!df$adjusted))
})

test_that("est: non-target-side adjustment adds 1 to adj_non_target_n", {
  # node3: non_target_n=0 → boost non-target → adj_non_target_n=1
  df <- cta_staging_table(.est_perfect_ep_tree())
  row <- df[df$endpoint_node_id == 3L, ]
  expect_equal(row$adjusted_non_target_n, 1.0)
  expect_equal(row$adjusted_target_n, 4.0)           # unchanged
  expect_equal(row$adjusted_denominator, 5.0)         # 4+1
  expect_equal(row$adjusted_target_proportion, 4/5)
  expect_equal(row$adjusted_odds, 4/1)
})

test_that("est: adjusted columns equal empirical for non-adjusted endpoints", {
  df <- cta_staging_table(.est_perfect_ep_tree())
  row <- df[df$endpoint_node_id == 2L, ]
  expect_equal(row$adjusted_target_n,          row$target_n)
  expect_equal(row$adjusted_denominator,        row$denominator)
  expect_equal(row$adjusted_target_proportion,  row$target_proportion)
  expect_equal(row$adjusted_non_target_n,       row$non_target_n)
  expect_equal(row$adjusted_odds,               row$odds)
})

test_that("est: adjust_perfect=FALSE leaves adjusted columns = empirical values", {
  df <- cta_staging_table(.est_perfect_ep_tree(), adjust_perfect = FALSE)
  row <- df[df$endpoint_node_id == 3L, ]
  expect_equal(row$adjusted_target_n,     row$target_n)
  expect_equal(row$adjusted_denominator,  row$denominator)
  expect_equal(row$adjusted_non_target_n, row$non_target_n)
  # adjusted_odds is NA (non_target_n=0), same as odds
  expect_true(is.na(row$adjusted_odds))
})

# ---- weighted parameter ------------------------------------------------------

test_that("est: weighted column stores the weighted= argument (FALSE)", {
  df <- cta_staging_table(.est_2leaf_known(), weighted = FALSE)
  expect_true(all(df$weighted == FALSE))
})

test_that("est: weighted column stores the weighted= argument (TRUE)", {
  df <- cta_staging_table(.est_weighted_known(), weighted = TRUE)
  expect_true(all(df$weighted == TRUE))
})

test_that("est: weighted=FALSE uses n_raw counts", {
  df <- cta_staging_table(.est_weighted_known(), weighted = FALSE)
  # node2: n_raw class1 = 2, denom_raw = 6
  expect_equal(df$target_n[df$endpoint_node_id == 2L], 2.0)
  expect_equal(df$denominator[df$endpoint_node_id == 2L], 6.0)
})

test_that("est: weighted=TRUE uses n_weighted counts", {
  df <- cta_staging_table(.est_weighted_known(), weighted = TRUE)
  # node2: w_class1 = 6.0, w_denom = 14.0
  expect_equal(df$target_n[df$endpoint_node_id == 2L], 6.0)
  expect_equal(df$denominator[df$endpoint_node_id == 2L], 14.0)
})

# ---- n_obs / n_weighted passthrough -----------------------------------------

test_that("est: n_obs matches cta_endpoint_summary()$n_obs", {
  tree <- .est_2leaf_known()
  df   <- cta_staging_table(tree)
  es   <- cta_endpoint_summary(tree)
  # Align by endpoint_node_id (staging table may be reordered)
  for (nid in es$endpoint_node_id) {
    expect_equal(df$n_obs[df$endpoint_node_id == nid],
                 es$n_obs[es$endpoint_node_id == nid])
  }
})

test_that("est: n_weighted matches cta_endpoint_summary()$n_weighted", {
  tree <- .est_2leaf_known()
  df   <- cta_staging_table(tree)
  es   <- cta_endpoint_summary(tree)
  for (nid in es$endpoint_node_id) {
    expect_equal(df$n_weighted[df$endpoint_node_id == nid],
                 es$n_weighted[es$endpoint_node_id == nid])
  }
})

# ---- Stump fit (real oda_cta_fit call) --------------------------------------

test_that("est: stump fit returns 2 rows with correct columns", {
  df <- cta_staging_table(.est_stump_fit())
  expect_equal(nrow(df), 2L)
  expect_equal(names(df), c(
    "stage", "endpoint_id", "endpoint_node_id", "path",
    "terminal_prediction", "target_class", "target_n", "denominator",
    "target_proportion", "non_target_n", "odds", "perfectly_predicted",
    "adjusted", "adjusted_target_n", "adjusted_denominator",
    "adjusted_target_proportion", "adjusted_non_target_n", "adjusted_odds",
    "weighted", "n_obs", "n_weighted"
  ))
})

test_that("est: stump fit denominator sums match n_obs from endpoint_summary", {
  tree <- .est_stump_fit()
  df   <- cta_staging_table(tree)
  es   <- cta_endpoint_summary(tree)
  for (nid in es$endpoint_node_id) {
    expect_equal(df$denominator[df$endpoint_node_id == nid],
                 as.numeric(es$n_obs[es$endpoint_node_id == nid]))
  }
})

test_that("est: stump fit target_proportion is within [0,1]", {
  df <- cta_staging_table(.est_stump_fit())
  expect_true(all(df$target_proportion >= 0 & df$target_proportion <= 1,
                  na.rm = TRUE))
})

test_that("est: weighted stump fit with weighted=TRUE: denominator = n_weighted", {
  tree <- .est_weighted_fit()
  df   <- cta_staging_table(tree, weighted = TRUE)
  es   <- cta_endpoint_summary(tree)
  for (nid in es$endpoint_node_id) {
    expect_equal(df$denominator[df$endpoint_node_id == nid],
                 es$n_weighted[es$endpoint_node_id == nid])
  }
})

# =============================================================================
# Slow fixture tests — myeloma canon
# =============================================================================

.est_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.est_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.est_myeloma_fit <- function(mindenom) {
  df <- .est_load_myeloma()
  suppressMessages(
    oda_cta_fit(
      X           = df[, .est_myeloma_attrs],
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
      attr_names  = .est_myeloma_attrs
    )
  )
}

test_that("est: myeloma MINDENOM=56 (no-tree) returns 0 rows with 21 columns", {
  skip_if_slow_tests_disabled("cta-staging-table")
  df <- cta_staging_table(.est_myeloma_fit(56L))
  expect_equal(nrow(df), 0L)
  expect_equal(length(names(df)), 21L)
})

test_that("est: myeloma MINDENOM=30 (stump) returns 2 rows", {
  skip_if_slow_tests_disabled("cta-staging-table")
  df <- cta_staging_table(.est_myeloma_fit(30L))
  expect_equal(nrow(df), 2L)
})

test_that("est: myeloma MINDENOM=30 stage is 1:2", {
  skip_if_slow_tests_disabled("cta-staging-table")
  df <- cta_staging_table(.est_myeloma_fit(30L))
  expect_equal(df$stage, c(1L, 2L))
})

test_that("est: myeloma MINDENOM=30 adjusted_target_proportion ascending", {
  skip_if_slow_tests_disabled("cta-staging-table")
  df <- cta_staging_table(.est_myeloma_fit(30L))
  expect_true(df$adjusted_target_proportion[1L] <=
                df$adjusted_target_proportion[2L])
})

test_that("est: myeloma MINDENOM=30 target_class = 1 by default", {
  skip_if_slow_tests_disabled("cta-staging-table")
  df <- cta_staging_table(.est_myeloma_fit(30L))
  expect_equal(unique(df$target_class), 1L)
})

test_that("est: myeloma MINDENOM=30 weighted=TRUE denominator matches n_weighted", {
  skip_if_slow_tests_disabled("cta-staging-table")
  tree <- .est_myeloma_fit(30L)
  df   <- cta_staging_table(tree, weighted = TRUE)
  es   <- cta_endpoint_summary(tree)
  for (nid in es$endpoint_node_id) {
    expect_equal(df$denominator[df$endpoint_node_id == nid],
                 es$n_weighted[es$endpoint_node_id == nid])
  }
})

test_that("est: myeloma MINDENOM=30 sum of target_n matches confusion table", {
  skip_if_slow_tests_disabled("cta-staging-table")
  tree <- .est_myeloma_fit(30L)
  df   <- cta_staging_table(tree)
  ct   <- cta_confusion_table(tree)
  expect_equal(sum(df$target_n), as.numeric(sum(ct$n[ct$actual == 1L])))
})

test_that("est: myeloma MINDENOM=1 (V14->V15 tree) returns 3 rows", {
  skip_if_slow_tests_disabled("cta-staging-table")
  df <- cta_staging_table(.est_myeloma_fit(1L))
  expect_equal(nrow(df), 3L)
})

test_that("est: myeloma MINDENOM=1 stage 1 has lowest target_proportion", {
  skip_if_slow_tests_disabled("cta-staging-table")
  df <- cta_staging_table(.est_myeloma_fit(1L))
  expect_true(df$target_proportion[1L] == min(df$target_proportion, na.rm = TRUE))
})
