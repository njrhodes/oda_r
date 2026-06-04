###############################################################################
# test-cta-propensity-weights.R - cta_propensity_weights()
#
# Endpoint x class stabilized propensity-style weights.
# Consumes cta_endpoint_counts(tree) only. Also covers lean-fit invariants.
###############################################################################

# ---- Helpers ----------------------------------------------------------------

.epw_no_tree_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 999L, mc_iter = 100L,
                mc_seed = 1L, loo = "off")
  )
}

.epw_stump_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off")
  )
}

# Synthetic 2-leaf tree with known counts.
#   node2: class0=6, class1=2 (n_s=8); node3: class0=1, class1=7 (n_s=8)
#   Marginals: Pr(0)=7/16, Pr(1)=9/16.  Known weights:
#   ep1/0: 8*(7/16)/6 = 7/12;  ep1/1: 8*(9/16)/2 = 9/4
#   ep2/0: 8*(7/16)/1 = 7/2;   ep2/1: 8*(9/16)/7 = 9/14
.epw_2leaf_known <- function() {
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
    list(nodes = list(node1, node2, node3), no_tree = FALSE,
         overall_ess = 65.0, n_nodes = 3L, root_id = 1L,
         has_weights = FALSE, mindenom = 2L, alpha_split = 0.05,
         prune_alpha = 1.0, loo = "off"),
    class = "cta_tree"
  )
}

# Synthetic tree with one perfectly predicted endpoint (class0 absent at node3).
.epw_perfect_ep_tree <- function() {
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
    list(nodes = list(node1, node2, node3), no_tree = FALSE,
         overall_ess = 60.0, n_nodes = 3L, root_id = 1L,
         has_weights = FALSE, mindenom = 2L, alpha_split = 0.05,
         prune_alpha = 1.0, loo = "off"),
    class = "cta_tree"
  )
}

# Synthetic 3-class tree for target_class=NULL error.
.epw_3class_tree <- function() {
  leaf <- list(
    node_id = 1L, parent_id = 0L, depth = 1L, leaf = TRUE,
    majority_class = 1L, n_obs = 9L, n_weighted = 9,
    class_counts_raw      = c("0" = 3L, "1" = 3L, "2" = 3L),
    class_counts_weighted = c("0" = 3.0, "1" = 3.0, "2" = 3.0),
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    rule = NULL, ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
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
# Contract tests
# =============================================================================

test_that("epw lean: no cached artifacts on tree; no row indices on leaves; immutability", {
  tree <- .epw_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  # Top-level: no reporting cache fields
  forbidden <- c("endpoint_summary", "endpoint_counts", "staging_table",
                 "propensity_weights", "propensity_table", "count_table",
                 "staging", "weights_cache")
  expect_equal(length(intersect(names(tree), forbidden)), 0L)
  # Leaves: no row-index fields
  leaves <- Filter(function(nd) !is.null(nd) && isTRUE(nd$leaf), tree$nodes)
  for (nd in leaves) {
    expect_null(nd$idx);          expect_null(nd$row_indices)
    expect_null(nd$obs_idx);      expect_null(nd$membership)
    # class_counts_raw: correct length, integer, named, not data.frame
    expect_equal(length(nd$class_counts_raw), 2L)
    expect_true(is.integer(nd$class_counts_raw))
    expect_true(!is.null(names(nd$class_counts_raw)))
    expect_false(is.data.frame(nd$class_counts_raw))
  }
  # Immutability: endpoint_counts, staging_table, propensity_weights do not mutate tree
  t2 <- .epw_2leaf_known(); t_before <- t2
  invisible(cta_endpoint_counts(t2))
  expect_identical(t_before, t2)
  invisible(cta_staging_table(t2))
  expect_identical(t_before, t2)
  invisible(cta_propensity_weights(t2))
  expect_identical(t_before, t2)
})

test_that("epw: schema - 0 rows no-tree, 21 columns, types, no forbidden columns", {
  df <- cta_propensity_weights(.epw_no_tree_fit())
  expect_equal(nrow(df), 0L)
  expected_cols <- c(
    "endpoint_id", "endpoint_node_id", "path", "terminal_prediction",
    "class", "target_class", "class_n", "endpoint_n",
    "marginal_class_n", "marginal_total_n", "marginal_class_probability",
    "propensity_weight", "undefined_empirical", "perfectly_predicted_endpoint",
    "adjusted", "adjusted_class_n", "adjusted_endpoint_n",
    "adjusted_marginal_class_n", "adjusted_marginal_total_n",
    "adjusted_marginal_class_probability", "adjusted_propensity_weight"
  )
  expect_equal(names(df), expected_cols)
  # Types on non-empty tree
  df2 <- cta_propensity_weights(.epw_2leaf_known())
  expect_true(is.integer(df2$endpoint_id))
  expect_true(is.character(df2$class))
  expect_true(is.integer(df2$class_n))
  expect_true(is.integer(df2$endpoint_n))
  expect_true(is.numeric(df2$propensity_weight))
  expect_true(is.logical(df2$undefined_empirical))
  expect_true(is.logical(df2$adjusted))
  # No model-selection columns
  forbidden <- c("ess", "wess", "p_mc", "loo_status", "loo_ess", "loo_p")
  expect_equal(length(intersect(names(df2), forbidden)), 0L)
  # Stump: 4 rows (2 endpoints x 2 classes)
  tree <- .epw_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  expect_equal(nrow(cta_propensity_weights(tree)), 4L)
})

test_that("epw: target_class resolution and errors", {
  expect_equal(unique(cta_propensity_weights(.epw_2leaf_known())$target_class), 1L)
  expect_equal(unique(cta_propensity_weights(.epw_2leaf_known(), target_class=0L)$target_class), 0L)
  expect_error(cta_propensity_weights(.epw_3class_tree()), regexp = "target_class must be specified")
  expect_error(cta_propensity_weights(.epw_2leaf_known(), target_class=9L), regexp = "not found")
  expect_equal(length(unique(cta_propensity_weights(.epw_2leaf_known())$target_class)), 1L)
})

test_that("epw: empirical formula - 4 exact weights, marginal sum, endpoint_n, within-ep mean", {
  df <- cta_propensity_weights(.epw_2leaf_known())
  # Formula: weight = endpoint_n * marginal_prob / class_n
  expect_equal(df$propensity_weight,
               df$endpoint_n * df$marginal_class_probability / df$class_n,
               tolerance = 1e-12)
  # Known exact weights
  expect_equal(df$propensity_weight[df$endpoint_node_id == 2L & df$class == "0"],  7/12, tolerance = 1e-12)
  expect_equal(df$propensity_weight[df$endpoint_node_id == 2L & df$class == "1"],  9/4,  tolerance = 1e-12)
  expect_equal(df$propensity_weight[df$endpoint_node_id == 3L & df$class == "0"],  7/2,  tolerance = 1e-12)
  expect_equal(df$propensity_weight[df$endpoint_node_id == 3L & df$class == "1"],  9/14, tolerance = 1e-12)
  # Marginal probs sum to 1
  by_cls <- tapply(df$marginal_class_probability, df$class, unique)
  expect_equal(sum(unlist(by_cls)), 1.0, tolerance = 1e-12)
  # marginal_total_n = N = 16
  expect_equal(unique(df$marginal_total_n), 16L)
  # endpoint_n = sum of class_n within endpoint
  for (eid in unique(df$endpoint_id)) {
    rows <- df[df$endpoint_id == eid, ]
    expect_equal(unique(rows$endpoint_n), as.integer(sum(rows$class_n)))
  }
  # Within-endpoint weighted mean of propensity_weight = 1.0 (algebraic identity)
  for (eid in unique(df$endpoint_id)) {
    rows <- df[df$endpoint_id == eid, ]
    wmean <- sum(rows$propensity_weight * rows$class_n) / unique(rows$endpoint_n)
    expect_equal(wmean, 1.0, tolerance = 1e-10)
  }
})

test_that("epw: perfect endpoints - undefined_empirical, Inf, adjustment arithmetic", {
  df <- cta_propensity_weights(.epw_perfect_ep_tree())
  # undefined_empirical <-> class_n == 0; propensity_weight = Inf there
  expect_true(all(df$undefined_empirical == (df$class_n == 0L)))
  expect_true(all(is.infinite(df$propensity_weight[df$class_n == 0L])))
  # perfectly_predicted_endpoint for both rows of node3
  rows3 <- df[df$endpoint_node_id == 3L, ]
  expect_true(all(rows3$perfectly_predicted_endpoint))
  expect_false(any(df$perfectly_predicted_endpoint[df$endpoint_node_id == 2L]))
  # adjusted = TRUE <-> class_n == 0
  expect_true(all(df$adjusted == (df$class_n == 0L)))
  # adjusted_class_n = class_n + 1 for adjusted rows; unchanged for others
  adj  <- df[ df$adjusted, ]; nadj <- df[!df$adjusted, ]
  expect_equal(adj$adjusted_class_n,  adj$class_n  + 1.0)
  expect_equal(nadj$adjusted_class_n, as.numeric(nadj$class_n))
  # adjusted_endpoint_n = endpoint_n + 1 for adjusted rows
  expect_equal(adj$adjusted_endpoint_n, adj$endpoint_n + 1.0)
  # adjusted_marginal_total_n = marginal_total_n + n_adjustments
  n_adj <- sum(df$adjusted)
  expect_equal(unique(df$adjusted_marginal_total_n), unique(df$marginal_total_n) + n_adj)
  # adjusted marginal probs sum to 1
  by_cls_adj <- tapply(df$adjusted_marginal_class_probability, df$class, unique)
  expect_equal(sum(unlist(by_cls_adj)), 1.0, tolerance = 1e-12)
  # adjusted_propensity_weight finite for all rows; exact for node3/class0
  expect_true(all(is.finite(df$adjusted_propensity_weight)))
  row <- df[df$endpoint_node_id == 3L & df$class == "0", ]
  # adj_class_n=1, adj_ep_n=5, adj_Pr(0)=7/13 -> 5*(7/13)/1 = 35/13
  expect_equal(row$adjusted_propensity_weight, 35/13, tolerance = 1e-12)
  # Within-endpoint adjusted weighted mean = 1.0
  for (eid in unique(df$endpoint_id)) {
    rows  <- df[df$endpoint_id == eid, ]
    wmean <- sum(rows$adjusted_propensity_weight * rows$adjusted_class_n) /
             unique(rows$adjusted_endpoint_n)
    expect_equal(wmean, 1.0, tolerance = 1e-10)
  }
  # adjusted=FALSE: all adjusted=FALSE; Inf preserved; adjusted cols = empirical for finite rows
  df_f <- cta_propensity_weights(.epw_perfect_ep_tree(), adjusted = FALSE)
  expect_true(all(!df_f$adjusted))
  expect_true(is.infinite(df_f$adjusted_propensity_weight[df_f$endpoint_node_id == 3L & df_f$class == "0"]))
  df_n <- cta_propensity_weights(.epw_2leaf_known(), adjusted = FALSE)
  expect_equal(df_n$adjusted_class_n,    as.numeric(df_n$class_n))
  expect_equal(df_n$adjusted_endpoint_n, as.numeric(df_n$endpoint_n))
  expect_equal(df_n$adjusted_propensity_weight, df_n$propensity_weight, tolerance = 1e-12)
})

# =============================================================================
# Smoke: myeloma canon
# =============================================================================

.epw_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.epw_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.epw_myeloma_fit <- local({
  fits <- list()
  function(mindenom) {
    key <- as.character(mindenom)
    if (is.null(fits[[key]])) {
      df <- .epw_load_myeloma()
      fits[[key]] <<- suppressMessages(oda_cta_fit(
        X = df[, .epw_myeloma_attrs], y = as.integer(df$V1), w = df$V2,
        priors_on = TRUE, miss_codes = -9, alpha_split = 0.05,
        mindenom = mindenom, prune_alpha = 0.05, max_depth = 20L,
        mc_iter = 5000L, mc_target = 0.05, mc_stop = 99.9,
        mc_seed = NULL, loo = "stable", attr_names = .epw_myeloma_attrs
      ))
    }
    fits[[key]]
  }
})

test_that("epw: myeloma - rows, marginal probs, adjusted weights for MINDENOM=56/30/1", {
  skip_if_slow_tests_disabled("cta-propensity-weights")

  # MINDENOM=56: no-tree -> 0 rows, 21 columns
  df56 <- cta_propensity_weights(.epw_myeloma_fit(56L))
  expect_equal(nrow(df56), 0L)
  expect_equal(length(names(df56)), 21L)

  # MINDENOM=30: stump -> 4 rows, no NA endpoint_id, marginal probs sum to 1
  df30 <- cta_propensity_weights(.epw_myeloma_fit(30L))
  expect_equal(nrow(df30), 4L)
  expect_false(any(is.na(df30$endpoint_id)))
  by30 <- tapply(df30$marginal_class_probability, df30$class, unique)
  expect_equal(sum(unlist(by30)), 1.0, tolerance = 1e-12)
  by30a <- tapply(df30$adjusted_marginal_class_probability, df30$class, unique)
  expect_equal(sum(unlist(by30a)), 1.0, tolerance = 1e-12)

  # MINDENOM=1: 3 endpoints x 2 classes = 6 rows; adjusted weights finite;
  # within-endpoint adjusted weighted mean = 1.0
  df1 <- cta_propensity_weights(.epw_myeloma_fit(1L))
  expect_equal(nrow(df1), 6L)
  expect_true(all(is.finite(df1$adjusted_propensity_weight)))
  for (eid in unique(df1$endpoint_id)) {
    rows  <- df1[df1$endpoint_id == eid, ]
    wmean <- sum(rows$adjusted_propensity_weight * rows$adjusted_class_n) /
             unique(rows$adjusted_endpoint_n)
    expect_equal(wmean, 1.0, tolerance = 1e-10)
  }
})
