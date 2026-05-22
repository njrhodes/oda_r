###############################################################################
# test-cta-propensity-weights.R — cta_propensity_weights()
#
# Endpoint × class stabilized propensity-style weights.
# Consumes cta_endpoint_counts(tree) only — no refitting, no prediction,
# no mutation of the tree object.
#
# Also contains lean-fit invariant tests: verifies that the cta_tree and
# its leaf nodes do not cache any reporting artifacts or training row
# indices (class_counts_raw/weighted are the only per-leaf additions).
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

# Synthetic 2-leaf tree: known exact counts for arithmetic verification.
#   node2: class0=6, class1=2  (endpoint 1, n_s=8)
#   node3: class0=1, class1=7  (endpoint 2, n_s=8)
# Marginals: class0 total=7, class1 total=9, N=16
# Pr(class0)=7/16, Pr(class1)=9/16
# Empirical weights:
#   ep1/class0: 8*(7/16)/6 = 7/12
#   ep1/class1: 8*(9/16)/2 = 9/4
#   ep2/class0: 8*(7/16)/1 = 7/2
#   ep2/class1: 8*(9/16)/7 = 9/14
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

# Synthetic 2-leaf tree: one perfectly predicted endpoint.
#   node2: class0=6, class1=2  (not perfect)
#   node3: class0=0, class1=4  (perfectly predicted — class0 absent)
# Empirical: N=12, Pr(0)=6/12=0.5, Pr(1)=6/12=0.5
# Adjusted (add 1 to class0 at node3): adj_N=13, adj_Pr(0)=7/13, adj_Pr(1)=6/13
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

# Synthetic 3-class tree (for target_class=NULL error).
.epw_3class_tree <- function() {
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
# Lean-fit invariant tests
# =============================================================================

test_that("epw lean: cta_tree top-level names do not include any cached reporting artifact", {
  tree <- .epw_stump_fit()
  forbidden <- c("endpoint_summary", "endpoint_counts", "staging_table",
                 "propensity_weights", "propensity_table", "count_table",
                 "staging", "weights_cache")
  expect_true(length(intersect(names(tree), forbidden)) == 0L,
              info = paste("found:", paste(intersect(names(tree), forbidden), collapse=", ")))
})

test_that("epw lean: terminal leaf nodes do not store row indices", {
  tree  <- .epw_stump_fit()
  leaves <- Filter(function(nd) !is.null(nd) && isTRUE(nd$leaf), tree$nodes)
  for (nd in leaves) {
    expect_null(nd$idx,         label = "leaf$idx must be NULL")
    expect_null(nd$row_indices, label = "leaf$row_indices must be NULL")
    expect_null(nd$obs_idx,     label = "leaf$obs_idx must be NULL")
    expect_null(nd$membership,  label = "leaf$membership must be NULL")
  }
})

test_that("epw lean: cta_endpoint_counts() does not mutate the tree object", {
  tree_before <- .epw_2leaf_known()
  tree_after  <- tree_before
  invisible(cta_endpoint_counts(tree_after))
  expect_identical(tree_before, tree_after)
})

test_that("epw lean: cta_staging_table() does not mutate the tree object", {
  tree_before <- .epw_2leaf_known()
  tree_after  <- tree_before
  invisible(cta_staging_table(tree_after))
  expect_identical(tree_before, tree_after)
})

test_that("epw lean: cta_propensity_weights() does not mutate the tree object", {
  tree_before <- .epw_2leaf_known()
  tree_after  <- tree_before
  invisible(cta_propensity_weights(tree_after))
  expect_identical(tree_before, tree_after)
})

test_that("epw lean: leaf class_counts_raw length equals number of classes", {
  tree   <- .epw_stump_fit()
  leaves <- Filter(function(nd) !is.null(nd) && isTRUE(nd$leaf), tree$nodes)
  for (nd in leaves) {
    expect_equal(length(nd$class_counts_raw), 2L)     # binary fit
    expect_equal(length(nd$class_counts_weighted), 2L)
  }
})

test_that("epw lean: leaf class_counts_raw contains only named integer counts, no data frames", {
  tree   <- .epw_2leaf_known()
  leaves <- Filter(function(nd) !is.null(nd) && isTRUE(nd$leaf), tree$nodes)
  for (nd in leaves) {
    expect_true(is.integer(nd$class_counts_raw))
    expect_true(is.numeric(nd$class_counts_weighted))
    expect_true(!is.null(names(nd$class_counts_raw)))
    expect_false(is.data.frame(nd$class_counts_raw))
  }
})

# =============================================================================
# Fast tests — schema and structure
# =============================================================================

test_that("epw: no-tree fit returns 0 rows with exact 21-column schema", {
  df <- cta_propensity_weights(.epw_no_tree_fit())
  expect_equal(nrow(df), 0L)
  expected <- c(
    "endpoint_id", "endpoint_node_id", "path", "terminal_prediction",
    "class", "target_class", "class_n", "endpoint_n",
    "marginal_class_n", "marginal_total_n", "marginal_class_probability",
    "propensity_weight", "undefined_empirical", "perfectly_predicted_endpoint",
    "adjusted", "adjusted_class_n", "adjusted_endpoint_n",
    "adjusted_marginal_class_n", "adjusted_marginal_total_n",
    "adjusted_marginal_class_probability", "adjusted_propensity_weight"
  )
  expect_equal(names(df), expected)
})

test_that("epw: column types are correct", {
  df <- cta_propensity_weights(.epw_2leaf_known())
  expect_true(is.integer(df$endpoint_id))
  expect_true(is.integer(df$endpoint_node_id))
  expect_true(is.character(df$path))
  expect_true(is.integer(df$terminal_prediction))
  expect_true(is.character(df$class))
  expect_true(is.integer(df$target_class))
  expect_true(is.integer(df$class_n))
  expect_true(is.integer(df$endpoint_n))
  expect_true(is.integer(df$marginal_class_n))
  expect_true(is.integer(df$marginal_total_n))
  expect_true(is.numeric(df$marginal_class_probability))
  expect_true(is.numeric(df$propensity_weight))
  expect_true(is.logical(df$undefined_empirical))
  expect_true(is.logical(df$perfectly_predicted_endpoint))
  expect_true(is.logical(df$adjusted))
  expect_true(is.numeric(df$adjusted_class_n))
  expect_true(is.numeric(df$adjusted_endpoint_n))
  expect_true(is.numeric(df$adjusted_marginal_class_n))
  expect_true(is.numeric(df$adjusted_marginal_total_n))
  expect_true(is.numeric(df$adjusted_marginal_class_probability))
  expect_true(is.numeric(df$adjusted_propensity_weight))
})

test_that("epw: no forbidden model-selection columns", {
  df <- cta_propensity_weights(.epw_2leaf_known())
  forbidden <- c("ess", "wess", "p_mc", "loo_status", "loo_ess",
                 "loo_p", "alpha_split", "prune_alpha")
  expect_equal(length(intersect(names(df), forbidden)), 0L)
})

test_that("epw: stump returns 4 rows (2 endpoints x 2 classes)", {
  tree <- .epw_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed split")
  df <- cta_propensity_weights(tree)
  expect_equal(nrow(df), 4L)
})

# =============================================================================
# Fast tests — target_class resolution
# =============================================================================

test_that("epw: NULL target_class defaults to 1 for binary tree", {
  df <- cta_propensity_weights(.epw_2leaf_known())
  expect_equal(unique(df$target_class), 1L)
})

test_that("epw: explicit target_class=0 uses class 0", {
  df <- cta_propensity_weights(.epw_2leaf_known(), target_class = 0L)
  expect_equal(unique(df$target_class), 0L)
})

test_that("epw: NULL target_class on 3-class tree stops with error", {
  expect_error(cta_propensity_weights(.epw_3class_tree()),
               regexp = "target_class must be specified")
})

test_that("epw: invalid target_class stops with informative error", {
  expect_error(cta_propensity_weights(.epw_2leaf_known(), target_class = 9L),
               regexp = "not found in tree classes")
})

test_that("epw: all rows have same target_class value", {
  df <- cta_propensity_weights(.epw_2leaf_known())
  expect_equal(length(unique(df$target_class)), 1L)
})

# =============================================================================
# Fast tests — empirical formula
# =============================================================================

test_that("epw: empirical formula: weight = endpoint_n * marginal_prob / class_n", {
  df <- cta_propensity_weights(.epw_2leaf_known())
  # No perfect endpoints, so all rows have class_n > 0.
  expected_w <- df$endpoint_n * df$marginal_class_probability / df$class_n
  expect_equal(df$propensity_weight, expected_w, tolerance = 1e-12)
})

test_that("epw: known exact weight for 2-leaf tree (ep1/class0)", {
  df  <- cta_propensity_weights(.epw_2leaf_known())
  row <- df[df$endpoint_node_id == 2L & df$class == "0", ]
  # ep_n=8, marginal_prob=7/16, class_n=6 → 8*(7/16)/6 = 7/12
  expect_equal(row$propensity_weight, 7/12, tolerance = 1e-12)
})

test_that("epw: known exact weight for 2-leaf tree (ep1/class1)", {
  df  <- cta_propensity_weights(.epw_2leaf_known())
  row <- df[df$endpoint_node_id == 2L & df$class == "1", ]
  # ep_n=8, marginal_prob=9/16, class_n=2 → 8*(9/16)/2 = 9/4
  expect_equal(row$propensity_weight, 9/4, tolerance = 1e-12)
})

test_that("epw: known exact weight for 2-leaf tree (ep2/class0)", {
  df  <- cta_propensity_weights(.epw_2leaf_known())
  row <- df[df$endpoint_node_id == 3L & df$class == "0", ]
  # ep_n=8, marginal_prob=7/16, class_n=1 → 8*(7/16)/1 = 7/2
  expect_equal(row$propensity_weight, 7/2, tolerance = 1e-12)
})

test_that("epw: known exact weight for 2-leaf tree (ep2/class1)", {
  df  <- cta_propensity_weights(.epw_2leaf_known())
  row <- df[df$endpoint_node_id == 3L & df$class == "1", ]
  # ep_n=8, marginal_prob=9/16, class_n=7 → 8*(9/16)/7 = 9/14
  expect_equal(row$propensity_weight, 9/14, tolerance = 1e-12)
})

test_that("epw: marginal_class_probability sums to 1 across unique classes", {
  df      <- cta_propensity_weights(.epw_2leaf_known())
  # One marginal_class_probability per class (same value across rows for that class).
  by_cls  <- tapply(df$marginal_class_probability, df$class, unique)
  expect_equal(sum(unlist(by_cls)), 1.0, tolerance = 1e-12)
})

test_that("epw: marginal_total_n equals total n_raw across all endpoints", {
  df <- cta_propensity_weights(.epw_2leaf_known())
  # marginal_total_n = sum over all endpoint × class rows of class_n,
  # counting each endpoint only once per class: for a binary tree with
  # 2 endpoints and 2 classes, sum(class_n) = 6+2+1+7 = 16.
  expect_equal(unique(df$marginal_total_n), 16L)
})

test_that("epw: endpoint_n equals sum of class_n within each endpoint", {
  df <- cta_propensity_weights(.epw_2leaf_known())
  for (eid in unique(df$endpoint_id)) {
    rows <- df[df$endpoint_id == eid, ]
    expect_equal(unique(rows$endpoint_n), as.integer(sum(rows$class_n)))
  }
})

# Algebraic identity: within-endpoint weighted mean of propensity_weight
# (weighted by class_n) equals 1.0 exactly for any endpoint with no absent class.
test_that("epw: within-endpoint weighted mean of propensity_weight = 1.0", {
  df <- cta_propensity_weights(.epw_2leaf_known())
  for (eid in unique(df$endpoint_id)) {
    rows  <- df[df$endpoint_id == eid, ]
    wmean <- sum(rows$propensity_weight * rows$class_n) / unique(rows$endpoint_n)
    expect_equal(wmean, 1.0, tolerance = 1e-10)
  }
})

# =============================================================================
# Fast tests — perfect endpoint and undefined_empirical
# =============================================================================

test_that("epw: undefined_empirical TRUE where class_n == 0", {
  df <- cta_propensity_weights(.epw_perfect_ep_tree())
  expect_true(all(df$undefined_empirical == (df$class_n == 0L)))
})

test_that("epw: propensity_weight is Inf where class_n == 0", {
  df <- cta_propensity_weights(.epw_perfect_ep_tree())
  expect_true(all(is.infinite(df$propensity_weight[df$class_n == 0L])))
})

test_that("epw: perfectly_predicted_endpoint TRUE for all rows of a pure endpoint", {
  df  <- cta_propensity_weights(.epw_perfect_ep_tree())
  # node3 (endpoint 2) has class0=0 → perfectly_predicted_endpoint for both rows
  rows3 <- df[df$endpoint_node_id == 3L, ]
  expect_true(all(rows3$perfectly_predicted_endpoint))
})

test_that("epw: perfectly_predicted_endpoint FALSE for all rows of a mixed endpoint", {
  df    <- cta_propensity_weights(.epw_perfect_ep_tree())
  rows2 <- df[df$endpoint_node_id == 2L, ]
  expect_false(any(rows2$perfectly_predicted_endpoint))
})

# =============================================================================
# Fast tests — adjustment (adjusted=TRUE default)
# =============================================================================

test_that("epw: adjusted TRUE only for absent-class rows, FALSE for present-class rows", {
  df <- cta_propensity_weights(.epw_perfect_ep_tree())
  expect_true(all(df$adjusted == (df$class_n == 0L)))
})

test_that("epw: adjusted_class_n = class_n + 1 for adjusted rows", {
  df <- cta_propensity_weights(.epw_perfect_ep_tree())
  adj_rows  <- df[df$adjusted, ]
  nadj_rows <- df[!df$adjusted, ]
  expect_equal(adj_rows$adjusted_class_n,  adj_rows$class_n  + 1.0)
  expect_equal(nadj_rows$adjusted_class_n, as.numeric(nadj_rows$class_n))
})

test_that("epw: adjusted_endpoint_n = endpoint_n + 1 for adjusted rows", {
  df <- cta_propensity_weights(.epw_perfect_ep_tree())
  adj_rows <- df[df$adjusted, ]
  expect_equal(adj_rows$adjusted_endpoint_n, adj_rows$endpoint_n + 1.0)
})

test_that("epw: adjusted_marginal_total_n = marginal_total_n + n_adjustments", {
  df  <- cta_propensity_weights(.epw_perfect_ep_tree())
  n_adj <- sum(df$adjusted)   # 1 row adjusted (class0 at node3)
  expect_equal(unique(df$adjusted_marginal_total_n),
               unique(df$marginal_total_n) + n_adj)
})

test_that("epw: adjusted marginal class probabilities sum to 1", {
  df     <- cta_propensity_weights(.epw_perfect_ep_tree())
  by_cls <- tapply(df$adjusted_marginal_class_probability, df$class, unique)
  expect_equal(sum(unlist(by_cls)), 1.0, tolerance = 1e-12)
})

test_that("epw: adjusted_propensity_weight finite for all rows when adjusted=TRUE", {
  df <- cta_propensity_weights(.epw_perfect_ep_tree(), adjusted = TRUE)
  expect_true(all(is.finite(df$adjusted_propensity_weight)))
})

test_that("epw: known exact adjusted weight for perfect endpoint (class0 at node3)", {
  df  <- cta_propensity_weights(.epw_perfect_ep_tree())
  row <- df[df$endpoint_node_id == 3L & df$class == "0", ]
  # adj_class_n=1, adj_ep_n=5, adj_Pr(0)=7/13 → 5*(7/13)/1 = 35/13
  expect_equal(row$adjusted_propensity_weight, 35/13, tolerance = 1e-12)
})

test_that("epw: within-endpoint weighted mean of adjusted_propensity_weight = 1.0", {
  df <- cta_propensity_weights(.epw_perfect_ep_tree())
  for (eid in unique(df$endpoint_id)) {
    rows  <- df[df$endpoint_id == eid, ]
    wmean <- sum(rows$adjusted_propensity_weight * rows$adjusted_class_n) /
             unique(rows$adjusted_endpoint_n)
    expect_equal(wmean, 1.0, tolerance = 1e-10,
                 label = paste("endpoint", eid, "adjusted weighted mean"))
  }
})

# =============================================================================
# Fast tests — adjusted=FALSE behavior
# =============================================================================

test_that("epw: adjusted=FALSE sets all adjusted column to FALSE", {
  df <- cta_propensity_weights(.epw_perfect_ep_tree(), adjusted = FALSE)
  expect_true(all(!df$adjusted))
})

test_that("epw: adjusted=FALSE keeps Inf in adjusted_propensity_weight for absent class", {
  df  <- cta_propensity_weights(.epw_perfect_ep_tree(), adjusted = FALSE)
  row <- df[df$endpoint_node_id == 3L & df$class == "0", ]
  expect_true(is.infinite(row$adjusted_propensity_weight))
})

test_that("epw: adjusted=FALSE adjusted columns equal empirical columns for finite rows", {
  df <- cta_propensity_weights(.epw_2leaf_known(), adjusted = FALSE)
  # No perfect endpoints: all rows finite; adjusted should equal empirical.
  expect_equal(df$adjusted_class_n,    as.numeric(df$class_n))
  expect_equal(df$adjusted_endpoint_n, as.numeric(df$endpoint_n))
  expect_equal(df$adjusted_marginal_class_probability,
               df$marginal_class_probability,
               tolerance = 1e-12)
  expect_equal(df$adjusted_propensity_weight, df$propensity_weight, tolerance = 1e-12)
})

# =============================================================================
# Slow fixture tests — myeloma canon
# =============================================================================

.epw_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.epw_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.epw_myeloma_fit <- function(mindenom) {
  df <- .epw_load_myeloma()
  suppressMessages(
    oda_cta_fit(
      X           = df[, .epw_myeloma_attrs],
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
      attr_names  = .epw_myeloma_attrs
    )
  )
}

test_that("epw: myeloma MINDENOM=56 (no-tree) returns 0 rows with 21 columns", {
  skip_if_slow_tests_disabled("cta-propensity-weights")
  df <- cta_propensity_weights(.epw_myeloma_fit(56L))
  expect_equal(nrow(df), 0L)
  expect_equal(length(names(df)), 21L)
})

test_that("epw: myeloma MINDENOM=30 (stump) returns 4 rows", {
  skip_if_slow_tests_disabled("cta-propensity-weights")
  df <- cta_propensity_weights(.epw_myeloma_fit(30L))
  expect_equal(nrow(df), 4L)
})

test_that("epw: myeloma MINDENOM=1 returns 6 rows (3 endpoints x 2 classes)", {
  skip_if_slow_tests_disabled("cta-propensity-weights")
  df <- cta_propensity_weights(.epw_myeloma_fit(1L))
  expect_equal(nrow(df), 6L)
})

test_that("epw: myeloma MINDENOM=30 empirical marginal probabilities sum to 1", {
  skip_if_slow_tests_disabled("cta-propensity-weights")
  df     <- cta_propensity_weights(.epw_myeloma_fit(30L))
  by_cls <- tapply(df$marginal_class_probability, df$class, unique)
  expect_equal(sum(unlist(by_cls)), 1.0, tolerance = 1e-12)
})

test_that("epw: myeloma MINDENOM=30 adjusted marginal probabilities sum to 1", {
  skip_if_slow_tests_disabled("cta-propensity-weights")
  df     <- cta_propensity_weights(.epw_myeloma_fit(30L))
  by_cls <- tapply(df$adjusted_marginal_class_probability, df$class, unique)
  expect_equal(sum(unlist(by_cls)), 1.0, tolerance = 1e-12)
})

test_that("epw: myeloma MINDENOM=30 no rows dropped", {
  skip_if_slow_tests_disabled("cta-propensity-weights")
  df <- cta_propensity_weights(.epw_myeloma_fit(30L))
  expect_equal(nrow(df), 4L)
  expect_false(any(is.na(df$endpoint_id)))
})

test_that("epw: myeloma MINDENOM=1 adjusted weighted mean per endpoint = 1.0", {
  skip_if_slow_tests_disabled("cta-propensity-weights")
  df <- cta_propensity_weights(.epw_myeloma_fit(1L))
  for (eid in unique(df$endpoint_id)) {
    rows  <- df[df$endpoint_id == eid, ]
    wmean <- sum(rows$adjusted_propensity_weight * rows$adjusted_class_n) /
             unique(rows$adjusted_endpoint_n)
    expect_equal(wmean, 1.0, tolerance = 1e-10,
                 label = paste("myeloma ep", eid, "adjusted weighted mean"))
  }
})

test_that("epw: myeloma MINDENOM=1 adjusted_propensity_weight all finite", {
  skip_if_slow_tests_disabled("cta-propensity-weights")
  df <- cta_propensity_weights(.epw_myeloma_fit(1L), adjusted = TRUE)
  expect_true(all(is.finite(df$adjusted_propensity_weight)))
})
