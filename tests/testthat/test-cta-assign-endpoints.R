###############################################################################
# test-cta-assign-endpoints.R - cta_assign_endpoints()
#
# On-demand observation-to-endpoint assignment.
# Traverses the fitted cta_tree for each row of newdata; returns
# row_id, endpoint_node_id, endpoint_id.
###############################################################################

# ---- Helpers ----------------------------------------------------------------

.ae_no_tree_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 999L, mc_iter = 100L,
                mc_seed = 1L, loo = "off")
  )
}

.ae_stump_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off")
  )
}

# Synthetic 2-leaf tree for merge-sketch tests.
.ae_2leaf_tree <- function() {
  node1 <- list(
    node_id = 1L, parent_id = 0L, depth = 1L, leaf = FALSE,
    attribute = "x", attr_col = 1L, attr_type = "ordered",
    n_obs = 8L, n_weighted = 8.0,
    rule = list(type = "ordered_cut", cut_value = 4.5, direction = "0->1"),
    ess = 60.0, ess_weighted = NA_real_, p_mc = 0.01,
    loo_status = "STABLE", loo_ess = 58.0, loo_p = NA_real_,
    confusion = NULL, split_labels = c(0L, 1L), child_ids = c(2L, 3L),
    majority_class = 0L
  )
  node2 <- list(
    node_id = 2L, parent_id = 1L, depth = 2L, leaf = TRUE,
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    n_obs = 4L, n_weighted = 4.0, majority_class = 0L, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    class_counts_raw      = c("0" = 4L, "1" = 0L),
    class_counts_weighted = c("0" = 4.0, "1" = 0.0),
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  node3 <- list(
    node_id = 3L, parent_id = 1L, depth = 2L, leaf = TRUE,
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    n_obs = 4L, n_weighted = 4.0, majority_class = 1L, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    class_counts_raw      = c("0" = 0L, "1" = 4L),
    class_counts_weighted = c("0" = 0.0, "1" = 4.0),
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  structure(
    list(nodes = list(node1, node2, node3),
         no_tree = FALSE, overall_ess = 100.0, n_nodes = 3L,
         root_id = 1L, has_weights = FALSE, mindenom = 2L,
         alpha_split = 0.05, prune_alpha = 1.0, max_depth = 10L,
         loo = "off", miss_codes = NULL, attr_names = "x", n_attrs = 1L,
         training_confusion = matrix(c(4L,0L,0L,4L), nrow=2L)),
    class = "cta_tree"
  )
}

# =============================================================================
# Contract tests
# =============================================================================

test_that("ae: schema - columns/types/row-count; no-tree returns all-NA endpoints", {
  # Schema on stump
  tree <- .ae_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  ep <- cta_assign_endpoints(tree, data.frame(x = 1:6))
  expect_equal(names(ep), c("row_id", "endpoint_node_id", "endpoint_id"))
  expect_type(ep$row_id,           "integer")
  expect_type(ep$endpoint_node_id, "integer")
  expect_type(ep$endpoint_id,      "integer")
  expect_equal(ep$row_id, 1:6)
  expect_equal(nrow(ep), 6L)
  # No-tree: correct nrow and all-NA endpoints
  nt <- .ae_no_tree_fit()
  ep_nt <- cta_assign_endpoints(nt, data.frame(x = 1:5))
  expect_equal(nrow(ep_nt), 5L)
  expect_true(all(is.na(ep_nt$endpoint_node_id)))
  expect_true(all(is.na(ep_nt$endpoint_id)))
})

test_that("ae: stump routing, missingness, immutability, and merge sketch", {
  tree <- .ae_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")

  # Routing: all endpoint_ids/node_ids match cta_endpoint_summary
  ep <- cta_assign_endpoints(tree, data.frame(x = 1:8))
  es <- cta_endpoint_summary(tree)
  expect_true(all(ep$endpoint_id      %in% es$endpoint_id))
  expect_true(all(ep$endpoint_node_id %in% es$endpoint_node_id))
  # Two distinct endpoints for perfectly separated data
  ep2 <- cta_assign_endpoints(tree, data.frame(x = c(1L, 8L)))
  expect_equal(length(unique(ep2$endpoint_id)), 2L)

  # Missingness: missing_action='na' -> NA for missing row; others still valid
  nd <- data.frame(x = c(1L, 2L, NA_integer_, 6L, 7L))
  ep_na <- cta_assign_endpoints(tree, nd, missing_action = "na")
  expect_true(is.na(ep_na$endpoint_id[3L]))
  expect_true(is.na(ep_na$endpoint_node_id[3L]))
  expect_true(all(ep_na$endpoint_id[-3L] %in% es$endpoint_id))
  # missing_action='majority' -> non-NA endpoint for missing row
  ep_maj <- cta_assign_endpoints(tree, nd, missing_action = "majority")
  expect_false(is.na(ep_maj$endpoint_id[3L]))
  expect_true(ep_maj$endpoint_id[3L] %in% es$endpoint_id)

  # Immutability: tree not modified by call
  tree2 <- .ae_2leaf_tree()
  before <- tree2
  invisible(cta_assign_endpoints(tree2, data.frame(x = 1:8)))
  expect_identical(before, tree2)

  # Merge sketch: join with propensity weights - 8 rows, all non-NA weights
  X <- data.frame(x = 1:8)
  y <- c(0L,0L,0L,0L,1L,1L,1L,1L)
  ep_m <- cta_assign_endpoints(.ae_2leaf_tree(), X, missing_action = "na")
  pw   <- cta_propensity_weights(.ae_2leaf_tree(), target_class = 1L, adjusted = TRUE)
  obs_df <- data.frame(row_id = seq_len(nrow(X)), class = as.character(y),
                       stringsAsFactors = FALSE)
  joined <- merge(obs_df,
                  merge(ep_m, pw[, c("endpoint_id","class","adjusted_propensity_weight")],
                        by = "endpoint_id"),
                  by = c("row_id","class"))
  expect_equal(nrow(joined), 8L)
  expect_true(all(!is.na(joined$adjusted_propensity_weight)))
  # NA row drops from merge naturally
  X_miss <- data.frame(x = c(1L, 2L, NA_integer_, 6L, 7L))
  ep_miss <- cta_assign_endpoints(.ae_2leaf_tree(), X_miss, missing_action = "na")
  obs2 <- data.frame(row_id = seq_len(nrow(X_miss)),
                     class = as.character(c(0L,0L,0L,1L,1L)),
                     stringsAsFactors = FALSE)
  joined2 <- merge(obs2,
                   merge(ep_miss, pw[, c("endpoint_id","class","adjusted_propensity_weight")],
                         by = "endpoint_id"),
                   by = c("row_id","class"))
  expect_equal(nrow(joined2), 4L)
})

test_that("ae: errors - non-cta_tree, bad missing_action, zero-column data", {
  tree <- .ae_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  expect_error(cta_assign_endpoints(list(), data.frame(x = 1)), regexp = "inherits")
  expect_error(cta_assign_endpoints(tree, data.frame(x = 1), missing_action = "xyz"))
  expect_error(cta_assign_endpoints(tree, data.frame()), regexp = "column")
})

# =============================================================================
# Smoke: myeloma canon
# =============================================================================

.ae_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.ae_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.ae_myeloma_fit <- local({
  fits <- list()
  function(mindenom) {
    key <- as.character(mindenom)
    if (is.null(fits[[key]])) {
      df <- .ae_load_myeloma()
      fits[[key]] <<- suppressMessages(oda_cta_fit(
        X = df[, .ae_myeloma_attrs], y = as.integer(df$V1), w = df$V2,
        priors_on = TRUE, miss_codes = -9, alpha_split = 0.05,
        mindenom = mindenom, prune_alpha = 0.05, max_depth = 20L,
        mc_iter = 5000L, mc_target = 0.05, mc_stop = 99.9,
        mc_seed = NULL, loo = "stable", attr_names = .ae_myeloma_attrs
      ))
    }
    fits[[key]]
  }
})

test_that("ae: myeloma - classified counts and merge sketch for MINDENOM=56/1/30", {
  skip_if_slow_tests_disabled("cta-assign-endpoints")
  df <- .ae_load_myeloma()
  X  <- df[, .ae_myeloma_attrs]
  y  <- as.integer(df$V1)

  # MINDENOM=56: no-tree -> all NA
  t56 <- .ae_myeloma_fit(56L)
  ep56 <- cta_assign_endpoints(t56, X, missing_action = "na")
  expect_equal(nrow(ep56), nrow(df))
  expect_true(all(is.na(ep56$endpoint_id)))

  # MINDENOM=1: classified = confusion total (255)
  t1 <- .ae_myeloma_fit(1L)
  ep1 <- cta_assign_endpoints(t1, X, missing_action = "na")
  expect_equal(sum(!is.na(ep1$endpoint_id)), sum(cta_confusion_table(t1)$n))

  # MINDENOM=30: classified = confusion total (186)
  t30 <- .ae_myeloma_fit(30L)
  ep30 <- cta_assign_endpoints(t30, X, missing_action = "na")
  expect_equal(sum(!is.na(ep30$endpoint_id)), sum(cta_confusion_table(t30)$n))

  # MINDENOM=1 merge sketch: classified rows with weights
  pw <- cta_propensity_weights(t1, target_class = 1L, adjusted = TRUE)
  obs_df <- data.frame(row_id = seq_len(nrow(X)), class = as.character(y),
                       stringsAsFactors = FALSE)
  joined <- merge(obs_df,
                  merge(ep1, pw[, c("endpoint_id","class","adjusted_propensity_weight")],
                        by = "endpoint_id"),
                  by = c("row_id","class"))
  expect_equal(nrow(joined), sum(cta_confusion_table(t1)$n))
})
