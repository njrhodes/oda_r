###############################################################################
# test-cta-assign-endpoints.R — cta_assign_endpoints()
#
# On-demand observation-to-endpoint assignment.
# Traverses the fitted cta_tree for each row of newdata; returns
# row_id, endpoint_node_id, endpoint_id.
#
# No endpoint membership is stored at fit time. This function is the
# bridge between a fitted tree and observation-level propensity weights
# (via user-side merge with cta_propensity_weights()).
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

# Stump fit with a known missing-value row.
# x=1..4 -> class 0, x=5..8 -> class 1. newdata row 3 gets x=NA.
.ae_stump_newdata_with_missing <- function() {
  data.frame(x = c(1L, 2L, NA_integer_, 6L, 7L))
}

# Synthetic 2-leaf tree for merge-sketch test.
#   node1: split on col 1 at 4.5
#   node2 (leaf): n_obs=4, majority=0
#   node3 (leaf): n_obs=4, majority=1
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
    list(nodes       = list(node1, node2, node3),
         no_tree     = FALSE,
         overall_ess = 100.0,
         n_nodes     = 3L,
         root_id     = 1L,
         has_weights = FALSE,
         mindenom    = 2L,
         alpha_split = 0.05,
         prune_alpha = 1.0,
         max_depth   = 10L,
         loo         = "off",
         miss_codes  = NULL,
         attr_names  = "x",
         n_attrs     = 1L,
         training_confusion = matrix(c(4L,0L,0L,4L), nrow=2L)),
    class = "cta_tree"
  )
}

# =============================================================================
# Fast tests — schema and types
# =============================================================================

test_that("ae: output has exactly columns row_id, endpoint_node_id, endpoint_id", {
  tree <- .ae_stump_fit()
  ep   <- cta_assign_endpoints(tree, data.frame(x = 1:4))
  expect_equal(names(ep), c("row_id", "endpoint_node_id", "endpoint_id"))
})

test_that("ae: row_id is integer 1..n", {
  tree <- .ae_stump_fit()
  ep   <- cta_assign_endpoints(tree, data.frame(x = 1:6))
  expect_equal(ep$row_id, 1:6)
  expect_type(ep$row_id, "integer")
})

test_that("ae: endpoint_node_id is integer", {
  tree <- .ae_stump_fit()
  ep   <- cta_assign_endpoints(tree, data.frame(x = 1:4))
  expect_type(ep$endpoint_node_id, "integer")
})

test_that("ae: endpoint_id is integer", {
  tree <- .ae_stump_fit()
  ep   <- cta_assign_endpoints(tree, data.frame(x = 1:4))
  expect_type(ep$endpoint_id, "integer")
})

test_that("ae: output has one row per newdata row", {
  tree <- .ae_stump_fit()
  for (n in c(1L, 5L, 8L)) {
    ep <- cta_assign_endpoints(tree, data.frame(x = seq_len(n)))
    expect_equal(nrow(ep), n, label = paste("n =", n))
  }
})

# =============================================================================
# Fast tests — no-tree behaviour
# =============================================================================

test_that("ae: no_tree fit returns nrow(newdata) rows", {
  tree <- .ae_no_tree_fit()
  ep   <- cta_assign_endpoints(tree, data.frame(x = 1:5))
  expect_equal(nrow(ep), 5L)
})

test_that("ae: no_tree fit has all NA endpoint_node_id", {
  tree <- .ae_no_tree_fit()
  ep   <- cta_assign_endpoints(tree, data.frame(x = 1:5))
  expect_true(all(is.na(ep$endpoint_node_id)))
})

test_that("ae: no_tree fit has all NA endpoint_id", {
  tree <- .ae_no_tree_fit()
  ep   <- cta_assign_endpoints(tree, data.frame(x = 1:5))
  expect_true(all(is.na(ep$endpoint_id)))
})

# =============================================================================
# Fast tests — stump routing
# =============================================================================

test_that("ae: stump fit assigns each row to one of two endpoint_ids", {
  tree  <- .ae_stump_fit()
  ep    <- cta_assign_endpoints(tree, data.frame(x = 1:8))
  es    <- cta_endpoint_summary(tree)
  expect_true(all(ep$endpoint_id %in% es$endpoint_id))
})

test_that("ae: stump endpoint_node_id values match cta_endpoint_summary()$endpoint_node_id", {
  tree <- .ae_stump_fit()
  ep   <- cta_assign_endpoints(tree, data.frame(x = 1:8))
  es   <- cta_endpoint_summary(tree)
  expect_true(all(ep$endpoint_node_id %in% es$endpoint_node_id))
})

test_that("ae: stump endpoint_id values match cta_endpoint_summary()$endpoint_id", {
  tree <- .ae_stump_fit()
  ep   <- cta_assign_endpoints(tree, data.frame(x = 1:8))
  es   <- cta_endpoint_summary(tree)
  expect_true(all(ep$endpoint_id %in% es$endpoint_id))
})

test_that("ae: stump produces two distinct endpoint_ids for separated data", {
  tree <- .ae_stump_fit()
  # x=1..4 predicts class 0, x=5..8 predicts class 1 -> should land on different endpoints
  ep <- cta_assign_endpoints(tree, data.frame(x = c(1L, 8L)))
  expect_equal(length(unique(ep$endpoint_id)), 2L)
})

# =============================================================================
# Fast tests — missingness
# =============================================================================

test_that("ae: missing_action='na' returns NA endpoint for missing root attribute", {
  tree    <- .ae_stump_fit()
  newdata <- .ae_stump_newdata_with_missing()
  ep      <- cta_assign_endpoints(tree, newdata, missing_action = "na")
  expect_true(is.na(ep$endpoint_id[3L]))
  expect_true(is.na(ep$endpoint_node_id[3L]))
})

test_that("ae: missing_action='na' non-missing rows still get valid endpoint", {
  tree    <- .ae_stump_fit()
  newdata <- .ae_stump_newdata_with_missing()
  ep      <- cta_assign_endpoints(tree, newdata, missing_action = "na")
  non_missing <- ep[!is.na(ep$endpoint_id), ]
  es          <- cta_endpoint_summary(tree)
  expect_true(all(non_missing$endpoint_id %in% es$endpoint_id))
})

test_that("ae: missing_action='majority' returns non-NA endpoint for missing root attribute", {
  tree    <- .ae_stump_fit()
  newdata <- .ae_stump_newdata_with_missing()
  ep      <- cta_assign_endpoints(tree, newdata, missing_action = "majority")
  expect_false(is.na(ep$endpoint_id[3L]))
})

test_that("ae: missing_action='majority' assigns missing row to a valid endpoint_id", {
  tree    <- .ae_stump_fit()
  newdata <- .ae_stump_newdata_with_missing()
  ep      <- cta_assign_endpoints(tree, newdata, missing_action = "majority")
  es      <- cta_endpoint_summary(tree)
  expect_true(ep$endpoint_id[3L] %in% es$endpoint_id)
})

# =============================================================================
# Fast tests — immutability
# =============================================================================

test_that("ae: cta_assign_endpoints() does not mutate the tree object", {
  tree_before <- .ae_2leaf_tree()
  tree_after  <- tree_before
  invisible(cta_assign_endpoints(tree_after, data.frame(x = 1:8)))
  expect_identical(tree_before, tree_after)
})

# =============================================================================
# Fast tests — error handling
# =============================================================================

test_that("ae: non-cta_tree input errors", {
  expect_error(cta_assign_endpoints(list(), data.frame(x = 1)),
               regexp = "inherits")
})

test_that("ae: invalid missing_action errors", {
  tree <- .ae_stump_fit()
  expect_error(cta_assign_endpoints(tree, data.frame(x = 1), missing_action = "xyz"))
})

test_that("ae: too-few newdata columns errors with informative message", {
  tree <- .ae_stump_fit()
  # stump uses attr_col=1; zero-column data.frame should error
  expect_error(cta_assign_endpoints(tree, data.frame()),
               regexp = "column")
})

# =============================================================================
# Fast tests — merge sketch: observation-level weight join works
# =============================================================================

test_that("ae: merge sketch with cta_propensity_weights yields one row per classified obs", {
  # Use the 2-leaf tree: 8 obs, all non-missing, no perfect endpoints.
  tree <- .ae_2leaf_tree()
  X    <- data.frame(x = 1:8)
  y    <- c(0L,0L,0L,0L,1L,1L,1L,1L)

  ep  <- cta_assign_endpoints(tree, X, missing_action = "na")
  pw  <- cta_propensity_weights(tree, target_class = 1L, adjusted = TRUE)

  obs_df <- data.frame(row_id = seq_len(nrow(X)),
                       class  = as.character(y),
                       stringsAsFactors = FALSE)
  joined <- merge(
    obs_df,
    merge(ep, pw[, c("endpoint_id","class","adjusted_propensity_weight")],
          by = "endpoint_id"),
    by = c("row_id","class")
  )
  # All 8 obs are classified and non-missing; expect 8 rows
  expect_equal(nrow(joined), 8L)
  expect_true(all(!is.na(joined$adjusted_propensity_weight)))
})

test_that("ae: merge sketch drops NA-endpoint rows naturally", {
  # Introduce a missing row: x=NA -> endpoint_id=NA -> drops from merge
  tree <- .ae_2leaf_tree()
  X    <- data.frame(x = c(1L, 2L, NA_integer_, 6L, 7L))
  y    <- c(0L, 0L, 0L, 1L, 1L)

  ep  <- cta_assign_endpoints(tree, X, missing_action = "na")
  pw  <- cta_propensity_weights(tree, target_class = 1L, adjusted = TRUE)

  obs_df <- data.frame(row_id = seq_len(nrow(X)),
                       class  = as.character(y),
                       stringsAsFactors = FALSE)
  joined <- merge(
    obs_df,
    merge(ep, pw[, c("endpoint_id","class","adjusted_propensity_weight")],
          by = "endpoint_id"),
    by = c("row_id","class")
  )
  # row 3 (NA) dropped; expect 4 rows
  expect_equal(nrow(joined), 4L)
  expect_true(all(!is.na(joined$adjusted_propensity_weight)))
})

# =============================================================================
# Slow fixture tests — myeloma canon
# =============================================================================

.ae_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.ae_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.ae_myeloma_fit <- function(mindenom) {
  df <- .ae_load_myeloma()
  suppressMessages(
    oda_cta_fit(
      X           = df[, .ae_myeloma_attrs],
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
      attr_names  = .ae_myeloma_attrs
    )
  )
}

test_that("ae: myeloma MINDENOM=56 (no-tree) — all endpoint_id NA", {
  skip_if_slow_tests_disabled("cta-assign-endpoints")
  df   <- .ae_load_myeloma()
  tree <- .ae_myeloma_fit(56L)
  ep   <- cta_assign_endpoints(tree, df[, .ae_myeloma_attrs], missing_action = "na")
  expect_true(all(is.na(ep$endpoint_id)))
})

test_that("ae: myeloma MINDENOM=56 (no-tree) — nrow equals nrow(newdata)", {
  skip_if_slow_tests_disabled("cta-assign-endpoints")
  df   <- .ae_load_myeloma()
  tree <- .ae_myeloma_fit(56L)
  ep   <- cta_assign_endpoints(tree, df[, .ae_myeloma_attrs])
  expect_equal(nrow(ep), nrow(df))
})

test_that("ae: myeloma MINDENOM=1 — classified endpoint count equals 255", {
  skip_if_slow_tests_disabled("cta-assign-endpoints")
  df   <- .ae_load_myeloma()
  tree <- .ae_myeloma_fit(1L)
  ep   <- cta_assign_endpoints(tree, df[, .ae_myeloma_attrs], missing_action = "na")
  classified_n <- sum(!is.na(ep$endpoint_id))
  ct_n         <- sum(cta_confusion_table(tree)$n)
  expect_equal(classified_n, ct_n)   # expected 255
})

test_that("ae: myeloma MINDENOM=30 — classified endpoint count equals 186", {
  skip_if_slow_tests_disabled("cta-assign-endpoints")
  df   <- .ae_load_myeloma()
  tree <- .ae_myeloma_fit(30L)
  ep   <- cta_assign_endpoints(tree, df[, .ae_myeloma_attrs], missing_action = "na")
  classified_n <- sum(!is.na(ep$endpoint_id))
  ct_n         <- sum(cta_confusion_table(tree)$n)
  expect_equal(classified_n, ct_n)   # expected 186
})

test_that("ae: myeloma MINDENOM=1 — merge sketch yields 255 classified rows", {
  skip_if_slow_tests_disabled("cta-assign-endpoints")
  df   <- .ae_load_myeloma()
  tree <- .ae_myeloma_fit(1L)
  X    <- df[, .ae_myeloma_attrs]
  y    <- as.integer(df$V1)

  ep  <- cta_assign_endpoints(tree, X, missing_action = "na")
  pw  <- cta_propensity_weights(tree, target_class = 1L, adjusted = TRUE)

  obs_df <- data.frame(row_id = seq_len(nrow(X)),
                       class  = as.character(y),
                       stringsAsFactors = FALSE)
  joined <- merge(
    obs_df,
    merge(ep, pw[, c("endpoint_id","class","adjusted_propensity_weight")],
          by = "endpoint_id"),
    by = c("row_id","class")
  )
  ct_n <- sum(cta_confusion_table(tree)$n)
  expect_equal(nrow(joined), ct_n)   # expected 255
})
