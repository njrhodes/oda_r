###############################################################################
# test-cta-endpoint-table.R — cta_endpoint_table() + cta_node_table()
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

# Synthetic 3-leaf tree — stable anchor for known structural values.
#   node1 (split, x<=4.5): left->node2(leaf,class0), right->node3(split)
#   node3 (split, x<=6.5): left->node4(leaf,class0), right->node5(leaf,class1)
.synthetic_3leaf_tree <- function() {
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
    n_obs = 4L, n_weighted = 4, majority_class = 0L,
    class_counts_raw      = c("0" = 4L, "1" = 0L),
    class_counts_weighted = c("0" = 4,  "1" = 0),
    rule = NULL, ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
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
    n_obs = 3L, n_weighted = 3, majority_class = 0L,
    class_counts_raw      = c("0" = 3L, "1" = 0L),
    class_counts_weighted = c("0" = 3,  "1" = 0),
    rule = NULL, ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  node5 <- list(
    node_id = 5L, parent_id = 3L, depth = 3L, leaf = TRUE,
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    n_obs = 3L, n_weighted = 3, majority_class = 1L,
    class_counts_raw      = c("0" = 0L, "1" = 3L),
    class_counts_weighted = c("0" = 0,  "1" = 3),
    rule = NULL, ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  structure(
    list(nodes = list(node1, node2, node3, node4, node5),
         no_tree = FALSE, overall_ess = 75.0, n_nodes = 5L,
         root_id = 1L, has_weights = FALSE, mindenom = 2L,
         alpha_split = 0.05, prune_alpha = 1.0, loo = "off"),
    class = "cta_tree"
  )
}

# =============================================================================
# Contract tests
# =============================================================================

test_that("eptable: schema — data.frame, zero rows for no-tree, required columns", {
  df <- cta_endpoint_table(.eptest_no_tree_fit())
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 0L)
  required <- c("endpoint_id", "leaf_node_id", "terminal_marker", "terminal",
                "depth", "parent_split_node_id", "path", "n",
                "class_counts_raw", "class_counts_weighted",
                "predicted_class", "target_n", "target_prop",
                "parent_split_attribute", "parent_split_ess",
                "parent_split_wess", "parent_split_loo_status",
                "parent_split_loo_ess", "parent_split_p_mc")
  expect_true(all(required %in% names(df)),
              info = paste("missing:", paste(setdiff(required, names(df)), collapse = ", ")))
})

test_that("eptable: stump — 2 rows, structural and count contracts", {
  tree <- .eptest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_table(tree)
  # Row count and sequential endpoint_id
  expect_equal(nrow(df), 2L)
  expect_equal(df$endpoint_id, c(1L, 2L))
  # n sums to training obs; leaf_node_ids are integers
  expect_equal(sum(df$n), 8L)
  expect_type(df$leaf_node_id, "integer")
  # Paths non-empty and contain attribute name
  expect_true(all(nzchar(df$path)))
  expect_true(all(grepl("x", df$path, fixed = TRUE)))
  # predicted_class in {0,1}; terminal flags
  expect_true(all(df$predicted_class %in% c(0L, 1L)))
  expect_true(all(df$terminal_marker == "*"))
  expect_true(all(df$terminal == TRUE))
  # class_counts_raw row sums equal n; total class-1 = 4
  for (i in seq_len(nrow(df))) {
    ccr <- df$class_counts_raw[[i]]
    if (!is.null(ccr))
      expect_equal(sum(as.integer(ccr)), df$n[i])
  }
  total_class1 <- sum(vapply(df$class_counts_raw,
                             function(ccr) if (!is.null(ccr)) as.integer(ccr["1"]) else 0L,
                             integer(1L)))
  expect_equal(total_class1, 4L)
  # target_prop = target_n / n; target_class=0 complements target_class=1
  ok <- !is.na(df$target_n) & df$n > 0L
  if (any(ok))
    expect_equal(df$target_prop[ok], df$target_n[ok] / df$n[ok])
  df0 <- cta_endpoint_table(tree, target_class = 0L)
  df1 <- cta_endpoint_table(tree, target_class = 1L)
  ok0 <- !is.na(df0$target_n); ok1 <- !is.na(df1$target_n)
  if (any(ok0) && any(ok1))
    expect_equal(df1$target_n[ok1] + df0$target_n[ok0], df1$n[ok1])
  # Parent split lineage: attribute non-NA, ESS > 0
  expect_true(all(!is.na(df$parent_split_attribute)))
  expect_true(all(df$parent_split_ess > 0))
})

test_that("eptable: synthetic 3-leaf — node IDs, paths, depths, predictions, lineage", {
  df <- cta_endpoint_table(.synthetic_3leaf_tree())
  expect_equal(nrow(df), 3L)
  expect_equal(sort(df$leaf_node_id), c(2L, 4L, 5L))
  # Paths
  expect_equal(df$path[df$leaf_node_id == 2L], "x<=4.5")
  expect_equal(df$path[df$leaf_node_id == 4L], "x>4.5 AND x<=6.5")
  expect_equal(df$path[df$leaf_node_id == 5L], "x>4.5 AND x>6.5")
  # Depths
  expect_equal(df$depth[df$leaf_node_id == 2L], 2L)
  expect_equal(df$depth[df$leaf_node_id == 4L], 3L)
  expect_equal(df$depth[df$leaf_node_id == 5L], 3L)
  # Predictions
  expect_equal(df$predicted_class[df$leaf_node_id == 2L], 0L)
  expect_equal(df$predicted_class[df$leaf_node_id == 4L], 0L)
  expect_equal(df$predicted_class[df$leaf_node_id == 5L], 1L)
  # n values
  expect_equal(df$n[df$leaf_node_id == 2L], 4L)
  expect_equal(df$n[df$leaf_node_id == 4L], 3L)
  expect_equal(df$n[df$leaf_node_id == 5L], 3L)
  # target_prop
  expect_equal(df$target_prop[df$leaf_node_id == 2L], 0.0)
  expect_equal(df$target_prop[df$leaf_node_id == 5L], 1.0)
  # Terminal markers
  expect_true(all(df$terminal_marker == "*"))
  # Parent split lineage
  expect_equal(df$parent_split_node_id[df$leaf_node_id == 2L], 1L)
  expect_equal(df$parent_split_node_id[df$leaf_node_id == 4L], 3L)
  expect_equal(df$parent_split_node_id[df$leaf_node_id == 5L], 3L)
  expect_true(all(df$parent_split_attribute == "x"))
})

test_that("nodetable: required columns, model strings, ess/p on split nodes", {
  tree <- .eptest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  nt <- cta_node_table(tree)
  required <- c("node_id", "level", "depth", "leaf", "attribute", "attr_type",
                "n_obs", "p_mc", "ess", "ess_weighted",
                "loo_status", "loo_ess", "loo_p", "model")
  expect_true(all(required %in% names(nt)),
              info = paste("missing:", paste(setdiff(required, names(nt)), collapse = ", ")))
  expect_equal(nt$level, nt$depth)
  splits <- nt[!nt$leaf, , drop = FALSE]
  leaves <- nt[ nt$leaf, , drop = FALSE]
  expect_true(all(!is.na(splits$model)))
  expect_true(all(grepl("-->", splits$model, fixed = TRUE)))
  expect_true(all(is.na(leaves$model)))
  expect_true(all(!is.na(splits$ess)))
  expect_true(all(!is.na(splits$p_mc)))
})

# =============================================================================
# Smoke: myeloma canon
# =============================================================================

.eptest_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.eptest_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.eptest_myeloma_fit <- local({
  fits <- list()
  function(mindenom) {
    key <- as.character(mindenom)
    if (is.null(fits[[key]])) {
      df <- .eptest_load_myeloma()
      fits[[key]] <<- suppressMessages(oda_cta_fit(
        X = df[, .eptest_myeloma_attrs], y = as.integer(df$V1), w = df$V2,
        priors_on = TRUE, miss_codes = -9, alpha_split = 0.05,
        mindenom = mindenom, prune_alpha = 0.05, max_depth = 20L,
        mc_iter = 5000L, mc_target = 0.05, mc_stop = 99.9,
        mc_seed = NULL, loo = "stable", attr_names = .eptest_myeloma_attrs
      ))
    }
    fits[[key]]
  }
})

test_that("eptable: myeloma — rows/paths/n for MINDENOM=1/30/56", {
  skip_if_slow_tests_disabled("cta-endpoint-table")

  # MINDENOM=1: 3 endpoints, V14 in all paths, n=255
  df1 <- cta_endpoint_table(.eptest_myeloma_fit(1L))
  expect_equal(nrow(df1), 3L)
  expect_true(all(grepl("V14", df1$path, fixed = TRUE)))
  expect_equal(sum(df1$n), 255L)

  # MINDENOM=30: 2 endpoints (stump), V17 in all paths
  df30 <- cta_endpoint_table(.eptest_myeloma_fit(30L))
  expect_equal(nrow(df30), 2L)
  expect_true(all(grepl("V17", df30$path, fixed = TRUE)))

  # MINDENOM=56: no-tree → zero rows with 'path' column
  df56 <- cta_endpoint_table(.eptest_myeloma_fit(56L))
  expect_equal(nrow(df56), 0L)
  expect_true("path" %in% names(df56))
})
