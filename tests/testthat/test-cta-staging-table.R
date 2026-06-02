###############################################################################
# test-cta-staging-table.R — cta_staging_table()
#
# Per-endpoint staging table ordered by ascending adjusted target-class
# propensity. Consumes cta_endpoint_summary() + cta_endpoint_counts() only.
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

# Synthetic 2-leaf tree with known counts for exact arithmetic.
#   node2 (x<=6.5): class0=6, class1=2 → target_n=2, denom=8, prop=0.25
#   node3 (x>6.5):  class0=1, class1=7 → target_n=7, denom=8, prop=0.875
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
    list(nodes = list(node1, node2, node3), no_tree = FALSE,
         overall_ess = 65.0, n_nodes = 3L, root_id = 1L,
         has_weights = FALSE, mindenom = 2L, alpha_split = 0.05,
         prune_alpha = 1.0, loo = "off"),
    class = "cta_tree"
  )
}

# Synthetic tree with one perfect endpoint: node3 has class0=0, class1=4.
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
    list(nodes = list(node1, node2, node3), no_tree = FALSE,
         overall_ess = 60.0, n_nodes = 3L, root_id = 1L,
         has_weights = FALSE, mindenom = 2L, alpha_split = 0.05,
         prune_alpha = 1.0, loo = "off"),
    class = "cta_tree"
  )
}

# Synthetic tree with weighted counts differing from raw.
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
    list(nodes = list(node1, node2, node3), no_tree = FALSE,
         overall_ess = 60.0, n_nodes = 3L, root_id = 1L,
         has_weights = TRUE, mindenom = 2L, alpha_split = 0.05,
         prune_alpha = 1.0, loo = "off"),
    class = "cta_tree"
  )
}

# Synthetic 3-class tree for target_class=NULL error check.
.est_3class_tree <- function() {
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

test_that("est: schema — zero rows, 21 columns, correct types for no-tree", {
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
  # Column types on a non-empty tree
  df2 <- cta_staging_table(.est_2leaf_known())
  expect_true(is.integer(df2$stage))
  expect_true(is.integer(df2$endpoint_node_id))
  expect_true(is.integer(df2$target_class))
  expect_true(is.logical(df2$perfectly_predicted))
  expect_true(is.logical(df2$adjusted))
  expect_true(is.logical(df2$weighted))
  expect_true(is.numeric(df2$target_proportion))
  expect_true(is.numeric(df2$odds))
})

test_that("est: arithmetic + ordering — known 2-leaf tree, target_class resolution/errors", {
  df <- cta_staging_table(.est_2leaf_known())
  expect_equal(nrow(df), 2L)
  expect_equal(df$stage, c(1L, 2L))
  # stage 1 = node2 (lower prop 0.25), stage 2 = node3 (higher prop 0.875)
  expect_equal(df$endpoint_node_id[1L], 2L)
  expect_equal(df$endpoint_node_id[2L], 3L)
  # target_n and denominator
  expect_equal(df$target_n[df$endpoint_node_id == 2L], 2.0)
  expect_equal(df$denominator[df$endpoint_node_id == 2L], 8.0)
  expect_equal(df$target_n[df$endpoint_node_id == 3L], 7.0)
  expect_equal(df$denominator[df$endpoint_node_id == 3L], 8.0)
  # Derived: proportion, non_target_n, odds
  expect_equal(df$target_proportion, df$target_n / df$denominator, tolerance = 1e-12)
  expect_equal(df$non_target_n, df$denominator - df$target_n)
  expect_equal(df$odds, df$target_n / df$non_target_n, tolerance = 1e-12)
  # Stage ordering by ascending adjusted_target_proportion
  expect_true(all(diff(df$adjusted_target_proportion) >= 0))
  # NULL target_class → binary default = 1
  expect_equal(unique(df$target_class), 1L)
  # Explicit target_class=0 uses class-0 counts
  df0 <- cta_staging_table(.est_2leaf_known(), target_class = 0L)
  expect_equal(unique(df0$target_class), 0L)
  expect_equal(sort(df0$target_n), c(1.0, 6.0))
  # Errors
  expect_error(cta_staging_table(.est_3class_tree()),   regexp = "target_class must be specified")
  expect_error(cta_staging_table(.est_2leaf_known(), target_class = 9L), regexp = "not found")
  # n_obs/n_weighted passthrough matches endpoint_summary
  es <- cta_endpoint_summary(.est_2leaf_known())
  for (nid in es$endpoint_node_id) {
    expect_equal(df$n_obs[df$endpoint_node_id == nid],
                 es$n_obs[es$endpoint_node_id == nid])
    expect_equal(df$n_weighted[df$endpoint_node_id == nid],
                 es$n_weighted[es$endpoint_node_id == nid])
  }
  # weighted=FALSE/TRUE stored in column
  expect_true(all(!cta_staging_table(.est_2leaf_known(), weighted = FALSE)$weighted))
  expect_true(all(cta_staging_table(.est_weighted_known(), weighted = TRUE)$weighted))
  # weighted=FALSE uses raw counts; weighted=TRUE uses n_weighted counts
  df_w0 <- cta_staging_table(.est_weighted_known(), weighted = FALSE)
  df_w1 <- cta_staging_table(.est_weighted_known(), weighted = TRUE)
  expect_equal(df_w0$target_n[df_w0$endpoint_node_id == 2L], 2.0)     # raw class1
  expect_equal(df_w0$denominator[df_w0$endpoint_node_id == 2L], 6.0)  # raw total
  expect_equal(df_w1$target_n[df_w1$endpoint_node_id == 2L], 6.0)     # wt class1
  expect_equal(df_w1$denominator[df_w1$endpoint_node_id == 2L], 14.0) # wt total
})

test_that("est: perfect endpoint — detection, odds NA, adjustment arithmetic", {
  df <- cta_staging_table(.est_perfect_ep_tree())
  # node3: class0=0, class1=4 → non_target_n=0 → perfectly_predicted
  expect_true(df$perfectly_predicted[df$endpoint_node_id == 3L])
  expect_false(df$perfectly_predicted[df$endpoint_node_id == 2L])
  expect_true(is.na(df$odds[df$endpoint_node_id == 3L]))
  expect_false(is.na(df$odds[df$endpoint_node_id == 2L]))
  # Adjustment: node3 adjusted; non_target_n boost by 1
  expect_true(df$adjusted[df$endpoint_node_id == 3L])
  expect_false(df$adjusted[df$endpoint_node_id == 2L])
  row3 <- df[df$endpoint_node_id == 3L, ]
  expect_equal(row3$adjusted_non_target_n,       1.0)
  expect_equal(row3$adjusted_target_n,           4.0)
  expect_equal(row3$adjusted_denominator,        5.0)
  expect_equal(row3$adjusted_target_proportion,  4/5)
  expect_equal(row3$adjusted_odds,               4/1)
  # Non-adjusted endpoint: adjusted columns = empirical
  row2 <- df[df$endpoint_node_id == 2L, ]
  expect_equal(row2$adjusted_target_n,          row2$target_n)
  expect_equal(row2$adjusted_denominator,       row2$denominator)
  expect_equal(row2$adjusted_target_proportion, row2$target_proportion)
  # adjust_perfect=FALSE: all adjusted=FALSE; perfect endpoint unadjusted
  df_no <- cta_staging_table(.est_perfect_ep_tree(), adjust_perfect = FALSE)
  expect_true(all(!df_no$adjusted))
  row3_no <- df_no[df_no$endpoint_node_id == 3L, ]
  expect_true(is.na(row3_no$adjusted_odds))
  expect_equal(row3_no$adjusted_target_n, row3_no$target_n)
})

test_that("est: stump fit — 2 rows, denominator sums, weighted variant", {
  tree <- .est_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_staging_table(tree)
  expect_equal(nrow(df), 2L)
  expect_equal(names(df)[1L], "stage")
  expect_true(all(df$target_proportion >= 0 & df$target_proportion <= 1, na.rm = TRUE))
  # Denominator sums match n_obs from endpoint_summary
  es <- cta_endpoint_summary(tree)
  for (nid in es$endpoint_node_id)
    expect_equal(df$denominator[df$endpoint_node_id == nid],
                 as.numeric(es$n_obs[es$endpoint_node_id == nid]))
  # Weighted stump: denominator = n_weighted
  wtree <- .est_weighted_fit()
  skip_if(isTRUE(wtree$no_tree), "mc sampling missed")
  dfw <- cta_staging_table(wtree, weighted = TRUE)
  esw <- cta_endpoint_summary(wtree)
  for (nid in esw$endpoint_node_id)
    expect_equal(dfw$denominator[dfw$endpoint_node_id == nid],
                 esw$n_weighted[esw$endpoint_node_id == nid])
})

# =============================================================================
# Smoke: myeloma canon
# =============================================================================

.est_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.est_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.est_myeloma_fit <- local({
  fits <- list()
  function(mindenom) {
    key <- as.character(mindenom)
    if (is.null(fits[[key]])) {
      df <- .est_load_myeloma()
      fits[[key]] <<- suppressMessages(oda_cta_fit(
        X = df[, .est_myeloma_attrs], y = as.integer(df$V1), w = df$V2,
        priors_on = TRUE, miss_codes = -9, alpha_split = 0.05,
        mindenom = mindenom, prune_alpha = 0.05, max_depth = 20L,
        mc_iter = 5000L, mc_target = 0.05, mc_stop = 99.9,
        mc_seed = NULL, loo = "stable", attr_names = .est_myeloma_attrs
      ))
    }
    fits[[key]]
  }
})

test_that("est: myeloma — rows, stage order, target reconcile for MINDENOM=56/30/1", {
  skip_if_slow_tests_disabled("cta-staging-table")

  # MINDENOM=56: no-tree → 0 rows, 21 columns
  df56 <- cta_staging_table(.est_myeloma_fit(56L))
  expect_equal(nrow(df56), 0L)
  expect_equal(length(names(df56)), 21L)

  # MINDENOM=30: stump → 2 rows, stage 1:2, ascending prop, target_class=1
  df30 <- cta_staging_table(.est_myeloma_fit(30L))
  expect_equal(nrow(df30), 2L)
  expect_equal(df30$stage, c(1L, 2L))
  expect_true(df30$adjusted_target_proportion[1L] <= df30$adjusted_target_proportion[2L])
  expect_equal(unique(df30$target_class), 1L)
  # sum(target_n) matches confusion table
  expect_equal(sum(df30$target_n),
               as.numeric(sum(cta_confusion_table(.est_myeloma_fit(30L))$n[
                 cta_confusion_table(.est_myeloma_fit(30L))$actual == 1L])))
  # weighted=TRUE denominator matches n_weighted
  t30 <- .est_myeloma_fit(30L)
  dfw30 <- cta_staging_table(t30, weighted = TRUE)
  es30  <- cta_endpoint_summary(t30)
  for (nid in es30$endpoint_node_id)
    expect_equal(dfw30$denominator[dfw30$endpoint_node_id == nid],
                 es30$n_weighted[es30$endpoint_node_id == nid])

  # MINDENOM=1: 3 rows, stage-1 has lowest target_proportion
  df1 <- cta_staging_table(.est_myeloma_fit(1L))
  expect_equal(nrow(df1), 3L)
  expect_true(df1$target_proportion[1L] == min(df1$target_proportion, na.rm = TRUE))
})
