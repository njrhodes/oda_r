###############################################################################
# test-cta-endpoint-count-table.R — cta_endpoint_counts()
#
# Per-endpoint × class count table.
# Reads stored leaf$class_counts_raw and leaf$class_counts_weighted only.
# Returns one row per terminal endpoint per actual class.
###############################################################################

# ---- Helpers ----------------------------------------------------------------

.ect_no_tree_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 999L, mc_iter = 100L,
                mc_seed = 1L, loo = "off")
  )
}

.ect_stump_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off")
  )
}

.ect_weighted_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  w <- c(2, 2, 2, 2, 3, 3, 3, 3)
  suppressMessages(
    oda_cta_fit(d$X, d$y, w = w, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off")
  )
}

# Synthetic 3-leaf tree with class_counts_raw/weighted on leaf nodes.
# node2: 4 obs all class-0; node4: 3 obs all class-0; node5: 3 obs all class-1.
.ect_3leaf_tree <- function() {
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
    class_counts_raw      = c("0" = 4L, "1" = 0L),
    class_counts_weighted = c("0" = 4.0, "1" = 0.0),
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
    class_counts_raw      = c("0" = 3L, "1" = 0L),
    class_counts_weighted = c("0" = 3.0, "1" = 0.0),
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  node5 <- list(
    node_id = 5L, parent_id = 3L, depth = 3L, leaf = TRUE,
    attribute = NA_character_, attr_col = NA_integer_, attr_type = NA_character_,
    n_obs = 3L, n_weighted = 3, majority_class = 1L, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    class_counts_raw      = c("0" = 0L, "1" = 3L),
    class_counts_weighted = c("0" = 0.0, "1" = 3.0),
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

# Synthetic cta_tree with a leaf missing class_counts_raw — for guard test.
# no_tree = FALSE so the guard is reached.
.ect_missing_counts_tree <- function() {
  leaf <- list(
    node_id = 1L, parent_id = 0L, depth = 1L, leaf = TRUE,
    majority_class = 0L, n_obs = 4L, n_weighted = 4,
    class_counts_raw = NULL, class_counts_weighted = NULL,
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
# Return type and column structure
# =============================================================================

test_that("ect: returns a data.frame", {
  expect_s3_class(cta_endpoint_counts(.ect_no_tree_fit()), "data.frame")
})

test_that("ect: no-tree returns zero rows", {
  expect_equal(nrow(cta_endpoint_counts(.ect_no_tree_fit())), 0L)
})

test_that("ect: no-tree has exactly the required columns", {
  expected <- c("endpoint_id", "endpoint_node_id", "path",
                "terminal_prediction", "class", "n_raw", "n_weighted")
  expect_equal(names(cta_endpoint_counts(.ect_no_tree_fit())), expected)
})

test_that("ect: no-tree column types are correct", {
  df <- cta_endpoint_counts(.ect_no_tree_fit())
  expect_type(df$endpoint_id,         "integer")
  expect_type(df$endpoint_node_id,    "integer")
  expect_type(df$path,                "character")
  expect_type(df$terminal_prediction, "integer")
  expect_type(df$class,               "character")
  expect_type(df$n_raw,               "integer")
  expect_type(df$n_weighted,          "double")
})

# =============================================================================
# Stump: 2 endpoints × 2 classes = 4 rows
# =============================================================================

test_that("ect: stump returns 4 rows", {
  tree <- .ect_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed — stump tests need a split")
  expect_equal(nrow(cta_endpoint_counts(tree)), 4L)
})

test_that("ect: stump endpoint_id is c(1L,1L,2L,2L)", {
  tree <- .ect_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_counts(tree)
  expect_equal(df$endpoint_id, c(1L, 1L, 2L, 2L))
})

test_that("ect: stump each endpoint has class rows '0' and '1'", {
  tree <- .ect_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_counts(tree)
  for (eid in c(1L, 2L)) {
    expect_equal(sort(df$class[df$endpoint_id == eid]), c("0", "1"))
  }
})

test_that("ect: stump per-endpoint sum(n_raw) equals cta_endpoint_summary() denominator", {
  tree  <- .ect_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df    <- cta_endpoint_counts(tree)
  epsum <- cta_endpoint_summary(tree)
  raw_totals <- as.integer(tapply(df$n_raw, df$endpoint_id, sum))
  expect_equal(raw_totals, epsum$denominator)
})

test_that("ect: stump per-endpoint sum(n_weighted) equals cta_endpoint_summary() n_weighted", {
  tree  <- .ect_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df    <- cta_endpoint_counts(tree)
  epsum <- cta_endpoint_summary(tree)
  wt_totals <- as.numeric(tapply(df$n_weighted, df$endpoint_id, sum))
  expect_equal(wt_totals, epsum$n_weighted)
})

test_that("ect: stump total n_raw equals sum(cta_confusion_table()$n)", {
  tree <- .ect_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_counts(tree)
  ct <- cta_confusion_table(tree)
  expect_equal(sum(df$n_raw), sum(ct$n))
})

test_that("ect: stump actual-class n_raw totals match cta_confusion_table() marginals", {
  tree <- .ect_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_counts(tree)
  ct <- cta_confusion_table(tree)
  expect_equal(sum(df$n_raw[df$class == "0"]), sum(ct$n[ct$actual == 0L]))
  expect_equal(sum(df$n_raw[df$class == "1"]), sum(ct$n[ct$actual == 1L]))
})

test_that("ect: stump zero-count class rows are preserved", {
  # y = c(0,0,0,0,1,1,1,1) with cut ~4.5 → left leaf all class-0, right all
  # class-1.  Both "0" and "1" rows must be present on each endpoint.
  tree <- .ect_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_counts(tree)
  for (eid in unique(df$endpoint_id)) {
    sub <- df[df$endpoint_id == eid, ]
    expect_true(setequal(sub$class, c("0", "1")))
  }
  expect_true(any(df$n_raw == 0L))
})

test_that("ect: unweighted stump n_weighted equals n_raw numerically", {
  tree <- .ect_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_counts(tree)
  expect_equal(unname(df$n_weighted), unname(as.numeric(df$n_raw)))
})

test_that("ect: weighted fit has at least one n_weighted != as.numeric(n_raw)", {
  tree <- .ect_weighted_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_counts(tree)
  expect_true(any(df$n_weighted != as.numeric(df$n_raw)))
})

# =============================================================================
# Guard: missing class_counts_raw on a valid tree
# =============================================================================

test_that("ect: valid tree missing class_counts_raw stops with clear message", {
  expect_error(cta_endpoint_counts(.ect_missing_counts_tree()),
               "endpoint class counts are unavailable")
})

# =============================================================================
# Synthetic 3-leaf tree
# =============================================================================

test_that("ect: synthetic 3-leaf returns 6 rows", {
  expect_equal(nrow(cta_endpoint_counts(.ect_3leaf_tree())), 6L)
})

test_that("ect: synthetic 3-leaf node2 n_raw[class=='0']==4 and n_raw[class=='1']==0", {
  df   <- cta_endpoint_counts(.ect_3leaf_tree())
  sub2 <- df[df$endpoint_node_id == 2L, ]
  expect_equal(sub2$n_raw[sub2$class == "0"], 4L)
  expect_equal(sub2$n_raw[sub2$class == "1"], 0L)
})

test_that("ect: synthetic 3-leaf node5 n_raw[class=='1']==3 (majority class)", {
  df   <- cta_endpoint_counts(.ect_3leaf_tree())
  sub5 <- df[df$endpoint_node_id == 5L, ]
  expect_equal(sub5$n_raw[sub5$class == "0"], 0L)
  expect_equal(sub5$n_raw[sub5$class == "1"], 3L)
})

# =============================================================================
# Forbidden columns
# =============================================================================

test_that("ect: no forbidden event-rate or staging columns present", {
  df       <- cta_endpoint_counts(.ect_3leaf_tree())
  forbidden <- c("target_class", "event_rate", "target_proportion",
                 "odds", "staging_order", "status", "staging_tier")
  expect_true(!any(forbidden %in% names(df)))
})

# =============================================================================
# Slow fixture tests — myeloma canon
# =============================================================================

.ect_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.ect_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.ect_myeloma_fit <- function(mindenom) {
  df <- .ect_load_myeloma()
  suppressMessages(
    oda_cta_fit(
      X           = df[, .ect_myeloma_attrs],
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
      attr_names  = .ect_myeloma_attrs
    )
  )
}

test_that("ect: myeloma MINDENOM=1 returns 6 rows (3 endpoints x 2 classes)", {
  skip_if_slow_tests_disabled("cta-endpoint-count-table")
  df <- cta_endpoint_counts(.ect_myeloma_fit(1L))
  expect_equal(nrow(df), 6L)
})

test_that("ect: myeloma MINDENOM=30 returns 4 rows (2 endpoints x 2 classes)", {
  skip_if_slow_tests_disabled("cta-endpoint-count-table")
  df <- cta_endpoint_counts(.ect_myeloma_fit(30L))
  expect_equal(nrow(df), 4L)
})

test_that("ect: myeloma MINDENOM=56 returns zero rows with exact columns", {
  skip_if_slow_tests_disabled("cta-endpoint-count-table")
  df <- cta_endpoint_counts(.ect_myeloma_fit(56L))
  expect_equal(nrow(df), 0L)
  expected <- c("endpoint_id", "endpoint_node_id", "path",
                "terminal_prediction", "class", "n_raw", "n_weighted")
  expect_equal(names(df), expected)
})

test_that("ect: myeloma MINDENOM=1 sum(n_raw) = 255", {
  skip_if_slow_tests_disabled("cta-endpoint-count-table")
  df <- cta_endpoint_counts(.ect_myeloma_fit(1L))
  expect_equal(sum(df$n_raw), 255L)
})

test_that("ect: myeloma MINDENOM=30 sum(n_raw) = 186", {
  skip_if_slow_tests_disabled("cta-endpoint-count-table")
  df <- cta_endpoint_counts(.ect_myeloma_fit(30L))
  expect_equal(sum(df$n_raw), 186L)
})

test_that("ect: myeloma MINDENOM=1 class totals reconcile with cta_confusion_table()", {
  skip_if_slow_tests_disabled("cta-endpoint-count-table")
  tree <- .ect_myeloma_fit(1L)
  df   <- cta_endpoint_counts(tree)
  ct   <- cta_confusion_table(tree)
  expect_equal(sum(df$n_raw[df$class == "0"]), sum(ct$n[ct$actual == 0L]))
  expect_equal(sum(df$n_raw[df$class == "1"]), sum(ct$n[ct$actual == 1L]))
})

test_that("ect: myeloma MINDENOM=30 class totals reconcile with cta_confusion_table()", {
  skip_if_slow_tests_disabled("cta-endpoint-count-table")
  tree <- .ect_myeloma_fit(30L)
  df   <- cta_endpoint_counts(tree)
  ct   <- cta_confusion_table(tree)
  expect_equal(sum(df$n_raw[df$class == "0"]), sum(ct$n[ct$actual == 0L]))
  expect_equal(sum(df$n_raw[df$class == "1"]), sum(ct$n[ct$actual == 1L]))
})

test_that("ect: myeloma MINDENOM=1 n_weighted all nonnegative", {
  skip_if_slow_tests_disabled("cta-endpoint-count-table")
  df <- cta_endpoint_counts(.ect_myeloma_fit(1L))
  expect_true(all(df$n_weighted >= 0))
})

test_that("ect: myeloma MINDENOM=30 n_weighted all nonnegative", {
  skip_if_slow_tests_disabled("cta-endpoint-count-table")
  df <- cta_endpoint_counts(.ect_myeloma_fit(30L))
  expect_true(all(df$n_weighted >= 0))
})

test_that("ect: myeloma MINDENOM=1 per-endpoint sum(n_weighted) matches cta_endpoint_summary()$n_weighted", {
  skip_if_slow_tests_disabled("cta-endpoint-count-table")
  tree  <- .ect_myeloma_fit(1L)
  df    <- cta_endpoint_counts(tree)
  epsum <- cta_endpoint_summary(tree)
  wt_totals <- vapply(sort(unique(df$endpoint_id)), function(eid) {
    sum(df$n_weighted[df$endpoint_id == eid])
  }, numeric(1L))
  expect_equal(wt_totals, epsum$n_weighted, tolerance = 1e-10)
})
