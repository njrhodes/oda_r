###############################################################################
# test-cta-endpoint-count-table.R — cta_endpoint_counts()
#
# Per-endpoint × class count table.
# Reads stored leaf$class_counts_raw and leaf$class_counts_weighted only.
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

# Synthetic 3-leaf tree (constructed) — stable anchor for known values.
# node2: 4 obs class-0 only; node4: 3 obs class-0; node5: 3 obs class-1.
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
    list(nodes = list(node1, node2, node3, node4, node5),
         no_tree = FALSE, overall_ess = 75.0, n_nodes = 5L,
         root_id = 1L, has_weights = FALSE, mindenom = 2L,
         alpha_split = 0.05, prune_alpha = 1.0, loo = "off"),
    class = "cta_tree"
  )
}

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
# Contract tests
# =============================================================================

test_that("ect: schema — data.frame, correct columns/types, no-tree zero rows", {
  df <- cta_endpoint_counts(.ect_no_tree_fit())
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 0L)
  expected_cols <- c("endpoint_id", "endpoint_node_id", "path",
                     "terminal_prediction", "class", "n_raw", "n_weighted")
  expect_equal(names(df), expected_cols)
  expect_type(df$endpoint_id, "integer")
  expect_type(df$path,        "character")
  expect_type(df$class,       "character")
  expect_type(df$n_raw,       "integer")
  expect_type(df$n_weighted,  "double")
  # No forbidden staging/event-rate columns
  forbidden <- c("target_class", "event_rate", "target_proportion",
                 "odds", "staging_order", "staging_tier")
  expect_true(!any(forbidden %in% names(df)))
})

test_that("ect: stump — 4 rows, class labels, zero-count preserved, totals reconcile", {
  tree <- .ect_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df   <- cta_endpoint_counts(tree)
  expect_equal(nrow(df), 4L)
  # Both class "0" and "1" on each endpoint
  for (eid in unique(df$endpoint_id)) {
    expect_true(setequal(df$class[df$endpoint_id == eid], c("0", "1")))
  }
  # Zero-count class preserved (perfectly separated data)
  expect_true(any(df$n_raw == 0L))
  # Totals reconcile with confusion table
  ct <- cta_confusion_table(tree)
  expect_equal(sum(df$n_raw[df$class == "0"]), sum(ct$n[ct$actual == 0L]))
  expect_equal(sum(df$n_raw[df$class == "1"]), sum(ct$n[ct$actual == 1L]))
  # Per-endpoint raw totals match endpoint_summary denominators
  epsum      <- cta_endpoint_summary(tree)
  raw_totals <- as.integer(tapply(df$n_raw, df$endpoint_id, sum))
  expect_equal(raw_totals, epsum$denominator)
})

test_that("ect: weighted stump — n_weighted differs from n_raw", {
  tree <- .ect_weighted_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_counts(tree)
  expect_true(any(df$n_weighted != as.numeric(df$n_raw)))
})

test_that("ect: missing class_counts_raw on valid tree errors with message", {
  expect_error(cta_endpoint_counts(.ect_missing_counts_tree()),
               "endpoint class counts are unavailable")
})

test_that("ect: synthetic 3-leaf — known node values correct", {
  df   <- cta_endpoint_counts(.ect_3leaf_tree())
  expect_equal(nrow(df), 6L)
  sub2 <- df[df$endpoint_node_id == 2L, ]
  expect_equal(sub2$n_raw[sub2$class == "0"], 4L)
  expect_equal(sub2$n_raw[sub2$class == "1"], 0L)
  sub5 <- df[df$endpoint_node_id == 5L, ]
  expect_equal(sub5$n_raw[sub5$class == "0"], 0L)
  expect_equal(sub5$n_raw[sub5$class == "1"], 3L)
})

# =============================================================================
# Smoke: myeloma canon integration
# =============================================================================

.ect_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.ect_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.ect_myeloma_fit <- local({
  fits <- list()
  function(mindenom) {
    key <- as.character(mindenom)
    if (is.null(fits[[key]])) {
      df <- .ect_load_myeloma()
      fits[[key]] <<- suppressMessages(oda_cta_fit(
        X = df[, .ect_myeloma_attrs], y = as.integer(df$V1), w = df$V2,
        priors_on = TRUE, miss_codes = -9, alpha_split = 0.05,
        mindenom = mindenom, prune_alpha = 0.05, max_depth = 20L,
        mc_iter = 5000L, mc_target = 0.05, mc_stop = 99.9, mc_stopup = 99.9,
        mc_seed = NULL, loo = "stable", attr_names = .ect_myeloma_attrs
      ))
    }
    fits[[key]]
  }
})

test_that("ect: myeloma — rows/totals/confusion reconciliation for MINDENOM=1/30/56", {
  skip_if_slow_tests_disabled("cta-endpoint-count-table")

  # MINDENOM=1: 3 endpoints × 2 classes = 6 rows, total 255 obs
  t1 <- .ect_myeloma_fit(1L)
  df1 <- cta_endpoint_counts(t1)
  expect_equal(nrow(df1), 6L)
  expect_equal(sum(df1$n_raw), 255L)
  ct1 <- cta_confusion_table(t1)
  expect_equal(sum(df1$n_raw[df1$class == "0"]), sum(ct1$n[ct1$actual == 0L]))
  expect_equal(sum(df1$n_raw[df1$class == "1"]), sum(ct1$n[ct1$actual == 1L]))

  # MINDENOM=30: 2 endpoints (stump) × 2 classes = 4 rows, total 186 obs
  t30 <- .ect_myeloma_fit(30L)
  df30 <- cta_endpoint_counts(t30)
  expect_equal(nrow(df30), 4L)
  expect_equal(sum(df30$n_raw), 186L)
  ct30 <- cta_confusion_table(t30)
  expect_equal(sum(df30$n_raw[df30$class == "0"]), sum(ct30$n[ct30$actual == 0L]))
  expect_equal(sum(df30$n_raw[df30$class == "1"]), sum(ct30$n[ct30$actual == 1L]))

  # MINDENOM=56: no-tree → zero rows with correct columns
  df56 <- cta_endpoint_counts(.ect_myeloma_fit(56L))
  expect_equal(nrow(df56), 0L)
  expect_equal(names(df56), c("endpoint_id", "endpoint_node_id", "path",
                               "terminal_prediction", "class", "n_raw", "n_weighted"))
})
