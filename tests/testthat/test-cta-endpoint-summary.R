###############################################################################
# test-cta-endpoint-summary.R - cta_endpoint_summary()
#
# One row per terminal leaf with structural fields only.
# No class counts, no target-class proportions, no staging order.
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

# Synthetic 3-leaf tree - stable anchor for path/structure known values.
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

test_that("epsum: schema - data.frame, correct columns/types, no-tree zero rows", {
  df <- cta_endpoint_summary(.epsum_no_tree_fit())
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 0L)
  expected <- c("endpoint_id", "endpoint_node_id", "path", "depth",
                "terminal_prediction", "n_obs", "n_weighted", "denominator")
  expect_equal(names(df), expected)
  expect_type(df$endpoint_id,         "integer")
  expect_type(df$path,                "character")
  expect_type(df$depth,               "integer")
  expect_type(df$n_obs,               "integer")
  expect_type(df$n_weighted,          "double")
  expect_type(df$denominator,         "integer")
  # No forbidden staging/event-rate columns
  expect_false(any(c("target_class","target_proportion","odds","staging_tier") %in% names(df)))
})

test_that("epsum: stump - 2 rows, denominator=n_obs, valid terminal_prediction", {
  tree <- .epsum_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_endpoint_summary(tree)
  expect_equal(nrow(df), 2L)
  expect_equal(df$endpoint_id, c(1L, 2L))
  expect_equal(df$denominator, df$n_obs)
  expect_true(all(df$terminal_prediction %in% c(0L, 1L)))
  expect_true(all(nzchar(df$path)))
})

test_that("epsum: synthetic 3-leaf - known node IDs, paths, predictions, depths", {
  df <- cta_endpoint_summary(.epsum_3leaf_tree())
  expect_equal(nrow(df), 3L)
  expect_equal(df$endpoint_node_id, c(2L, 4L, 5L))
  expect_equal(df$path[df$endpoint_node_id == 2L], "x<=4.5")
  expect_equal(df$terminal_prediction[df$endpoint_node_id == 2L], 0L)
  expect_equal(df$terminal_prediction[df$endpoint_node_id == 5L], 1L)
  expect_equal(df$depth[df$endpoint_node_id == 2L], 2L)
  expect_equal(df$depth[df$endpoint_node_id == 5L], 3L)
  expect_equal(df$denominator, df$n_obs)
})

# =============================================================================
# Smoke: myeloma canon
# =============================================================================

.epsum_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.epsum_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.epsum_myeloma_fit <- local({
  fits <- list()
  function(mindenom) {
    key <- as.character(mindenom)
    if (is.null(fits[[key]])) {
      df <- .epsum_load_myeloma()
      fits[[key]] <<- suppressMessages(oda_cta_fit(
        X = df[, .epsum_myeloma_attrs], y = as.integer(df$V1), w = df$V2,
        priors_on = TRUE, miss_codes = -9, alpha_split = 0.05,
        mindenom = mindenom, prune_alpha = 0.05, max_depth = 20L,
        mc_iter = 5000L, mc_target = 0.05, mc_stop = 99.9,
        mc_seed = NULL, loo = "stable", attr_names = .epsum_myeloma_attrs
      ))
    }
    fits[[key]]
  }
})

test_that("epsum: myeloma - rows/paths for MINDENOM=1/30/56; denominators reconcile", {
  skip_if_slow_tests_disabled("cta-endpoint-summary")

  # MINDENOM=1: 3 endpoints, paths all mention V14 (root)
  df1 <- cta_endpoint_summary(.epsum_myeloma_fit(1L))
  expect_equal(nrow(df1), 3L)
  expect_true(all(grepl("V14", df1$path, fixed = TRUE)))
  expect_equal(sort(df1$denominator),
               sort(unname(cta_endpoint_denominators(.epsum_myeloma_fit(1L)))))

  # MINDENOM=30: 2 endpoints (stump), paths all mention V17
  df30 <- cta_endpoint_summary(.epsum_myeloma_fit(30L))
  expect_equal(nrow(df30), 2L)
  expect_true(all(grepl("V17", df30$path, fixed = TRUE)))

  # MINDENOM=56: no-tree -> zero rows with correct columns
  df56 <- cta_endpoint_summary(.epsum_myeloma_fit(56L))
  expect_equal(nrow(df56), 0L)
  expect_equal(names(df56), c("endpoint_id", "endpoint_node_id", "path", "depth",
                               "terminal_prediction", "n_obs", "n_weighted", "denominator"))
})
