###############################################################################
# test-cta-confusion-table.R — cta_confusion_table() (Phase 2B)
#
# Semantic contract:
#   cta_confusion_table(tree) returns FINAL SELECTED TREE training confusion,
#   stored on tree$training_confusion at fit time.
#   It does NOT return split-node local confusion.
###############################################################################

# ---- Helpers ----------------------------------------------------------------

.ctbl_no_tree_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 999L, mc_iter = 100L,
                mc_seed = 1L, loo = "off")
  )
}

.ctbl_stump_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off")
  )
}

# Manually constructed cta_tree with TWO different confusion values:
#   - split node 1 has $confusion = [[3,1],[1,3]]  (local split confusion)
#   - $training_confusion = [[10,2],[3,9]]          (stored final-tree confusion)
# cta_confusion_table() must return the final-tree value, not the split confusion.
# This proves final-tree semantics without a slow fit.
.constructed_semantic_tree <- function() {
  # Split-node local confusion (intentionally different from final-tree confusion)
  split_conf <- matrix(c(3L, 1L, 1L, 3L), nrow = 2L,
                       dimnames = list(actual = c("0","1"),
                                       predicted = c("0","1")))
  # Final-tree training confusion (what cta_confusion_table must return)
  final_conf  <- matrix(c(10L, 3L, 2L, 9L), nrow = 2L,
                        dimnames = list(actual = c("0","1"),
                                        predicted = c("0","1")))
  node1 <- list(
    node_id = 1L, parent_id = 0L, depth = 1L, leaf = FALSE,
    attribute = "x", attr_col = 1L, attr_type = "ordered",
    n_obs = 24L, n_weighted = 24,
    rule = list(type = "ordered_cut", cut_value = 4.5, direction = "0->1"),
    ess = 60.0, ess_weighted = NA_real_, p_mc = 0.01,
    loo_status = "STABLE", loo_ess = 58.0, loo_p = NA_real_,
    confusion = split_conf,            # local split-node confusion
    split_labels = c(0L, 1L), child_ids = c(2L, 3L)
  )
  node2 <- list(
    node_id = 2L, parent_id = 1L, depth = 2L, leaf = TRUE,
    n_obs = 12L, n_weighted = 12, majority_class = 0L,
    attribute = NA_character_, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  node3 <- list(
    node_id = 3L, parent_id = 1L, depth = 2L, leaf = TRUE,
    n_obs = 12L, n_weighted = 12, majority_class = 1L,
    attribute = NA_character_, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  structure(
    list(nodes              = list(node1, node2, node3),
         no_tree            = FALSE,
         overall_ess        = 62.5,
         n_nodes            = 3L,
         root_id            = 1L,
         has_weights        = FALSE,
         mindenom           = 2L,
         alpha_split        = 0.05,
         prune_alpha        = 1.0,
         loo                = "off",
         training_confusion = final_conf),   # final-tree confusion stored here
    class = "cta_tree"
  )
}

# No-tree with training_confusion = NULL (as produced by oda_cta_fit no-tree path)
.constructed_no_tree <- function() {
  node1 <- list(
    node_id = 1L, parent_id = 0L, depth = 1L, leaf = TRUE,
    n_obs = 8L, n_weighted = 8, majority_class = 0L,
    attribute = NA_character_, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  structure(
    list(nodes              = list(node1),
         no_tree            = TRUE,
         overall_ess        = NA_real_,
         n_nodes            = 1L,
         root_id            = 1L,
         has_weights        = FALSE,
         mindenom           = 999L,
         alpha_split        = 0.05,
         prune_alpha        = 1.0,
         loo                = "off",
         training_confusion = NULL),
    class = "cta_tree"
  )
}

# =============================================================================
# Return type and column structure
# =============================================================================

test_that("cta_confusion_table: returns a data.frame", {
  expect_s3_class(cta_confusion_table(.ctbl_no_tree_fit()), "data.frame")
})

test_that("cta_confusion_table: no-tree fit returns zero rows", {
  expect_equal(nrow(cta_confusion_table(.ctbl_no_tree_fit())), 0L)
})

test_that("cta_confusion_table: constructed no-tree returns zero rows", {
  expect_equal(nrow(cta_confusion_table(.constructed_no_tree())), 0L)
})

test_that("cta_confusion_table: columns are exactly actual/predicted/n", {
  df <- cta_confusion_table(.ctbl_no_tree_fit())
  expect_equal(names(df), c("actual", "predicted", "n"))
})

test_that("cta_confusion_table: column types correct for empty df", {
  df <- cta_confusion_table(.ctbl_no_tree_fit())
  expect_type(df$actual,    "integer")
  expect_type(df$predicted, "integer")
  expect_type(df$n,         "integer")
})

# =============================================================================
# Semantic test: final-tree confusion, NOT split-node confusion
# =============================================================================

test_that("cta_confusion_table: returns training_confusion, not split-node confusion", {
  tree <- .constructed_semantic_tree()

  # Verify the tree has intentionally different split-node and final-tree confusion
  root_split_conf <- tree$nodes[[1L]]$confusion
  final_conf      <- tree$training_confusion
  expect_false(identical(root_split_conf, final_conf),
               label = "test setup: split confusion must differ from final confusion")

  df <- cta_confusion_table(tree)

  # Must match final_conf = [[10,2],[3,9]], not split_conf = [[3,1],[1,3]]
  expect_equal(df$n[df$actual == 0L & df$predicted == 0L], 10L)
  expect_equal(df$n[df$actual == 0L & df$predicted == 1L],  2L)
  expect_equal(df$n[df$actual == 1L & df$predicted == 0L],  3L)
  expect_equal(df$n[df$actual == 1L & df$predicted == 1L],  9L)
})

test_that("cta_confusion_table: does not return split-node local confusion values", {
  tree <- .constructed_semantic_tree()
  df   <- cta_confusion_table(tree)
  # split_conf values are 3,1,1,3 — none of these should appear as n values
  # (final_conf values are 10,2,3,9 which are all distinct from 3,1,1,3 except n=3)
  # The key differentiator: actual=0,predicted=0 must be 10 not 3.
  n_00 <- df$n[df$actual == 0L & df$predicted == 0L]
  expect_equal(n_00, 10L,
               label = "TN must be 10 (final-tree) not 3 (split-node)")
})

test_that("cta_confusion_table: semantic tree total n = 24", {
  df <- cta_confusion_table(.constructed_semantic_tree())
  expect_equal(sum(df$n), 24L)
})

# Class-level stability: training_confusion spans fit-level class universe.
# A tree where the final confusion might only observe one class in predictions
# should still have rows/columns for all fit-level classes.
test_that("cta_confusion_table: class-level stable — both classes present in matrix", {
  # constructed_semantic_tree has training_confusion with dimnames "0"/"1"
  tree <- .constructed_semantic_tree()
  conf <- tree$training_confusion
  expect_true(is.matrix(conf))
  expect_equal(rownames(conf), c("0", "1"))
  expect_equal(colnames(conf), c("0", "1"))
  # The tidy output covers all four cells, including any zero cells
  df <- cta_confusion_table(tree)
  expect_equal(nrow(df), 4L)
  expect_true(all(c(0L, 1L) %in% df$actual))
  expect_true(all(c(0L, 1L) %in% df$predicted))
})

# =============================================================================
# Stump fit: final-tree confusion stored at fit time
# =============================================================================

test_that("cta_confusion_table: stump fit has 4 rows (2x2)", {
  tree <- .ctbl_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed — stump needs a split")
  df <- cta_confusion_table(tree)
  expect_equal(nrow(df), 4L)
})

test_that("cta_confusion_table: stump fit has training_confusion stored", {
  tree <- .ctbl_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  expect_true(is.matrix(tree$training_confusion))
})

test_that("cta_confusion_table: stump n sums to total classified obs", {
  tree <- .ctbl_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_confusion_table(tree)
  # Perfectly separable 8-obs dataset: all 8 obs classified
  expect_equal(sum(df$n), 8L)
})

test_that("cta_confusion_table: stump actual and predicted are 0 or 1", {
  tree <- .ctbl_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_confusion_table(tree)
  expect_true(all(df$actual    %in% c(0L, 1L)))
  expect_true(all(df$predicted %in% c(0L, 1L)))
})

test_that("cta_confusion_table: stump n values are non-negative integers", {
  tree <- .ctbl_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_confusion_table(tree)
  expect_type(df$n, "integer")
  expect_true(all(df$n >= 0L))
})

test_that("cta_confusion_table: stump perfect separation has zero off-diagonal", {
  tree <- .ctbl_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_confusion_table(tree)
  fp <- df$n[df$actual == 0L & df$predicted == 1L]
  fn <- df$n[df$actual == 1L & df$predicted == 0L]
  expect_equal(fp, 0L)
  expect_equal(fn, 0L)
})

test_that("cta_confusion_table: stump rows sorted actual then predicted", {
  tree <- .ctbl_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  df <- cta_confusion_table(tree)
  expect_equal(df$actual,    c(0L, 0L, 1L, 1L))
  expect_equal(df$predicted, c(0L, 1L, 0L, 1L))
})

# =============================================================================
# Slow fixture tests — canon anchor validation
# =============================================================================

.ctbl_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.ctbl_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.ctbl_myeloma_fit <- function(mindenom) {
  df <- .ctbl_load_myeloma()
  suppressMessages(
    oda_cta_fit(
      X           = df[, .ctbl_myeloma_attrs],
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
      attr_names  = .ctbl_myeloma_attrs
    )
  )
}

test_that("cta_confusion_table: myeloma MINDENOM=56 no-tree returns zero rows", {
  skip_if_slow_tests_disabled("cta-confusion-fixture")
  tree <- .ctbl_myeloma_fit(56L)
  expect_true(isTRUE(tree$no_tree))
  df <- cta_confusion_table(tree)
  expect_equal(nrow(df), 0L)
  expect_equal(names(df), c("actual", "predicted", "n"))
})

test_that("cta_confusion_table: myeloma MINDENOM=30 stump = [[101,34],[30,21]], n=186", {
  skip_if_slow_tests_disabled("cta-confusion-fixture")
  tree <- .ctbl_myeloma_fit(30L)
  # LOO STABLE selects V17 stump, path-local (186 obs classified)
  expect_equal(cta_strata(tree), 2L)
  df <- cta_confusion_table(tree)
  expect_equal(nrow(df), 4L)
  expect_equal(sum(df$n), 186L)
  expect_equal(df$n[df$actual == 0L & df$predicted == 0L], 101L)
  expect_equal(df$n[df$actual == 0L & df$predicted == 1L],  34L)
  expect_equal(df$n[df$actual == 1L & df$predicted == 0L],  30L)
  expect_equal(df$n[df$actual == 1L & df$predicted == 1L],  21L)
})

test_that("cta_confusion_table: myeloma MINDENOM=1 tree = [[146,40],[36,33]], n=255", {
  skip_if_slow_tests_disabled("cta-confusion-fixture")
  tree <- .ctbl_myeloma_fit(1L)
  # LOO STABLE selects V14->V15 (3 leaves), n_classified=255
  expect_equal(cta_strata(tree), 3L)
  df <- cta_confusion_table(tree)
  expect_equal(nrow(df), 4L)
  expect_equal(sum(df$n), 255L)
  expect_equal(df$n[df$actual == 0L & df$predicted == 0L], 146L)
  expect_equal(df$n[df$actual == 0L & df$predicted == 1L],  40L)
  expect_equal(df$n[df$actual == 1L & df$predicted == 0L],  36L)
  expect_equal(df$n[df$actual == 1L & df$predicted == 1L],  33L)
})
