###############################################################################
# test-cta-confusion-table.R — cta_confusion_table()
#
# Reads training_confusion stored at fit time (final selected tree).
# Must NOT return split-node local confusion.
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

# Constructed tree with intentionally different split-node ($confusion) and
# final-tree ($training_confusion) values. Proves final-tree semantics.
#   split_conf = [[3,1],[1,3]];  final_conf = [[10,2],[3,9]]
.constructed_semantic_tree <- function() {
  split_conf <- matrix(c(3L, 1L, 1L, 3L), nrow = 2L,
                       dimnames = list(actual = c("0","1"), predicted = c("0","1")))
  final_conf  <- matrix(c(10L, 3L, 2L, 9L), nrow = 2L,
                        dimnames = list(actual = c("0","1"), predicted = c("0","1")))
  node1 <- list(
    node_id = 1L, parent_id = 0L, depth = 1L, leaf = FALSE,
    attribute = "x", attr_col = 1L, attr_type = "ordered",
    n_obs = 24L, n_weighted = 24,
    rule = list(type = "ordered_cut", cut_value = 4.5, direction = "0->1"),
    ess = 60.0, ess_weighted = NA_real_, p_mc = 0.01,
    loo_status = "STABLE", loo_ess = 58.0, loo_p = NA_real_,
    confusion = split_conf, split_labels = c(0L, 1L), child_ids = c(2L, 3L)
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
    list(nodes = list(node1, node2, node3),
         no_tree = FALSE, overall_ess = 62.5, n_nodes = 3L,
         root_id = 1L, has_weights = FALSE, mindenom = 2L,
         alpha_split = 0.05, prune_alpha = 1.0, loo = "off",
         training_confusion = final_conf),
    class = "cta_tree"
  )
}

# =============================================================================
# Contract tests
# =============================================================================

test_that("ctbl: schema — data.frame, columns, types, zero rows for no-tree", {
  df <- cta_confusion_table(.ctbl_no_tree_fit())
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 0L)
  expect_equal(names(df), c("actual", "predicted", "n"))
  expect_type(df$actual,    "integer")
  expect_type(df$predicted, "integer")
  expect_type(df$n,         "integer")
  # Constructed no-tree also returns 0 rows
  nt <- structure(
    list(nodes = list(list(node_id=1L,leaf=TRUE,n_obs=8L,n_weighted=8,
                           majority_class=0L,attribute=NA_character_,rule=NULL,
                           ess=NA_real_,ess_weighted=NA_real_,p_mc=NA_real_,
                           loo_status=NA_character_,loo_ess=NA_real_,loo_p=NA_real_,
                           confusion=NULL,split_labels=integer(0),child_ids=integer(0),
                           parent_id=0L,depth=1L)),
         no_tree=TRUE, overall_ess=NA_real_, n_nodes=1L, root_id=1L,
         has_weights=FALSE, mindenom=999L, alpha_split=0.05, prune_alpha=1.0,
         loo="off", training_confusion=NULL),
    class = "cta_tree"
  )
  expect_equal(nrow(cta_confusion_table(nt)), 0L)
})

test_that("ctbl: semantic — returns training_confusion, not split-node confusion", {
  tree <- .constructed_semantic_tree()
  # Verify intentionally different split vs final confusion
  expect_false(identical(tree$nodes[[1L]]$confusion, tree$training_confusion))
  df <- cta_confusion_table(tree)
  # Must match final_conf = [[10,2],[3,9]]
  expect_equal(df$n[df$actual == 0L & df$predicted == 0L], 10L)
  expect_equal(df$n[df$actual == 0L & df$predicted == 1L],  2L)
  expect_equal(df$n[df$actual == 1L & df$predicted == 0L],  3L)
  expect_equal(df$n[df$actual == 1L & df$predicted == 1L],  9L)
  expect_equal(sum(df$n), 24L)
  # Both classes present; 4 cells
  expect_equal(nrow(df), 4L)
  expect_true(all(c(0L, 1L) %in% df$actual))
  expect_true(all(c(0L, 1L) %in% df$predicted))
})

test_that("ctbl: stump — 4 rows, n=8, zero off-diagonal, canonical row order", {
  tree <- .ctbl_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  expect_true(is.matrix(tree$training_confusion))
  df <- cta_confusion_table(tree)
  expect_equal(nrow(df), 4L)
  expect_equal(sum(df$n), 8L)
  expect_type(df$n, "integer")
  expect_true(all(df$n >= 0L))
  expect_true(all(df$actual    %in% c(0L, 1L)))
  expect_true(all(df$predicted %in% c(0L, 1L)))
  # Perfect separation: zero off-diagonal
  expect_equal(df$n[df$actual == 0L & df$predicted == 1L], 0L)
  expect_equal(df$n[df$actual == 1L & df$predicted == 0L], 0L)
  # Canonical row order: sorted by actual then predicted
  expect_equal(df$actual,    c(0L, 0L, 1L, 1L))
  expect_equal(df$predicted, c(0L, 1L, 0L, 1L))
})

# =============================================================================
# Smoke: myeloma canon
# =============================================================================

.ctbl_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.ctbl_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.ctbl_myeloma_fit <- local({
  fits <- list()
  function(mindenom) {
    key <- as.character(mindenom)
    if (is.null(fits[[key]])) {
      df <- .ctbl_load_myeloma()
      fits[[key]] <<- suppressMessages(oda_cta_fit(
        X = df[, .ctbl_myeloma_attrs], y = as.integer(df$V1), w = df$V2,
        priors_on = TRUE, miss_codes = -9, alpha_split = 0.05,
        mindenom = mindenom, prune_alpha = 0.05, max_depth = 20L,
        mc_iter = 5000L, mc_target = 0.05, mc_stop = 99.9,
        mc_seed = NULL, loo = "stable", attr_names = .ctbl_myeloma_attrs
      ))
    }
    fits[[key]]
  }
})

test_that("ctbl: myeloma — exact confusion values for MINDENOM=56/30/1", {
  skip_if_slow_tests_disabled("cta-confusion-fixture")

  # MINDENOM=56: no-tree → zero rows
  df56 <- cta_confusion_table(.ctbl_myeloma_fit(56L))
  expect_equal(nrow(df56), 0L)
  expect_equal(names(df56), c("actual", "predicted", "n"))

  # MINDENOM=30: V17 stump → [[101,34],[30,21]], n=186
  t30 <- .ctbl_myeloma_fit(30L)
  expect_equal(cta_strata(t30), 2L)
  df30 <- cta_confusion_table(t30)
  expect_equal(sum(df30$n), 186L)
  expect_equal(df30$n[df30$actual == 0L & df30$predicted == 0L], 101L)
  expect_equal(df30$n[df30$actual == 0L & df30$predicted == 1L],  34L)
  expect_equal(df30$n[df30$actual == 1L & df30$predicted == 0L],  30L)
  expect_equal(df30$n[df30$actual == 1L & df30$predicted == 1L],  21L)

  # MINDENOM=1: V14→V15 tree → [[146,40],[36,33]], n=255
  t1 <- .ctbl_myeloma_fit(1L)
  expect_equal(cta_strata(t1), 3L)
  df1 <- cta_confusion_table(t1)
  expect_equal(sum(df1$n), 255L)
  expect_equal(df1$n[df1$actual == 0L & df1$predicted == 0L], 146L)
  expect_equal(df1$n[df1$actual == 0L & df1$predicted == 1L],  40L)
  expect_equal(df1$n[df1$actual == 1L & df1$predicted == 0L],  36L)
  expect_equal(df1$n[df1$actual == 1L & df1$predicted == 1L],  33L)
})
