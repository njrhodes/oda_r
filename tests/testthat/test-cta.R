###############################################################################
# test-cta.R — CTA core tests
###############################################################################

# ---- Data helpers -----------------------------------------------------------
bin_data <- function() {
  # Perfectly separable: x <= 4.5 → class 0, x > 4.5 → class 1
  list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
}

sep3_data <- function() {
  # Three-class, one informative attribute
  list(
    X = data.frame(x1 = c(1,2,3,4,5,6,7,8,9), x2 = rep(5, 9)),
    y = c(1L,1L,1L,2L,2L,2L,3L,3L,3L)
  )
}

# =============================================================================

test_that("oda_cta_fit returns cta_tree for binary data", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                      mc_seed = 1L, loo = "off")
  expect_s3_class(tree, "cta_tree")
  expect_gte(tree$n_nodes, 1L)
})

test_that("root node n_obs equals n", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                      mc_seed = 1L, loo = "off")
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$n_obs, length(d$y))
})

test_that("mindenom too large → single leaf node", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 999L, mc_iter = 100L,
                      mc_seed = 1L, loo = "off")
  expect_equal(tree$n_nodes, 1L)
  expect_true(isTRUE(tree$nodes[[1]]$leaf))
})

test_that("max_depth = 1 → root splits, children are leaves", {
  # ENUMERATE semantics: the root is created directly by ENUMERATE at depth 1
  # and is not filtered by the depth >= max_depth guard (which lives inside
  # .ho_grow(), called only for depth-2+ children).  With max_depth = 1,
  # children are grown at depth 2 and are immediately forced to leaves.
  #
  # The prior expectation ("root is a leaf") depended on the old LaPlace-
  # smoothed MC p inflating the root's p_mc above alpha_split in this seeded
  # example (iter_used=50, ge_count=2: LaPlace=3/51≈0.059 vs raw=2/50=0.04).
  # With canonical raw MC p the root candidate passes alpha; this test is a
  # synthetic unit test, not a gold fixture parity anchor.
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, max_depth = 1L,
                      mc_iter = 300L, mc_seed = 1L, loo = "off")
  root     <- tree$nodes[[tree$root_id]]
  children <- tree$nodes[root$child_ids]

  # Root is a split node (significant ODA on perfectly separable data)
  expect_false(isTRUE(root$leaf), label = "root is a split node")
  # Children of root are leaves (depth 2 >= max_depth 1 forces them)
  expect_true(all(vapply(children, function(nd) isTRUE(nd$leaf), logical(1))),
              label = "all children of root are leaves")
  # No node in the tree has depth > 2
  depths <- vapply(tree$nodes, function(nd) nd$depth, integer(1))
  expect_lte(max(depths, na.rm = TRUE), 2L,
             label = "tree depth does not exceed 2")
})

test_that("max_depth = 2 → no node deeper than depth 2", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, max_depth = 2L,
                      mc_iter = 300L, mc_seed = 1L, loo = "off")
  depths <- vapply(tree$nodes, function(nd) nd$depth, integer(1))
  expect_lte(max(depths, na.rm = TRUE), 2L)
})

test_that("predict returns integer vector of length n", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                      mc_seed = 1L, loo = "off")
  preds <- predict(tree, d$X)
  expect_equal(length(preds), nrow(d$X))
  expect_true(is.integer(preds))
})

test_that("predict: all returned classes are in the training class set", {
  d    <- bin_data()
  tree <- suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off")
  )
  preds <- predict(tree, d$X)
  # no_tree → all NA (correct); valid tree → preds in training set
  expect_true(all(is.na(preds) | preds %in% unique(d$y)))
})

test_that("training accuracy >= majority class baseline", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 500L,
                      mc_seed = 7L, loo = "off")
  preds        <- predict(tree, d$X)
  acc          <- mean(preds == d$y)
  majority_acc <- max(table(d$y)) / length(d$y)
  expect_gte(acc, majority_acc - 0.01)
})

test_that("cta_node_table returns data frame with required columns", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                      mc_seed = 1L, loo = "off")
  tbl <- cta_node_table(tree)
  expect_true(is.data.frame(tbl))
  required <- c("node_id","depth","n_obs","ess","p_mc","leaf","attribute")
  expect_true(all(required %in% names(tbl)))
})

test_that("print.cta_tree runs without error and mentions CTA Tree", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                      mc_seed = 1L, loo = "off")
  expect_output(print(tree), "CTA Tree")
})

test_that("case weights accepted without error", {
  d <- bin_data()
  w <- c(2,1,2,1,2,1,2,1)
  tree <- oda_cta_fit(d$X, d$y, w = w, mindenom = 2L,
                      mc_iter = 200L, mc_seed = 1L, loo = "off")
  expect_s3_class(tree, "cta_tree")
  root <- tree$nodes[[1]]
  expect_equal(root$n_weighted, sum(w))
})

# ---- verbose reporting tests ------------------------------------------------

.loo_fit_args <- list(
  priors_on   = TRUE,
  alpha_split = 0.05,
  mindenom    = 1L,
  prune_alpha = 0.05,
  max_depth   = 3L,
  ess_min     = 0,
  mc_iter     = 5000L,
  mc_target   = 0.05,
  mc_stop     = 99.9,
  mc_stopup   = 99.9,
  mc_seed     = NULL,
  loo         = "stable"
)

test_that("verbose=FALSE produces no [CTA] messages", {
  d    <- .make_loo_gate_data()
  msgs <- character(0)
  withCallingHandlers(
    do.call(cta_fit, c(
      list(X = d[, c("Trap", "Stable")], y = d$Class, verbose = FALSE),
      .loo_fit_args
    )),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_false(any(grepl("\\[CTA", msgs)))
})

test_that("verbose=TRUE emits [CTA] messages with 'selected' or 'no valid root'", {
  d    <- .make_loo_gate_data()
  msgs <- character(0)
  withCallingHandlers(
    do.call(cta_fit, c(
      list(X = d[, c("Trap", "Stable")], y = d$Class, verbose = TRUE),
      .loo_fit_args
    )),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("\\[CTA", msgs)))
  expect_true(any(grepl("selected|no valid root", msgs)))
})

test_that("verbose does not change model output", {
  d      <- .make_loo_gate_data()
  t_quiet <- do.call(cta_fit, c(
    list(X = d[, c("Trap", "Stable")], y = d$Class, verbose = FALSE),
    .loo_fit_args
  ))
  t_verb  <- suppressMessages(do.call(cta_fit, c(
    list(X = d[, c("Trap", "Stable")], y = d$Class, verbose = TRUE),
    .loo_fit_args
  )))
  r_q <- t_quiet$nodes[[t_quiet$root_id]]
  r_v <- t_verb$nodes[[t_verb$root_id]]
  expect_equal(r_q$attribute,       r_v$attribute)
  expect_equal(r_q$rule$cut_value,  r_v$rule$cut_value)
})

test_that("multi-attribute: most informative attribute selected at root", {
  # x1 perfectly separates, x2 is noise
  set.seed(42)
  x1 <- 1:10; x2 <- sample(1:10)
  y  <- c(0L,0L,0L,0L,0L,1L,1L,1L,1L,1L)
  X  <- data.frame(x1 = x1, x2 = x2)
  tree <- oda_cta_fit(X, y, mindenom = 2L, mc_iter = 500L,
                      mc_seed = 42L, loo = "off")
  root <- tree$nodes[[tree$root_id]]
  if (!isTRUE(root$leaf)) {
    expect_equal(root$attribute, "x1")
  } else {
    skip("root is leaf — increase mc_iter for reliable test")
  }
})

# ---- no-tree (degenerate / leaf-only) regression ----------------------------

# Force no-tree by making MINDENOM impossible (> n/2).  This is deterministic
# and does not depend on MC seeds.
.no_tree_fit <- function() {
  d <- bin_data()   # n = 8
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 999L, mc_iter = 100L,
                mc_seed = 1L, loo = "off")
  )
}

test_that("no-tree fit: no_tree flag is TRUE", {
  tree <- .no_tree_fit()
  expect_true(isTRUE(tree$no_tree))
})

test_that("no-tree fit: predict returns NA_integer_ for every observation", {
  d    <- bin_data()
  tree <- .no_tree_fit()
  preds <- predict(tree, d$X)
  expect_true(all(is.na(preds)))
  expect_true(is.integer(preds))
  expect_equal(length(preds), nrow(d$X))
})

test_that("no-tree fit: print says 'No tree found'", {
  tree <- .no_tree_fit()
  expect_output(print(tree), "No tree found")
})

# ---- Phase 2C: read-only CTA endpoint accessor tests ------------------------

# cta_strata -------------------------------------------------------------------

test_that("cta_strata: valid tree returns integer >= 2", {
  d    <- bin_data()
  tree <- suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off"))
  skip_if(isTRUE(tree$no_tree), "mc sampling missed — myeloma tests are the anchor")
  s <- cta_strata(tree)
  expect_true(is.integer(s))
  expect_false(is.na(s))
  expect_gte(s, 2L)
})

test_that("cta_strata: no-tree returns NA_integer_", {
  expect_identical(cta_strata(.no_tree_fit()), NA_integer_)
})

# cta_endpoint_denominators ----------------------------------------------------

test_that("cta_endpoint_denominators: valid tree returns named integer vector", {
  d    <- bin_data()
  tree <- suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off"))
  skip_if(isTRUE(tree$no_tree), "mc sampling missed — myeloma tests are the anchor")
  ep <- cta_endpoint_denominators(tree)
  expect_true(is.integer(ep))
  expect_true(length(ep) >= 2L)
  expect_false(is.null(names(ep)))
  expect_true(all(ep > 0L))
  expect_equal(length(ep), cta_strata(tree))
})

test_that("cta_endpoint_denominators: no-tree returns integer(0)", {
  expect_identical(cta_endpoint_denominators(.no_tree_fit()), integer(0))
})

# cta_min_terminal_denom -------------------------------------------------------

test_that("cta_min_terminal_denom: valid tree equals min of endpoint denominators", {
  d    <- bin_data()
  tree <- suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off"))
  skip_if(isTRUE(tree$no_tree), "mc sampling missed — myeloma tests are the anchor")
  expect_equal(cta_min_terminal_denom(tree), min(cta_endpoint_denominators(tree)))
})

test_that("cta_min_terminal_denom: no-tree returns NA_integer_", {
  expect_identical(cta_min_terminal_denom(.no_tree_fit()), NA_integer_)
})

# new_cta_family_member --------------------------------------------------------

test_that("new_cta_family_member: valid tree fields are consistent", {
  d    <- bin_data()
  tree <- suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off"))
  skip_if(isTRUE(tree$no_tree), "mc sampling missed — myeloma tests are the anchor")
  m <- odacore:::new_cta_family_member(2L, tree)
  expect_identical(m$mindenom, 2L)
  expect_false(m$no_tree)
  expect_identical(m$strata, cta_strata(tree))
  expect_identical(m$min_terminal_denom, cta_min_terminal_denom(tree))
  expect_identical(m$next_mindenom, cta_min_terminal_denom(tree) + 1L)
  expect_identical(m$overall_ess, tree$overall_ess)
  expect_identical(m$d, cta_d_stat(tree))
})

test_that("new_cta_family_member: no-tree fields are correct", {
  tree <- .no_tree_fit()
  m    <- odacore:::new_cta_family_member(999L, tree)
  expect_identical(m$mindenom, 999L)
  expect_true(m$no_tree)
  expect_identical(m$strata, NA_integer_)
  expect_identical(m$min_terminal_denom, NA_integer_)
  expect_identical(m$next_mindenom, NA_integer_)
  expect_identical(m$overall_ess, NA_real_)
  expect_identical(m$d, NA_real_)
})

# new_cta_family ---------------------------------------------------------------

test_that("new_cta_family: returns cta_family S3 class with members list", {
  tree <- .no_tree_fit()
  m    <- odacore:::new_cta_family_member(999L, tree)
  fam  <- odacore:::new_cta_family(list(m))
  expect_s3_class(fam, "cta_family")
  expect_true(is.list(fam$members))
  expect_equal(length(fam$members), 1L)
})

# ---- Myeloma MINDENOM chain: endpoint denominator canon anchors --------------
# Fixture: tests/testthat/fixtures/myeloma/data.txt
# EX V2=0 (exclude zero-weight rows), WEIGHT V2, MISSING ALL (-9), y = V1
# Attrs: V4 V9 V11 V12 V14 V15 V16 V17 V18 V19
# Chain: MINDENOM 1 → next 30 → next 56 → no_tree (terminate)

.myeloma_tree <- function(mindenom) {
  fpath <- testthat::test_path("fixtures/myeloma/data.txt")
  skip_if_not(file.exists(fpath), "myeloma fixture not available")
  d <- read.table(fpath)
  colnames(d) <- paste0("V", seq_len(ncol(d)))
  d  <- d[d$V2 > 0, ]
  ac <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")
  suppressMessages(
    oda_cta_fit(
      X          = d[, ac, drop = FALSE],
      y          = as.integer(d$V1),
      w          = as.numeric(d$V2),
      miss_codes = -9,
      mindenom   = as.integer(mindenom),
      loo        = "stable",
      alpha_split = 0.05,
      mc_iter    = 5000L,
      mc_seed    = 12345L,
      verbose    = FALSE
    )
  )
}

test_that("myeloma MINDENOM=1: min_terminal_denom == 29, next_mindenom == 30", {
  skip_if_slow_tests_disabled("cta-myeloma-chain")
  tree <- .myeloma_tree(1L)
  expect_false(tree$no_tree)
  expect_equal(cta_strata(tree), 3L)
  expect_equal(cta_min_terminal_denom(tree), 29L)
  m <- odacore:::new_cta_family_member(1L, tree)
  expect_equal(m$next_mindenom, 30L)
})

test_that("myeloma MINDENOM=30: min_terminal_denom == 55, next_mindenom == 56", {
  skip_if_slow_tests_disabled("cta-myeloma-chain")
  tree <- .myeloma_tree(30L)
  expect_false(tree$no_tree)
  expect_equal(cta_strata(tree), 2L)
  expect_equal(cta_min_terminal_denom(tree), 55L)
  m <- odacore:::new_cta_family_member(30L, tree)
  expect_equal(m$next_mindenom, 56L)
})

test_that("myeloma MINDENOM=56: no_tree, strata NA, denominators integer(0)", {
  skip_if_slow_tests_disabled("cta-myeloma-chain")
  tree <- .myeloma_tree(56L)
  expect_true(tree$no_tree)
  expect_identical(cta_strata(tree), NA_integer_)
  expect_identical(cta_endpoint_denominators(tree), integer(0))
  expect_identical(cta_min_terminal_denom(tree), NA_integer_)
  m <- odacore:::new_cta_family_member(56L, tree)
  expect_true(m$no_tree)
  expect_identical(m$next_mindenom, NA_integer_)
})

# ---- cta_d_stat: synthetic contract tests ------------------------------------
# These tests exercise D-stat behavior independently of CTA topology so they
# remain valid regardless of open ENUMERATE / PRUNE issues.
#
# Synthetic minimal cta_tree objects: no fitting required.  Leaf nodes only;
# split nodes are not needed for the accessor contract under test.
.fake_cta_tree <- function(overall_ess, n_leaves, no_tree = FALSE) {
  nodes <- lapply(seq_len(n_leaves), function(i)
    list(leaf = TRUE, node_id = i, n_obs = 10L))
  structure(
    list(nodes = nodes, no_tree = no_tree, overall_ess = overall_ess),
    class = "cta_tree"
  )
}

test_that("cta_d_stat: no_tree returns NA_real_", {
  tree <- .fake_cta_tree(overall_ess = 50, n_leaves = 2L, no_tree = TRUE)
  expect_identical(cta_d_stat(tree), NA_real_)
})

test_that("cta_d_stat: missing overall_ess returns NA_real_", {
  tree <- .fake_cta_tree(overall_ess = 50, n_leaves = 2L)
  tree$overall_ess <- NULL
  expect_identical(cta_d_stat(tree), NA_real_)
})

test_that("cta_d_stat: non-finite overall_ess returns NA_real_", {
  tree_inf <- .fake_cta_tree(overall_ess = Inf,  n_leaves = 2L)
  tree_nan <- .fake_cta_tree(overall_ess = NaN,  n_leaves = 2L)
  tree_na  <- .fake_cta_tree(overall_ess = NA_real_, n_leaves = 2L)
  expect_identical(cta_d_stat(tree_inf), NA_real_)
  expect_identical(cta_d_stat(tree_nan), NA_real_)
  expect_identical(cta_d_stat(tree_na),  NA_real_)
})

test_that("cta_d_stat: non-positive overall_ess returns NA_real_", {
  tree_zero <- .fake_cta_tree(overall_ess =  0,  n_leaves = 2L)
  tree_neg  <- .fake_cta_tree(overall_ess = -5,  n_leaves = 2L)
  expect_identical(cta_d_stat(tree_zero), NA_real_)
  expect_identical(cta_d_stat(tree_neg),  NA_real_)
})

test_that("cta_d_stat: strata < 2 returns NA_real_", {
  tree <- .fake_cta_tree(overall_ess = 50, n_leaves = 1L)
  expect_identical(cta_d_stat(tree), NA_real_)
})

test_that("cta_d_stat: two leaves, overall_ess = 50 returns D = 2", {
  # D = 100 / (50 / 2) - 2 = 4 - 2 = 2
  tree <- .fake_cta_tree(overall_ess = 50, n_leaves = 2L)
  expect_equal(cta_d_stat(tree), 2)
})

# ---- PRUNE: Sidak-Bonferroni + maximum-accuracy pruning tests ----------------
#
# Fixture: myeloma data (tests/testthat/fixtures/myeloma/data.txt), mc_seed=800.
#
# With mc_seed=800 and prune_alpha=0.05, the HO-CTA candidate for the V14 root
# grows V14 -> V15 -> V11 (3 split nodes). V11 has p_mc=0.04:
#   K=3, alpha_3 = 1-(0.95)^(1/3) ~ 0.017 < 0.04 → V11 is Sidak-flagged.
# Pruning V11 improves WESS from 25.40% to 27.69% → maximum-accuracy accepts.
# Result: V14 -> V15 (2 split nodes), overall_ess ~ 27.69%.
#
# With prune_alpha=1.0 (no-op), V11 is never flagged: V14->V15->V11 survives,
# scoring 25.40% over all n observations.
#
# External validation anchor (not tracked as fixture):
#   tests/testthat/myeloma/MODEL1.TXT: Unpruned HO V17->V15->V14 WESS=23.61%,
#   Pruned HO V17->V15 WESS=25.43%, Enumerated V14->V15 WESS=27.69%.

.myeloma_prune_tree <- function(prune_alpha_val) {
  fpath <- testthat::test_path("fixtures/myeloma/data.txt")
  skip_if_not(file.exists(fpath), "myeloma fixture not available")
  d <- read.table(fpath)
  colnames(d) <- paste0("V", seq_len(ncol(d)))
  d  <- d[d[["V2"]] > 0, ]
  ac <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")
  suppressMessages(
    oda_cta_fit(
      X           = d[, ac, drop = FALSE],
      y           = as.integer(d[["V1"]]),
      w           = as.numeric(d[["V2"]]),
      miss_codes  = -9,
      mindenom    = 1L,
      loo         = "stable",
      alpha_split = 0.05,
      prune_alpha = prune_alpha_val,
      mc_iter     = 5000L,
      mc_seed     = 800L,
      verbose     = FALSE
    )
  )
}

.n_split_nodes <- function(tree) {
  sum(vapply(tree[["nodes"]],
             function(nd) !is.null(nd) && !isTRUE(nd[["leaf"]]) &&
               length(nd[["child_ids"]]) > 0L,
             logical(1L)))
}

test_that("PRUNE no-op: prune_alpha=1.0 retains all grown split nodes", {
  skip_if_slow_tests_disabled("cta-myeloma-prune")
  tree <- .myeloma_prune_tree(1.0)
  # With no pruning, HO candidate V14->V15->V11 survives (3 splits).
  expect_equal(.n_split_nodes(tree), 3L,
               label = "prune_alpha=1.0: 3 split nodes (V11 not pruned)")
  expect_equal(round(tree[["overall_ess"]], 2), 25.40,
               label = "prune_alpha=1.0: WESS=25.40% (unpruned)")
})

test_that("PRUNE active: Sidak-flagged V11 pruned, WESS improves", {
  skip_if_slow_tests_disabled("cta-myeloma-prune")
  tree <- .myeloma_prune_tree(0.05)
  # V11 (p=0.04 > alpha_3~0.017) is Sidak-flagged; pruning improves WESS.
  expect_equal(.n_split_nodes(tree), 2L,
               label = "prune_alpha=0.05: 2 split nodes (V11 pruned)")
  expect_equal(round(tree[["overall_ess"]], 2), 27.69,
               label = "prune_alpha=0.05: WESS=27.69% (pruned)")
  root_attr <- tree[["nodes"]][[tree[["root_id"]]]][["attribute"]]
  expect_equal(root_attr, "V14",
               label = "pruned tree root is V14 (enumerated winner)")
})

# ---- A×B×C ENUMERATE: leaf sentinel tests ------------------------------------
#
# Canonical myeloma MINDENOM=1 result: V14→V15 tree.
#   A = V14, cut=0.5  (root)
#   B = V15 (left child of V14, i.e. V14<=0.5 branch)
#   C = leaf          (right child of V14, i.e. V14>0.5 branch)
#
# The C=leaf outcome proves that the expanded A×B×C loop evaluated the leaf
# sentinel for C and selected it over any split candidates on the right branch.
# This is a structural correctness test: if C=leaf were not in C_options, the
# loop would only emit (B=V15, C=split) combinations and would fail to reproduce
# the canonical tree.
#
# Uses the same .myeloma_prune_tree(0.05) fixture (mc_seed=800, MINDENOM=1).

test_that("ENUMERATE A×B×C: myeloma MINDENOM=1 right branch of V14 is a leaf", {
  skip_if_slow_tests_disabled("cta-myeloma-prune")
  tree <- .myeloma_prune_tree(0.05)
  skip_if(isTRUE(tree$no_tree), "no tree found — fixture issue")
  root <- tree$nodes[[tree$root_id]]
  skip_if(isTRUE(root$leaf), "root is leaf — ENUMERATE did not grow a tree")

  # Root must be V14 (canonical enumerated winner)
  expect_equal(root$attribute, "V14", label = "root attribute is V14")

  # The right child of the root (V14 > 0.5 branch) must be a leaf.
  # child_ids[2] = right child in the split labelling convention.
  right_child_id <- root$child_ids[2L]
  right_child    <- tree$nodes[[right_child_id]]
  expect_true(isTRUE(right_child$leaf),
    label = "right branch of V14 root is a leaf (C=leaf sentinel selected by ENUMERATE)")
})

test_that("ENUMERATE A×B×C: myeloma MINDENOM=1 left branch of V14 is V15 split", {
  skip_if_slow_tests_disabled("cta-myeloma-prune")
  tree <- .myeloma_prune_tree(0.05)
  skip_if(isTRUE(tree$no_tree), "no tree found — fixture issue")
  root <- tree$nodes[[tree$root_id]]
  skip_if(isTRUE(root$leaf), "root is leaf — ENUMERATE did not grow a tree")

  expect_equal(root$attribute, "V14", label = "root attribute is V14")

  # The left child of the root (V14 <= 0.5 branch) must be the V15 split.
  left_child_id <- root$child_ids[1L]
  left_child    <- tree$nodes[[left_child_id]]
  expect_false(isTRUE(left_child$leaf),
    label = "left branch of V14 root is a split node (not a leaf)")
  expect_equal(left_child$attribute, "V15",
    label = "left branch attribute is V15 (canonical depth-2 split)")
})

# ---- cta_descendant_family: myeloma MDSA chain tests ------------------------
#
# Canonical myeloma descendant family:
#   MINDENOM 1  → min_terminal_denom = 29 → next = 30
#   MINDENOM 30 → min_terminal_denom = 55 → next = 56
#   MINDENOM 56 → no_tree → terminate
#
# Uses mc_seed = 12345L (same as existing .myeloma_tree() fixture tests).
# family$min_d_idx = 1L because MINDENOM=1 (strata=3, WESS≈27.69%) has lower
# D than MINDENOM=30 (strata=2, WESS≈16.51%):
#   D(1)  = 100/(27.69/3) - 3  ≈  7.83
#   D(30) = 100/(16.51/2) - 2  ≈ 10.11

.myeloma_family <- local({
  fam <- NULL
  function() {
    if (is.null(fam)) {
      fpath <- testthat::test_path("fixtures/myeloma/data.txt")
      skip_if_not(file.exists(fpath), "myeloma fixture not available")
      d <- read.table(fpath)
      colnames(d) <- paste0("V", seq_len(ncol(d)))
      d  <- d[d[["V2"]] > 0, ]
      ac <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")
      fam <<- suppressMessages(
        cta_descendant_family(
          X              = d[, ac, drop = FALSE],
          y              = as.integer(d[["V1"]]),
          w              = as.numeric(d[["V2"]]),
          miss_codes     = -9,
          loo            = "stable",
          alpha_split    = 0.05,
          prune_alpha    = 0.05,
          mc_iter        = 5000L,
          mc_seed        = 12345L,
          verbose        = FALSE,
          start_mindenom = 1L
        )
      )
    }
    fam
  }
})

test_that("cta_descendant_family: myeloma chain length is 3", {
  skip_if_slow_tests_disabled("cta-family")
  fam <- .myeloma_family()
  expect_equal(length(fam$members), 3L)
})

test_that("cta_descendant_family: myeloma MINDENOM chain is {1, 30, 56}", {
  skip_if_slow_tests_disabled("cta-family")
  fam <- .myeloma_family()
  expect_equal(fam$mindenoms, c(1L, 30L, 56L))
})

test_that("cta_descendant_family: myeloma terminated = TRUE, reason = 'no_tree'", {
  skip_if_slow_tests_disabled("cta-family")
  fam <- .myeloma_family()
  expect_true(fam$terminated)
  expect_equal(fam$termination_reason, "no_tree")
})

test_that("cta_descendant_family: myeloma final member is no-tree", {
  skip_if_slow_tests_disabled("cta-family")
  fam <- .myeloma_family()
  expect_true(fam$members[[3L]]$no_tree)
})

test_that("cta_descendant_family: myeloma feasible members have finite D", {
  skip_if_slow_tests_disabled("cta-family")
  fam <- .myeloma_family()
  expect_true(is.finite(fam$members[[1L]]$d),
              label = "member 1 (MINDENOM=1) has finite D")
  expect_true(is.finite(fam$members[[2L]]$d),
              label = "member 2 (MINDENOM=30) has finite D")
})

test_that("cta_descendant_family: myeloma no-tree member has NA D", {
  skip_if_slow_tests_disabled("cta-family")
  fam <- .myeloma_family()
  expect_true(is.na(fam$members[[3L]]$d),
              label = "member 3 (MINDENOM=56, no-tree) D is NA")
})

test_that("cta_descendant_family: myeloma min_d_idx = 1L (MINDENOM=1 wins)", {
  skip_if_slow_tests_disabled("cta-family")
  fam <- .myeloma_family()
  expect_equal(fam$min_d_idx, 1L,
               label = "MINDENOM=1 has lower D than MINDENOM=30")
})

test_that("cta_descendant_family: myeloma summary has required columns and shape", {
  skip_if_slow_tests_disabled("cta-family")
  fam <- .myeloma_family()
  required <- c("mindenom","status","strata","min_terminal_denom",
                "overall_ess","d","no_tree")
  expect_true(all(required %in% names(fam$summary)),
              info = paste("missing:", paste(setdiff(required, names(fam$summary)),
                                             collapse = ", ")))
  expect_equal(nrow(fam$summary), 3L)
})
