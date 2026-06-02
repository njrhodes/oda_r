###############################################################################
# test-fixture-cta-demo-model8.R
#
# Canon coverage for CTA_DEMO MINDENOM=8 — the only active-PRUNE fixture
# across all current EXE canon files.
#
# EXE command (fixtures/cta_demo/cta_8.pgm):
#   OPEN CTA_DEMO.CSV;
#   VARS V1 TO V6; CLASS V1; ATTRIBUTE V2 V3 V4 V5 V6;
#   MC ITER 25000 CUTOFF 0.05 STOP 99.9;
#   MINDENOM 8; PRUNE 0.05; ENUMERATE; LOO STABLE;
#
# Canon (from MODEL8.TXT):
#
#   Unpruned tree (7 split nodes, OVERALL ESS = 60.65%):
#     Node  1 V2  lev=1  n=200
#     Node  2 V6  lev=2  n=110
#     Node  3 V2  lev=2  n=90
#     Node  4 V3  lev=3  n=29
#     Node  5 V5  lev=3  n=81   <- PRUNE target
#     Node  6 V6  lev=3  n=47
#     Node 13 V4  lev=4  n=37
#
#   Pruned tree (6 split nodes, OVERALL ESS = 68.08%):
#     Node  1 V2  lev=1  n=200
#     Node  2 V6  lev=2  n=110
#     Node  3 V2  lev=2  n=90
#     Node  4 V3  lev=3  n=29
#     Node  6 V6  lev=3  n=47
#     Node 13 V4  lev=4  n=37
#     (V5 / Node 5 removed; global ESS improved by 7.43pp)
#
#   Enumerated tree == Pruned tree (OVERALL ESS = 68.08%)
#   Confusion (Pruned/Enumerated): [[104,20],[12,64]], n=200
#
# Test structure:
#   1. Static canon fixture — hard-coded EXE anchors; no R fitting.
#      Verifies the fixture file records the active-PRUNE event.
#   2. R end-to-end MODEL8 — fits with mc_iter=25000 and asserts
#      prune_info and tree structure against EXE canon.
#      Gated: full tier (mc_iter=25000; EXE-grade parity test).
#
# Never rewrite these expected values from R output.  They derive solely
# from MODEL8.TXT (EXE-produced).
###############################################################################

# ---- Static canon fixture test -------------------------------------------
# These assertions encode what MODEL8.TXT says, not what R produces.
# If any of these fail the fixture file has been altered.

test_that("MODEL8 static canon: unpruned tree contains V5 node 5 as a split", {
  # Unpruned split attributes from MODEL8.TXT (in node-ID order)
  EXE_UNPRUNED_SPLIT_ATTRS <- c("V2", "V6", "V2", "V3", "V5", "V6", "V4")
  EXE_UNPRUNED_NODE_IDS    <- c(1L,   2L,   3L,   4L,   5L,   6L,  13L)
  EXE_UNPRUNED_ESS         <- 60.65

  expect_length(EXE_UNPRUNED_SPLIT_ATTRS, 7L)
  expect_true("V5" %in% EXE_UNPRUNED_SPLIT_ATTRS, label = "V5 in unpruned")
  expect_equal(EXE_UNPRUNED_NODE_IDS[EXE_UNPRUNED_SPLIT_ATTRS == "V5"], 5L,
               label = "V5 is node 5 in unpruned")
  expect_equal(EXE_UNPRUNED_ESS, 60.65)
})

test_that("MODEL8 static canon: pruned tree removes V5 node 5; ESS improves to 68.08", {
  EXE_UNPRUNED_SPLIT_ATTRS <- c("V2", "V6", "V2", "V3", "V5", "V6", "V4")
  EXE_PRUNED_SPLIT_ATTRS   <- c("V2", "V6", "V2", "V3",       "V6", "V4")
  EXE_UNPRUNED_ESS         <- 60.65
  EXE_PRUNED_ESS           <- 68.08

  expect_length(EXE_PRUNED_SPLIT_ATTRS, 6L)
  expect_false("V5" %in% EXE_PRUNED_SPLIT_ATTRS, label = "V5 absent from pruned")
  expect_equal(setdiff(EXE_UNPRUNED_SPLIT_ATTRS, EXE_PRUNED_SPLIT_ATTRS), "V5",
               label = "only V5 removed by pruning")
  expect_gt(EXE_PRUNED_ESS, EXE_UNPRUNED_ESS,
            label  = "pruning improves ESS")
  expect_equal(EXE_PRUNED_ESS - EXE_UNPRUNED_ESS, 7.43, tolerance = 0.01,
               label = "pruning gain = 7.43pp")
})

test_that("MODEL8 static canon: enumerated tree equals pruned tree", {
  EXE_PRUNED_SPLIT_ATTRS     <- c("V2", "V6", "V2", "V3", "V6", "V4")
  EXE_ENUMERATED_SPLIT_ATTRS <- c("V2", "V6", "V2", "V3", "V6", "V4")
  EXE_PRUNED_ESS             <- 68.08
  EXE_ENUMERATED_ESS         <- 68.08

  expect_identical(EXE_PRUNED_SPLIT_ATTRS, EXE_ENUMERATED_SPLIT_ATTRS)
  expect_equal(EXE_PRUNED_ESS, EXE_ENUMERATED_ESS)
})

test_that("MODEL8 static canon: pruned confusion = [[104,20],[12,64]], n=200", {
  # Confusion matrix from MODEL8.TXT Pruned section
  EXE_CONF <- matrix(c(104L, 12L, 20L, 64L), nrow = 2L, ncol = 2L)
  expect_equal(sum(EXE_CONF), 200L, label = "n=200")
  expect_equal(EXE_CONF[1, 1], 104L); expect_equal(EXE_CONF[1, 2], 20L)
  expect_equal(EXE_CONF[2, 1],  12L); expect_equal(EXE_CONF[2, 2], 64L)
})

# ---- R end-to-end MODEL8 test -----------------------------------------------
# Fits CTA_DEMO with MINDENOM=8 and verifies prune_info + tree structure.
# Gated at full tier: mc_iter=25000 matches EXE canon.

.load_cta_demo_model8 <- function() {
  f <- tryCatch(testthat::test_path("fixtures/cta_demo/CTA_DEMO.CSV"),
                error = function(e) "")
  if (!nzchar(f) || !file.exists(f)) testthat::skip("CTA_DEMO.CSV not found")
  read.csv(f, header = FALSE,
           col.names = c("V1", "V2", "V3", "V4", "V5", "V6"))
}

.cta_demo_model8_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      d <- .load_cta_demo_model8()
      fit <<- oda_cta_fit(
        X           = d[, c("V2", "V3", "V4", "V5", "V6")],
        y           = d$V1,
        priors_on   = TRUE,
        alpha_split = 0.05,
        mindenom    = 8L,
        prune_alpha = 0.05,
        max_depth   = 20L,
        ess_min     = 0,
        mc_iter     = 25000L,
        mc_target   = 0.05,
        mc_stop     = 99.9,
        mc_stopup   = 99.9,
        mc_seed     = 42L,    # fixed seed: deterministic pruning for canon tests
        loo         = "stable",
        attr_names  = c("V2", "V3", "V4", "V5", "V6")
      )
    }
    fit
  }
})

test_that("MODEL8 R e2e: prune_info is non-NULL and records V5 removal", {
  skip_if_not_full("fixture-cta-demo-model8")
  tree <- .cta_demo_model8_fit()

  expect_false(isTRUE(tree$no_tree), label = "tree found")
  expect_false(is.null(tree$prune_info),
               label = "prune_info is attached to cta_tree")

  pi <- tree$prune_info
  expect_equal(pi$unpruned_n_splits, 7L,
               label = "7 split nodes before pruning")
  expect_equal(pi$pruned_n_splits,   6L,
               label = "6 split nodes after pruning")
  expect_equal(length(pi$removed_node_ids), 1L,
               label = "exactly one node removed")
  expect_equal(pi$removed_node_ids, 5L,
               label = "removed node ID is 5 (canonical EXE node 5 = V5)")
  expect_true("V5" %in% pi$removed_attrs,
              label = "removed attr is V5")
})

test_that("MODEL8 R e2e: unpruned ESS ≈ 60.65, pruned ESS ≈ 68.08", {
  skip_if_not_full("fixture-cta-demo-model8")
  tree <- .cta_demo_model8_fit()
  pi   <- tree$prune_info

  expect_equal(round(pi$unpruned_ess, 2), 60.65, tolerance = 0.5,
               label = "unpruned ESS ≈ 60.65%")
  expect_equal(round(pi$pruned_ess,   2), 68.08, tolerance = 0.5,
               label = "pruned ESS ≈ 68.08%")
  expect_equal(pi$pruned_ess, tree$overall_ess,
               label = "pruned_ess == overall_ess")
})

test_that("MODEL8 R e2e: pruned tree has EXE-canonical node_id:attr (V2 V6 V2 V3 V6 V4)", {
  # EXE canonical pruned tree (MODEL8.TXT Pruned section):
  #   node 1=V2, node 2=V6, node 3=V2, node 4=V3, node 6=V6, node 13=V4
  #   node 5 removed (was V5); all other positions are leaves.
  # Canonical geometry (MPE Appendix C, Figure C.1): left child = 2x, right = 2x+1.
  skip_if_not_full("fixture-cta-demo-model8")
  tree  <- .cta_demo_model8_fit()
  nodes <- tree$nodes

  .is_split <- function(nd) !is.null(nd) && !isTRUE(nd$leaf) && length(nd$child_ids) > 0L

  # Exact node_id:attribute for every split in the pruned tree
  expect_true(.is_split(nodes[[1L]]),           label = "node 1 is a split node")
  expect_equal(nodes[[1L]]$attribute, "V2",     label = "node 1 attr = V2")
  expect_true(.is_split(nodes[[2L]]),           label = "node 2 is a split node")
  expect_equal(nodes[[2L]]$attribute, "V6",     label = "node 2 attr = V6")
  expect_true(.is_split(nodes[[3L]]),           label = "node 3 is a split node")
  expect_equal(nodes[[3L]]$attribute, "V2",     label = "node 3 attr = V2")
  expect_true(.is_split(nodes[[4L]]),           label = "node 4 is a split node")
  expect_equal(nodes[[4L]]$attribute, "V3",     label = "node 4 attr = V3")
  # Node 5 was V5 — pruned to leaf
  expect_false(is.null(nodes[[5L]]),            label = "node 5 exists in node list")
  expect_true(isTRUE(nodes[[5L]]$leaf),         label = "node 5 is leaf (V5 pruned)")
  expect_true(.is_split(nodes[[6L]]),           label = "node 6 is a split node")
  expect_equal(nodes[[6L]]$attribute, "V6",     label = "node 6 attr = V6")
  expect_true(.is_split(nodes[[13L]]),          label = "node 13 is a split node")
  expect_equal(nodes[[13L]]$attribute, "V4",    label = "node 13 attr = V4")

  # Exactly 6 split nodes total
  split_ids <- sort(which(vapply(nodes, .is_split, logical(1L))))
  expect_equal(split_ids, c(1L, 2L, 3L, 4L, 6L, 13L),
               label = "split node IDs are exactly {1,2,3,4,6,13}")
})

test_that("MODEL8 R e2e: overall ESS ≈ 68.08% (tolerance 0.5pp)", {
  skip_if_not_full("fixture-cta-demo-model8")
  tree <- .cta_demo_model8_fit()

  expect_equal(round(tree$overall_ess, 2), 68.08, tolerance = 0.5,
               label = "overall ESS ≈ 68.08%")
})

test_that("MODEL8 R e2e: pruned confusion = [[104,20],[12,64]] (tolerance 2)", {
  skip_if_not_full("fixture-cta-demo-model8")
  tree <- .cta_demo_model8_fit()

  conf <- tree$training_confusion

  expect_equal(conf[1, 1], 104L, tolerance = 2L, label = "TP class 1")
  expect_equal(conf[1, 2],  20L, tolerance = 2L, label = "FN class 1")
  expect_equal(conf[2, 1],  12L, tolerance = 2L, label = "FP class 2")
  expect_equal(conf[2, 2],  64L, tolerance = 2L, label = "TP class 2")
})

test_that("MODEL8 R e2e: prune_info removed_attrs is exactly 'V5' (exact equality, parallel to removed_node_ids)", {
  # Existing test checks "V5" %in% pi$removed_attrs but not exact equality.
  # This test locks: removed_attrs == c("V5"), length == 1, parallel to removed_node_ids.
  skip_if_not_full("fixture-cta-demo-model8")
  tree <- .cta_demo_model8_fit()
  pi   <- tree$prune_info

  expect_equal(length(pi$removed_attrs), 1L,
               label = "exactly one removed_attrs entry")
  expect_equal(pi$removed_attrs, "V5",
               label = "removed_attrs is exactly 'V5'")
  expect_equal(length(pi$removed_attrs), length(pi$removed_node_ids),
               label = "removed_attrs and removed_node_ids are parallel (same length)")
})

# ---- Canonical node geometry invariant tests ---------------------------------
# These tests lock down the Appendix C node-ID geometry and ordered-cut branch
# direction (LEFT = x<=cut) introduced in the fix.  All expected values are
# derived directly from MODEL8.TXT (EXE-produced); do not update from R output.
#
# Split nodes and their EXE-canonical branch sizes (MODEL8.TXT Pruned section):
#   Node  1  V2 <=4.5 : left(node  2) n=110, right(node  3) n= 90
#   Node  2  V6 <=7.5 : left(node  4) n= 29, right(node  5) n= 81
#   Node  3  V2 <=6.5 : left(node  6) n= 47, right(node  7) n= 43
#   Node  4  V3 <=4.5 : left(node  8) n= 16, right(node  9) n= 13
#   Node  6  V6 <=5.5 : left(node 12) n= 10, right(node 13) n= 37
#   Node 13  V4 <=0.5 : left(node 26) n= 19, right(node 27) n= 18

test_that("MODEL8 R e2e: every split node has child_ids = c(2*nid, 2*nid+1) and correct parent_id", {
  skip_if_not_full("fixture-cta-demo-model8")
  tree  <- .cta_demo_model8_fit()
  nodes <- tree$nodes

  .is_split <- function(nd) !is.null(nd) && !isTRUE(nd$leaf) && length(nd$child_ids) > 0L

  split_ids <- which(vapply(nodes, .is_split, logical(1L)))
  for (nid in split_ids) {
    nd        <- nodes[[nid]]
    left_id   <- 2L * nid
    right_id  <- 2L * nid + 1L
    expected  <- c(left_id, right_id)
    actual    <- sort(nd$child_ids)
    expect_equal(actual, sort(expected),
                 label = paste0("node ", nid, " child_ids = {", left_id, ",", right_id, "}"))

    for (cid in nd$child_ids) {
      child_nd <- nodes[[cid]]
      expect_false(is.null(child_nd),
                   label = paste0("node ", cid, " exists (child of ", nid, ")"))
      expect_equal(child_nd$parent_id, nid,
                   label = paste0("node ", cid, " parent_id = ", nid))
    }
  }
})

test_that("MODEL8 R e2e: ordered_cut left branch = x<=cut_value (EXE canonical obs counts)", {
  # Branch sizes are exact obs counts from MODEL8.TXT — no tolerance.
  # left child  = x <= cut_value (Appendix C left = 2*nid).
  # right child = x >  cut_value (Appendix C right = 2*nid+1).
  skip_if_not_full("fixture-cta-demo-model8")
  tree  <- .cta_demo_model8_fit()
  nodes <- tree$nodes

  .n_obs <- function(nid) {
    nd <- nodes[[nid]]
    if (is.null(nd)) return(NA_integer_)
    nd$n_obs
  }

  # Node 1 (V2 <=4.5): left=node2, right=node3
  expect_equal(.n_obs(2L),  110L, label = "node 2 n_obs = 110 (V2<=4.5 side)")
  expect_equal(.n_obs(3L),   90L, label = "node 3 n_obs = 90  (V2>4.5 side)")

  # Node 2 (V6 <=7.5): left=node4, right=node5
  expect_equal(.n_obs(4L),   29L, label = "node 4 n_obs = 29  (V6<=7.5 side)")
  expect_equal(.n_obs(5L),   81L, label = "node 5 n_obs = 81  (V6>7.5 side)")

  # Node 3 (V2 <=6.5): left=node6, right=node7
  expect_equal(.n_obs(6L),   47L, label = "node 6 n_obs = 47  (V2<=6.5 side)")
  expect_equal(.n_obs(7L),   43L, label = "node 7 n_obs = 43  (V2>6.5 side)")

  # Node 4 (V3 <=4.5): left=node8, right=node9
  expect_equal(.n_obs(8L),   16L, label = "node 8 n_obs = 16  (V3<=4.5 side)")
  expect_equal(.n_obs(9L),   13L, label = "node 9 n_obs = 13  (V3>4.5 side)")

  # Node 6 (V6 <=5.5): left=node12, right=node13
  expect_equal(.n_obs(12L),  10L, label = "node 12 n_obs = 10 (V6<=5.5 side)")
  expect_equal(.n_obs(13L),  37L, label = "node 13 n_obs = 37 (V6>5.5 side)")

  # Node 13 (V4 <=0.5): left=node26, right=node27
  expect_equal(.n_obs(26L),  19L, label = "node 26 n_obs = 19 (V4<=0.5 side)")
  expect_equal(.n_obs(27L),  18L, label = "node 27 n_obs = 18 (V4>0.5 side)")
})
