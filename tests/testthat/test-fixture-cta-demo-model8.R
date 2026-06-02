###############################################################################
# test-fixture-cta-demo-model8.R — CTA_DEMO MINDENOM=8 canon coverage
#
# EXE: MINDENOM 8; PRUNE 0.05; ENUMERATE; LOO STABLE; MC ITER 25000 STOP 99.9.
# Canon from MODEL8.TXT:
#   Unpruned: 7 splits (V2 V6 V2 V3 V5 V6 V4), ESS = 60.65%
#   Pruned:   6 splits (V2 V6 V2 V3 V6 V4), node5/V5 removed, ESS = 68.08%
#   Enumerated == Pruned. Confusion: [[104,20],[12,64]], n=200.
#
# Block 1 (CRAN):  static canon — encodes MODEL8.TXT, no R fitting.
# Block 2 (smoke): R e2e — prune_info fields and ESS values.
# Block 3 (smoke): R e2e — tree geometry (node IDs, branch sizes, confusion).
#
# Never rewrite expected values from R output — all derive from MODEL8.TXT.
###############################################################################

# ---- R e2e helpers ----------------------------------------------------------

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
        mc_seed     = 42L,
        loo         = "stable",
        attr_names  = c("V2", "V3", "V4", "V5", "V6")
      )
    }
    fit
  }
})

# =============================================================================
# 1. Static canon — MODEL8.TXT anchors (CRAN-safe, no R fitting)
# =============================================================================

test_that("MODEL8 static canon: unpruned/pruned/enumerated structure and confusion", {
  EXE_UNPRUNED_SPLIT_ATTRS <- c("V2", "V6", "V2", "V3", "V5", "V6", "V4")
  EXE_PRUNED_SPLIT_ATTRS   <- c("V2", "V6", "V2", "V3",       "V6", "V4")
  EXE_ENUM_SPLIT_ATTRS     <- c("V2", "V6", "V2", "V3",       "V6", "V4")
  EXE_UNPRUNED_ESS         <- 60.65
  EXE_PRUNED_ESS           <- 68.08
  EXE_CONF                 <- matrix(c(104L, 12L, 20L, 64L), nrow = 2L)

  # unpruned: 7 splits, V5 at position 5
  expect_length(EXE_UNPRUNED_SPLIT_ATTRS, 7L)
  expect_true("V5" %in% EXE_UNPRUNED_SPLIT_ATTRS)
  expect_equal(which(EXE_UNPRUNED_SPLIT_ATTRS == "V5"), 5L)

  # pruning removes exactly V5 → 6 splits
  expect_length(EXE_PRUNED_SPLIT_ATTRS, 6L)
  expect_false("V5" %in% EXE_PRUNED_SPLIT_ATTRS)
  expect_equal(setdiff(EXE_UNPRUNED_SPLIT_ATTRS, EXE_PRUNED_SPLIT_ATTRS), "V5")

  # ESS improves 7.43pp; pruned == enumerated
  expect_gt(EXE_PRUNED_ESS, EXE_UNPRUNED_ESS)
  expect_equal(EXE_PRUNED_ESS - EXE_UNPRUNED_ESS, 7.43, tolerance = 0.01)
  expect_identical(EXE_ENUM_SPLIT_ATTRS, EXE_PRUNED_SPLIT_ATTRS)

  # pruned confusion [[104,20],[12,64]], n=200
  expect_equal(sum(EXE_CONF), 200L)
  expect_equal(EXE_CONF[1,1], 104L); expect_equal(EXE_CONF[1,2], 20L)
  expect_equal(EXE_CONF[2,1],  12L); expect_equal(EXE_CONF[2,2], 64L)
})

# =============================================================================
# 2. R e2e — prune_info fields and ESS values
# =============================================================================

test_that("MODEL8 R e2e: prune_info records V5/node5 removal; ESS values match EXE", {
  skip_if_not_smoke("fixture-cta-demo-model8")
  tree <- .cta_demo_model8_fit()
  expect_false(isTRUE(tree$no_tree))

  pi <- tree$prune_info
  expect_false(is.null(pi))
  expect_equal(pi$unpruned_n_splits,        7L)
  expect_equal(pi$pruned_n_splits,          6L)
  expect_equal(length(pi$removed_node_ids), 1L)
  expect_equal(pi$removed_node_ids,         5L)
  expect_equal(pi$removed_attrs,           "V5")
  expect_equal(length(pi$removed_attrs), length(pi$removed_node_ids))

  expect_equal(round(pi$unpruned_ess, 2), 60.65, tolerance = 0.5,
               label = "unpruned ESS ≈ 60.65%")
  expect_equal(round(pi$pruned_ess,   2), 68.08, tolerance = 0.5,
               label = "pruned ESS ≈ 68.08%")
  expect_equal(pi$pruned_ess, tree$overall_ess,
               label = "pruned_ess == overall_ess")
})

# =============================================================================
# 3. R e2e — tree geometry and training performance
# =============================================================================

test_that("MODEL8 R e2e: split IDs, node attrs, child geometry, branch sizes, ESS, confusion", {
  # Split nodes and EXE-canonical branch sizes (MODEL8.TXT Pruned section):
  #   Node  1 V2 <=4.5: left(2) n=110, right(3) n= 90
  #   Node  2 V6 <=7.5: left(4) n= 29, right(5) n= 81
  #   Node  3 V2 <=6.5: left(6) n= 47, right(7) n= 43
  #   Node  4 V3 <=4.5: left(8) n= 16, right(9) n= 13
  #   Node  6 V6 <=5.5: left(12)n= 10, right(13)n= 37
  #   Node 13 V4 <=0.5: left(26)n= 19, right(27)n= 18
  skip_if_not_smoke("fixture-cta-demo-model8")
  tree  <- .cta_demo_model8_fit()
  nodes <- tree$nodes

  .is_split <- function(nd) !is.null(nd) && !isTRUE(nd$leaf) && length(nd$child_ids) > 0L
  .n        <- function(nid) nodes[[nid]]$n_obs

  # split IDs and attributes
  split_ids <- sort(which(vapply(nodes, .is_split, logical(1L))))
  expect_equal(split_ids, c(1L, 2L, 3L, 4L, 6L, 13L))
  expect_equal(nodes[[1L]]$attribute, "V2")
  expect_equal(nodes[[2L]]$attribute, "V6")
  expect_equal(nodes[[3L]]$attribute, "V2")
  expect_equal(nodes[[4L]]$attribute, "V3")
  expect_true(isTRUE(nodes[[5L]]$leaf))    # V5 pruned to leaf
  expect_equal(nodes[[6L]]$attribute, "V6")
  expect_equal(nodes[[13L]]$attribute, "V4")

  # child_ids: left=2*nid, right=2*nid+1; parent_id round-trips
  for (nid in split_ids) {
    expect_equal(sort(nodes[[nid]]$child_ids), sort(c(2L*nid, 2L*nid+1L)),
                 label = paste0("node ", nid, " child_ids"))
    for (cid in nodes[[nid]]$child_ids)
      expect_equal(nodes[[cid]]$parent_id, nid,
                   label = paste0("node ", cid, " parent_id"))
  }

  # EXE-canonical n_obs per branch (exact, no tolerance)
  expect_equal(.n(2L), 110L); expect_equal(.n(3L),  90L)
  expect_equal(.n(4L),  29L); expect_equal(.n(5L),  81L)
  expect_equal(.n(6L),  47L); expect_equal(.n(7L),  43L)
  expect_equal(.n(8L),  16L); expect_equal(.n(9L),  13L)
  expect_equal(.n(12L), 10L); expect_equal(.n(13L), 37L)
  expect_equal(.n(26L), 19L); expect_equal(.n(27L), 18L)

  # overall ESS and stored confusion
  expect_equal(round(tree$overall_ess, 2), 68.08, tolerance = 0.5)
  conf <- tree$training_confusion
  expect_equal(conf[1,1], 104L, tolerance = 2L)
  expect_equal(conf[1,2],  20L, tolerance = 2L)
  expect_equal(conf[2,1],  12L, tolerance = 2L)
  expect_equal(conf[2,2],  64L, tolerance = 2L)
})
