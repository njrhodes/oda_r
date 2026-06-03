###############################################################################
# test-cta.R — CTA core tests
###############################################################################

# ---- Data helpers -----------------------------------------------------------
bin_data <- function() {
  list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
}

# ---- No-tree helper ---------------------------------------------------------
.no_tree_fit <- function() {
  d <- bin_data()
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 999L, mc_iter = 100L,
                mc_seed = 1L, loo = "off")
  )
}

# =============================================================================
# 1. Basic fit contract
# =============================================================================

test_that("basic fit contract: class, root n_obs, mindenom, max_depth, case weights", {
  d <- bin_data()

  # cta_tree class and root n_obs
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L, mc_seed = 1L, loo = "off")
  expect_s3_class(tree, "cta_tree")
  expect_gte(tree$n_nodes, 1L)
  expect_equal(tree$nodes[[tree$root_id]]$n_obs, length(d$y))

  # mindenom too large → single leaf
  nt <- oda_cta_fit(d$X, d$y, mindenom = 999L, mc_iter = 100L, mc_seed = 1L, loo = "off")
  expect_equal(nt$n_nodes, 1L)
  expect_true(isTRUE(nt$nodes[[1]]$leaf))

  # case weights stored on root
  w    <- c(2,1,2,1,2,1,2,1)
  wtree <- oda_cta_fit(d$X, d$y, w = w, mindenom = 2L, mc_iter = 200L, mc_seed = 1L, loo = "off")
  expect_s3_class(wtree, "cta_tree")
  expect_equal(wtree$nodes[[1]]$n_weighted, sum(w))
})

# =============================================================================
# 2. max_depth enforcement + node_table columns
# =============================================================================

test_that("max_depth: depth-1 root splits → children are leaves; depth-2 limit enforced; node_table schema", {
  d <- bin_data()

  # max_depth=1: root is split, all children are leaves, no node deeper than 2
  t1   <- oda_cta_fit(d$X, d$y, mindenom = 2L, max_depth = 1L, mc_iter = 300L, mc_seed = 1L, loo = "off")
  root <- t1$nodes[[t1$root_id]]
  expect_false(isTRUE(root$leaf), label = "root is a split node")
  expect_true(all(vapply(t1$nodes[root$child_ids], function(nd) isTRUE(nd$leaf), logical(1))),
              label = "all children of root are leaves")
  expect_lte(max(vapply(t1$nodes, function(nd) nd$depth, integer(1))), 2L)

  # max_depth=2: no node deeper than 2
  t2     <- oda_cta_fit(d$X, d$y, mindenom = 2L, max_depth = 2L, mc_iter = 300L, mc_seed = 1L, loo = "off")
  depths <- vapply(t2$nodes, function(nd) nd$depth, integer(1))
  expect_lte(max(depths), 2L)

  # node_table: data.frame with required columns
  tbl <- cta_node_table(t2)
  expect_true(is.data.frame(tbl))
  expect_true(all(c("node_id","depth","n_obs","ess","p_mc","leaf","attribute") %in% names(tbl)))
})

# =============================================================================
# 3. predict contract
# =============================================================================

test_that("predict: length/type/classes in training set/accuracy >= majority", {
  d    <- bin_data()
  tree <- oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 500L, mc_seed = 7L, loo = "off")
  preds <- predict(tree, d$X)

  expect_equal(length(preds), nrow(d$X))
  expect_true(is.integer(preds))
  expect_true(all(is.na(preds) | preds %in% unique(d$y)))

  acc          <- mean(preds == d$y, na.rm = TRUE)
  majority_acc <- max(table(d$y)) / length(d$y)
  expect_gte(acc, majority_acc - 0.01)
})

# =============================================================================
# 4. print.cta_tree + cta_fit() public wrapper
# =============================================================================

test_that("print.cta_tree mentions 'CTA Tree'; cta_fit() matches oda_cta_fit() root", {
  d <- bin_data()
  expect_output(
    print(suppressMessages(
      oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L, mc_seed = 1L, loo = "off")
    )),
    "CTA Tree"
  )

  X <- data.frame(x1 = 1:8, x2 = c(0L,0L,1L,0L,1L,1L,0L,1L))
  y <- c(1L,1L,1L,1L,2L,2L,2L,2L)
  args <- list(X = X, y = y, priors_on = TRUE, mindenom = 1L,
               mc_iter = 300L, mc_seed = 7L, loo = "off",
               attr_names = c("x1","x2"))
  t1 <- do.call(cta_fit,     args)
  t2 <- do.call(oda_cta_fit, args)
  expect_s3_class(t1, "cta_tree")
  r1 <- t1$nodes[[t1$root_id]]; r2 <- t2$nodes[[t2$root_id]]
  expect_equal(r1$attribute,      r2$attribute)
  expect_equal(r1$rule$cut_value, r2$rule$cut_value)
  expect_equal(r1$n_obs,          r2$n_obs)
  expect_true(is.data.frame(cta_node_table(t1)) && nrow(cta_node_table(t1)) >= 1L)
})

# =============================================================================
# 5. Verbose reporting
# =============================================================================

.loo_fit_args <- list(
  priors_on = TRUE, alpha_split = 0.05, mindenom = 1L, prune_alpha = 0.05,
  max_depth = 3L, ess_min = 0, mc_iter = 5000L, mc_target = 0.05,
  mc_stop = 99.9, mc_seed = NULL, loo = "stable"
)

test_that("verbose=FALSE no [CTA] messages; verbose=TRUE emits them; no model change", {
  d <- .make_loo_gate_data()
  X <- d[, c("Trap", "Stable")]
  y <- d$Class

  collect_msgs <- function(...) {
    msgs <- character(0)
    withCallingHandlers(
      do.call(cta_fit, c(list(X = X, y = y, ...), .loo_fit_args)),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    msgs
  }

  msgs_quiet <- collect_msgs(verbose = FALSE)
  msgs_verb  <- collect_msgs(verbose = TRUE)
  expect_false(any(grepl("\\[CTA", msgs_quiet)))
  expect_true(any(grepl("\\[CTA", msgs_verb)))
  expect_true(any(grepl("selected|no valid root", msgs_verb)))

  # Model output unchanged
  t_q <- do.call(cta_fit, c(list(X = X, y = y, verbose = FALSE), .loo_fit_args))
  t_v <- suppressMessages(do.call(cta_fit, c(list(X = X, y = y, verbose = TRUE), .loo_fit_args)))
  expect_equal(t_q$nodes[[t_q$root_id]]$attribute,      t_v$nodes[[t_v$root_id]]$attribute)
  expect_equal(t_q$nodes[[t_q$root_id]]$rule$cut_value, t_v$nodes[[t_v$root_id]]$rule$cut_value)
})

# =============================================================================
# 6. Multi-attribute: most informative attribute selected at root
# =============================================================================

test_that("multi-attribute: x1 (perfect separator) selected over x2 (noise)", {
  set.seed(42)
  x1 <- 1:10; x2 <- sample(1:10)
  y  <- c(0L,0L,0L,0L,0L,1L,1L,1L,1L,1L)
  tree <- oda_cta_fit(data.frame(x1=x1, x2=x2), y, mindenom = 2L, mc_iter = 500L,
                      mc_seed = 42L, loo = "off")
  root <- tree$nodes[[tree$root_id]]
  if (!isTRUE(root$leaf))
    expect_equal(root$attribute, "x1")
  else
    skip("root is leaf — increase mc_iter for reliable test")
})

# =============================================================================
# 7. No-tree regression
# =============================================================================

test_that("no-tree: flag TRUE, predict all NA_integer_, print 'No tree found'", {
  d    <- bin_data()
  tree <- .no_tree_fit()
  expect_true(isTRUE(tree$no_tree))
  preds <- predict(tree, d$X)
  expect_true(all(is.na(preds)))
  expect_true(is.integer(preds))
  expect_equal(length(preds), nrow(d$X))
  expect_output(print(tree), "No tree found")
})

# =============================================================================
# 8. predict: 1->0 rule routes correctly (training_confusion consistency)
# =============================================================================

test_that("predict 1->0 rule: public predict consistent with training_confusion; ESS > 50%", {
  X <- data.frame(V1 = 1:10)
  y <- c(1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L)
  tree <- oda_cta_fit(X, y, mindenom = 2L, mc_iter = 500L, mc_seed = 1L, loo = "off")
  skip_if(isTRUE(tree$no_tree), "no tree found; increase mc_iter")
  skip_if(is.null(tree$nodes[[tree$root_id]]$rule), "root has no rule")

  preds <- predict(tree, X)
  conf_pred <- matrix(0L, 2L, 2L, dimnames = list(c("0","1"), c("0","1")))
  for (i in seq_along(y)) {
    a <- as.character(y[i]); p <- as.character(preds[i])
    if (!is.na(p)) conf_pred[a, p] <- conf_pred[a, p] + 1L
  }
  tc <- tree$training_confusion
  expect_equal(conf_pred["0","0"], tc[rownames(tc)=="0", colnames(tc)=="0"][[1L]])
  expect_equal(conf_pred["1","1"], tc[rownames(tc)=="1", colnames(tc)=="1"][[1L]])

  pac0 <- conf_pred["0","0"] / sum(y == 0L)
  pac1 <- conf_pred["1","1"] / sum(y == 1L)
  expect_gt((pac0 + pac1 - 1) * 100, 50)
})

# =============================================================================
# 9. mc_stopup bypass: same root and ESS as default
# =============================================================================

test_that("mc_stopup=NA bypass: same root attr, cut, ESS as mc_stopup=NULL", {
  d    <- bin_data()
  args <- list(X = d$X, y = d$y, mindenom = 2L, mc_iter = 500L,
               mc_stop = 99.9, mc_seed = 1L, loo = "off")
  t_def <- do.call(oda_cta_fit, c(args, list(mc_stopup = NULL)))
  t_byp <- do.call(oda_cta_fit, c(args, list(mc_stopup = NA)))
  r_d <- t_def$nodes[[t_def$root_id]]
  r_b <- t_byp$nodes[[t_byp$root_id]]
  expect_equal(isTRUE(r_d$leaf),   isTRUE(r_b$leaf))
  expect_equal(r_d$attribute,      r_b$attribute)
  expect_equal(r_d$rule$cut_value, r_b$rule$cut_value)
  expect_equal(t_def$overall_ess,  t_byp$overall_ess)
})

# =============================================================================
# 10. Class-domain guard and canonical recoding
# =============================================================================

test_that("class-domain guard: 3-class errors; non-canonical {1,2} recoding; y_levels stored", {
  # 3-class must error
  expect_error(
    oda_cta_fit(data.frame(x = 1:9), c(0L,0L,0L,1L,1L,1L,2L,2L,2L),
                mindenom = 2L, mc_iter = 100L, mc_seed = 1L, loo = "off"),
    "CTA currently supports exactly two class levels; got 3"
  )

  # Non-canonical {1,2}: y_levels stored, predict returns {1,2}
  X <- data.frame(V1 = 1:8)
  y <- c(1L,1L,1L,1L,2L,2L,2L,2L)
  tree <- suppressMessages(
    oda_cta_fit(X, y, mindenom = 2L, mc_iter = 500L, mc_seed = 1L, loo = "off"))
  expect_equal(tree$y_levels, c(1L, 2L))
  if (!isTRUE(tree$no_tree))
    expect_true(all(predict(tree, X) %in% c(1L, 2L)))

  # No-tree with non-canonical still stores y_levels
  nt <- suppressMessages(
    oda_cta_fit(X, y, mindenom = 999L, mc_iter = 100L, mc_seed = 1L, loo = "off"))
  expect_true(isTRUE(nt$no_tree))
  expect_equal(nt$y_levels, c(1L, 2L))

  # Canonical {0,1}: y_levels is c(0L,1L)
  d <- bin_data()
  t0 <- suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L, mc_seed = 1L, loo = "off"))
  expect_equal(t0$y_levels, c(0L, 1L))
  if (!isTRUE(t0$no_tree))
    expect_true(all(predict(t0, d$X) %in% c(0L, 1L)))
})

# =============================================================================
# 11. Read-only CTA accessors: strata, denominators, min_terminal_denom
# =============================================================================

test_that("cta_strata / cta_endpoint_denominators / cta_min_terminal_denom: valid + no-tree", {
  d    <- bin_data()
  tree <- suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L, mc_seed = 1L, loo = "off"))

  if (!isTRUE(tree$no_tree)) {
    s  <- cta_strata(tree)
    ep <- cta_endpoint_denominators(tree)
    expect_true(is.integer(s) && !is.na(s) && s >= 2L)
    expect_true(is.integer(ep) && length(ep) >= 2L && !is.null(names(ep)) && all(ep > 0L))
    expect_equal(length(ep), s)
    expect_equal(cta_min_terminal_denom(tree), min(ep))
  }

  expect_identical(cta_strata(.no_tree_fit()),              NA_integer_)
  expect_identical(cta_endpoint_denominators(.no_tree_fit()), integer(0))
  expect_identical(cta_min_terminal_denom(.no_tree_fit()),  NA_integer_)
})

# =============================================================================
# 12. new_cta_family_member + new_cta_family
# =============================================================================

test_that("new_cta_family_member: valid + no-tree fields; new_cta_family S3 class", {
  d    <- bin_data()
  tree <- suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L, mc_seed = 1L, loo = "off"))

  if (!isTRUE(tree$no_tree)) {
    m <- odacore:::new_cta_family_member(2L, tree)
    expect_identical(m$mindenom,          2L)
    expect_false(m$no_tree)
    expect_identical(m$strata,            cta_strata(tree))
    expect_identical(m$min_terminal_denom, cta_min_terminal_denom(tree))
    expect_identical(m$next_mindenom,     cta_min_terminal_denom(tree) + 1L)
    expect_identical(m$overall_ess,       tree$overall_ess)
    expect_identical(m$d,                 cta_d_stat(tree))
  }

  nt <- .no_tree_fit()
  mn <- odacore:::new_cta_family_member(999L, nt)
  expect_identical(mn$mindenom,          999L)
  expect_true(mn$no_tree)
  expect_identical(mn$strata,            NA_integer_)
  expect_identical(mn$min_terminal_denom, NA_integer_)
  expect_identical(mn$next_mindenom,     NA_integer_)
  expect_identical(mn$overall_ess,       NA_real_)
  expect_identical(mn$d,                 NA_real_)

  fam <- odacore:::new_cta_family(list(mn))
  expect_s3_class(fam, "cta_family")
  expect_true(is.list(fam$members))
  expect_equal(length(fam$members), 1L)
})

# =============================================================================
# 13. Myeloma endpoint chain (SLOW): MINDENOM=1/30/56 denominator anchors
# =============================================================================

.myeloma_tree <- function(mindenom) {
  fpath <- testthat::test_path("fixtures/myeloma/data.txt")
  skip_if_not(file.exists(fpath), "myeloma fixture not available")
  d <- read.table(fpath); colnames(d) <- paste0("V", seq_len(ncol(d)))
  d  <- d[d$V2 > 0, ]
  ac <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")
  suppressMessages(oda_cta_fit(
    X = d[, ac, drop = FALSE], y = as.integer(d$V1), w = as.numeric(d$V2),
    miss_codes = -9, mindenom = as.integer(mindenom), loo = "stable",
    alpha_split = 0.05, mc_iter = 5000L, mc_seed = 12345L, verbose = FALSE
  ))
}

test_that("myeloma endpoint chain: MINDENOM=1→30→56 denominators and strata (SLOW)", {
  skip_if_slow_tests_disabled("cta-myeloma-chain")

  t1 <- .myeloma_tree(1L)
  expect_false(t1$no_tree)
  expect_equal(cta_strata(t1), 3L)
  expect_equal(cta_min_terminal_denom(t1), 29L)
  m1 <- odacore:::new_cta_family_member(1L, t1)
  expect_equal(m1$next_mindenom, 30L)

  t30 <- .myeloma_tree(30L)
  expect_false(t30$no_tree)
  expect_equal(cta_strata(t30), 2L)
  expect_equal(cta_min_terminal_denom(t30), 55L)
  m30 <- odacore:::new_cta_family_member(30L, t30)
  expect_equal(m30$next_mindenom, 56L)

  t56 <- .myeloma_tree(56L)
  expect_true(t56$no_tree)
  expect_identical(cta_strata(t56),              NA_integer_)
  expect_identical(cta_endpoint_denominators(t56), integer(0))
  m56 <- odacore:::new_cta_family_member(56L, t56)
  expect_true(m56$no_tree)
  expect_identical(m56$next_mindenom, NA_integer_)
})

# =============================================================================
# 14. cta_d_stat: NA conditions and correct value
# =============================================================================

.synthetic_cta_tree <- function(overall_ess, n_leaves, no_tree = FALSE) {
  nodes <- lapply(seq_len(n_leaves), function(i) list(leaf = TRUE, node_id = i, n_obs = 10L))
  structure(list(nodes = nodes, no_tree = no_tree, overall_ess = overall_ess), class = "cta_tree")
}

test_that("cta_d_stat: NA for no_tree/missing/non-finite/non-positive ess/strata<2; correct D=2", {
  expect_identical(cta_d_stat(.synthetic_cta_tree(50, 2L, no_tree=TRUE)), NA_real_)
  t_null <- .synthetic_cta_tree(50, 2L); t_null$overall_ess <- NULL
  expect_identical(cta_d_stat(t_null), NA_real_)
  expect_identical(cta_d_stat(.synthetic_cta_tree(Inf,       2L)), NA_real_)
  expect_identical(cta_d_stat(.synthetic_cta_tree(NA_real_,  2L)), NA_real_)
  expect_identical(cta_d_stat(.synthetic_cta_tree(0,         2L)), NA_real_)
  expect_identical(cta_d_stat(.synthetic_cta_tree(-5,        2L)), NA_real_)
  expect_identical(cta_d_stat(.synthetic_cta_tree(50,        1L)), NA_real_)
  # D = 100 / (50 / 2) - 2 = 4 - 2 = 2
  expect_equal(cta_d_stat(.synthetic_cta_tree(50, 2L)), 2)
})

# =============================================================================
# 15. PRUNE: Sidak-Bonferroni + A×B×C ENUMERATE sentinel (SLOW)
# =============================================================================

.myeloma_prune_tree <- function(prune_alpha_val) {
  fpath <- testthat::test_path("fixtures/myeloma/data.txt")
  skip_if_not(file.exists(fpath), "myeloma fixture not available")
  d <- read.table(fpath); colnames(d) <- paste0("V", seq_len(ncol(d)))
  d  <- d[d[["V2"]] > 0, ]
  ac <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")
  suppressMessages(oda_cta_fit(
    X = d[, ac, drop = FALSE], y = as.integer(d[["V1"]]), w = as.numeric(d[["V2"]]),
    miss_codes = -9, mindenom = 1L, loo = "stable", alpha_split = 0.05,
    prune_alpha = prune_alpha_val, mc_iter = 5000L, mc_seed = 800L, verbose = FALSE
  ))
}

.n_split_nodes <- function(tree) {
  sum(vapply(tree[["nodes"]],
             function(nd) !is.null(nd) && !isTRUE(nd[["leaf"]]) &&
               length(nd[["child_ids"]]) > 0L, logical(1L)))
}

test_that("PRUNE + A×B×C: prune no-op and active both select V14->V15 WESS=27.69%; right=leaf, left=V15 (SLOW)", {
  skip_if_slow_tests_disabled("cta-myeloma-prune")

  # prune_alpha=1.0: pruning disabled; enumerated winner V14->V15 returned unchanged
  t_noop <- .myeloma_prune_tree(1.0)
  expect_equal(.n_split_nodes(t_noop), 2L,
               label = "prune no-op: 2 split nodes (V14->V15)")
  expect_equal(round(t_noop[["overall_ess"]], 2), 27.69,
               label = "prune no-op: WESS=27.69%")

  # prune_alpha=0.05: same enumerated winner (V11 pruned or not, ENUMERATE wins V14->V15)
  t_act <- .myeloma_prune_tree(0.05)
  expect_equal(.n_split_nodes(t_act), 2L,
               label = "prune active: 2 split nodes")
  expect_equal(round(t_act[["overall_ess"]], 2), 27.69,
               label = "prune active: WESS=27.69%")
  expect_equal(t_act[["nodes"]][[t_act[["root_id"]]]][["attribute"]], "V14")

  # A×B×C sentinel: root=V14; right child=leaf (C=leaf); left child=V15 split
  skip_if(isTRUE(t_act$no_tree) || isTRUE(t_act$nodes[[t_act$root_id]]$leaf),
          "ENUMERATE did not grow a tree")
  root       <- t_act$nodes[[t_act$root_id]]
  right_child <- t_act$nodes[[root$child_ids[2L]]]
  left_child  <- t_act$nodes[[root$child_ids[1L]]]
  expect_true(isTRUE(right_child$leaf),
              label = "right branch of V14 is a leaf (C=leaf sentinel)")
  expect_false(isTRUE(left_child$leaf),
               label = "left branch of V14 is a split node")
  expect_equal(left_child$attribute, "V15",
               label = "left branch attribute is V15")
})

# =============================================================================
# 16. .cta_predictions_degenerate unit tests
# =============================================================================

test_that(".cta_predictions_degenerate: all degenerate cases true; non-degenerate false; required_classes", {
  expect_true(odacore:::.cta_predictions_degenerate(c(0L, 0L, NA_integer_, 0L)))
  expect_true(odacore:::.cta_predictions_degenerate(c(NA_integer_, 1L, 1L, 1L)))
  expect_true(odacore:::.cta_predictions_degenerate(c(NA_integer_, NA_integer_)))
  expect_false(odacore:::.cta_predictions_degenerate(c(0L, 1L, 0L, NA_integer_, 1L)))
  expect_false(odacore:::.cta_predictions_degenerate(c(0L, NA_integer_, 1L)))
  expect_true(odacore:::.cta_predictions_degenerate(c(1L, 2L, NA_integer_),
                                                     required_classes = 3L))
  expect_false(odacore:::.cta_predictions_degenerate(c(1L, 2L, 3L, NA_integer_),
                                                      required_classes = 3L))
})

# =============================================================================
# 17. Degeneracy denom-admissibility (B1/B2/B4/B5) and arg contracts (B3)
# =============================================================================

test_that("denom-admissibility: B1 no_tree when MINDENOM>child; B2 valid non-degenerate; B4 structure; B5 output invariant", {
  X_nd <- data.frame(V1 = 1:10)
  y_nd <- c(rep(0L, 5L), rep(1L, 5L))

  # B1: MINDENOM=6 > n_child=5 → no_tree
  fit6 <- oda_cta_fit(X_nd, y_nd, mindenom = 6L, mc_iter = 500L, mc_seed = 1L, loo = "off")
  expect_true(isTRUE(fit6$no_tree))

  # B2: MINDENOM=4 ≤ n_child=5 → valid, predicts both classes
  fit4 <- oda_cta_fit(X_nd, y_nd, mindenom = 4L, mc_iter = 500L, mc_seed = 1L, loo = "off")
  expect_false(isTRUE(fit4$no_tree))
  expect_gte(length(unique(predict(fit4, X_nd)[!is.na(predict(fit4, X_nd))])), 2L)

  # B4: no_tree structure is well-formed; predict() doesn't crash
  expect_true(is.list(fit6$nodes))
  expect_equal(fit6$root_id, 1L)
  expect_equal(fit6$n, 10L)
  expect_equal(fit6$C, 2L)
  expect_equal(length(predict(fit6, X_nd)), 10L)

  # B5: sweep md=1,3,6 — valid trees predict both classes; md=6 no_tree confirmed
  for (md in c(1L, 3L)) {
    fit <- oda_cta_fit(X_nd, y_nd, mindenom = md, mc_iter = 500L, mc_seed = 1L, loo = "off")
    if (!isTRUE(fit$no_tree))
      expect_gte(length(unique(predict(fit, X_nd)[!is.na(predict(fit, X_nd))])), 2L,
                 label = sprintf("MINDENOM=%d valid tree predicts both classes", md))
  }
})

test_that("arg contracts: B3a mindenom rejected for recursive=TRUE; B3b min_n guard fires", {
  X_nd <- data.frame(V1 = 1:10)
  y_nd <- c(rep(0L, 5L), rep(1L, 5L))

  # B3a: recursive=TRUE rejects explicit mindenom
  expect_error(
    cta_fit(X_nd, y_nd, recursive = TRUE, mindenom = 4L,
            mc_iter = 100L, mc_seed = 1L, loo = "off"),
    regexp = "mindenom"
  )

  # B3b: min_n guard → root is terminal, stop_reason='min_n'
  ort  <- cta_fit(X_nd, y_nd, recursive = TRUE, min_n = 20L,
                  mc_iter = 100L, mc_seed = 1L, loo = "off")
  root <- ort$ort_nodes[["1"]]
  expect_true(isTRUE(root$is_terminal))
  expect_equal(root$stop_reason, "min_n")
})

# =============================================================================
# 18. Degeneracy integration: stochastic regression + family contracts
# =============================================================================

test_that("degeneracy: non-no_tree CTA predicts both classes; aggressive pruning safe; family contracts", {
  # Integration: any valid tree on balanced separable data predicts both classes
  d   <- bin_data()
  fit <- oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L, mc_seed = 1L, loo = "off")
  if (!isTRUE(fit$no_tree))
    expect_gte(length(unique(predict(fit, d$X)[!is.na(predict(fit, d$X))])), 2L)

  # Stochastic regression: severe imbalance + aggressive pruning
  set.seed(7)
  n0 <- 97; n1 <- 3
  y_imb <- c(rep(0L, n0), rep(1L, n1))
  X_imb <- data.frame(V1 = c(runif(n0,0,4), runif(n1,7,10)), V2 = runif(n0+n1))
  fit_imb <- oda_cta_fit(X_imb, y_imb, mindenom = 1L, mc_iter = 500L,
                          mc_seed = 17L, prune_alpha = 0.05, loo = "off")
  if (!isTRUE(fit_imb$no_tree)) {
    preds_imb <- predict(fit_imb, X_imb)
    expect_gte(length(unique(preds_imb[!is.na(preds_imb)])), 2L,
               label = "pruned imbalanced CTA must not be all-same-class")
    conf_df <- cta_confusion_table(fit_imb)
    expect_gte(length(unique(conf_df$predicted[conf_df$n > 0L])), 2L)
  }

  # Family: every non-no_tree member predicts both classes
  fam <- cta_descendant_family(d$X, d$y, mc_iter = 300L, mc_seed = 1L,
                                loo = "off", start_mindenom = 2L)
  for (i in seq_along(fam$members)) {
    m <- fam$members[[i]]
    if (!isTRUE(m$no_tree) && !is.null(m$tree))
      expect_gte(length(unique(predict(m$tree, d$X)[!is.na(predict(m$tree, d$X))])), 2L,
                 label = sprintf("family member %d predicts both classes", i))
  }

  # Family: min_d_idx never points to all-same-class tree
  mid <- fam$min_d_idx
  if (!is.na(mid) && !is.null(fam$members[[mid]]$tree) && !isTRUE(fam$members[[mid]]$tree$no_tree))
    expect_gte(length(unique(predict(fam$members[[mid]]$tree, d$X)[
      !is.na(predict(fam$members[[mid]]$tree, d$X))])), 2L)

  # Family: min_d_idx is NA when all members are no_tree
  fam_nt <- cta_descendant_family(d$X, d$y, mc_iter = 50L, mc_seed = 1L,
                                   loo = "off", start_mindenom = 9999L)
  expect_true(all(fam_nt$summary$no_tree))
  expect_true(is.na(fam_nt$min_d_idx))

  # LORT: node models non-degenerate
  ort <- cta_fit(d$X, d$y, recursive = TRUE, mc_iter = 300L, mc_seed = 1L, loo = "off")
  expect_false(is.null(ort$ort_nodes))
  for (nd_key in names(ort$ort_nodes)) {
    nd <- ort$ort_nodes[[nd_key]]
    if (is.null(nd) || isTRUE(nd$is_terminal) || is.null(nd$model) || isTRUE(nd$model$no_tree)) next
    preds_nd <- predict(nd$model, d$X)
    expect_gte(length(unique(preds_nd[!is.na(preds_nd)])), 2L,
               label = sprintf("ORT node '%s' model predicts both classes", nd_key))
  }
})

# (cta_descendant_family myeloma family block removed: unique assertions absorbed
# into test-cta-family-reporting.R B7. Independent family fit eliminated ~15s.)
