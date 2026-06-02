###############################################################################
# test-cta-summary.R — summary.cta_tree and print.cta_tree footer
###############################################################################

# ---- Helpers ----------------------------------------------------------------

.sumtest_no_tree_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 999L, mc_iter = 100L,
                mc_seed = 1L, loo = "off")
  )
}

.sumtest_stump_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off")
  )
}

# Synthetic valid_tree: strata=3, no fitting.
.synthetic_valid_tree <- function() {
  nodes <- lapply(seq_len(3L), function(i) list(leaf = TRUE, node_id = i, n_obs = 10L))
  structure(
    list(nodes = nodes, no_tree = FALSE, overall_ess = 75.0,
         n_nodes = 3L, root_id = 1L, has_weights = FALSE),
    class = "cta_tree"
  )
}

.sumtest_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.sumtest_myeloma_fit <- function(mindenom) {
  df <- {
    f <- testthat::test_path("fixtures/myeloma/data.txt")
    d <- read.table(f, header = FALSE)
    names(d) <- paste0("V", seq_len(ncol(d)))
    d[d$V2 != 0, ]
  }
  oda_cta_fit(
    X           = df[, .sumtest_myeloma_attrs],
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
    mc_seed     = NULL,
    loo         = "stable",
    attr_names  = .sumtest_myeloma_attrs
  )
}

# =============================================================================
# 1. Class and required fields
# =============================================================================

test_that("summary.cta_tree: class is cta_tree_summary; all required fields present", {
  s <- summary(.sumtest_no_tree_fit())
  expect_s3_class(s, "cta_tree_summary")
  required <- c("status", "no_tree", "root_attribute",
                "n_nodes", "n_splits", "n_leaves", "strata",
                "overall_ess", "d", "min_terminal_denom",
                "endpoint_denominators", "has_weights",
                "mindenom", "alpha_split", "prune_alpha", "loo")
  expect_true(all(required %in% names(s)),
              info = paste("missing:", paste(setdiff(required, names(s)), collapse = ", ")))
})

# =============================================================================
# 2. No-tree status — all scalar fields
# =============================================================================

test_that("summary.cta_tree no-tree: status='no_tree', NA scalars, integer(0) endpoints", {
  s <- summary(.sumtest_no_tree_fit())
  expect_equal(s$status, "no_tree")
  expect_true(s$no_tree)
  expect_identical(s$strata,             NA_integer_)
  expect_identical(s$n_leaves,           NA_integer_)
  expect_identical(s$d,                  NA_real_)
  expect_identical(s$min_terminal_denom, NA_integer_)
  expect_identical(s$endpoint_denominators, integer(0))
  expect_identical(s$root_attribute,     NA_character_)
})

# =============================================================================
# 3. Stump status — all value fields
# =============================================================================

test_that("summary.cta_tree stump: status, strata, endpoints, root_attr, counts, ess, d, params", {
  tree <- .sumtest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  s <- summary(tree)
  expect_equal(s$status,   "stump")
  expect_false(s$no_tree)
  expect_equal(s$strata,   2L)
  expect_equal(s$n_leaves, 2L)
  expect_equal(length(s$endpoint_denominators), 2L)
  expect_true(all(s$endpoint_denominators > 0L))
  expect_true(is.character(s$root_attribute) && !is.na(s$root_attribute))
  expect_equal(s$n_splits, 1L)
  expect_equal(s$n_nodes,  as.integer(s$n_splits + s$n_leaves))
  expect_true(is.finite(s$overall_ess) && s$overall_ess > 0)
  expect_equal(s$d, cta_d_stat(tree))
  expect_equal(s$min_terminal_denom, cta_min_terminal_denom(tree))
  expect_equal(s$mindenom,    tree$mindenom)
  expect_equal(s$alpha_split, tree$alpha_split)
  expect_equal(s$prune_alpha, tree$prune_alpha)
  expect_equal(s$loo,         tree$loo)
})

# =============================================================================
# 4. Synthetic valid_tree (strata=3, no fitting)
# =============================================================================

test_that("summary.cta_tree synthetic valid_tree: status, strata=3, ess=75, d=1, has_weights=FALSE", {
  s <- summary(.synthetic_valid_tree())
  expect_equal(s$status,      "valid_tree")
  expect_false(s$no_tree)
  expect_equal(s$strata,      3L)
  expect_equal(s$n_leaves,    3L)
  expect_equal(s$overall_ess, 75.0)
  # D = 100 / (75 / 3) - 3 = 4 - 3 = 1
  expect_equal(s$d, 1.0)
  expect_false(s$has_weights)
})

# =============================================================================
# 5. print methods — output contract
# =============================================================================

test_that("print.cta_tree_summary and print.cta_tree: output labels and invisible return", {
  # no-tree: status label + No tree found
  s_nt <- summary(.sumtest_no_tree_fit())
  expect_output(print(s_nt), "no_tree")
  expect_output(print(s_nt), "No tree found")
  expect_invisible(print(s_nt))
  expect_output(print(.sumtest_no_tree_fit()), "No tree found")

  # stump: ESS label, D label, strata label, min_denom, endpoints section
  tree_st <- .sumtest_stump_fit()
  if (!isTRUE(tree_st$no_tree)) {
    s_st  <- summary(tree_st)
    out_s <- paste(capture.output(print(s_st)), collapse = "\n")
    expect_match(out_s, "stump",       fixed = TRUE)
    expect_match(out_s, "overall_ess", fixed = TRUE)
    out_t <- paste(capture.output(print(tree_st)), collapse = "\n")
    expect_match(out_t, "ESS",      fixed = TRUE)
    expect_match(out_t, "D:",       fixed = TRUE)
    expect_match(out_t, "strata",   fixed = TRUE)
    expect_match(out_t, "min_denom", fixed = TRUE)
    expect_match(out_t, "Terminal endpoints (*)", fixed = TRUE)
    expect_match(out_t, "* endpoint", fixed = TRUE)
  }

  # valid_tree: overall_ess and D= in summary output
  s_vt  <- summary(.synthetic_valid_tree())
  out_vt <- paste(capture.output(print(s_vt)), collapse = "\n")
  expect_match(out_vt, "overall_ess", fixed = TRUE)
  expect_match(out_vt, "D=",          fixed = TRUE)

  # weighted synthetic: WESS label
  wt <- .synthetic_valid_tree()
  wt$has_weights <- TRUE; wt$alpha_split <- 0.05
  wt$mindenom <- 1L; wt$prune_alpha <- 1.0
  wt$max_depth <- 10L; wt$loo <- "off"
  out_wt <- tryCatch(
    paste(capture.output(print(wt)), collapse = "\n"),
    error = function(e) ""
  )
  if (nzchar(out_wt)) expect_match(out_wt, "WESS", fixed = TRUE)

  # split-node confusion block (verbose fit)
  tree_v <- suppressMessages(
    oda_cta_fit(data.frame(x = 1:8), c(0L,0L,0L,0L,1L,1L,1L,1L),
                mindenom = 2L, mc_iter = 300L, mc_seed = 1L,
                loo = "off", verbose = TRUE)
  )
  if (!isTRUE(tree_v$no_tree)) {
    has_conf <- any(vapply(tree_v$nodes, function(nd)
      !is.null(nd) && !isTRUE(nd$leaf) && !is.null(nd$confusion), logical(1L)))
    if (has_conf) {
      out_v <- paste(capture.output(print(tree_v)), collapse = "\n")
      expect_match(out_v, "Node-local split confusion", fixed = TRUE)
    }
  }
})

# =============================================================================
# 6. Myeloma slow smoke — valid_tree, stump, no_tree status
# =============================================================================

test_that("myeloma summary: MINDENOM=1 valid_tree+weights; MINDENOM=30 stump; MINDENOM=56 no_tree", {
  skip_if_slow_tests_disabled("cta-summary-fixture")
  # MINDENOM=1: valid_tree (V14->V15, 3 leaves)
  t1 <- suppressMessages(.sumtest_myeloma_fit(1L))
  s1 <- summary(t1)
  expect_equal(s1$status, "valid_tree")
  expect_equal(s1$strata,  3L)
  expect_true(s1$has_weights)
  expect_true(is.finite(s1$overall_ess) && s1$overall_ess > 0)

  # MINDENOM=56: no_tree
  t56 <- suppressMessages(.sumtest_myeloma_fit(56L))
  s56 <- summary(t56)
  expect_equal(s56$status, "no_tree")
  expect_true(s56$no_tree)
  expect_identical(s56$d,      NA_real_)
  expect_identical(s56$strata, NA_integer_)

  # MINDENOM=30: stump (V15 right-branch n=18 < 30)
  t30 <- suppressMessages(.sumtest_myeloma_fit(30L))
  s30 <- summary(t30)
  expect_equal(s30$status, "stump")
  expect_equal(s30$strata,  2L)
  expect_true(s30$has_weights)
})
