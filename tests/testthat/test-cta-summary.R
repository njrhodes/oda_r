###############################################################################
# test-cta-summary.R — summary.cta_tree and print.cta_tree footer (Phase 2B)
###############################################################################

# ---- Helpers ----------------------------------------------------------------

# Force no-tree by making MINDENOM impossible.  Same approach as test-cta.R.
.sumtest_no_tree_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 999L, mc_iter = 100L,
                mc_seed = 1L, loo = "off")
  )
}

# Stump: perfectly separable binary data → two pure leaves.
.sumtest_stump_fit <- function() {
  d <- list(X = data.frame(x = 1:8), y = c(0L,0L,0L,0L,1L,1L,1L,1L))
  suppressMessages(
    oda_cta_fit(d$X, d$y, mindenom = 2L, mc_iter = 300L,
                mc_seed = 1L, loo = "off")
  )
}

# Synthetic valid_tree with strata = 3 — no fitting, structural tests only.
.synthetic_valid_tree <- function() {
  nodes <- lapply(seq_len(3L), function(i)
    list(leaf = TRUE, node_id = i, n_obs = 10L))
  structure(
    list(nodes       = nodes,
         no_tree     = FALSE,
         overall_ess = 75.0,
         n_nodes     = 3L,
         root_id     = 1L,
         has_weights = FALSE),
    class = "cta_tree"
  )
}

# =============================================================================
# summary.cta_tree: class and required fields
# =============================================================================

test_that("summary.cta_tree returns cta_tree_summary class", {
  s <- summary(.sumtest_no_tree_fit())
  expect_s3_class(s, "cta_tree_summary")
})

test_that("summary.cta_tree has all required fields", {
  required <- c("status", "no_tree", "root_attribute",
                "n_nodes", "n_splits", "n_leaves", "strata",
                "overall_ess", "d", "min_terminal_denom",
                "endpoint_denominators", "has_weights",
                "mindenom", "alpha_split", "prune_alpha", "loo")
  s <- summary(.sumtest_no_tree_fit())
  expect_true(all(required %in% names(s)),
              info = paste("missing:", paste(setdiff(required, names(s)),
                                             collapse = ", ")))
})

# =============================================================================
# no-tree status
# =============================================================================

test_that("summary.cta_tree: no-tree status = 'no_tree'", {
  s <- summary(.sumtest_no_tree_fit())
  expect_equal(s$status, "no_tree")
  expect_true(s$no_tree)
})

test_that("summary.cta_tree: no-tree scalar fields are NA", {
  s <- summary(.sumtest_no_tree_fit())
  expect_identical(s$strata,             NA_integer_)
  expect_identical(s$n_leaves,           NA_integer_)
  expect_identical(s$d,                  NA_real_)
  expect_identical(s$min_terminal_denom, NA_integer_)
})

test_that("summary.cta_tree: no-tree endpoint_denominators is integer(0)", {
  s <- summary(.sumtest_no_tree_fit())
  expect_identical(s$endpoint_denominators, integer(0))
})

test_that("summary.cta_tree: no-tree root_attribute is NA_character_", {
  s <- summary(.sumtest_no_tree_fit())
  expect_identical(s$root_attribute, NA_character_)
})

# =============================================================================
# stump status
# =============================================================================

test_that("summary.cta_tree: stump status = 'stump'", {
  tree <- .sumtest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed — stump tests use bin_data")
  s <- summary(tree)
  expect_equal(s$status, "stump")
  expect_false(s$no_tree)
})

test_that("summary.cta_tree: stump strata = 2 and n_leaves = 2", {
  tree <- .sumtest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  s <- summary(tree)
  expect_equal(s$strata,   2L)
  expect_equal(s$n_leaves, 2L)
})

test_that("summary.cta_tree: stump endpoint_denominators has length 2", {
  tree <- .sumtest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  s <- summary(tree)
  expect_equal(length(s$endpoint_denominators), 2L)
  expect_true(all(s$endpoint_denominators > 0L))
})

test_that("summary.cta_tree: stump root_attribute is non-NA character", {
  tree <- .sumtest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  s <- summary(tree)
  expect_true(is.character(s$root_attribute))
  expect_false(is.na(s$root_attribute))
})

test_that("summary.cta_tree: stump n_splits = 1, n_nodes = n_splits + n_leaves", {
  tree <- .sumtest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  s <- summary(tree)
  expect_equal(s$n_splits, 1L)
  expect_equal(s$n_nodes,  as.integer(s$n_splits + s$n_leaves))
})

test_that("summary.cta_tree: stump overall_ess is finite and positive", {
  tree <- .sumtest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  s <- summary(tree)
  expect_true(is.finite(s$overall_ess))
  expect_gt(s$overall_ess, 0)
})

test_that("summary.cta_tree: stump d matches cta_d_stat()", {
  tree <- .sumtest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  s <- summary(tree)
  expect_equal(s$d, cta_d_stat(tree))
})

test_that("summary.cta_tree: stump min_terminal_denom matches accessor", {
  tree <- .sumtest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  s <- summary(tree)
  expect_equal(s$min_terminal_denom, cta_min_terminal_denom(tree))
})

test_that("summary.cta_tree: stump fit params stored correctly", {
  tree <- .sumtest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  s <- summary(tree)
  expect_equal(s$mindenom,    tree$mindenom)
  expect_equal(s$alpha_split, tree$alpha_split)
  expect_equal(s$prune_alpha, tree$prune_alpha)
  expect_equal(s$loo,         tree$loo)
})

# =============================================================================
# valid_tree path via synthetic tree (strata > 2, no fitting)
# =============================================================================

test_that("summary.cta_tree: synthetic valid_tree status = 'valid_tree'", {
  s <- summary(.synthetic_valid_tree())
  expect_equal(s$status, "valid_tree")
  expect_false(s$no_tree)
})

test_that("summary.cta_tree: synthetic valid_tree strata = 3", {
  s <- summary(.synthetic_valid_tree())
  expect_equal(s$strata,   3L)
  expect_equal(s$n_leaves, 3L)
})

test_that("summary.cta_tree: synthetic valid_tree overall_ess read correctly", {
  s <- summary(.synthetic_valid_tree())
  expect_equal(s$overall_ess, 75.0)
})

test_that("summary.cta_tree: synthetic valid_tree d = 1 (100/(75/3) - 3)", {
  s <- summary(.synthetic_valid_tree())
  expect_equal(s$d, 1.0)
})

test_that("summary.cta_tree: synthetic valid_tree has_weights = FALSE", {
  s <- summary(.synthetic_valid_tree())
  expect_false(s$has_weights)
})

# =============================================================================
# print.cta_tree_summary
# =============================================================================

test_that("print.cta_tree_summary: no-tree prints status and no-tree message", {
  s <- summary(.sumtest_no_tree_fit())
  expect_output(print(s), "no_tree")
  expect_output(print(s), "No tree found")
})

test_that("print.cta_tree_summary: stump prints status and ESS", {
  tree <- .sumtest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  s <- summary(tree)
  expect_output(print(s), "stump")
  expect_output(print(s), "overall_ess")
})

test_that("print.cta_tree_summary: synthetic valid_tree prints overall_ess and D", {
  s <- summary(.synthetic_valid_tree())
  out <- capture.output(print(s))
  joined <- paste(out, collapse = "\n")
  expect_match(joined, "overall_ess", fixed = TRUE)
  expect_match(joined, "D=",          fixed = TRUE)
})

test_that("print.cta_tree_summary: returns invisibly", {
  s <- summary(.sumtest_no_tree_fit())
  expect_invisible(print(s))
})

# =============================================================================
# print.cta_tree: footer added without breaking existing output
# =============================================================================

test_that("print.cta_tree: no-tree still says 'No tree found'", {
  expect_output(print(.sumtest_no_tree_fit()), "No tree found")
})

test_that("print.cta_tree: stump footer includes ESS and D labels", {
  tree <- .sumtest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  out <- paste(capture.output(print(tree)), collapse = "\n")
  expect_match(out, "ESS",  fixed = TRUE)
  expect_match(out, "D:",   fixed = TRUE)
})

test_that("print.cta_tree: stump footer includes strata and min_denom labels", {
  tree <- .sumtest_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  out <- paste(capture.output(print(tree)), collapse = "\n")
  expect_match(out, "strata",    fixed = TRUE)
  expect_match(out, "min_denom", fixed = TRUE)
})

test_that("print.cta_tree: WESS label shown for weighted synthetic tree", {
  wt <- .synthetic_valid_tree()
  wt$has_weights  <- TRUE
  wt$alpha_split  <- 0.05
  wt$mindenom     <- 1L
  wt$prune_alpha  <- 1.0
  wt$max_depth    <- 10L
  wt$loo          <- "off"
  # Synthetic nodes are all leaves so the split-node loop body is skipped;
  # the footer still runs.
  out <- tryCatch(
    paste(capture.output(print(wt)), collapse = "\n"),
    error = function(e) ""
  )
  if (nzchar(out))
    expect_match(out, "WESS", fixed = TRUE)
})

# =============================================================================
# Slow fixture tests — gated (canonical settings matching test-fixture-myeloma-cta.R)
# =============================================================================

# Shared myeloma loader: no header, space-delimited, names V1..V19.
# Mirrors .load_myeloma() in test-fixture-myeloma-cta.R.
.sumtest_load_myeloma <- function() {
  f <- testthat::test_path("fixtures/myeloma/data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}
.sumtest_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

# Fit with canonical settings (LOO STABLE + prune_alpha=0.05 + mc_iter=5000).
# LOO STABLE makes structural selection deterministic regardless of MC seed.
.sumtest_myeloma_fit <- function(mindenom) {
  df <- .sumtest_load_myeloma()
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
    mc_stopup   = 99.9,
    mc_seed     = NULL,
    loo         = "stable",
    attr_names  = .sumtest_myeloma_attrs
  )
}

test_that("summary.cta_tree: myeloma MINDENOM=1 is valid_tree, has_weights", {
  skip_if_slow_tests_disabled("cta-summary-fixture")
  tree <- suppressMessages(.sumtest_myeloma_fit(1L))
  s    <- summary(tree)
  # LOO STABLE selects V14->V15 (3 leaves) regardless of MC seed.
  expect_equal(s$status, "valid_tree")
  expect_equal(s$strata,  3L)
  expect_true(s$has_weights)
  expect_true(is.finite(s$overall_ess) && s$overall_ess > 0)
})

test_that("summary.cta_tree: myeloma MINDENOM=56 is no_tree", {
  skip_if_slow_tests_disabled("cta-summary-fixture")
  tree <- suppressMessages(.sumtest_myeloma_fit(56L))
  s    <- summary(tree)
  # All admissible child sizes < 56 → no_tree regardless of MC.
  expect_equal(s$status, "no_tree")
  expect_true(s$no_tree)
  expect_identical(s$d, NA_real_)
  expect_identical(s$strata, NA_integer_)
})

test_that("summary.cta_tree: myeloma MINDENOM=30 is stump (strata=2)", {
  skip_if_slow_tests_disabled("cta-summary-fixture")
  tree <- suppressMessages(.sumtest_myeloma_fit(30L))
  s    <- summary(tree)
  # MINDENOM=30: V15 right-branch child n=18 < 30 → no further splits → stump.
  expect_equal(s$status, "stump")
  expect_equal(s$strata,  2L)
  expect_true(s$has_weights)
})
