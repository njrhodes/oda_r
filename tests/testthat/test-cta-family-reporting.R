###############################################################################
# test-cta-family-reporting.R — cta_family_table(), summary.cta_family(),
#                                print.cta_family_summary(), print.cta_family()
###############################################################################

# ---- Constructed helpers (no fitting) ----------------------------------------

# Feasible stump: two leaf nodes, overall_ess = 60, no weights.
# min_terminal_denom = 10, next_mindenom = 11, strata = 2.
# D = 100 / (60 / 2) - 2 = 100/30 - 2 ≈ 1.333
.report_feasible_tree <- function() {
  node1 <- list(
    node_id = 1L, parent_id = 0L, depth = 1L, leaf = FALSE,
    attribute = "x", attr_col = 1L, attr_type = "ordered",
    n_obs = 20L, n_weighted = 20,
    rule = list(type = "ordered_cut", cut_value = 3.5, direction = "0->1"),
    ess = 60.0, ess_weighted = NA_real_, p_mc = 0.02,
    loo_status = "STABLE", loo_ess = 58.0, loo_p = NA_real_,
    confusion = NULL, split_labels = c(0L, 1L), child_ids = c(2L, 3L)
  )
  node2 <- list(
    node_id = 2L, parent_id = 1L, depth = 2L, leaf = TRUE,
    n_obs = 10L, n_weighted = 10, majority_class = 0L,
    attribute = NA_character_, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  node3 <- list(
    node_id = 3L, parent_id = 1L, depth = 2L, leaf = TRUE,
    n_obs = 10L, n_weighted = 10, majority_class = 1L,
    attribute = NA_character_, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  structure(
    list(nodes = list(node1, node2, node3), no_tree = FALSE,
         overall_ess = 60.0, n_nodes = 3L, root_id = 1L,
         has_weights = FALSE, mindenom = 2L, alpha_split = 0.05,
         prune_alpha = 1.0, loo = "off", training_confusion = NULL),
    class = "cta_tree"
  )
}

# No-tree leaf: single leaf, no_tree = TRUE.
.report_no_tree <- function() {
  node1 <- list(
    node_id = 1L, parent_id = 0L, depth = 1L, leaf = TRUE,
    n_obs = 20L, n_weighted = 20, majority_class = 0L,
    attribute = NA_character_, rule = NULL,
    ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
    loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
    confusion = NULL, split_labels = integer(0), child_ids = integer(0)
  )
  structure(
    list(nodes = list(node1), no_tree = TRUE, overall_ess = NA_real_,
         n_nodes = 1L, root_id = 1L, has_weights = FALSE,
         mindenom = 11L, alpha_split = 0.05, prune_alpha = 1.0,
         loo = "off", training_confusion = NULL),
    class = "cta_tree"
  )
}

# Two-member family: member 1 feasible (MINDENOM=2), member 2 no-tree (MINDENOM=11).
# min_d_idx = 1L.
.report_two_member_family <- function() {
  tree1 <- .report_feasible_tree()
  tree2 <- .report_no_tree()
  m1 <- odacore:::new_cta_family_member(2L,  tree1)
  m2 <- odacore:::new_cta_family_member(11L, tree2)
  members   <- list(m1, m2)
  mindenoms <- c(2L, 11L)
  d_vals    <- c(m1$d, m2$d)
  no_trees  <- c(FALSE, TRUE)
  strata_v  <- c(m1$strata,             NA_integer_)
  mintd_v   <- c(m1$min_terminal_denom, NA_integer_)
  ess_v     <- c(m1$overall_ess,        NA_real_)
  status_v  <- c("stump", "no_tree")
  summary_df <- data.frame(
    mindenom           = mindenoms,
    status             = status_v,
    strata             = strata_v,
    min_terminal_denom = mintd_v,
    overall_ess        = ess_v,
    d                  = d_vals,
    no_tree            = no_trees,
    stringsAsFactors   = FALSE
  )
  odacore:::new_cta_family(
    members            = members,
    mindenoms          = mindenoms,
    summary            = summary_df,
    min_d_idx          = 1L,
    terminated         = TRUE,
    termination_reason = "no_tree"
  )
}

# No-tree-only family: single no-tree member, min_d_idx = NA.
.report_no_tree_only_family <- function() {
  tree <- .report_no_tree()
  m    <- odacore:::new_cta_family_member(2L, tree)
  members   <- list(m)
  mindenoms <- 2L
  summary_df <- data.frame(
    mindenom           = 2L,
    status             = "no_tree",
    strata             = NA_integer_,
    min_terminal_denom = NA_integer_,
    overall_ess        = NA_real_,
    d                  = NA_real_,
    no_tree            = TRUE,
    stringsAsFactors   = FALSE
  )
  odacore:::new_cta_family(
    members            = members,
    mindenoms          = mindenoms,
    summary            = summary_df,
    min_d_idx          = NA_integer_,
    terminated         = TRUE,
    termination_reason = "no_tree"
  )
}

# =============================================================================
# cta_family_table: return type and columns
# =============================================================================

test_that("cta_family_table: returns a data.frame", {
  expect_s3_class(cta_family_table(.report_two_member_family()), "data.frame")
})

test_that("cta_family_table: column names are exactly as specified", {
  expected_cols <- c("index", "mindenom", "status", "no_tree", "strata",
                     "min_terminal_denom", "next_mindenom", "overall_ess",
                     "has_weights", "d", "selected_min_d")
  df <- cta_family_table(.report_two_member_family())
  expect_equal(names(df), expected_cols)
})

test_that("cta_family_table: rows equal n_members", {
  df <- cta_family_table(.report_two_member_family())
  expect_equal(nrow(df), 2L)
})

test_that("cta_family_table: index is 1..n_members", {
  df <- cta_family_table(.report_two_member_family())
  expect_equal(df$index, c(1L, 2L))
})

test_that("cta_family_table: mindenom column matches family mindenoms", {
  df <- cta_family_table(.report_two_member_family())
  expect_equal(df$mindenom, c(2L, 11L))
})

test_that("cta_family_table: selected_min_d is logical", {
  df <- cta_family_table(.report_two_member_family())
  expect_type(df$selected_min_d, "logical")
})

test_that("cta_family_table: exactly one TRUE in selected_min_d when feasible exists", {
  df <- cta_family_table(.report_two_member_family())
  expect_equal(sum(df$selected_min_d), 1L)
  expect_true(df$selected_min_d[1L])
  expect_false(df$selected_min_d[2L])
})

test_that("cta_family_table: no-tree-only family has no selected_min_d TRUE", {
  df <- cta_family_table(.report_no_tree_only_family())
  expect_false(any(df$selected_min_d))
})

test_that("cta_family_table: no-tree row has status 'no_tree' and NA d/strata/min_terminal_denom", {
  df <- cta_family_table(.report_two_member_family())
  nt_row <- df[df$no_tree, , drop = FALSE]
  expect_equal(nt_row$status, "no_tree")
  expect_true(is.na(nt_row$d))
  expect_true(is.na(nt_row$strata))
  expect_true(is.na(nt_row$min_terminal_denom))
})

test_that("cta_family_table: next_mindenom for no-tree row is NA", {
  df <- cta_family_table(.report_two_member_family())
  expect_true(is.na(df$next_mindenom[df$no_tree]))
})

test_that("cta_family_table: feasible row has correct next_mindenom", {
  df <- cta_family_table(.report_two_member_family())
  # min_terminal_denom of .report_feasible_tree() = min(10, 10) = 10; next = 11
  expect_equal(df$next_mindenom[!df$no_tree], 11L)
})

test_that("cta_family_table: has_weights FALSE for constructed unweighted family", {
  df <- cta_family_table(.report_two_member_family())
  expect_true(all(!df$has_weights))
})

test_that("cta_family_table: empty when zero members", {
  fam <- odacore:::new_cta_family()
  df  <- cta_family_table(fam)
  expect_equal(nrow(df), 0L)
  expect_true("selected_min_d" %in% names(df))
})

# =============================================================================
# summary.cta_family: class and required fields
# =============================================================================

test_that("summary.cta_family: returns cta_family_summary class", {
  s <- summary(.report_two_member_family())
  expect_s3_class(s, "cta_family_summary")
})

test_that("summary.cta_family: has all required fields", {
  required <- c("n_members", "min_d_idx", "terminated",
                "termination_reason", "has_weights", "table")
  s <- summary(.report_two_member_family())
  expect_true(all(required %in% names(s)),
              info = paste("missing:", paste(setdiff(required, names(s)),
                                             collapse = ", ")))
})

test_that("summary.cta_family: n_members matches family length", {
  s <- summary(.report_two_member_family())
  expect_equal(s$n_members, 2L)
})

test_that("summary.cta_family: min_d_idx = 1L for two-member family", {
  s <- summary(.report_two_member_family())
  expect_equal(s$min_d_idx, 1L)
})

test_that("summary.cta_family: min_d_idx = NA for no-tree-only family", {
  s <- summary(.report_no_tree_only_family())
  expect_identical(s$min_d_idx, NA_integer_)
})

test_that("summary.cta_family: termination_reason = 'no_tree'", {
  s <- summary(.report_two_member_family())
  expect_equal(s$termination_reason, "no_tree")
})

test_that("summary.cta_family: terminated = TRUE", {
  s <- summary(.report_two_member_family())
  expect_true(s$terminated)
})

test_that("summary.cta_family: table is a data.frame with required columns", {
  s <- summary(.report_two_member_family())
  expect_s3_class(s$table, "data.frame")
  expected_cols <- c("index", "mindenom", "status", "no_tree", "strata",
                     "min_terminal_denom", "next_mindenom", "overall_ess",
                     "has_weights", "d", "selected_min_d")
  expect_equal(names(s$table), expected_cols)
})

test_that("summary.cta_family: has_weights = FALSE for unweighted family", {
  s <- summary(.report_two_member_family())
  expect_false(s$has_weights)
})

# =============================================================================
# print methods: smoke tests
# =============================================================================

test_that("print.cta_family_summary: output contains 'MINDENOM'", {
  s   <- summary(.report_two_member_family())
  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "MINDENOM", fixed = TRUE)
})

test_that("print.cta_family_summary: output contains 'D'", {
  s   <- summary(.report_two_member_family())
  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "D", fixed = TRUE)
})

test_that("print.cta_family_summary: output contains 'no_tree'", {
  s   <- summary(.report_two_member_family())
  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "no_tree", fixed = TRUE)
})

test_that("print.cta_family_summary: returns invisibly", {
  s <- summary(.report_two_member_family())
  expect_invisible(print(s))
})

test_that("print.cta_family: output contains 'MINDENOM'", {
  out <- paste(capture.output(print(.report_two_member_family())), collapse = "\n")
  expect_match(out, "MINDENOM", fixed = TRUE)
})

test_that("print.cta_family: output contains 'no_tree'", {
  out <- paste(capture.output(print(.report_two_member_family())), collapse = "\n")
  expect_match(out, "no_tree", fixed = TRUE)
})

test_that("print.cta_family: returns invisibly", {
  expect_invisible(print(.report_two_member_family()))
})

test_that("print.cta_family_summary: ESS label shown for unweighted family", {
  s   <- summary(.report_two_member_family())
  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "ESS", fixed = TRUE)
})

# =============================================================================
# Slow fixture tests — myeloma MDSA chain
# =============================================================================

.report_myeloma_family <- local({
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

test_that("cta_family_table: myeloma chain has 3 rows", {
  skip_if_slow_tests_disabled("cta-family-reporting")
  df <- cta_family_table(.report_myeloma_family())
  expect_equal(nrow(df), 3L)
})

test_that("cta_family_table: myeloma mindenoms are c(1, 30, 56)", {
  skip_if_slow_tests_disabled("cta-family-reporting")
  df <- cta_family_table(.report_myeloma_family())
  expect_equal(df$mindenom, c(1L, 30L, 56L))
})

test_that("cta_family_table: myeloma next_mindenom chain is c(30, 56, NA)", {
  skip_if_slow_tests_disabled("cta-family-reporting")
  df <- cta_family_table(.report_myeloma_family())
  expect_equal(df$next_mindenom[1L], 30L)
  expect_equal(df$next_mindenom[2L], 56L)
  expect_true(is.na(df$next_mindenom[3L]))
})

test_that("cta_family_table: myeloma terminal row is no_tree", {
  skip_if_slow_tests_disabled("cta-family-reporting")
  df <- cta_family_table(.report_myeloma_family())
  expect_true(df$no_tree[3L])
  expect_equal(df$status[3L], "no_tree")
})

test_that("cta_family_table: myeloma selected_min_d TRUE at index 1", {
  skip_if_slow_tests_disabled("cta-family-reporting")
  df <- cta_family_table(.report_myeloma_family())
  expect_true(df$selected_min_d[1L])
  expect_false(df$selected_min_d[2L])
  expect_false(df$selected_min_d[3L])
})

test_that("cta_family_table: myeloma has_weights TRUE for all members", {
  skip_if_slow_tests_disabled("cta-family-reporting")
  df <- cta_family_table(.report_myeloma_family())
  expect_true(all(df$has_weights))
})

test_that("summary.cta_family: myeloma n_members=3, min_d_idx=1, reason='no_tree'", {
  skip_if_slow_tests_disabled("cta-family-reporting")
  s <- summary(.report_myeloma_family())
  expect_equal(s$n_members,          3L)
  expect_equal(s$min_d_idx,          1L)
  expect_equal(s$termination_reason, "no_tree")
  expect_true(s$has_weights)
})

test_that("print.cta_family: myeloma output contains 'MINDENOM' and 'no_tree'", {
  skip_if_slow_tests_disabled("cta-family-reporting")
  out <- paste(capture.output(print(.report_myeloma_family())), collapse = "\n")
  expect_match(out, "MINDENOM", fixed = TRUE)
  expect_match(out, "no_tree",  fixed = TRUE)
})

test_that("print.cta_family: myeloma WESS label shown (weighted family)", {
  skip_if_slow_tests_disabled("cta-family-reporting")
  out <- paste(capture.output(print(.report_myeloma_family())), collapse = "\n")
  expect_match(out, "WESS", fixed = TRUE)
})
