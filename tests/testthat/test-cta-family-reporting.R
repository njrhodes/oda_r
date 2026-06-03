###############################################################################
# test-cta-family-reporting.R — cta_family_table(), summary.cta_family(),
#                                print.cta_family_summary(), print.cta_family()
###############################################################################

# ---- Constructed helpers (no fitting) ----------------------------------------

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

.report_two_member_family <- function() {
  tree1 <- .report_feasible_tree()
  tree2 <- .report_no_tree()
  m1 <- odacore:::new_cta_family_member(2L,  tree1)
  m2 <- odacore:::new_cta_family_member(11L, tree2)
  members   <- list(m1, m2)
  mindenoms <- c(2L, 11L)
  summary_df <- data.frame(
    mindenom           = mindenoms,
    status             = c("stump", "no_tree"),
    strata             = c(m1$strata,             NA_integer_),
    min_terminal_denom = c(m1$min_terminal_denom, NA_integer_),
    overall_ess        = c(m1$overall_ess,        NA_real_),
    d                  = c(m1$d, m2$d),
    no_tree            = c(FALSE, TRUE),
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

.report_no_tree_only_family <- function() {
  tree <- .report_no_tree()
  m    <- odacore:::new_cta_family_member(2L, tree)
  summary_df <- data.frame(
    mindenom = 2L, status = "no_tree", strata = NA_integer_,
    min_terminal_denom = NA_integer_, overall_ess = NA_real_,
    d = NA_real_, no_tree = TRUE, stringsAsFactors = FALSE
  )
  odacore:::new_cta_family(
    members = list(m), mindenoms = 2L, summary = summary_df,
    min_d_idx = NA_integer_, terminated = TRUE, termination_reason = "no_tree"
  )
}

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
          X = d[, ac, drop = FALSE], y = as.integer(d[["V1"]]),
          w = as.numeric(d[["V2"]]), miss_codes = -9, loo = "stable",
          alpha_split = 0.05, prune_alpha = 0.05, mc_iter = 5000L,
          mc_seed = 12345L, verbose = FALSE, start_mindenom = 1L
        )
      )
    }
    fam
  }
})

# =============================================================================
# 1. cta_family_table — schema and column contract
# =============================================================================

test_that("cta_family_table schema: data.frame with exact columns", {
  df <- cta_family_table(.report_two_member_family())
  expect_s3_class(df, "data.frame")
  expected_cols <- c("index", "mindenom", "status", "no_tree", "strata",
                     "min_terminal_denom", "next_mindenom", "overall_ess",
                     "has_weights", "d", "selected_min_d")
  expect_equal(names(df), expected_cols)
})

# =============================================================================
# 2. cta_family_table — value contract
# =============================================================================

test_that("cta_family_table values: rows, index, mindenoms, selected_min_d, no-tree fields, next_mindenom, empty", {
  fam <- .report_two_member_family()
  df  <- cta_family_table(fam)

  expect_equal(nrow(df), 2L)
  expect_equal(df$index, c(1L, 2L))
  expect_equal(df$mindenom, c(2L, 11L))

  # selected_min_d: logical, exactly one TRUE at index 1
  expect_type(df$selected_min_d, "logical")
  expect_equal(sum(df$selected_min_d), 1L)
  expect_true(df$selected_min_d[1L])
  expect_false(df$selected_min_d[2L])

  # no-tree row fields
  nt <- df[df$no_tree, , drop = FALSE]
  expect_equal(nt$status, "no_tree")
  expect_true(is.na(nt$d))
  expect_true(is.na(nt$strata))
  expect_true(is.na(nt$min_terminal_denom))
  expect_true(is.na(nt$next_mindenom))

  # feasible row next_mindenom: min_terminal_denom + 1 = 10+1 = 11
  expect_equal(df$next_mindenom[!df$no_tree], 11L)

  # has_weights FALSE for unweighted family
  expect_true(all(!df$has_weights))

  # no-tree-only family: no selected_min_d TRUE
  df_nt <- cta_family_table(.report_no_tree_only_family())
  expect_false(any(df_nt$selected_min_d))

  # empty family: 0 rows with required columns
  df_empty <- cta_family_table(odacore:::new_cta_family())
  expect_equal(nrow(df_empty), 0L)
  expect_true("selected_min_d" %in% names(df_empty))
})

# =============================================================================
# 3. summary.cta_family — class, fields, values
# =============================================================================

test_that("summary.cta_family: class, required fields, and correct values", {
  s <- summary(.report_two_member_family())
  expect_s3_class(s, "cta_family_summary")

  required <- c("n_members", "min_d_idx", "terminated",
                "termination_reason", "has_weights", "table")
  expect_true(all(required %in% names(s)),
              info = paste("missing:", paste(setdiff(required, names(s)), collapse = ", ")))

  expect_equal(s$n_members,          2L)
  expect_equal(s$min_d_idx,          1L)
  expect_equal(s$termination_reason, "no_tree")
  expect_true(s$terminated)
  expect_false(s$has_weights)

  # table is cta_family_table output
  expect_s3_class(s$table, "data.frame")
  expect_equal(names(s$table), c("index", "mindenom", "status", "no_tree", "strata",
                                  "min_terminal_denom", "next_mindenom", "overall_ess",
                                  "has_weights", "d", "selected_min_d"))

  # no-tree-only: min_d_idx = NA
  s_nt <- summary(.report_no_tree_only_family())
  expect_identical(s_nt$min_d_idx, NA_integer_)
})

# =============================================================================
# 4. print methods — output contract (CRAN-safe)
# =============================================================================

test_that("print.cta_family and print.cta_family_summary: output contract and invisible return", {
  fam <- .report_two_member_family()
  s   <- summary(fam)

  # print.cta_family_summary: MINDENOM, D, no_tree, ESS label present; returns invisible
  out_s <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out_s, "MINDENOM", fixed = TRUE)
  expect_match(out_s, "D",        fixed = TRUE)
  expect_match(out_s, "no_tree",  fixed = TRUE)
  expect_match(out_s, "ESS",      fixed = TRUE)
  expect_invisible(print(s))

  # print.cta_family: MINDENOM, no_tree present; returns invisible
  out_f <- paste(capture.output(print(fam)), collapse = "\n")
  expect_match(out_f, "MINDENOM", fixed = TRUE)
  expect_match(out_f, "no_tree",  fixed = TRUE)
  expect_invisible(print(fam))
})

# =============================================================================
# 5. Myeloma slow smoke — cta_family_table values
# =============================================================================

test_that("myeloma family table: 3 rows, mindenoms {1,30,56}, chain, terminal, selected_min_d, has_weights; summary shape", {
  skip_if_slow_tests_disabled("cta-family-reporting")
  fam <- .report_myeloma_family()
  df  <- cta_family_table(fam)

  expect_equal(nrow(df), 3L)
  expect_equal(df$mindenom, c(1L, 30L, 56L))
  expect_equal(df$next_mindenom[1L], 30L)
  expect_equal(df$next_mindenom[2L], 56L)
  expect_true(is.na(df$next_mindenom[3L]))
  expect_true(df$no_tree[3L])
  expect_equal(df$status[3L], "no_tree")
  expect_true(df$selected_min_d[1L])
  expect_false(df$selected_min_d[2L])
  expect_false(df$selected_min_d[3L])
  expect_true(all(df$has_weights))

  # member D-values: md=1 and md=30 finite; md=56 no_tree -> NA
  expect_true(is.finite(fam$members[[1L]]$d))
  expect_true(is.finite(fam$members[[2L]]$d))
  expect_true(is.na(fam$members[[3L]]$d))

  # summary shape (absorbed from cta-family block in test-cta.R)
  required <- c("mindenom","status","strata","min_terminal_denom",
                "overall_ess","d","no_tree")
  expect_true(all(required %in% names(fam$summary)))
  expect_equal(nrow(fam$summary), 3L)
})

# =============================================================================
# 6. Myeloma slow smoke — summary and print
# =============================================================================

test_that("myeloma summary: n_members=3, min_d_idx=1, has_weights=TRUE; print shows WESS", {
  skip_if_slow_tests_disabled("cta-family-reporting")
  s <- summary(.report_myeloma_family())
  expect_equal(s$n_members,          3L)
  expect_equal(s$min_d_idx,          1L)
  expect_equal(s$termination_reason, "no_tree")
  expect_true(s$has_weights)

  out <- paste(capture.output(print(.report_myeloma_family())), collapse = "\n")
  expect_match(out, "MINDENOM", fixed = TRUE)
  expect_match(out, "no_tree",  fixed = TRUE)
  expect_match(out, "WESS",     fixed = TRUE)
})
