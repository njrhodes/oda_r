###############################################################################
# test-cta-ort.R — Tests for cta_fit(recursive = TRUE) / LORT
#
# Tier: CRAN-safe. All fits use mc_iter = 100L, mc_seed = 42L, loo = "off".
###############################################################################

# ---------------------------------------------------------------------------
# Synthetic two-level dataset (n=60, 2 attributes A and B)
#   B perfectly separates: B>0.5 → all y=1; B<=0.5 → all y=0
#   LORT produces depth=1 ORT: root at 0, terminal children at 1 → 3 strata
# ---------------------------------------------------------------------------
syn_X <- data.frame(
  A = c(rep(0, 20), rep(1, 20), rep(1, 20)),
  B = c(rep(0, 20), rep(0, 20), rep(1, 20))
)
syn_y <- c(rep(0L, 20), rep(0L, 20), rep(1L, 20))

syn_ort_args <- list(
  X = syn_X, y = syn_y,
  recursive = TRUE, mc_iter = 100L, mc_seed = 42L, loo = "off", min_n = 5L
)

# Module-level fixture: computed once and reused
syn_ort <- do.call(cta_fit, syn_ort_args)

# Non-recursive tree (needed for assign_endpoints tests T38-T43)
syn_tree_named <- cta_fit(syn_X, syn_y, mc_iter = 100L, mc_seed = 42L,
                           loo = "off", min_n = 5L)

syn_X_wide <- data.frame(A = syn_X$A, B = syn_X$B, C_extra = seq_len(60L))
syn_X_wide_reordered <- data.frame(C_extra = seq_len(60L), B = syn_X$B, A = syn_X$A)

# =============================================================================
# Core contract tests
# =============================================================================

test_that("cta_ort: errors — mindenom with recursive, non-recursive is cta_tree not cta_ort", {
  expect_error(cta_fit(syn_X, syn_y, recursive = TRUE, mindenom = 5L),
               regexp = "mindenom")
  tree <- cta_fit(syn_X, syn_y, recursive = FALSE, mindenom = 1L,
                  mc_iter = 100L, mc_seed = 42L, loo = "off")
  expect_false(inherits(tree, "cta_ort"))
  expect_true(inherits(tree, "cta_tree"))
})

test_that("cta_ort: termination guards — min_n, max_depth, max_nodes", {
  # min_n guard: 1 stratum with stop_reason "min_n"
  ort_mn <- cta_fit(syn_X, syn_y, recursive = TRUE,
                    min_n = nrow(syn_X) + 1L,
                    mc_iter = 100L, mc_seed = 42L, loo = "off")
  expect_equal(ort_mn$n_strata, 1L)
  expect_equal(ort_mn$strata$stop_reason[1L], "min_n")
  expect_equal(ort_mn$strata$n[1L], nrow(syn_X))
  # max_depth guard: all terminal nodes at depth <= 1
  ort_md <- cta_fit(syn_X, syn_y, recursive = TRUE, max_depth = 1L,
                    mc_iter = 100L, mc_seed = 42L, loo = "off", min_n = 5L)
  depths <- vapply(Filter(function(nd) isTRUE(nd$is_terminal), ort_md$ort_nodes),
                   `[[`, integer(1L), "depth")
  expect_true(all(depths <= 1L))
  # max_nodes guard: stop_reason "max_nodes"
  ort_mx <- cta_fit(syn_X, syn_y, recursive = TRUE, max_nodes = 1L,
                    mc_iter = 100L, mc_seed = 42L, loo = "off", min_n = 5L)
  stop_reasons <- vapply(Filter(function(nd) isTRUE(nd$is_terminal), ort_mx$ort_nodes),
                         `[[`, character(1L), "stop_reason")
  expect_true(all(stop_reasons == "max_nodes"))
})

test_that("cta_ort: object identity, strata structure, settings, reproducibility", {
  ort <- syn_ort
  expect_s3_class(ort, "cta_ort")
  expect_s3_class(ort, "cta_tree")
  expect_true(isTRUE(ort$recursive))
  # Strata structure
  expect_equal(ort$n_strata, 3L)
  expect_equal(nrow(ort$strata), 3L)
  expect_equal(sum(ort$strata$n), nrow(syn_X))
  expect_true(isTRUE(ort$strata_check_passed))
  expect_true(all(diff(ort$strata$prop_class1) >= 0))
  # Seed stored
  expect_equal(ort$ort_settings$mc_seed, 42L)
  # Reproducibility: same seed → identical strata and predictions
  ort2 <- do.call(cta_fit, syn_ort_args)
  expect_equal(ort$strata$n,           ort2$strata$n)
  expect_equal(ort$strata$path,        ort2$strata$path)
  expect_equal(ort$strata$prop_class1, ort2$strata$prop_class1)
  expect_identical(predict(ort, syn_X, type="all")$stratum_id,
                   predict(ort2, syn_X, type="all")$stratum_id)
  # Different seed stored differently
  ort_b <- do.call(cta_fit, modifyList(syn_ort_args, list(mc_seed = 99L)))
  expect_equal(ort_b$ort_settings$mc_seed, 99L)
})

test_that("cta_ort: predict types — class/stratum/path/all; stratum counts match strata$n", {
  ort <- syn_ort
  # type="class"
  pc <- predict(ort, syn_X, type = "class")
  expect_true(is.integer(pc)); expect_equal(length(pc), nrow(syn_X))
  # type="stratum"
  ps <- predict(ort, syn_X, type = "stratum")
  expect_true(is.integer(ps)); expect_equal(length(ps), nrow(syn_X))
  # type="path"
  pp <- predict(ort, syn_X, type = "path")
  expect_true(is.character(pp)); expect_equal(length(pp), nrow(syn_X))
  # type="all"
  df <- predict(ort, syn_X, type = "all")
  expect_s3_class(df, "data.frame"); expect_equal(nrow(df), nrow(syn_X))
  expect_true(all(c("predicted_class","stratum_id","path","prop_class1","stop_reason") %in% names(df)))
  # Stratum row counts match strata$n
  for (sid in ort$strata$stratum_id) {
    predicted_n <- sum(!is.na(df$stratum_id) & df$stratum_id == sid)
    expect_equal(predicted_n, ort$strata$n[ort$strata$stratum_id == sid])
  }
})

test_that("cta_ort: predict edge cases — NA root attribute, no-tree root majority class", {
  ort <- syn_ort
  # NA in root-split attribute → NA with missing_action='na'
  new_row <- syn_X[1L, ]; new_row$A <- NA_integer_
  pred_na <- predict(ort, new_row, type = "all", missing_action = "na")
  expect_true(is.na(pred_na$predicted_class[1L]))
  expect_true(is.na(pred_na$stratum_id[1L]))
  # No-tree root (min_n guard): predict assigns majority class, no NAs
  ort_nt <- cta_fit(syn_X, syn_y, recursive = TRUE,
                    min_n = nrow(syn_X) + 1L, mc_iter = 100L, mc_seed = 42L, loo = "off")
  pred_nt <- predict(ort_nt, syn_X, type = "all")
  expect_true(all(!is.na(pred_nt$stratum_id)))
  expect_equal(length(unique(pred_nt$stratum_id)), 1L)
  expect_true(all(pred_nt$predicted_class == 0L, na.rm = TRUE))
})

test_that("cta_ort: print, summary, plot, ort_plot_data work without error", {
  ort <- syn_ort
  expect_silent(capture.output(print(ort)))
  sm <- summary(ort)
  expect_s3_class(sm, "cta_ort_summary")
  expect_silent(capture.output(print(sm)))
  pd <- ort_plot_data(ort)
  expect_type(pd, "list")
  expect_true(all(c("nodes","edges","strata") %in% names(pd)))
  expect_equal(nrow(pd$nodes), length(ort$ort_nodes))
  pdf(nullfile()); on.exit(dev.off())
  expect_silent(plot(ort))
  expect_silent(plot(ort, target_class = 1L))
})

# =============================================================================
# MC threading and family_max_steps
# =============================================================================

test_that("cta_ort: MC threading — mc_stop/mc_stopup stored; defaults set", {
  # Explicit values stored
  ort_mc <- do.call(cta_fit, modifyList(syn_ort_args,
                                        list(mc_stop = 99.0, mc_stopup = 10)))
  expect_equal(ort_mc$ort_settings$mc_stop,   99.0)
  expect_equal(ort_mc$ort_settings$mc_stopup, 10)
  # Defaults when not supplied
  expect_equal(syn_ort$ort_settings$mc_stop,   99.9)
  expect_true(is.na(syn_ort$ort_settings$mc_stopup))
})

test_that("cta_ort: family_max_steps — stored, default=20, invalid errors, non-recursive errors", {
  # Explicit stored
  ort_fms <- do.call(cta_fit, modifyList(syn_ort_args, list(family_max_steps = 5L)))
  expect_equal(ort_fms$ort_settings$family_max_steps, 5L)
  # Default
  expect_equal(syn_ort$ort_settings$family_max_steps, 20L)
  # Invalid errors
  expect_error(do.call(cta_fit, modifyList(syn_ort_args, list(family_max_steps = 0L))),
               regexp = "family_max_steps")
  # Non-recursive with explicit family_max_steps errors
  expect_error(cta_fit(syn_X, syn_y, recursive = FALSE, family_max_steps = 5L,
                       mindenom = 1L, mc_iter = 100L, mc_seed = 42L, loo = "off"),
               regexp = "family_max_steps")
  # tiny value=1 runs without error
  ort_tiny <- cta_fit(syn_X, syn_y, recursive = TRUE,
                      family_max_steps = 1L, mc_iter = 100L, mc_seed = 42L,
                      loo = "off", min_n = 5L)
  expect_true(inherits(ort_tiny, "cta_ort"))
  expect_true(ort_tiny$n_strata >= 1L)
})

# =============================================================================
# Wide-newdata column routing (T38-T44)
# =============================================================================

test_that("cta_ort routing: positional (narrow) and name-based (wide/reordered) give same result", {
  skip_if_not(inherits(syn_tree_named, "cta_tree"))
  ep_ref      <- cta_assign_endpoints(syn_tree_named, syn_X)
  ep_wide     <- cta_assign_endpoints(syn_tree_named, syn_X_wide)
  ep_reordered <- cta_assign_endpoints(syn_tree_named, syn_X_wide_reordered)
  expect_identical(ep_ref$endpoint_id, ep_wide$endpoint_id)
  expect_identical(ep_ref$endpoint_id, ep_reordered$endpoint_id)
  # Narrow unnamed (positional routing)
  nd_unnamed <- as.data.frame(matrix(c(syn_X$A, syn_X$B), ncol = 2L))
  expect_identical(ep_ref$endpoint_id,
                   cta_assign_endpoints(syn_tree_named, nd_unnamed)$endpoint_id)
  # predict.cta_ort with wide newdata same as narrow
  expect_identical(predict(syn_ort, syn_X), predict(syn_ort, syn_X_wide))
})

test_that("cta_ort routing: missing required variable errors; auto-named wide errors", {
  skip_if_not(inherits(syn_tree_named, "cta_tree"))
  nd_missing_B <- data.frame(A = syn_X$A, C_extra = seq_len(60L), D_extra = seq_len(60L))
  expect_error(cta_assign_endpoints(syn_tree_named, nd_missing_B),
               regexp = "missing required split variable")
  nd_autonamed_wide <- as.data.frame(matrix(c(syn_X$A, syn_X$B, seq_len(60L)), ncol = 3L))
  expect_error(cta_assign_endpoints(syn_tree_named, nd_autonamed_wide),
               regexp = "missing required split variable")
})

# =============================================================================
# LORT method taxonomy metadata (T45-T55)
# =============================================================================

test_that("LORT taxonomy: ort_settings fields all correct", {
  ort <- syn_ort
  expect_equal(ort$ort_settings$method,              "lort")
  expect_equal(ort$ort_settings$method_label,        "Locally Optimal Recursive Tree")
  expect_false(isTRUE(ort$ort_settings$global_optimization))
  expect_false(isTRUE(ort$ort_settings$global_lookahead))
  expect_false(isTRUE(ort$ort_settings$sda_anchored))
  expect_equal(ort$ort_settings$recursive_selection, "greedy_local_min_d")
})

test_that("LORT taxonomy: print output and summary metadata fields", {
  ort <- syn_ort
  out <- capture.output(print(ort))
  expect_true(any(grepl("Locally Optimal Recursive Tree", out, fixed = TRUE)))
  expect_true(any(grepl("greedy local min-D", out, fixed = TRUE)))
  expect_true(any(grepl("global optimization: no", out, fixed = TRUE)))
  expect_true(any(grepl("SDA anchored: no", out, fixed = TRUE)))
  sm <- summary(ort)
  expect_equal(sm$method,              "lort")
  expect_equal(sm$method_label,        "Locally Optimal Recursive Tree")
  expect_equal(sm$recursive_selection, "greedy_local_min_d")
  expect_false(isTRUE(sm$global_optimization))
  expect_false(isTRUE(sm$global_lookahead))
  expect_false(isTRUE(sm$sda_anchored))
})

# =============================================================================
# cta_ort_node_table (T56-T63)
# =============================================================================

test_that("cta_ort_node_table: schema — rows, required columns, method/flags all rows", {
  ort <- syn_ort
  tbl <- cta_ort_node_table(ort)
  expect_s3_class(tbl, "data.frame")
  expect_equal(nrow(tbl), length(ort$ort_nodes))
  required <- c("ort_node_id", "parent_ort_node_id", "depth", "n",
                "class_counts", "terminal", "stop_reason",
                "selected_mindenom", "selected_ess", "selected_d",
                "selected_root_attribute", "selected_tree_nodes",
                "selected_tree_leaves", "selected_endpoint_count",
                "child_ids", "method", "selection_scope",
                "global_optimization", "sda_anchored")
  expect_true(all(required %in% names(tbl)),
              info = paste("missing:", paste(setdiff(required, names(tbl)), collapse = ", ")))
  expect_true(all(tbl$method == "lort"))
  expect_true(all(!tbl$global_optimization))
  expect_true(all(!tbl$sda_anchored))
})

test_that("cta_ort_node_table: structure — root NA parent, non-root int parent, split nodes, depth>=1", {
  ort <- syn_ort
  tbl <- cta_ort_node_table(ort)
  # Root has NA parent; others have integer parent
  root_row <- tbl[tbl$ort_node_id == 1L, ]
  expect_true(is.na(root_row$parent_ort_node_id))
  expect_true(all(!is.na(tbl$parent_ort_node_id[tbl$ort_node_id != 1L])))
  # Non-terminal nodes expose selected_ess and root_attr
  non_term <- tbl[!tbl$terminal, ]
  expect_true(nrow(non_term) >= 1L)
  expect_true(all(!is.na(non_term$selected_ess)))
  expect_true(all(!is.na(non_term$selected_root_attribute)))
  # Max depth >= 1 (multi-level structure)
  expect_true(max(tbl$depth) >= 1L)
})

# =============================================================================
# Canonical CTA node geometry in LORT local trees (T64, T65)
# =============================================================================

test_that("lort_local_tree: every split node has canonical child_ids and parent_id (T64)", {
  tr <- lort_local_tree(syn_ort, 1L)
  expect_s3_class(tr, "cta_tree")
  expect_false(isTRUE(tr$no_tree))

  nodes <- tr$nodes
  expect_equal(tr$root_id, 1L)

  .is_split <- function(nd) !is.null(nd) && !isTRUE(nd$leaf) && length(nd$child_ids) > 0L
  split_ids <- which(vapply(nodes, .is_split, logical(1L)))
  expect_gt(length(split_ids), 0L)

  for (nid in split_ids) {
    nd       <- nodes[[nid]]
    left_id  <- 2L * nid
    right_id <- 2L * nid + 1L
    expect_equal(sort(nd$child_ids), c(left_id, right_id))
    expect_false(is.null(nodes[[left_id]]))
    expect_false(is.null(nodes[[right_id]]))
    expect_equal(nodes[[left_id]]$parent_id,  nid)
    expect_equal(nodes[[right_id]]$parent_id, nid)
  }
})

test_that("lort_local_tree: ordered_cut left branch = x<=cut, verified from fixture data (T65)", {
  tr    <- lort_local_tree(syn_ort, 1L)
  nodes <- tr$nodes
  .is_split <- function(nd) !is.null(nd) && !isTRUE(nd$leaf) && length(nd$child_ids) > 0L
  obs_at <- vector("list", length(nodes) + 1L)
  obs_at[[tr$root_id]] <- seq_len(nrow(syn_X))
  split_ids <- sort(which(vapply(nodes, .is_split, logical(1L))))
  checked <- 0L

  for (nid in split_ids) {
    nd  <- nodes[[nid]]
    idx <- obs_at[[nid]]
    if (is.null(idx) || length(idx) == 0L) next
    rule <- nd$rule
    if (!identical(rule$type, "ordered_cut")) next

    attr_nm  <- nd$attribute
    left_id  <- 2L * nid;  right_id <- 2L * nid + 1L
    x_vals   <- syn_X[[attr_nm]][idx]
    sides    <- oda_rule_side(x_vals, rule)
    n_left   <- sum(sides == 0L);  n_right <- sum(sides == 1L)

    expect_equal(nodes[[left_id]]$n_obs,  n_left)
    expect_equal(nodes[[right_id]]$n_obs, n_right)

    obs_at[[left_id]]  <- idx[sides == 0L]
    obs_at[[right_id]] <- idx[sides == 1L]
    checked <- checked + 1L
  }
  expect_gt(checked, 0L)
})
