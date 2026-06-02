###############################################################################
# test-cta-ort.R — Tests for cta_fit(recursive = TRUE) / ORT
#
# Tier: CRAN-safe (all tests run at default tier).
# All fits use mc_iter = 100L, mc_seed = 42L, loo = "off" for speed.
###############################################################################

# ---------------------------------------------------------------------------
# Synthetic two-level dataset
#
# 60 observations, 2 attributes (A, B):
#   A = 0, B = 0 (n=20): y = 0
#   A = 1, B = 0 (n=20): y = 0
#   A = 1, B = 1 (n=20): y = 1
#
# B perfectly separates the full dataset (B>0.5 → all y=1; B<=0.5 → all y=0),
# so MDSA selects B as the root split with D=0. LORT produces a 2-level
# structure (ORT depth 0-indexed: root at 0, terminal children at 1):
#
#   Root LORT node 1 (depth=0): MDSA finds B stump; is_terminal = FALSE
#   Right child LORT node 2 (B > 0.5, depth=1): no_tree (pure y=1)
#   Left  child LORT node 3 (B <= 0.5, depth=1): no_tree (pure y=0)
#
# 2 terminal strata (nodes 2, 3). Max ORT depth = 1.
# ---------------------------------------------------------------------------
syn_X <- data.frame(
  A = c(rep(0, 20), rep(1, 20), rep(1, 20)),
  B = c(rep(0, 20), rep(0, 20), rep(1, 20))
)
syn_y <- c(rep(0L, 20), rep(0L, 20), rep(1L, 20))

syn_ort_args <- list(
  X           = syn_X,
  y           = syn_y,
  recursive   = TRUE,
  mc_iter     = 100L,
  mc_seed     = 42L,
  loo         = "off",
  min_n       = 5L
)

# Module-level fixture: computed once and reused across all tests that share
# these args.  Eliminates repeated identical fits (37 calls → 1) to keep
# the file CRAN-safe (< ~30 s instead of ~96 s).
# Tests that need variant args (max_depth, max_nodes, family_max_steps,
# different seed, non-recursive, etc.) continue to call cta_fit() directly.
syn_ort <- do.call(cta_fit, syn_ort_args)

# ---------------------------------------------------------------------------
# Test 1: mindenom error when recursive = TRUE
# ---------------------------------------------------------------------------
test_that("cta_fit: mindenom error when recursive = TRUE", {
  expect_error(
    cta_fit(syn_X, syn_y, recursive = TRUE, mindenom = 5L),
    regexp = "mindenom"
  )
})

# ---------------------------------------------------------------------------
# Test 2: min_n guard terminates at root
# ---------------------------------------------------------------------------
test_that("cta_fit recursive: min_n guard at root gives 1 stratum", {
  ort <- cta_fit(syn_X, syn_y, recursive = TRUE,
                 min_n   = nrow(syn_X) + 1L,
                 mc_iter = 100L, mc_seed = 42L, loo = "off")
  expect_true(inherits(ort, "cta_ort"))
  expect_true(inherits(ort, "cta_tree"))
  expect_equal(ort$n_strata, 1L)
  expect_equal(ort$strata$stop_reason[1L], "min_n")
  expect_equal(ort$strata$n[1L], nrow(syn_X))
})

# ---------------------------------------------------------------------------
# Test 3: class guard — inherits both cta_ort and cta_tree
# ---------------------------------------------------------------------------
test_that("cta_fit recursive: object inherits cta_ort and cta_tree", {
  ort <- syn_ort
  expect_s3_class(ort, "cta_ort")
  expect_s3_class(ort, "cta_tree")
  expect_true(isTRUE(ort$recursive))
})

# ---------------------------------------------------------------------------
# Test 4: two-level synthetic tree — 3 terminal strata
# ---------------------------------------------------------------------------
test_that("cta_fit recursive: two-level synthetic tree has 3 strata", {
  ort <- syn_ort
  expect_equal(ort$n_strata, 3L)
  expect_equal(nrow(ort$strata), 3L)
  # Strata sorted ascending by prop_class1
  expect_true(all(diff(ort$strata$prop_class1) >= 0))
})

# ---------------------------------------------------------------------------
# Test 5: strata n sums to total N
# ---------------------------------------------------------------------------
test_that("cta_fit recursive: strata n values sum to N", {
  ort <- syn_ort
  expect_equal(sum(ort$strata$n), nrow(syn_X))
})

# ---------------------------------------------------------------------------
# Test 6: strata_check_passed is TRUE
# ---------------------------------------------------------------------------
test_that("cta_fit recursive: strata_check_passed is TRUE", {
  ort <- syn_ort
  expect_true(isTRUE(ort$strata_check_passed))
})

# ---------------------------------------------------------------------------
# Test 7: predict type = "class" returns integer vector of length N
# ---------------------------------------------------------------------------
test_that("cta_fit recursive: predict type='class' has correct length and type", {
  ort  <- syn_ort
  pred <- predict(ort, syn_X, type = "class")
  expect_true(is.integer(pred))
  expect_equal(length(pred), nrow(syn_X))
})

# ---------------------------------------------------------------------------
# Test 8: predict type = "all" returns data.frame with correct columns
# ---------------------------------------------------------------------------
test_that("cta_fit recursive: predict type='all' returns correct data.frame", {
  ort <- syn_ort
  df  <- predict(ort, syn_X, type = "all")
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), nrow(syn_X))
  expect_true(all(c("predicted_class","stratum_id","path",
                    "prop_class1","stop_reason") %in% names(df)))
})

# ---------------------------------------------------------------------------
# Test 9: predict stratum counts match strata$n
# ---------------------------------------------------------------------------
test_that("cta_fit recursive: predict stratum row counts match strata$n", {
  ort <- syn_ort
  df  <- predict(ort, syn_X, type = "all")
  for (sid in ort$strata$stratum_id) {
    predicted_n <- sum(!is.na(df$stratum_id) & df$stratum_id == sid)
    stored_n    <- ort$strata$n[ort$strata$stratum_id == sid]
    expect_equal(predicted_n, stored_n,
                 info = sprintf("stratum_id %d", sid))
  }
})

# ---------------------------------------------------------------------------
# Test 10: max_depth guard
# ---------------------------------------------------------------------------
test_that("cta_fit recursive: max_depth = 1 produces depth-1 terminals only", {
  ort <- cta_fit(syn_X, syn_y, recursive = TRUE,
                 max_depth = 1L,
                 mc_iter   = 100L, mc_seed = 42L, loo = "off",
                 min_n     = 5L)
  # All terminal nodes must be at depth <= 1
  term_nodes <- Filter(function(nd) isTRUE(nd$is_terminal), ort$ort_nodes)
  depths <- vapply(term_nodes, `[[`, integer(1L), "depth")
  expect_true(all(depths <= 1L))
})

# ---------------------------------------------------------------------------
# Test 11: print and summary work without error
# ---------------------------------------------------------------------------
test_that("cta_fit recursive: print and summary run without error", {
  ort <- syn_ort
  expect_silent(capture.output(print(ort)))
  sm  <- summary(ort)
  expect_s3_class(sm, "cta_ort_summary")
  expect_silent(capture.output(print(sm)))
})

# ---------------------------------------------------------------------------
# Test 12: ort_plot_data returns valid structure
# ---------------------------------------------------------------------------
test_that("ort_plot_data: returns list with nodes, edges, strata", {
  ort <- syn_ort
  pd  <- ort_plot_data(ort)
  expect_type(pd, "list")
  expect_true(all(c("nodes","edges","strata") %in% names(pd)))
  expect_s3_class(pd$nodes, "data.frame")
  expect_s3_class(pd$edges, "data.frame")
  # All ort nodes should appear in nodes data frame
  expect_equal(nrow(pd$nodes), length(ort$ort_nodes))
  # Edge count: one per parent → child link
  n_split_nodes <- sum(!pd$nodes$is_terminal)
  # Each split node contributes one edge per child
  expect_true(nrow(pd$edges) >= n_split_nodes)
})

# ---------------------------------------------------------------------------
# Test 13: plot.cta_ort produces no error
# ---------------------------------------------------------------------------
test_that("plot.cta_ort: runs without error (no file output)", {
  ort <- syn_ort
  pdf(nullfile())
  on.exit(dev.off())
  expect_silent(plot(ort))
  expect_silent(plot(ort, target_class = 1L))
})

# ---------------------------------------------------------------------------
# Seed policy tests (Tests 14-16)
# ---------------------------------------------------------------------------

# Test 14: top-level mc_seed stored in ort_settings
test_that("cta_fit recursive: ort_settings$mc_seed equals supplied seed", {
  ort <- syn_ort   # syn_ort_args uses mc_seed = 42L
  expect_equal(ort$ort_settings$mc_seed, 42L)
})

# Test 15: reproducibility -- same seed gives identical strata and predictions
test_that("cta_fit recursive: same mc_seed gives identical strata and predictions", {
  ort1 <- syn_ort
  ort2 <- do.call(cta_fit, syn_ort_args)  # independent refit to verify reproducibility

  # Strata table identical (n, path, prop_class1, stop_reason)
  expect_equal(ort1$strata$n,           ort2$strata$n)
  expect_equal(ort1$strata$path,        ort2$strata$path)
  expect_equal(ort1$strata$prop_class1, ort2$strata$prop_class1)
  expect_equal(ort1$strata$stop_reason, ort2$strata$stop_reason)

  # Predictions identical
  pred1 <- predict(ort1, syn_X, type = "all")
  pred2 <- predict(ort2, syn_X, type = "all")
  expect_equal(pred1$stratum_id, pred2$stratum_id)
})

# Test 16: different mc_seed stores a different top-level seed
test_that("cta_fit recursive: different mc_seed stored in ort_settings", {
  ort_a <- syn_ort                                                           # seed 42 (module fixture)
  ort_b <- do.call(cta_fit, modifyList(syn_ort_args, list(mc_seed = 99L)))  # seed 99
  expect_equal(ort_a$ort_settings$mc_seed, 42L)
  expect_equal(ort_b$ort_settings$mc_seed, 99L)
})

# ---------------------------------------------------------------------------
# Hardening tests (Tests 17-21)
# ---------------------------------------------------------------------------

# Test 17: predict type = "stratum" and type = "path" return correct types
test_that("predict.cta_ort: type='stratum' is integer, type='path' is character", {
  ort <- syn_ort
  expect_true(is.integer(predict(ort, syn_X, type = "stratum")))
  expect_equal(length(predict(ort, syn_X, type = "stratum")), nrow(syn_X))
  expect_true(is.character(predict(ort, syn_X, type = "path")))
  expect_equal(length(predict(ort, syn_X, type = "path")), nrow(syn_X))
})

# Test 18: predict with NA in split attribute returns NA (missing_action = "na")
test_that("predict.cta_ort: NA in root-split attribute returns NA with missing_action='na'", {
  ort     <- syn_ort
  new_row <- syn_X[1L, ]
  new_row$A <- NA_integer_   # A is the root split attribute
  pred <- predict(ort, new_row, type = "all", missing_action = "na")
  expect_true(is.na(pred$predicted_class[1L]))
  expect_true(is.na(pred$stratum_id[1L]))
})

# Test 19: no-tree root (min_n guard) -- predict assigns majority class to all rows
test_that("predict.cta_ort: no-tree root assigns majority class, no NAs", {
  ort <- cta_fit(syn_X, syn_y, recursive = TRUE,
                 min_n   = nrow(syn_X) + 1L,
                 mc_iter = 100L, mc_seed = 42L, loo = "off")
  pred <- predict(ort, syn_X, type = "all")
  # All rows assigned; no NA stratum_ids
  expect_true(all(!is.na(pred$stratum_id)))
  expect_equal(length(unique(pred$stratum_id)), 1L)
  # Majority class of syn_y is 0 (40 zeros, 20 ones)
  expect_true(all(pred$predicted_class == 0L, na.rm = TRUE))
})

# Test 20: max_nodes guard produces stop_reason "max_nodes" and does not error
test_that("cta_fit recursive: max_nodes = 1 produces max_nodes stop_reason", {
  ort <- cta_fit(syn_X, syn_y, recursive = TRUE,
                 max_nodes = 1L,
                 mc_iter   = 100L, mc_seed = 42L, loo = "off",
                 min_n     = 5L)
  term_nodes   <- Filter(function(nd) isTRUE(nd$is_terminal), ort$ort_nodes)
  stop_reasons <- vapply(term_nodes, `[[`, character(1L), "stop_reason")
  expect_true(all(stop_reasons == "max_nodes"))
  # predict must not error
  pred <- predict(ort, syn_X, type = "all")
  expect_equal(nrow(pred), nrow(syn_X))
  expect_true(all(!is.na(pred$stratum_id)))
})

# Test 21: non-recursive cta_fit result is not cta_ort
test_that("cta_fit non-recursive: result is cta_tree but not cta_ort", {
  tree <- cta_fit(syn_X, syn_y,
                  recursive = FALSE,
                  mindenom  = 1L,
                  mc_iter   = 100L, mc_seed = 42L, loo = "off")
  expect_false(inherits(tree, "cta_ort"))
  expect_true(inherits(tree, "cta_tree"))
})

# ---------------------------------------------------------------------------
# MC threading contract tests (Tests 31-32)
# Confirm mc_stop / mc_stopup are accepted and stored.
# Do NOT assert wall-clock timing.
# ---------------------------------------------------------------------------

# Test 31: mc_stop / mc_stopup accepted and stored in ort_settings
test_that("cta_fit recursive: mc_stop and mc_stopup stored in ort_settings", {
  ort <- do.call(cta_fit, modifyList(syn_ort_args,
                                     list(mc_stop = 99.0, mc_stopup = 10)))
  expect_equal(ort$ort_settings$mc_stop,   99.0)
  expect_equal(ort$ort_settings$mc_stopup, 10)
})

# Test 32: mc_stop / mc_stopup default values stored when not supplied
test_that("cta_fit recursive: mc_stop/mc_stopup default to oda_cta_fit canonical values", {
  ort <- syn_ort   # no mc_stop / mc_stopup supplied
  expect_equal(ort$ort_settings$mc_stop,   99.9)
  expect_true(is.na(ort$ort_settings$mc_stopup))
})

# ---------------------------------------------------------------------------
# family_max_steps tests (Tests 33-36)
# ---------------------------------------------------------------------------

# Test 33: family_max_steps stored in ort_settings
test_that("cta_fit recursive: family_max_steps stored in ort_settings", {
  ort <- do.call(cta_fit, modifyList(syn_ort_args, list(family_max_steps = 5L)))
  expect_equal(ort$ort_settings$family_max_steps, 5L)
})

# Test 34: default family_max_steps = 20L stored when not supplied
test_that("cta_fit recursive: family_max_steps defaults to 20L in ort_settings", {
  ort <- syn_ort   # no family_max_steps supplied
  expect_equal(ort$ort_settings$family_max_steps, 20L)
})

# Test 35: invalid family_max_steps errors
test_that("cta_fit recursive: family_max_steps = 0 errors", {
  expect_error(
    do.call(cta_fit, modifyList(syn_ort_args, list(family_max_steps = 0L))),
    regexp = "family_max_steps"
  )
})

# Test 36: non-recursive explicit family_max_steps errors
test_that("cta_fit non-recursive: explicit family_max_steps errors", {
  expect_error(
    cta_fit(syn_X, syn_y, recursive = FALSE, family_max_steps = 5L,
            mindenom = 1L, mc_iter = 100L, mc_seed = 42L, loo = "off"),
    regexp = "family_max_steps"
  )
})

# Test 37: tiny family_max_steps = 1L does not error on synthetic ORT
test_that("cta_fit recursive: family_max_steps = 1L runs without error", {
  ort <- cta_fit(syn_X, syn_y, recursive = TRUE,
                 family_max_steps = 1L,
                 mc_iter = 100L, mc_seed = 42L, loo = "off",
                 min_n   = 5L)
  expect_true(inherits(ort, "cta_ort"))
  expect_true(ort$n_strata >= 1L)
  expect_true(isTRUE(ort$strata_check_passed))
})

# ---------------------------------------------------------------------------
# Wide-newdata column routing tests (T38-T44)
#
# Routing contract:
#   * ncol(newdata) == n_attrs  → positional routing (identity map, unchanged)
#   * ncol(newdata) != n_attrs  → name-based routing; no silent positional fallback
#
# Tests use a single non-recursive cta_tree fitted on syn_X (named, 2 cols)
# for direct cta_assign_endpoints() tests, plus one predict.cta_ort() test.
# ---------------------------------------------------------------------------

syn_tree_named <- cta_fit(syn_X, syn_y,
                           mc_iter = 100L, mc_seed = 42L, loo = "off",
                           min_n   = 5L)

# Wide newdata: extra unused column appended (3 cols, named)
syn_X_wide <- data.frame(
  A       = syn_X$A,
  B       = syn_X$B,
  C_extra = seq_len(60L)
)

# Wide newdata: split variables present but in different column order
syn_X_wide_reordered <- data.frame(
  C_extra = seq_len(60L),
  B       = syn_X$B,
  A       = syn_X$A
)

test_that("cta_assign_endpoints: named narrow newdata uses positional routing (T38)", {
  ep_ref  <- cta_assign_endpoints(syn_tree_named, syn_X)
  ep_same <- cta_assign_endpoints(syn_tree_named, syn_X)
  expect_identical(ep_ref$endpoint_id, ep_same$endpoint_id)
})

test_that("cta_assign_endpoints: wide named newdata with extra column gives same result (T39)", {
  skip_if_not(inherits(syn_tree_named, "cta_tree"))
  ep_narrow <- cta_assign_endpoints(syn_tree_named, syn_X)
  ep_wide   <- cta_assign_endpoints(syn_tree_named, syn_X_wide)
  expect_identical(ep_narrow$endpoint_id, ep_wide$endpoint_id)
})

test_that("cta_assign_endpoints: wide named newdata with reordered columns gives same result (T40)", {
  skip_if_not(inherits(syn_tree_named, "cta_tree"))
  ep_narrow   <- cta_assign_endpoints(syn_tree_named, syn_X)
  ep_reordered <- cta_assign_endpoints(syn_tree_named, syn_X_wide_reordered)
  expect_identical(ep_narrow$endpoint_id, ep_reordered$endpoint_id)
})

test_that("cta_assign_endpoints: wide newdata missing required split variable errors (T41)", {
  skip_if_not(inherits(syn_tree_named, "cta_tree"))
  # 3 named columns but B is absent — triggers wide routing and errors
  nd_missing_B <- data.frame(A = syn_X$A, C_extra = seq_len(60L), D_extra = seq_len(60L))
  expect_error(
    cta_assign_endpoints(syn_tree_named, nd_missing_B),
    regexp = "missing required split variable"
  )
})

test_that("cta_assign_endpoints: unnamed same-width newdata uses positional routing (T42)", {
  skip_if_not(inherits(syn_tree_named, "cta_tree"))
  ep_ref     <- cta_assign_endpoints(syn_tree_named, syn_X)
  nd_unnamed <- as.data.frame(matrix(c(syn_X$A, syn_X$B), ncol = 2L))
  # V1=A, V2=B — same column positions as training X
  expect_identical(ep_ref$endpoint_id,
                   cta_assign_endpoints(syn_tree_named, nd_unnamed)$endpoint_id)
})

test_that("cta_assign_endpoints: wide newdata with non-matching auto-names errors (T43)", {
  skip_if_not(inherits(syn_tree_named, "cta_tree"))
  # matrix() -> as.data.frame() auto-names columns V1/V2/V3 (not A/B) — wide routing
  # then errors because required split variable names A and B are absent
  nd_autonamed_wide <- as.data.frame(matrix(c(syn_X$A, syn_X$B, seq_len(60L)), ncol = 3L))
  expect_error(
    cta_assign_endpoints(syn_tree_named, nd_autonamed_wide),
    regexp = "missing required split variable"
  )
})

test_that("predict.cta_ort: wide named newdata gives same result as narrow newdata (T44)", {
  ort <- syn_ort  # same as syn_ort_args
  # default type = "class" returns integer vector
  pred_narrow <- predict(ort, syn_X)
  pred_wide   <- predict(ort, syn_X_wide)
  expect_identical(pred_narrow, pred_wide)
})

# ---------------------------------------------------------------------------
# LORT method taxonomy tests (T45-T55)
# Confirm recursive=TRUE fits identify as LORT, not generic ORT/SORT/GORT.
# ---------------------------------------------------------------------------

test_that("LORT metadata: method == 'lort' (T45)", {
  ort <- syn_ort
  expect_equal(ort$ort_settings$method, "lort")
})

test_that("LORT metadata: method_label correct (T46)", {
  ort <- syn_ort
  expect_equal(ort$ort_settings$method_label, "Locally Optimal Recursive Tree")
})

test_that("LORT metadata: global_optimization is FALSE (T47)", {
  ort <- syn_ort
  expect_false(isTRUE(ort$ort_settings$global_optimization))
})

test_that("LORT metadata: global_lookahead is FALSE (T48)", {
  ort <- syn_ort
  expect_false(isTRUE(ort$ort_settings$global_lookahead))
})

test_that("LORT metadata: sda_anchored is FALSE (T49)", {
  ort <- syn_ort
  expect_false(isTRUE(ort$ort_settings$sda_anchored))
})

test_that("LORT metadata: recursive_selection == 'greedy_local_min_d' (T50)", {
  ort <- syn_ort
  expect_equal(ort$ort_settings$recursive_selection, "greedy_local_min_d")
})

test_that("LORT print: header contains 'Locally Optimal Recursive Tree' (T51)", {
  ort <- syn_ort
  out <- capture.output(print(ort))
  expect_true(any(grepl("Locally Optimal Recursive Tree", out, fixed = TRUE)))
})

test_that("LORT print: contains 'greedy local min-D' (T52)", {
  ort <- syn_ort
  out <- capture.output(print(ort))
  expect_true(any(grepl("greedy local min-D", out, fixed = TRUE)))
})

test_that("LORT print: contains 'global optimization: no' (T53)", {
  ort <- syn_ort
  out <- capture.output(print(ort))
  expect_true(any(grepl("global optimization: no", out, fixed = TRUE)))
})

test_that("LORT print: contains 'SDA anchored: no' (T54)", {
  ort <- syn_ort
  out <- capture.output(print(ort))
  expect_true(any(grepl("SDA anchored: no", out, fixed = TRUE)))
})

test_that("LORT summary: method metadata fields present (T55)", {
  ort <- syn_ort
  sm  <- summary(ort)
  expect_equal(sm$method,              "lort")
  expect_equal(sm$method_label,        "Locally Optimal Recursive Tree")
  expect_equal(sm$recursive_selection, "greedy_local_min_d")
  expect_false(isTRUE(sm$global_optimization))
  expect_false(isTRUE(sm$global_lookahead))
  expect_false(isTRUE(sm$sda_anchored))
})

# ---------------------------------------------------------------------------
# cta_ort_node_table tests (T56-T63)
# ---------------------------------------------------------------------------

test_that("cta_ort_node_table: returns data.frame with correct row count (T56)", {
  ort <- syn_ort
  tbl <- cta_ort_node_table(ort)
  expect_s3_class(tbl, "data.frame")
  expect_equal(nrow(tbl), length(ort$ort_nodes))
})

test_that("cta_ort_node_table: required columns present (T57)", {
  ort <- syn_ort
  tbl <- cta_ort_node_table(ort)
  required <- c("ort_node_id", "parent_ort_node_id", "depth", "n",
                "class_counts", "terminal", "stop_reason",
                "selected_mindenom", "selected_ess", "selected_d",
                "selected_root_attribute", "selected_tree_nodes",
                "selected_tree_leaves", "selected_endpoint_count",
                "child_ids", "method", "selection_scope",
                "global_optimization", "sda_anchored")
  expect_true(all(required %in% names(tbl)),
              info = paste("missing:", paste(setdiff(required, names(tbl)), collapse = ", ")))
})

test_that("cta_ort_node_table: method column is always 'lort' (T58)", {
  ort <- syn_ort
  tbl <- cta_ort_node_table(ort)
  expect_true(all(tbl$method == "lort"))
})

test_that("cta_ort_node_table: global_optimization is FALSE for all rows (T59)", {
  ort <- syn_ort
  tbl <- cta_ort_node_table(ort)
  expect_true(all(!tbl$global_optimization))
})

test_that("cta_ort_node_table: sda_anchored is FALSE for all rows (T60)", {
  ort <- syn_ort
  tbl <- cta_ort_node_table(ort)
  expect_true(all(!tbl$sda_anchored))
})

test_that("cta_ort_node_table: root has NA parent, non-root has integer parent (T61)", {
  ort <- syn_ort
  tbl <- cta_ort_node_table(ort)
  root_row <- tbl[tbl$ort_node_id == 1L, ]
  expect_true(is.na(root_row$parent_ort_node_id))
  non_root <- tbl[tbl$ort_node_id != 1L, ]
  expect_true(all(!is.na(non_root$parent_ort_node_id)))
})

test_that("cta_ort_node_table: non-terminal nodes expose selected_ess and root_attr (T62)", {
  ort <- syn_ort
  tbl <- cta_ort_node_table(ort)
  non_term <- tbl[!tbl$terminal, ]
  expect_true(nrow(non_term) >= 1L)
  expect_true(all(!is.na(non_term$selected_ess)))
  expect_true(all(!is.na(non_term$selected_root_attribute)))
})

test_that("cta_ort_node_table: multi-level recursion visible — max depth >= 1 (T63)", {
  # syn_X produces a 2-level ORT: root MDSA selects the best attribute (depth=0),
  # leaving endpoint children at depth=1.  ORT depth is 0-indexed, so a
  # root-plus-children structure has max depth = 1.  Asserting >= 1L proves that
  # recursion actually ran (non-trivial structure exists beyond the root node).
  ort <- syn_ort
  tbl <- cta_ort_node_table(ort)
  expect_true(max(tbl$depth) >= 1L,
              info = sprintf("max depth was %d; expected >= 1 for two-level synthetic (root=0, children=1)", max(tbl$depth)))
  # At least one non-terminal node (the root ORT node) must exist
  non_term <- tbl[!tbl$terminal, ]
  expect_true(nrow(non_term) >= 1L)
})

# ---------------------------------------------------------------------------
# Canonical CTA node geometry in LORT local trees (T64, T65)
#
# Local trees from lort_local_tree() go through the same oda_cta_fit() code
# path as standalone CTA fits and must obey Appendix C geometry:
#   left child  = 2*nid  (x <= cut for ordered_cut rules)
#   right child = 2*nid+1 (x >  cut for ordered_cut rules)
#
# T64 — structural invariants: child_ids and parent_id links.
# T65 — branch-size invariants: for every ordered_cut split node, recompute
#        sides from the fixture data (no assumed root attribute) and verify
#        that n_obs of each child matches the computed side count.
# ---------------------------------------------------------------------------
test_that("lort_local_tree: every split node has canonical child_ids and parent_id (T64)", {
  tr <- lort_local_tree(syn_ort, 1L)
  expect_s3_class(tr, "cta_tree")
  expect_false(isTRUE(tr$no_tree), label = "local tree is not no_tree")

  nodes <- tr$nodes
  expect_equal(tr$root_id, 1L, label = "root_id = 1")

  .is_split <- function(nd) !is.null(nd) && !isTRUE(nd$leaf) && length(nd$child_ids) > 0L
  split_ids <- which(vapply(nodes, .is_split, logical(1L)))
  expect_gt(length(split_ids), 0L, label = "at least one split node exists")

  for (nid in split_ids) {
    nd       <- nodes[[nid]]
    left_id  <- 2L * nid
    right_id <- 2L * nid + 1L

    expect_equal(sort(nd$child_ids), c(left_id, right_id),
                 label = paste0("node ", nid, " child_ids = {", left_id, ",", right_id, "}"))

    left_nd  <- nodes[[left_id]]
    right_nd <- nodes[[right_id]]
    expect_false(is.null(left_nd),
                 label = paste0("node ", left_id, " exists (left child of ", nid, ")"))
    expect_false(is.null(right_nd),
                 label = paste0("node ", right_id, " exists (right child of ", nid, ")"))
    expect_equal(left_nd$parent_id,  nid,
                 label = paste0("node ", left_id,  " parent_id = ", nid))
    expect_equal(right_nd$parent_id, nid,
                 label = paste0("node ", right_id, " parent_id = ", nid))
  }
})

test_that("lort_local_tree: ordered_cut left branch = x<=cut, verified from fixture data (T65)", {
  # For every ordered_cut split node, determine which syn_X observations reach
  # that node by traversing from root.  Then recompute left/right using
  # oda_rule_side() and assert n_obs of each child matches the computed count.
  # No assumption about which attribute is root or what the cut value is.
  tr    <- lort_local_tree(syn_ort, 1L)
  nodes <- tr$nodes

  .is_split <- function(nd) !is.null(nd) && !isTRUE(nd$leaf) && length(nd$child_ids) > 0L

  # obs_at[[nid]] = integer vector of syn_X row indices at node nid.
  obs_at <- vector("list", length(nodes) + 1L)   # +1 for sparse safety
  obs_at[[tr$root_id]] <- seq_len(nrow(syn_X))

  split_ids <- sort(which(vapply(nodes, .is_split, logical(1L))))
  checked   <- 0L

  for (nid in split_ids) {
    nd  <- nodes[[nid]]
    idx <- obs_at[[nid]]
    if (is.null(idx) || length(idx) == 0L) next

    rule <- nd$rule
    if (!identical(rule$type, "ordered_cut")) next   # only test ordered_cut splits

    attr_nm  <- nd$attribute
    left_id  <- 2L * nid
    right_id <- 2L * nid + 1L

    x_vals <- syn_X[[attr_nm]][idx]
    sides  <- oda_rule_side(x_vals, rule)   # 0 = x<=cut, 1 = x>cut

    n_left  <- sum(sides == 0L)
    n_right <- sum(sides == 1L)

    left_nd  <- nodes[[left_id]]
    right_nd <- nodes[[right_id]]

    expect_equal(left_nd$n_obs,  n_left,
                 label = paste0("node ", left_id,  " n_obs = ", n_left,
                                " (", attr_nm, "<=", rule$cut_value, " side, canonical left)"))
    expect_equal(right_nd$n_obs, n_right,
                 label = paste0("node ", right_id, " n_obs = ", n_right,
                                " (", attr_nm, ">",  rule$cut_value, " side, canonical right)"))

    # Propagate observation indices to children for subsequent nodes
    obs_at[[left_id]]  <- idx[sides == 0L]
    obs_at[[right_id]] <- idx[sides == 1L]
    checked <- checked + 1L
  }

  expect_gt(checked, 0L, label = "at least one ordered_cut split was verified")
})
