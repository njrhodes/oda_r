###############################################################################
# test-graphics-v3.R
#
# Tests for Graphics v3C1: plot_cta_tree() and plot_lort_tree()
#         Graphics v3C2: plot_oda_balance(), plot_smd_balance(),
#                        plot_balance_love(), plot_cta_balance()
#
# All tests skip if ggplot2 is not installed.
#
# Test plan (v3C1):
#   T1   plot_cta_tree() returns ggplot for a fitted cta_tree
#   T2   plot_cta_tree() returns ggplot from cta_plot_data() output
#   T3   plot_lort_tree() returns ggplot for a fitted cta_ort (LORT)
#   T4   plot_lort_tree() returns ggplot from ort_plot_data() output
#   T5   plot_cta_tree(): color_by "target_rate", "prediction", "none" all work
#   T6   plot_lort_tree(): show_n, show_percent, wrap_width accepted without error
#   T7   no_tree CTA input returns a ggplot message panel, not an error
#   T8   wrong input class errors clearly
#   T9   No plot function calls fitting functions (checked by contract: functions
#        accept only pre-fitted objects or pre-computed plot-data)
#  T10   plot_lort_tree() color_by "target_rate", "prediction", "none"
#
# Test plan (v3C2):
#  B1   plot_oda_balance() returns ggplot from oda_balance_plot_data
#  B2   plot_oda_balance() accepts oda_balance_table (auto-coerces)
#  B3   plot_smd_balance() returns ggplot from smd_balance_table
#  B4   plot_balance_love() returns ggplot from smd_balance_table (wrapper)
#  B5   plot_cta_balance() no_tree returns ggplot message panel
#  B6   plot_cta_balance() valid tree/stump returns ggplot tree
#  B7   wrong input class errors clearly for all four balance renderers
#  B8   balance render functions have no fitting args in their signatures
###############################################################################

gg_skip <- function() {
  skip_if_not_installed("ggplot2")
}

# ---- Module-level fixtures ------------------------------------------------- #

# Small two-class tree (stump: B alone separates groups)
gv3_X <- data.frame(
  A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
  B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
)
gv3_y <- c(rep(0L, 40), rep(1L, 20))

# CTA fit — non-recursive (used for T1, T2, T5, T7, T8)
gv3_cta <- cta_fit(gv3_X, gv3_y, mindenom = 5L,
                    mc_iter = 200L, mc_seed = 42L, loo = "off")

# cta_plot_data output — used for T2
gv3_cpd <- cta_plot_data(gv3_cta, target_class = 1L)

# no_tree CTA — alpha_split = 0 forces no splits accepted
gv3_notree <- cta_fit(gv3_X, gv3_y, mindenom = 5L,
                       mc_iter = 50L, mc_seed = 42L, loo = "off",
                       alpha_split = 0)

# LORT fit — used for T3, T4, T6, T10
gv3_lort <- lort_fit(gv3_X, gv3_y,
                      mc_iter = 200L, mc_seed = 42L, loo = "off", min_n = 5L)

# ort_plot_data output — used for T4
gv3_opd <- ort_plot_data(gv3_lort, target_class = 1L)

# ---------------------------------------------------------------------------
# T1: plot_cta_tree() returns ggplot for a fitted cta_tree
# ---------------------------------------------------------------------------
test_that("plot_cta_tree: returns ggplot for fitted cta_tree (T1)", {
  gg_skip()
  p <- plot_cta_tree(gv3_cta)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# T2: plot_cta_tree() returns ggplot from cta_plot_data() output
# ---------------------------------------------------------------------------
test_that("plot_cta_tree: returns ggplot from cta_plot_data() output (T2)", {
  gg_skip()
  p <- plot_cta_tree(gv3_cpd)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# T3: plot_lort_tree() returns ggplot for a fitted cta_ort
# ---------------------------------------------------------------------------
test_that("plot_lort_tree: returns ggplot for fitted cta_ort (T3)", {
  gg_skip()
  p <- plot_lort_tree(gv3_lort)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# T4: plot_lort_tree() index-based: index=1L renders root sub-tree (v4 API)
# ---------------------------------------------------------------------------
test_that("plot_lort_tree: index=1L renders root sub-tree (T4)", {
  gg_skip()
  p <- plot_lort_tree(gv3_lort, index = 1L)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# T5: plot_cta_tree(): all color_by values work without error
# ---------------------------------------------------------------------------
test_that("plot_cta_tree: color_by target_rate / prediction / none (T5)", {
  gg_skip()
  expect_s3_class(plot_cta_tree(gv3_cta, color_by = "target_rate"),  "ggplot")
  expect_s3_class(plot_cta_tree(gv3_cta, color_by = "prediction"),   "ggplot")
  expect_s3_class(plot_cta_tree(gv3_cta, color_by = "none"),         "ggplot")
})

# ---------------------------------------------------------------------------
# T6: plot_lort_tree(): show_n, show_percent, wrap_width accepted
# ---------------------------------------------------------------------------
test_that("plot_lort_tree: show_n, show_percent, wrap_width accepted (T6)", {
  gg_skip()
  p <- plot_lort_tree(gv3_lort,
                       show_n = FALSE, show_percent = FALSE, wrap_width = 20L)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# T7: no_tree CTA returns a ggplot message panel, not an error
# ---------------------------------------------------------------------------
test_that("plot_cta_tree: no_tree returns ggplot message panel (T7)", {
  gg_skip()
  p <- plot_cta_tree(gv3_notree)
  expect_s3_class(p, "ggplot")
  # The build should succeed (no error during rendering)
  expect_silent(ggplot2::ggplot_build(p))
})

# ---------------------------------------------------------------------------
# T8: wrong input class errors clearly
# ---------------------------------------------------------------------------
test_that("plot_cta_tree / plot_lort_tree: wrong input class errors (T8)", {
  gg_skip()
  expect_error(plot_cta_tree(list(a = 1)), "cta_tree")
  expect_error(plot_lort_tree(list(a = 1)), "cta_ort")
  # Passing a cta_tree to plot_lort_tree should also error
  expect_error(plot_lort_tree(gv3_cta), "cta_ort")
})

# ---------------------------------------------------------------------------
# T9: plot functions do not accept group/X/y — no fitting path exposed
# ---------------------------------------------------------------------------
test_that("plot_cta_tree / plot_lort_tree: no fitting args in signature (T9)", {
  gg_skip()
  # Formal args must not include y, group, w, mc_iter, mindenom
  cta_formals  <- names(formals(plot_cta_tree))
  lort_formals <- names(formals(plot_lort_tree))
  fitting_args <- c("y", "group", "mc_iter", "mindenom", "mc_seed")
  expect_true(!any(fitting_args %in% cta_formals),
              info = paste("Fitting args in plot_cta_tree:",
                           paste(intersect(fitting_args, cta_formals), collapse=", ")))
  expect_true(!any(fitting_args %in% lort_formals),
              info = paste("Fitting args in plot_lort_tree:",
                           paste(intersect(fitting_args, lort_formals), collapse=", ")))
})

# ---------------------------------------------------------------------------
# T10: plot_lort_tree(): all color_by values work without error
# ---------------------------------------------------------------------------
test_that("plot_lort_tree: color_by target_rate / prediction / none (T10)", {
  gg_skip()
  expect_s3_class(plot_lort_tree(gv3_lort, color_by = "target_rate"),  "ggplot")
  expect_s3_class(plot_lort_tree(gv3_lort, color_by = "prediction"),   "ggplot")
  expect_s3_class(plot_lort_tree(gv3_lort, color_by = "none"),         "ggplot")
})

###############################################################################
# v3C2 — Balance renderer fixtures and tests
###############################################################################

# ---- Module-level balance fixtures ----------------------------------------- #

# Reuse gv3_X / gv3_y from v3C1 above (same data, group = gv3_y).
# ODA balance table — mcarlo=FALSE for speed (testing renderer, not p-values)
gv3_bt  <- oda_balance_table(gv3_y, gv3_X,
                               mcarlo  = FALSE,
                               mc_iter = 50L)

# oda_balance_plot_data (no SMD join — keeps it simple and self-contained)
gv3_obal_pd <- oda_balance_plot_data(gv3_bt)

# SMD table
gv3_smd_tbl <- smd_balance_table(gv3_y, gv3_X)

# CTA balance table — valid tree (B perfectly separates groups)
gv3_ctb_disc <- cta_balance_table(gv3_y, gv3_X,
                                    mindenom = 5L,
                                    mc_iter  = 200L,
                                    mc_seed  = 42L)

# CTA balance plot_data — valid tree
gv3_cbal_pd <- cta_balance_plot_data(gv3_ctb_disc)

# CTA balance plot_data — no_tree (alpha_split=0 forces rejection of all splits)
gv3_ctb_notree <- cta_balance_table(gv3_y, gv3_X,
                                      mindenom    = 5L,
                                      mc_iter     = 50L,
                                      mc_seed     = 42L,
                                      alpha_split = 0)
gv3_cbal_notree <- cta_balance_plot_data(gv3_ctb_notree)

# ---------------------------------------------------------------------------
# B1: plot_oda_balance() returns ggplot from oda_balance_plot_data
# ---------------------------------------------------------------------------
test_that("plot_oda_balance: returns ggplot from oda_balance_plot_data (B1)", {
  gg_skip()
  p <- plot_oda_balance(gv3_obal_pd)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# B2: plot_oda_balance() accepts oda_balance_table (auto-coerces, no fitting)
# ---------------------------------------------------------------------------
test_that("plot_oda_balance: accepts oda_balance_table and auto-coerces (B2)", {
  gg_skip()
  p <- plot_oda_balance(gv3_bt)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# B3: plot_smd_balance() returns ggplot from smd_balance_table
# ---------------------------------------------------------------------------
test_that("plot_smd_balance: returns ggplot from smd_balance_table (B3)", {
  gg_skip()
  p <- plot_smd_balance(gv3_smd_tbl)
  expect_s3_class(p, "ggplot")
  # ref_020 variant also works
  expect_s3_class(plot_smd_balance(gv3_smd_tbl, ref_020 = TRUE), "ggplot")
})

# ---------------------------------------------------------------------------
# B4: plot_balance_love() returns ggplot (wrapper of plot_smd_balance)
# ---------------------------------------------------------------------------
test_that("plot_balance_love: returns ggplot (wrapper) (B4)", {
  gg_skip()
  p <- plot_balance_love(gv3_smd_tbl)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# B5: plot_cta_balance() no_tree returns ggplot message panel
# ---------------------------------------------------------------------------
test_that("plot_cta_balance: no_tree returns ggplot message panel (B5)", {
  gg_skip()
  p <- plot_cta_balance(gv3_cbal_notree)
  expect_s3_class(p, "ggplot")
  # Build succeeds (no error during render)
  expect_silent(ggplot2::ggplot_build(p))
})

# ---------------------------------------------------------------------------
# B6: plot_cta_balance() valid tree/stump returns ggplot tree
# ---------------------------------------------------------------------------
test_that("plot_cta_balance: valid tree/stump returns ggplot (B6)", {
  gg_skip()
  p <- plot_cta_balance(gv3_cbal_pd)
  expect_s3_class(p, "ggplot")
  # Also accepts cta_balance_table directly (auto-coerces, no fitting)
  p2 <- plot_cta_balance(gv3_ctb_disc)
  expect_s3_class(p2, "ggplot")
})

# ---------------------------------------------------------------------------
# B7: wrong input class errors clearly for all four balance renderers
# ---------------------------------------------------------------------------
test_that("balance renderers: wrong input class errors clearly (B7)", {
  gg_skip()
  expect_error(plot_oda_balance(list()),      "oda_balance")
  expect_error(plot_smd_balance(list()),      "smd_balance_table")
  expect_error(plot_balance_love(list()),     "smd_balance_table")
  expect_error(plot_cta_balance(list()),      "cta_balance")
  # Passing wrong balance type to each function also errors
  expect_error(plot_oda_balance(gv3_smd_tbl), "oda_balance")
  expect_error(plot_smd_balance(gv3_obal_pd), "smd_balance_table")
  expect_error(plot_cta_balance(gv3_obal_pd), "cta_balance")
})

# ---------------------------------------------------------------------------
# B8: balance render functions have no fitting args in their signatures
# ---------------------------------------------------------------------------
test_that("balance renderers: no fitting args in formal signatures (B8)", {
  gg_skip()
  fitting_args <- c("y", "group", "mc_iter", "mindenom", "mc_seed")
  fns <- list(
    plot_oda_balance  = names(formals(plot_oda_balance)),
    plot_smd_balance  = names(formals(plot_smd_balance)),
    plot_balance_love = names(formals(plot_balance_love)),
    plot_cta_balance  = names(formals(plot_cta_balance))
  )
  for (fn_nm in names(fns)) {
    found <- intersect(fitting_args, fns[[fn_nm]])
    expect_true(length(found) == 0L,
                info = paste("Fitting args in", fn_nm, ":",
                             paste(found, collapse = ", ")))
  }
})

###############################################################################
# v3C3 - Evidence-interval balance builder and renderer tests
#
# Test plan:
#  E1  oda_balance_effect_table() returns correct class and row structure
#  E2  oda_balance_effect_table() compare_weights=TRUE produces unweighted
#      and weighted rows
#  E3  oda_balance_effect_table() balanced_by_interval / residual_imbalance
#      are logical
#  E4  plot_oda_balance_effects() returns ggplot from
#      oda_balance_effect_table
#  E5  plot_oda_balance_effects() wrong class errors clearly
#  E6  cta_balance_effect_summary() returns correct class and row structure
#  E7  cta_balance_effect_summary() no_tree result has estimate = 0
#  E8  plot_cta_balance_effects() returns ggplot from
#      cta_balance_effect_summary
#  E9  plot_cta_balance_effects() wrong class errors clearly
#  E10 plot_oda_balance_effects() and plot_cta_balance_effects() have no
#      fitting args in their signatures
###############################################################################

# ---- Module-level fixtures (v3C3) ------------------------------------------ #

# Reuse gv3_X / gv3_y from v3C1 (same deterministic data).
# Use small nboot/chance_iter so tests run in < 5s each.
gv3_eff_tbl <- oda_balance_effect_table(
  gv3_y, gv3_X,
  nboot       = 50L,
  chance_iter = 50L,
  mc_iter     = 200L,
  mc_seed     = 42L
)

gv3_cta_eff <- cta_balance_effect_summary(
  gv3_y, gv3_X,
  mindenom    = 5L,
  nboot       = 20L,
  chance_iter = 20L,
  mc_iter     = 200L,
  mc_seed     = 42L
)

# No-tree summary: alpha_split = 0 forces no splits accepted
gv3_cta_eff_nt <- cta_balance_effect_summary(
  gv3_y, gv3_X,
  mindenom    = 5L,
  nboot       = 10L,
  chance_iter = 10L,
  mc_iter     = 50L,
  mc_seed     = 42L,
  alpha_split = 0
)

# ---------------------------------------------------------------------------
# E1: oda_balance_effect_table() class and structure
# ---------------------------------------------------------------------------
test_that("oda_balance_effect_table: class and structure (E1)", {
  expect_s3_class(gv3_eff_tbl, "oda_balance_effect_table")
  expect_true(is.data.frame(gv3_eff_tbl$rows))
  # One row per covariate (no compare_weights)
  expect_equal(nrow(gv3_eff_tbl$rows), ncol(gv3_X))
  req_cols <- c("attribute", "analysis", "metric", "estimate",
                "boot_lo", "boot_hi", "chance_lo", "chance_hi",
                "p_mc", "p_sidak", "p_bonferroni",
                "rule_summary", "sensitivity", "specificity", "n_total",
                "balanced_by_interval", "residual_imbalance")
  for (col in req_cols)
    expect_true(col %in% names(gv3_eff_tbl$rows),
                info = paste("Missing column:", col))
  # All analysis values are "unweighted"
  expect_true(all(gv3_eff_tbl$rows$analysis == "unweighted"))
})

# ---------------------------------------------------------------------------
# E2: compare_weights = TRUE produces both unweighted and weighted rows
# ---------------------------------------------------------------------------
test_that("oda_balance_effect_table: compare_weights=TRUE adds weighted rows (E2)", {
  et2 <- oda_balance_effect_table(
    gv3_y, gv3_X,
    w           = rep(1, length(gv3_y)),
    compare_weights = TRUE,
    nboot       = 10L,
    chance_iter = 10L,
    mc_iter     = 50L,
    mc_seed     = 42L
  )
  analyses <- unique(et2$rows$analysis)
  expect_true("unweighted" %in% analyses)
  expect_true("weighted"   %in% analyses)
  expect_equal(nrow(et2$rows), 2L * ncol(gv3_X))
})

# ---------------------------------------------------------------------------
# E3: interval flag columns are logical
# ---------------------------------------------------------------------------
test_that("oda_balance_effect_table: interval flags are logical (E3)", {
  expect_true(is.logical(gv3_eff_tbl$rows$balanced_by_interval))
  expect_true(is.logical(gv3_eff_tbl$rows$residual_imbalance))
})

# ---------------------------------------------------------------------------
# E4: plot_oda_balance_effects() returns ggplot
# ---------------------------------------------------------------------------
test_that("plot_oda_balance_effects: returns ggplot (E4)", {
  gg_skip()
  p <- plot_oda_balance_effects(gv3_eff_tbl)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# E5: plot_oda_balance_effects() wrong class errors clearly
# ---------------------------------------------------------------------------
test_that("plot_oda_balance_effects: wrong class errors (E5)", {
  gg_skip()
  expect_error(plot_oda_balance_effects(list()), "oda_balance_effect_table")
  expect_error(plot_oda_balance_effects(gv3_bt),  "oda_balance_effect_table")
})

# ---------------------------------------------------------------------------
# E6: cta_balance_effect_summary() class and structure
# ---------------------------------------------------------------------------
test_that("cta_balance_effect_summary: class and structure (E6)", {
  expect_s3_class(gv3_cta_eff, "cta_balance_effect_summary")
  expect_true(is.data.frame(gv3_cta_eff$rows))
  # Tree-level schema: one row per analysis variant (unweighted / weighted)
  req_cols <- c("analysis", "metric", "estimate",
                "boot_lo", "boot_hi", "chance_lo", "chance_hi",
                "d_stat", "n_endpoints", "root_attribute",
                "status", "balance_interpretation")
  for (col in req_cols)
    expect_true(col %in% names(gv3_cta_eff$rows),
                info = paste("Missing column:", col))
  # At least one row returned
  expect_gte(nrow(gv3_cta_eff$rows), 1L)
})

# ---------------------------------------------------------------------------
# E7: cta_balance_effect_summary() no_tree result has estimate = 0
# ---------------------------------------------------------------------------
test_that("cta_balance_effect_summary: no_tree estimate = 0 (E7)", {
  expect_true(all(gv3_cta_eff_nt$rows$estimate == 0))
})

# ---------------------------------------------------------------------------
# E8: plot_cta_balance_effects() returns ggplot
# ---------------------------------------------------------------------------
test_that("plot_cta_balance_effects: returns ggplot (E8)", {
  gg_skip()
  p <- plot_cta_balance_effects(gv3_cta_eff)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# E9: plot_cta_balance_effects() wrong class / no_tree handled
# ---------------------------------------------------------------------------
test_that("plot_cta_balance_effects: wrong class errors; no_tree is ggplot (E9)", {
  gg_skip()
  expect_error(plot_cta_balance_effects(list()), "cta_balance_effect_summary")
  # no_tree summary should still render (message panel or similar)
  p_nt <- plot_cta_balance_effects(gv3_cta_eff_nt)
  expect_s3_class(p_nt, "ggplot")
})

# ---------------------------------------------------------------------------
# E10: evidence-interval renderers have no fitting args in signatures
# ---------------------------------------------------------------------------
test_that("evidence renderers: no fitting args in signatures (E10)", {
  gg_skip()
  fitting_args <- c("X", "y", "w", "mindenom", "mc_iter", "mc_seed",
                    "alpha_split", "prune_alpha", "loo", "miss_codes",
                    "nboot", "chance_iter")
  fns <- list(
    plot_oda_balance_effects = plot_oda_balance_effects,
    plot_cta_balance_effects = plot_cta_balance_effects
  )
  for (fn_nm in names(fns)) {
    found <- intersect(names(formals(fns[[fn_nm]])), fitting_args)
    expect_equal(length(found), 0L,
                 info = paste("Fitting args in", fn_nm, ":",
                              paste(found, collapse = ", ")))
  }
})

###############################################################################
# v3C5 — LORT path navigation API (L-series)
#  L1   lort_index_path() returns correct path data.frame
#  L2   lort_local_tree() returns cta_tree for non-terminal nodes
#  L3   at least one local CTA on path has >2 leaves
#  L4a  plot_lort_path(layout="list") returns named ggplot list
#  L4b  plot_lort_path(layout="multipanel") returns single patchwork object
#  L4c  terminal/no_tree panels are styled result cards, not black voids
#  L4d  internal split nodes rendered as ellipses (GeomPolygon layer)
#  L5   lort_path_table() prints and returns path df invisibly
###############################################################################

# Module-level LORT fixture: build_block data (n=140), deterministic seed.
# Root CTA: V1 stump (2 endpoints).
# Node 2 (V1=1 arm): V2+V3 3-endpoint tree.
# Path to node 5: root (1) -> V1=1 child (2) -> terminal (5).
.lort_build_block <- function(v1, v2, v3, n0, n1) {
  data.frame(
    V1    = rep(v1, n0 + n1),
    V2    = rep(v2, n0 + n1),
    V3    = rep(v3, n0 + n1),
    group = c(rep(0L, n0), rep(1L, n1)),
    stringsAsFactors = FALSE
  )
}
.lort_dat <- rbind(
  .lort_build_block(0L, 0L, 0L, 38L, 2L),
  .lort_build_block(0L, 1L, 0L, 22L, 18L),
  .lort_build_block(1L, 0L, 0L, 18L, 2L),
  .lort_build_block(1L, 0L, 1L, 2L,  18L),
  .lort_build_block(1L, 1L, 0L, 0L,  20L)
)
.lort_group   <- .lort_dat$group
.lort_X       <- .lort_dat[, c("V1", "V2", "V3")]
gv3_lort4 <- lort_fit(.lort_X, .lort_group,
                       mc_iter  = 1000L,
                       mc_seed  = 42L,
                       loo      = "off",
                       min_n    = 5L)

# ---------------------------------------------------------------------------
# L1: lort_index_path() returns path data.frame with correct structure
# ---------------------------------------------------------------------------
test_that("lort_index_path: path to index 5 has correct structure (L1)", {
  path5 <- lort_index_path(gv3_lort4, 5L)
  expect_s3_class(path5, "data.frame")
  expect_equal(nrow(path5), 3L)        # root, child, terminal
  expect_equal(path5$lort_index, c(1L, 2L, 5L))
  expect_equal(path5$depth, c(0L, 1L, 2L))
  expect_true(path5$is_terminal[[3L]])    # node 5 is terminal
  expect_false(path5$is_terminal[[1L]])   # root is not terminal
  expect_false(path5$is_terminal[[2L]])   # node 2 is not terminal
  # incoming path condition for last hop mentions V2
  expect_true(grepl("V2", path5$incoming_path_condition[[3L]]))
})

# ---------------------------------------------------------------------------
# L2: lort_local_tree() returns cta_tree for non-terminal path nodes
# ---------------------------------------------------------------------------
test_that("lort_local_tree: non-terminal nodes return cta_tree (L2)", {
  path5 <- lort_index_path(gv3_lort4, 5L)
  non_term_ids <- path5$lort_index[!path5$is_terminal]
  for (k in non_term_ids) {
    tr <- lort_local_tree(gv3_lort4, k)
    expect_true(inherits(tr, "cta_tree"),
                label = paste("Node", k, "returns cta_tree"))
  }
})

# ---------------------------------------------------------------------------
# L3: at least one local CTA on path has >2 leaves
# ---------------------------------------------------------------------------
test_that("lort_local_tree: at least one path node has >2 leaves (L3)", {
  path5  <- lort_index_path(gv3_lort4, 5L)
  trees  <- lapply(path5$lort_index, function(k) lort_local_tree(gv3_lort4, k))
  has_multi <- any(vapply(trees, function(tr) {
    !is.null(tr) && !isTRUE(tr$no_tree) &&
      sum(vapply(tr$nodes, function(nd) isTRUE(nd$leaf), logical(1L))) > 2L
  }, logical(1L)))
  expect_true(has_multi)
})

# ---------------------------------------------------------------------------
# L4a: plot_lort_path(layout="list") returns named list of ggplots
# ---------------------------------------------------------------------------
test_that("plot_lort_path: layout='list' returns named ggplot list (L4a)", {
  gg_skip()
  path5  <- lort_index_path(gv3_lort4, 5L)
  plots  <- plot_lort_path(gv3_lort4, index = 5L, layout = "list")
  expect_type(plots, "list")
  expect_equal(length(plots), nrow(path5))
  for (nm in names(plots))
    expect_true(inherits(plots[[nm]], "ggplot"),
                label = paste("Element", nm, "is ggplot"))
  expect_true(all(grepl("^index_[0-9]+$", names(plots))))
})

# ---------------------------------------------------------------------------
# L4b: plot_lort_path(layout="multipanel") returns single patchwork object
# ---------------------------------------------------------------------------
test_that("plot_lort_path: layout='multipanel' returns single patchwork object (L4b)", {
  gg_skip()
  if (!requireNamespace("patchwork", quietly = TRUE))
    skip("patchwork not installed")
  p <- plot_lort_path(gv3_lort4, index = 5L, layout = "multipanel")
  # patchwork objects inherit from ggplot or patchwork class
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"),
              label = "multipanel result is a gg/patchwork object")
})

# ---------------------------------------------------------------------------
# L4c: terminal/no_tree panels are result cards, not black voids
# ---------------------------------------------------------------------------
test_that("plot_lort_path: terminal panel is a result card not a void (L4c)", {
  gg_skip()
  # Node 5 is terminal (no_tree). Its panel should use .gg_result_card,
  # which has a border rect layer. Proxy: the plot has >1 layer.
  plots <- plot_lort_path(gv3_lort4, index = 5L, layout = "list")
  term_p <- plots[["index_5"]]
  expect_true(inherits(term_p, "ggplot"))
  expect_gt(length(term_p$layers), 1L)
})

# ---------------------------------------------------------------------------
# L4d: CTA plot-data for a non-trivial tree uses ellipse polygon for split nodes
# ---------------------------------------------------------------------------
test_that("plot_cta_tree: internal nodes rendered as ellipses (L4d)", {
  gg_skip()
  # Node 2 has a 3-endpoint V2+V3 tree -- has internal split nodes
  tr2 <- lort_local_tree(gv3_lort4, 2L)
  p   <- plot_cta_tree(tr2)
  # geom_polygon is used for ellipses; at least one polygon layer exists
  layer_classes <- vapply(p$layers,
                          function(l) class(l$geom)[1L], character(1L))
  expect_true("GeomPolygon" %in% layer_classes,
              label = "Internal nodes rendered via geom_polygon (ellipse)")
})

# ---------------------------------------------------------------------------
# L5: lort_path_table() prints and returns path df invisibly  (was L6)
# ---------------------------------------------------------------------------
test_that("lort_path_table: prints and returns path df (L5)", {
  expect_output(result <- lort_path_table(gv3_lort4, 5L), "LORT path")
  expect_s3_class(result, "data.frame")
  expect_equal(result$lort_index, c(1L, 2L, 5L))
})

###############################################################################
# v3C6 -- plot_cta_family() API (F-series)
#  F1   plot_cta_family(index=1L) returns a single ggplot
#  F2   plot_cta_family(min_d=TRUE) returns a single ggplot (min-D member)
#  F3   plot_cta_family(show_all=TRUE, layout="list") returns named list
#  F4   plot_cta_family(show_all=TRUE, layout="multipanel") returns patchwork
#  F5   cta_plot_data() nodes have p_mc and loo_p columns
###############################################################################

# Module-level CTA family fixture (synthetic, small, fast)
.fam_X <- data.frame(
  x1 = c(rep(0L, 30), rep(1L, 30)),
  x2 = c(rep(0L, 15), rep(1L, 15), rep(0L, 15), rep(1L, 15))
)
.fam_y <- c(rep(0L, 45), rep(1L, 15))
gv3_fam <- cta_descendant_family(.fam_X, .fam_y,
  start_mindenom = 5L,
  max_steps      = 3L,
  mc_iter        = 200L,
  mc_seed        = 42L,
  loo            = "off",
  attr_names     = c("x1", "x2"))

# CTA fit with LOO for F5 loo_p column test
gv3_cta_loo <- cta_fit(gv3_X, gv3_y, mindenom = 5L,
                        mc_iter = 500L, mc_seed = 42L, loo = "pvalue")

# ---------------------------------------------------------------------------
# F1: plot_cta_family(index=N) returns single ggplot
# ---------------------------------------------------------------------------
test_that("plot_cta_family: index=1L returns ggplot (F1)", {
  gg_skip()
  p <- plot_cta_family(gv3_fam, index = 1L)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# F2: plot_cta_family(min_d=TRUE) returns single ggplot for min-D member
# ---------------------------------------------------------------------------
test_that("plot_cta_family: min_d=TRUE returns ggplot (F2)", {
  gg_skip()
  p <- plot_cta_family(gv3_fam, min_d = TRUE)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# F3: show_all + layout="list" returns named list
# ---------------------------------------------------------------------------
test_that("plot_cta_family: show_all + layout='list' returns named list (F3)", {
  gg_skip()
  pl <- plot_cta_family(gv3_fam, show_all = TRUE, layout = "list")
  expect_true(is.list(pl))
  expect_true(length(pl) == length(gv3_fam$members))
  expect_true(all(vapply(pl, function(p) inherits(p, "ggplot"), logical(1L))))
})

# ---------------------------------------------------------------------------
# F4: show_all + layout="multipanel" returns patchwork object
# ---------------------------------------------------------------------------
test_that("plot_cta_family: show_all + layout='multipanel' returns patchwork (F4)", {
  gg_skip()
  skip_if_not_installed("patchwork")
  pp <- plot_cta_family(gv3_fam, show_all = TRUE, layout = "multipanel", ncol = 1L)
  expect_true(inherits(pp, "patchwork") || inherits(pp, "ggplot"))
})

# ---------------------------------------------------------------------------
# F5: cta_plot_data() nodes have p_mc and loo_p columns
# ---------------------------------------------------------------------------
test_that("cta_plot_data: nodes have p_mc, loo_p, loo_status columns (F5)", {
  pd_base <- cta_plot_data(gv3_cta)
  expect_true("p_mc"       %in% names(pd_base$nodes))
  expect_true("loo_p"      %in% names(pd_base$nodes))
  expect_true("loo_status" %in% names(pd_base$nodes))
  # Split nodes have numeric p_mc; leaves are NA
  split_rows <- pd_base$nodes[!pd_base$nodes$leaf, ]
  leaf_rows  <- pd_base$nodes[ pd_base$nodes$leaf, ]
  if (nrow(split_rows) > 0L)
    expect_true(all(!is.na(split_rows$p_mc)))
  if (nrow(leaf_rows) > 0L) {
    expect_true(all(is.na(leaf_rows$p_mc)))
    expect_true(all(is.na(leaf_rows$loo_status)))
  }
  # With LOO pvalue, split nodes have loo_p and loo_status == "PVALUE"
  pd_loo <- cta_plot_data(gv3_cta_loo)
  split_loo <- pd_loo$nodes[!pd_loo$nodes$leaf, ]
  if (nrow(split_loo) > 0L) {
    expect_true(all(!is.na(split_loo$loo_p)))
    expect_true(all(split_loo$loo_status == "PVALUE"))
  }
  # loo = "off" -> loo_status == "OFF" on split nodes
  split_off <- pd_base$nodes[!pd_base$nodes$leaf, ]
  if (nrow(split_off) > 0L)
    expect_true(all(split_off$loo_status == "OFF"))
})

###############################################################################
# v3C7 -- Canon wording and shape tests (W-series)
#  W1   CTA no-tree card header says "No CTA Tree Found" (not LORT wording)
#  W2   LORT terminal panel header says "No Local CTA Tree" with LORT index footer
#  W3   Endpoint/leaf nodes are rendered as GeomRect rectangles
###############################################################################

# Helper: extract all annotation text labels from a ggplot object.
# In ggplot2 >= 3.4 annotate() stores non-positional aesthetics (e.g. label)
# in lyr$aes_params, NOT in lyr$data.
.anno_labels <- function(p) {
  out <- character(0L)
  for (lyr in p$layers) {
    if (!is.null(lyr$aes_params) && "label" %in% names(lyr$aes_params))
      out <- c(out, as.character(lyr$aes_params[["label"]]))
  }
  out
}

# ---------------------------------------------------------------------------
# W1: CTA no-tree card header says "No CTA Tree Found"
# ---------------------------------------------------------------------------
test_that("plot_cta_tree: no-tree card header is 'No CTA Tree Found' (W1)", {
  gg_skip()
  p      <- plot_cta_tree(gv3_notree)
  labels <- .anno_labels(p)
  expect_true(any(grepl("No CTA Tree Found", labels, fixed = TRUE)),
              label = "Header annotation must be 'No CTA Tree Found'")
  expect_false(any(grepl("No Local CTA Tree", labels, fixed = TRUE)),
               label = "CTA no-tree card must not say 'No Local CTA Tree'")
})

# ---------------------------------------------------------------------------
# W2: LORT terminal panel header says "No Local CTA Tree" with LORT index footer
# ---------------------------------------------------------------------------
test_that("plot_lort_path: terminal panel uses LORT-specific wording (W2)", {
  gg_skip()
  plots   <- plot_lort_path(gv3_lort4, index = 5L, layout = "list")
  term_p  <- plots[["index_5"]]
  labels  <- .anno_labels(term_p)
  expect_true(any(grepl("No Local CTA Tree", labels, fixed = TRUE)),
              label = "LORT terminal panel header must be 'No Local CTA Tree'")
  # Footer must mention the LORT index number
  expect_true(any(grepl("LORT index", labels, fixed = TRUE)),
              label = "LORT terminal panel footer must mention 'LORT index'")
})

# ---------------------------------------------------------------------------
# W3: Leaf/endpoint nodes are rendered as GeomRect rectangles
# ---------------------------------------------------------------------------
test_that("plot_cta_tree: endpoint nodes rendered as rectangles (W3)", {
  gg_skip()
  p           <- plot_cta_tree(gv3_cta, color_by = "none")
  layer_geoms <- vapply(p$layers,
                        function(lyr) class(lyr$geom)[1L], character(1L))
  expect_true(any(layer_geoms == "GeomRect"),
              label = "At least one GeomRect layer must exist for endpoints")
})
