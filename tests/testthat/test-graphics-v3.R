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
# T4: plot_lort_tree() returns ggplot from ort_plot_data() output
# ---------------------------------------------------------------------------
test_that("plot_lort_tree: returns ggplot from ort_plot_data() output (T4)", {
  gg_skip()
  p <- plot_lort_tree(gv3_opd)
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
