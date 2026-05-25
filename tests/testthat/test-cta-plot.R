###############################################################################
# test-cta-plot.R — unit tests for cta_plot_data() and plot.cta_tree()
###############################################################################

# --- shared fixtures ---------------------------------------------------------

make_tree <- function() {
  data(mtcars, envir = environment())
  X <- mtcars[, c("cyl", "disp", "hp", "wt")]
  y <- as.integer(mtcars$am)
  suppressMessages(
    oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L,
                loo = "off")
  )
}

make_no_tree <- function() {
  # mindenom large enough to prevent any split
  data(mtcars, envir = environment())
  X <- mtcars[, c("cyl", "disp", "hp", "wt")]
  y <- as.integer(mtcars$am)
  suppressMessages(
    oda_cta_fit(X, y, mindenom = 9999L, mc_iter = 50L, mc_seed = 42L,
                loo = "off")
  )
}

# =============================================================================
# cta_plot_data() — structural contract (target_class = NULL)
# =============================================================================

test_that("cta_plot_data returns a list with correct names", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  expect_type(pd, "list")
  expect_setequal(names(pd), c("nodes", "edges", "no_tree", "has_weights"))
})

test_that("nodes data.frame has required columns", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  required <- c("node_id", "parent_id", "depth", "x", "y", "leaf",
                "attribute", "n_obs", "majority_class", "ess", "label")
  expect_true(all(required %in% names(pd$nodes)))
})

test_that("edges data.frame has required columns", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  required <- c("from_node_id", "to_node_id", "x0", "y0", "x1", "y1", "label")
  expect_true(all(required %in% names(pd$edges)))
})

test_that("nodes types are correct", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  nd   <- pd$nodes
  expect_type(nd$node_id,        "integer")
  expect_type(nd$parent_id,      "integer")
  expect_type(nd$depth,          "integer")
  expect_type(nd$x,              "double")
  expect_type(nd$y,              "double")
  expect_type(nd$leaf,           "logical")
  expect_type(nd$attribute,      "character")
  expect_type(nd$n_obs,          "integer")
  expect_type(nd$majority_class, "integer")
  expect_type(nd$ess,            "double")
  expect_type(nd$label,          "character")
})

test_that("edges types are correct", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  ed   <- pd$edges
  expect_type(ed$from_node_id, "integer")
  expect_type(ed$to_node_id,   "integer")
  expect_type(ed$x0,           "double")
  expect_type(ed$y0,           "double")
  expect_type(ed$x1,           "double")
  expect_type(ed$y1,           "double")
  expect_type(ed$label,        "character")
})

test_that("has_weights is FALSE for unit-weight fit", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  expect_false(pd$has_weights)
})

test_that("no_tree is FALSE for a real tree", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  expect_false(pd$no_tree)
})

test_that("nodes has at least 3 rows (root + 2 leaves) for a non-trivial tree", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  expect_gte(nrow(pd$nodes), 3L)
})

test_that("edges has at least 2 rows for a non-trivial tree", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  expect_gte(nrow(pd$edges), 2L)
})

test_that("leaf nodes have NA ess and NA attribute", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  nd   <- pd$nodes
  leaves <- nd[nd$leaf, ]
  expect_true(all(is.na(leaves$ess)))
  expect_true(all(is.na(leaves$attribute)))
})

test_that("split nodes have non-NA attribute", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  nd   <- pd$nodes
  splits <- nd[!nd$leaf, ]
  if (nrow(splits) > 0L)
    expect_true(all(!is.na(splits$attribute)))
})

test_that("leaf y = -depth and x positions are distinct numeric", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  nd   <- pd$nodes
  # y must equal -depth for every node
  expect_equal(nd$y, -nd$depth)
  # leaves have distinct x positions
  leaf_x <- nd$x[nd$leaf]
  expect_equal(length(leaf_x), length(unique(leaf_x)))
})

test_that("edge from/to node_ids match nodes in the tree", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  all_ids <- pd$nodes$node_id
  expect_true(all(pd$edges$from_node_id %in% all_ids))
  expect_true(all(pd$edges$to_node_id   %in% all_ids))
})

test_that("edge labels are non-empty strings", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  expect_true(all(nchar(pd$edges$label) > 0L))
})

test_that("leaf label contains 'class=' and 'n='", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  nd   <- pd$nodes
  leaves <- nd[nd$leaf, ]
  expect_true(all(grepl("class=", leaves$label)))
  expect_true(all(grepl("n=",     leaves$label)))
})

test_that("split label contains ESS or WESS", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  nd   <- pd$nodes
  splits <- nd[!nd$leaf, ]
  if (nrow(splits) > 0L)
    expect_true(all(grepl("ESS", splits$label)))
})

# =============================================================================
# cta_plot_data() — no-tree / stump edge cases
# =============================================================================

test_that("no-tree fit returns empty nodes and edges", {
  tree <- make_no_tree()
  pd   <- cta_plot_data(tree)
  expect_true(pd$no_tree)
  expect_equal(nrow(pd$nodes), 0L)
  expect_equal(nrow(pd$edges), 0L)
})

test_that("stump has exactly 3 nodes and 2 edges", {
  # Use mindenom just below half the data to force a stump
  data(mtcars, envir = environment())
  X <- mtcars[, c("cyl", "disp", "hp", "wt")]
  y <- as.integer(mtcars$am)
  tree <- suppressMessages(
    oda_cta_fit(X, y, mindenom = 14L, mc_iter = 500L, mc_seed = 42L,
                loo = "off")
  )
  pd <- cta_plot_data(tree)
  status <- .cta_tree_status(tree)
  if (identical(status, "stump")) {
    expect_equal(nrow(pd$nodes), 3L)
    expect_equal(nrow(pd$edges), 2L)
  } else {
    skip("tree is not a stump at this mindenom — skip stump-specific counts")
  }
})

# =============================================================================
# plot.cta_tree() — smoke tests (no errors, returns invisible pd list)
# =============================================================================

test_that("plot.cta_tree does not error on a valid tree", {
  tree <- make_tree()
  tmp  <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp, width = 600, height = 400)
  expect_no_error(plot(tree))
  grDevices::dev.off()
})

test_that("plot.cta_tree returns invisible cta_plot_data list", {
  tree <- make_tree()
  tmp  <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp, width = 600, height = 400)
  result <- plot(tree)
  grDevices::dev.off()
  expect_type(result, "list")
  expect_true(all(c("nodes", "edges", "no_tree", "has_weights") %in% names(result)))
})

test_that("plot.cta_tree does not error on no-tree fit", {
  tree <- make_no_tree()
  tmp  <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp, width = 400, height = 300)
  expect_no_error(plot(tree))
  grDevices::dev.off()
})

test_that("plot.cta_tree accepts main argument", {
  tree <- make_tree()
  tmp  <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp, width = 600, height = 400)
  expect_no_error(plot(tree, main = "Custom Title"))
  grDevices::dev.off()
})

# =============================================================================
# cta_plot_data() v2 — target-class enrichment fast-tier tests
# =============================================================================

test_that("cta_plot_data with target_class adds enrichment columns to nodes", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree, target_class = 1L)
  nd   <- pd$nodes
  enrichment_cols <- c("endpoint_id", "stage", "target_class", "target_n",
                        "denominator", "target_proportion", "target_rank",
                        "endpoint_fill_color", "predicted_label",
                        "target_label", "endpoint_label")
  expect_true(all(enrichment_cols %in% names(nd)))
})

test_that("cta_plot_data with target_class returns endpoints data.frame", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree, target_class = 1L)
  expect_true("endpoints" %in% names(pd))
  expect_s3_class(pd$endpoints, "data.frame")
})

test_that("cta_plot_data with target_class sets target_class_used", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree, target_class = 1L)
  expect_equal(pd$target_class_used, 1L)
})

test_that("endpoint_fill_color is character hex on leaf nodes", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree, target_class = 1L)
  nd   <- pd$nodes
  leaf_colors <- nd$endpoint_fill_color[nd$leaf]
  expect_type(leaf_colors, "character")
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", leaf_colors)))
})

test_that("target_rank is integer and assigns each leaf a unique rank", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree, target_class = 1L)
  nd   <- pd$nodes
  leaf_ranks <- nd$target_rank[nd$leaf]
  expect_type(leaf_ranks, "integer")
  n_ep <- sum(nd$leaf)
  expect_setequal(leaf_ranks, seq_len(n_ep))
})

test_that("target_proportion is numeric and in [0, 1] on leaves", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree, target_class = 1L)
  nd   <- pd$nodes
  props <- nd$target_proportion[nd$leaf]
  expect_type(props, "double")
  expect_true(all(props >= 0 & props <= 1))
})

test_that("no column named fill_key exists in nodes", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree, target_class = 1L)
  expect_false("fill_key" %in% names(pd$nodes))
})

test_that("class_labels named vector is respected in predicted_label", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree, target_class = 1L,
                         class_labels = c("0" = "Manual", "1" = "Auto"))
  nd   <- pd$nodes
  pred_labels <- nd$predicted_label[nd$leaf & !is.na(nd$predicted_label)]
  expect_true(all(pred_labels %in% c("Manual", "Auto")))
})

test_that("palette function(n) is accepted without error", {
  tree  <- make_tree()
  my_pal <- function(n) grDevices::colorRampPalette(c("#000000", "#ffffff"))(n)
  expect_no_error(cta_plot_data(tree, target_class = 1L, endpoint_palette = my_pal))
})

test_that("palette character vector changes endpoint_fill_color", {
  tree      <- make_tree()
  pd_def    <- cta_plot_data(tree, target_class = 1L)
  pd_custom <- cta_plot_data(tree, target_class = 1L,
                              endpoint_palette = c("#000000", "#ffffff"))
  nd_def    <- pd_def$nodes[pd_def$nodes$leaf, ]
  nd_cust   <- pd_custom$nodes[pd_custom$nodes$leaf, ]
  # Custom palette must differ from default on at least one leaf
  expect_false(all(nd_def$endpoint_fill_color == nd_cust$endpoint_fill_color))
})

test_that("no-tree fit with target_class does not error", {
  tree <- make_no_tree()
  expect_no_error(cta_plot_data(tree, target_class = 1L))
})

test_that("legacy node_col_split/node_col_leaf args accepted without error", {
  tree <- make_tree()
  tmp  <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp, width = 600, height = 400)
  expect_no_error(
    plot(tree, node_col_split = "#aabbcc", node_col_leaf = "#ddeeff")
  )
  grDevices::dev.off()
})

test_that("plot.cta_tree with target_class returns pd with target_class_used", {
  tree <- make_tree()
  tmp  <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp, width = 700, height = 500)
  result <- plot(tree, target_class = 1L)
  grDevices::dev.off()
  expect_equal(result$target_class_used, 1L)
})

# =============================================================================
# cta_plot_data() v2 — myeloma slow-tier tests
# =============================================================================

# ---- Myeloma fixture helpers ------------------------------------------------

.plot_myeloma_fixture_path <- function(f) {
  tryCatch(testthat::test_path(file.path("fixtures", "myeloma", f)),
           error = function(e) "")
}

.plot_myeloma_fixtures_ok <- function() {
  p <- .plot_myeloma_fixture_path("data.txt")
  nzchar(p) && file.exists(p)
}

.plot_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.plot_myeloma_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      if (!.plot_myeloma_fixtures_ok()) return(NULL)
      f  <- .plot_myeloma_fixture_path("data.txt")
      d  <- read.table(f, header = FALSE)
      names(d) <- paste0("V", seq_len(ncol(d)))
      d  <- d[d$V2 != 0, ]
      fit <<- suppressMessages(
        oda_cta_fit(
          X           = d[, .plot_myeloma_attrs],
          y           = d$V1,
          w           = d$V2,
          priors_on   = TRUE,
          miss_codes  = -9,
          alpha_split = 0.05,
          mindenom    = 1L,
          prune_alpha = 0.05,
          max_depth   = 20L,
          ess_min     = 0,
          mc_iter     = 5000L,
          mc_target   = 0.05,
          mc_stop     = 99.9,
          mc_stopup   = 99.9,
          mc_seed     = NULL,
          loo         = "stable",
          attr_names  = .plot_myeloma_attrs
        )
      )
    }
    fit
  }
})

test_that("myeloma plot: cta_plot_data with target_class=1 returns 3 endpoints", {
  skip_if_slow_tests_disabled("cta-plot-myeloma")
  if (!.plot_myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .plot_myeloma_fit()
  if (is.null(tree)) skip("myeloma fit returned NULL")
  pd <- cta_plot_data(tree, target_class = 1L)
  expect_equal(nrow(pd$endpoints), 3L)
})

test_that("myeloma plot: endpoint proportions match staging table", {
  skip_if_slow_tests_disabled("cta-plot-myeloma")
  if (!.plot_myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .plot_myeloma_fit()
  if (is.null(tree)) skip("myeloma fit returned NULL")
  pd <- cta_plot_data(tree, target_class = 1L)
  props <- sort(pd$endpoints$target_proportion)
  st    <- cta_staging_table(tree, target_class = 1L)
  st_props <- sort(st$target_proportion)
  expect_equal(props, st_props)
})

test_that("myeloma plot: endpoint stages are 1, 2, 3", {
  skip_if_slow_tests_disabled("cta-plot-myeloma")
  if (!.plot_myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .plot_myeloma_fit()
  if (is.null(tree)) skip("myeloma fit returned NULL")
  pd <- cta_plot_data(tree, target_class = 1L)
  expect_setequal(pd$endpoints$stage, 1:3)
})

test_that("myeloma plot: endpoint_label contains percentage string", {
  skip_if_slow_tests_disabled("cta-plot-myeloma")
  if (!.plot_myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .plot_myeloma_fit()
  if (is.null(tree)) skip("myeloma fit returned NULL")
  pd    <- cta_plot_data(tree, target_class = 1L)
  nd    <- pd$nodes
  ep_lbls <- nd$endpoint_label[nd$leaf & !is.na(nd$endpoint_label)]
  expect_true(all(grepl("%", ep_lbls)))
})

test_that("myeloma plot: all endpoint_fill_color values are distinct hex strings", {
  skip_if_slow_tests_disabled("cta-plot-myeloma")
  if (!.plot_myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .plot_myeloma_fit()
  if (is.null(tree)) skip("myeloma fit returned NULL")
  pd     <- cta_plot_data(tree, target_class = 1L)
  ep     <- pd$endpoints
  colors <- ep$endpoint_fill_color
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors)))
  expect_equal(length(unique(colors)), nrow(ep))
})

test_that("myeloma plot: plot.cta_tree with target_class=1 does not error", {
  skip_if_slow_tests_disabled("cta-plot-myeloma")
  if (!.plot_myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .plot_myeloma_fit()
  if (is.null(tree)) skip("myeloma fit returned NULL")
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp, width = 800, height = 600)
  expect_no_error(
    plot(tree, target_class = 1L,
         class_labels = c("0" = "Alive", "1" = "Deceased"))
  )
  grDevices::dev.off()
})

test_that("myeloma plot: custom palette changes endpoint colors", {
  skip_if_slow_tests_disabled("cta-plot-myeloma")
  if (!.plot_myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .plot_myeloma_fit()
  if (is.null(tree)) skip("myeloma fit returned NULL")
  pd_def    <- cta_plot_data(tree, target_class = 1L)
  pd_custom <- cta_plot_data(tree, target_class = 1L,
                              endpoint_palette = c("#ffffff", "#c62828"))
  def_cols    <- sort(pd_def$endpoints$endpoint_fill_color)
  custom_cols <- sort(pd_custom$endpoints$endpoint_fill_color)
  expect_false(all(def_cols == custom_cols))
})

# =============================================================================
# plot.cta_tree() v2 — new parameter smoke tests (border_col, text_col,
# arrow_col, show_caption)
# =============================================================================

test_that("plot.cta_tree accepts border_col and text_col without error", {
  tree <- make_tree()
  tmp  <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp, width = 600, height = 400)
  expect_no_error(
    plot(tree, border_col = "#333333", text_col = "#111111")
  )
  grDevices::dev.off()
})

test_that("plot.cta_tree accepts arrow_col without error", {
  tree <- make_tree()
  tmp  <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp, width = 600, height = 400)
  expect_no_error(
    plot(tree, arrow_col = "#c62828")
  )
  grDevices::dev.off()
})

test_that("show_caption = TRUE with target_class does not error", {
  tree <- make_tree()
  tmp  <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp, width = 700, height = 500)
  expect_no_error(
    plot(tree, target_class = 1L, show_caption = TRUE)
  )
  grDevices::dev.off()
})

test_that("show_caption = TRUE without target_class does not error", {
  tree <- make_tree()
  tmp  <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp, width = 600, height = 400)
  expect_no_error(
    plot(tree, show_caption = TRUE)
  )
  grDevices::dev.off()
})
