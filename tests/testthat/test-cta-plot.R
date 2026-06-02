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
  data(mtcars, envir = environment())
  X <- mtcars[, c("cyl", "disp", "hp", "wt")]
  y <- as.integer(mtcars$am)
  suppressMessages(
    oda_cta_fit(X, y, mindenom = 9999L, mc_iter = 50L, mc_seed = 42L,
                loo = "off")
  )
}

# =============================================================================
# 1. cta_plot_data() schema — list names and data.frame columns
# =============================================================================

test_that("cta_plot_data schema: correct top-level names and column sets", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  expect_type(pd, "list")
  expect_true(all(c("nodes", "edges", "no_tree", "has_weights",
                    "overall_ess", "ess_label", "d", "model_label",
                    "training_n") %in% names(pd)))
  expect_true(all(c("node_id", "parent_id", "depth", "x", "y", "leaf",
                    "attribute", "n_obs", "majority_class", "ess",
                    "label") %in% names(pd$nodes)))
  expect_true(all(c("from_node_id", "to_node_id", "x0", "y0",
                    "x1", "y1", "label") %in% names(pd$edges)))
})

# =============================================================================
# 2. cta_plot_data() — column types and data properties
# =============================================================================

test_that("cta_plot_data data properties: types, geometry, labels, flags", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree)
  nd   <- pd$nodes
  ed   <- pd$edges

  # scalar flags
  expect_false(pd$has_weights)
  expect_false(pd$no_tree)

  # node column types
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

  # edge column types
  expect_type(ed$from_node_id, "integer")
  expect_type(ed$to_node_id,   "integer")
  expect_type(ed$x0,           "double")
  expect_type(ed$y0,           "double")
  expect_type(ed$x1,           "double")
  expect_type(ed$y1,           "double")
  expect_type(ed$label,        "character")

  # size: root + 2 leaves minimum
  expect_gte(nrow(nd), 3L)
  expect_gte(nrow(ed), 2L)

  # leaf geometry
  leaves <- nd[nd$leaf, ]
  expect_true(all(is.na(leaves$ess)))
  expect_true(all(is.na(leaves$attribute)))
  expect_equal(nd$y, -nd$depth)
  expect_equal(length(leaves$x), length(unique(leaves$x)))

  # split nodes have non-NA attribute
  splits <- nd[!nd$leaf, ]
  if (nrow(splits) > 0L)
    expect_true(all(!is.na(splits$attribute)))

  # edge connectivity and labels
  all_ids <- nd$node_id
  expect_true(all(ed$from_node_id %in% all_ids))
  expect_true(all(ed$to_node_id   %in% all_ids))
  expect_true(all(nchar(ed$label) > 0L))

  # node label content
  expect_true(all(grepl("class=", leaves$label)))
  expect_true(all(grepl("n=",     leaves$label)))
  if (nrow(splits) > 0L)
    expect_true(all(grepl("ESS", splits$label)))
})

# =============================================================================
# 3. cta_plot_data() — no-tree and stump edge cases
# =============================================================================

test_that("cta_plot_data edge cases: no-tree is empty; stump has 3 nodes/2 edges", {
  # no-tree
  nt <- make_no_tree()
  pd <- cta_plot_data(nt)
  expect_true(pd$no_tree)
  expect_equal(nrow(pd$nodes), 0L)
  expect_equal(nrow(pd$edges), 0L)

  # stump (conditional on mindenom=14 producing a stump)
  data(mtcars, envir = environment())
  X <- mtcars[, c("cyl", "disp", "hp", "wt")]
  y <- as.integer(mtcars$am)
  st <- suppressMessages(
    oda_cta_fit(X, y, mindenom = 14L, mc_iter = 500L, mc_seed = 42L,
                loo = "off")
  )
  pd_st <- cta_plot_data(st)
  if (identical(.cta_tree_status(st), "stump")) {
    expect_equal(nrow(pd_st$nodes), 3L)
    expect_equal(nrow(pd_st$edges), 2L)
  } else {
    skip("tree is not a stump at mindenom=14 — skip stump-specific counts")
  }
})

# =============================================================================
# 4. plot.cta_tree() — smoke + parameter acceptance
# =============================================================================

test_that("plot.cta_tree: no errors across all parameter variants", {
  tree <- make_tree()
  nt   <- make_no_tree()

  .plot <- function(obj, ...) {
    tmp <- tempfile(fileext = ".png")
    on.exit(unlink(tmp), add = TRUE)
    grDevices::png(tmp, width = 700, height = 500)
    on.exit(grDevices::dev.off(), add = TRUE)
    plot(obj, ...)
  }

  # basic smoke
  expect_no_error(.plot(tree))
  expect_no_error(.plot(nt))

  # return value is a list with expected names
  result <- .plot(tree)
  expect_type(result, "list")
  expect_true(all(c("nodes", "edges", "no_tree", "has_weights") %in% names(result)))

  # parameter variants
  expect_no_error(.plot(tree, main = "Custom Title"))
  expect_no_error(.plot(tree, border_col = "#333333", text_col = "#111111"))
  expect_no_error(.plot(tree, arrow_col = "#c62828"))
  expect_no_error(.plot(tree, show_caption = TRUE))
  expect_no_error(.plot(tree, target_class = 1L, show_caption = TRUE))
  expect_no_error(.plot(tree, node_col_split = "#aabbcc", node_col_leaf = "#ddeeff"))

  # with target_class: return value carries target_class_used
  result2 <- .plot(tree, target_class = 1L)
  expect_equal(result2$target_class_used, 1L)
})

# =============================================================================
# 5. cta_plot_data() — target-class enrichment
# =============================================================================

test_that("cta_plot_data target-class enrichment: columns, types, contracts", {
  tree <- make_tree()
  pd   <- cta_plot_data(tree, target_class = 1L)
  nd   <- pd$nodes

  # enrichment columns present
  enrichment_cols <- c("endpoint_id", "stage", "target_class", "target_n",
                        "denominator", "target_proportion", "target_rank",
                        "endpoint_fill_color", "predicted_label",
                        "target_label", "endpoint_label")
  expect_true(all(enrichment_cols %in% names(nd)))

  # endpoints data.frame
  expect_true("endpoints" %in% names(pd))
  expect_s3_class(pd$endpoints, "data.frame")

  # scalar
  expect_equal(pd$target_class_used, 1L)

  # leaf fill_color: hex strings
  leaf_colors <- nd$endpoint_fill_color[nd$leaf]
  expect_type(leaf_colors, "character")
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", leaf_colors)))

  # target_rank: integer, each leaf gets unique rank 1..n_ep
  leaf_ranks <- nd$target_rank[nd$leaf]
  expect_type(leaf_ranks, "integer")
  expect_setequal(leaf_ranks, seq_len(sum(nd$leaf)))

  # target_proportion: double in [0, 1]
  props <- nd$target_proportion[nd$leaf]
  expect_type(props, "double")
  expect_true(all(props >= 0 & props <= 1))

  # no stale fill_key column
  expect_false("fill_key" %in% names(nd))

  # class_labels respected in predicted_label
  pd2 <- cta_plot_data(tree, target_class = 1L,
                        class_labels = c("0" = "Manual", "1" = "Auto"))
  pred_labels <- pd2$nodes$predicted_label[nd$leaf & !is.na(nd$predicted_label)]
  expect_true(all(pred_labels %in% c("Manual", "Auto")))

  # palette function accepted
  my_pal <- function(n) grDevices::colorRampPalette(c("#000000", "#ffffff"))(n)
  expect_no_error(cta_plot_data(tree, target_class = 1L, endpoint_palette = my_pal))

  # palette vector changes colors
  pd_def    <- cta_plot_data(tree, target_class = 1L)
  pd_custom <- cta_plot_data(tree, target_class = 1L,
                              endpoint_palette = c("#000000", "#ffffff"))
  nd_def  <- pd_def$nodes[pd_def$nodes$leaf, ]
  nd_cust <- pd_custom$nodes[pd_custom$nodes$leaf, ]
  expect_false(all(nd_def$endpoint_fill_color == nd_cust$endpoint_fill_color))

  # no-tree with target_class does not error
  expect_no_error(cta_plot_data(make_no_tree(), target_class = 1L))
})

# =============================================================================
# 6. Myeloma slow-tier smoke
# =============================================================================

# ---- helpers ----------------------------------------------------------------

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
          mc_seed     = NULL,
          loo         = "stable",
          attr_names  = .plot_myeloma_attrs
        )
      )
    }
    fit
  }
})

# ---- smoke block ------------------------------------------------------------

test_that("myeloma plot smoke: 3 endpoints, proportions, stages, labels, colors", {
  skip_if_slow_tests_disabled("cta-plot-myeloma")
  if (!.plot_myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .plot_myeloma_fit()
  if (is.null(tree)) skip("myeloma fit returned NULL")

  pd <- cta_plot_data(tree, target_class = 1L)
  ep <- pd$endpoints
  nd <- pd$nodes

  # 3 endpoints
  expect_equal(nrow(ep), 3L)

  # proportions match staging table
  props    <- sort(ep$target_proportion)
  st_props <- sort(cta_staging_table(tree, target_class = 1L)$target_proportion)
  expect_equal(props, st_props)

  # stages 1–3
  expect_setequal(ep$stage, 1:3)

  # endpoint_label contains "%" on leaf nodes
  ep_lbls <- nd$endpoint_label[nd$leaf & !is.na(nd$endpoint_label)]
  expect_true(all(grepl("%", ep_lbls)))

  # all endpoint colors are distinct hex strings
  colors <- ep$endpoint_fill_color
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors)))
  expect_equal(length(unique(colors)), nrow(ep))

  # plot with class_labels does not error
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp, width = 800, height = 600)
  expect_no_error(
    plot(tree, target_class = 1L,
         class_labels = c("0" = "Alive", "1" = "Deceased"))
  )
  grDevices::dev.off()

  # custom palette changes colors
  pd_def    <- cta_plot_data(tree, target_class = 1L)
  pd_custom <- cta_plot_data(tree, target_class = 1L,
                              endpoint_palette = c("#ffffff", "#c62828"))
  expect_false(all(sort(pd_def$endpoints$endpoint_fill_color) ==
                   sort(pd_custom$endpoints$endpoint_fill_color)))
})
