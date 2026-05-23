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
# cta_plot_data() — structural contract
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
# plot.cta_tree() — smoke tests (no errors, returns invisible)
# =============================================================================

test_that("plot.cta_tree does not error on a valid tree", {
  tree <- make_tree()
  tmp  <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp, width = 600, height = 400)
  expect_no_error(plot(tree))
  grDevices::dev.off()
})

test_that("plot.cta_tree returns invisible(x)", {
  tree <- make_tree()
  tmp  <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp, width = 600, height = 400)
  result <- plot(tree)
  grDevices::dev.off()
  expect_identical(result, tree)
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
