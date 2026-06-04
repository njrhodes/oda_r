# tests/testthat/test-fit-evidence.R
# Slice S - Fit Object Evidence Wiring Audit tests
#
# Covers:
#  E1. novo_boot_ci S3 dispatch
#  E2. cta_plot_data() evidence fields
#  E3. ort_plot_data() evidence fields
#  E4. plot_cta_tree(show_metrics=TRUE) subtitle
#  E5. plot_lort_tree(show_metrics=TRUE) subtitle
#  E6. No-fabrication / no-refitting invariants
#  E7. Reserved-name checks

# ---------------------------------------------------------------------------
# Shared small fixtures
# ---------------------------------------------------------------------------

local({
  # Binary ODA fit (no weights, simple ordered attribute)
  set.seed(1L)
  x_bin <<- c(rep(1, 30), rep(3, 30))
  y_bin <<- c(rep(0L, 20), rep(1L, 10), rep(0L, 5), rep(1L, 25))

  # Multiclass ODA fit (3 classes)
  x_mc <<- c(1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12)
  y_mc <<- c(0L,0L,0L,0L,1L,1L,1L,1L,2L,2L,2L,2L)
})

# ---------------------------------------------------------------------------
# E1: novo_boot_ci S3 dispatch
# ---------------------------------------------------------------------------

test_that("E1-1: novo_boot_ci.default still accepts 2x2 matrix", {
  m <- matrix(c(146L, 40L, 36L, 33L), nrow = 2L, byrow = TRUE)
  ci <- novo_boot_ci(m, nboot = 50L, seed = 42L)
  expect_s3_class(ci, "novo_boot_ci")
  expect_equal(ci$n, 255L)
})

test_that("E1-2: novo_boot_ci.oda_fit dispatches on binary oda_fit", {
  fit <- oda_fit(x_bin, y_bin, attr_type = "ordered",
                 mcarlo = FALSE, loo = "off")
  expect_true(fit$ok)
  ci <- novo_boot_ci(fit, nboot = 50L, seed = 42L)
  expect_s3_class(ci, "novo_boot_ci")
  expect_equal(ci$n, sum(x_bin != 0L | !is.na(x_bin)))  # all obs classified
  expect_equal(ci$n, 60L)
})

test_that("E1-3: novo_boot_ci.oda_fit builds correct 2x2 from TP/TN/FP/FN", {
  fit <- oda_fit(x_bin, y_bin, attr_type = "ordered",
                 mcarlo = FALSE, loo = "off")
  conf <- fit$confusion
  expected_m <- matrix(c(conf$TN, conf$FP, conf$FN, conf$TP),
                       nrow = 2L, byrow = TRUE)
  ci <- novo_boot_ci(fit, nboot = 50L, seed = 1L)
  expect_equal(ci$confusion, expected_m)
})

test_that("E1-4: novo_boot_ci.oda_fit errors on multiclass", {
  fit <- oda_fit(x_mc, y_mc, attr_type = "ordered",
                 mcarlo = FALSE, loo = "off")
  expect_s3_class(fit, "oda_fit_multiclass")
  expect_error(novo_boot_ci(fit, nboot = 50L), regexp = "binary ODA fits")
})

test_that("E1-5: novo_boot_ci.cta_tree dispatches on cta_tree", {
  skip_if(Sys.getenv("ODA_TEST_TIER", unset = Sys.getenv("ODACORE_TEST_TIER")) %in% c("", "cran"),
          "Skipping cta_tree novo_boot_ci dispatch at cran tier")
  X <- data.frame(
    A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
    B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
  )
  y <- c(rep(0L, 40), rep(1L, 20))
  tree <- cta_fit(X, y, mc_iter = 500L, mc_seed = 1L, loo = "off",
                  mindenom = 5L)
  if (isTRUE(tree$no_tree)) skip("no_tree - skip dispatch test")
  ci <- novo_boot_ci(tree, nboot = 50L, seed = 42L)
  expect_s3_class(ci, "novo_boot_ci")
  expect_equal(ci$confusion, tree$training_confusion)
})

test_that("E1-6: novo_boot_ci.cta_tree errors on no_tree", {
  # Manufacture a minimal no_tree object
  no_tree_obj <- structure(
    list(no_tree = TRUE, training_confusion = NULL, n = 10L,
         C = 2L, has_weights = FALSE, overall_ess = NA_real_,
         nodes = list(), root_id = 1L),
    class = "cta_tree"
  )
  expect_error(novo_boot_ci(no_tree_obj, nboot = 50L), regexp = "no_tree")
})

test_that("E1-7: novo_boot_ci.cta_ort dispatches on cta_ort", {
  skip_if(Sys.getenv("ODA_TEST_TIER", unset = Sys.getenv("ODACORE_TEST_TIER")) %in% c("", "cran"),
          "Skipping cta_ort novo_boot_ci dispatch at cran tier")
  X <- data.frame(
    A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
    B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
  )
  y <- c(rep(0L, 40), rep(1L, 20))
  lort <- lort_fit(X, y, mc_iter = 500L, mc_seed = 1L, loo = "off",
                   min_n = 5L)
  ci <- novo_boot_ci(lort, nboot = 50L, seed = 42L)
  expect_s3_class(ci, "novo_boot_ci")
  # total n should match classified obs in strata
  st <- lort$strata
  total_from_strata <- sum(vapply(st$class_counts, sum, integer(1L)))
  expect_equal(ci$n, total_from_strata)
})

test_that("E1-8: novo_boot_ci.cta_ort accumulates confusion correctly", {
  skip_if(Sys.getenv("ODA_TEST_TIER", unset = Sys.getenv("ODACORE_TEST_TIER")) %in% c("", "cran"),
          "Skipping LORT confusion accumulation at cran tier")
  # Construct a minimal strata data frame manually
  # Two strata: class-0 terminal (pred=0), class-1 terminal (pred=1)
  strata_df <- data.frame(
    stratum_id     = c(1L, 2L),
    node_id        = c(2L, 3L),
    terminal_class = c(0L, 1L),
    prop_class1    = c(0.1, 0.8),
    n              = c(20L, 20L),
    stringsAsFactors = FALSE
  )
  strata_df$class_counts <- list(c("0" = 18L, "1" = 2L),
                                  c("0" = 4L,  "1" = 16L))

  fake_ort <- structure(
    list(
      no_tree      = FALSE,
      strata       = strata_df,
      ort_nodes    = list(),
      n_strata     = 2L,
      overall_ess  = NA_real_,
      has_weights  = FALSE,
      n            = 40L,
      C            = 2L,
      ort_settings = list(method = "lort", global_optimization = FALSE,
                          sda_anchored = FALSE)
    ),
    class = c("cta_ort", "cta_tree")
  )
  ci <- novo_boot_ci(fake_ort, nboot = 50L, seed = 1L)
  # Expected confusion: actual x predicted
  # stratum 1 (pred=0): TN=18, FN=2
  # stratum 2 (pred=1): FP=4,  TP=16
  # conf[1,1]=TN=18, conf[1,2]=FP=4, conf[2,1]=FN=2, conf[2,2]=TP=16
  expected <- matrix(c(18L, 4L, 2L, 16L), nrow = 2L, byrow = TRUE)
  expect_equal(ci$confusion, expected)
  expect_equal(ci$n, 40L)
})

test_that("E1-9: novo_boot_ci.cta_ort errors when strata absent", {
  fake_ort <- structure(
    list(no_tree = FALSE, strata = data.frame(),
         ort_nodes = list(), n_strata = 0L,
         overall_ess = NA_real_, has_weights = FALSE, n = 0L, C = 2L),
    class = c("cta_ort", "cta_tree")
  )
  expect_error(novo_boot_ci(fake_ort), regexp = "strata absent or empty")
})

# ---------------------------------------------------------------------------
# E2: cta_plot_data() evidence fields
# ---------------------------------------------------------------------------

test_that("E2-1: cta_plot_data() returns evidence fields for a valid tree", {
  skip_if(Sys.getenv("ODA_TEST_TIER", unset = Sys.getenv("ODACORE_TEST_TIER")) %in% c("", "cran"),
          "Skipping cta_plot_data evidence fields at cran tier")
  X <- data.frame(
    A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20))
  )
  y <- c(rep(0L, 40), rep(1L, 20))
  tree <- cta_fit(X, y, mc_iter = 500L, mc_seed = 1L, loo = "off",
                  mindenom = 5L)
  if (isTRUE(tree$no_tree)) skip("no_tree")
  pd <- cta_plot_data(tree)
  expect_true("overall_ess" %in% names(pd))
  expect_true("ess_label"   %in% names(pd))
  expect_true("d"           %in% names(pd))
  expect_true("model_label" %in% names(pd))
  expect_true("training_n"  %in% names(pd))
  expect_equal(pd$model_label, "CTA")
  expect_true(is.numeric(pd$overall_ess))
  expect_true(is.character(pd$ess_label))
  expect_equal(pd$ess_label, if (isTRUE(tree$has_weights)) "WESS" else "ESS")
  expect_equal(pd$training_n, tree$n)
})

test_that("E2-2: cta_plot_data() returns NA evidence for no_tree", {
  no_tree_obj <- structure(
    list(no_tree = TRUE, nodes = list(), root_id = 1L,
         training_confusion = NULL, overall_ess = NA_real_,
         has_weights = FALSE, n = 30L, C = 2L,
         attr_names = "A", miss_codes = NULL, priors_on = FALSE,
         alpha_split = 0.05, mindenom = 5L, prune_alpha = 0.05,
         max_depth = 8L, loo = "off"),
    class = "cta_tree"
  )
  pd <- cta_plot_data(no_tree_obj)
  expect_true(pd$no_tree)
  expect_true(is.na(pd$overall_ess))
  expect_true(is.na(pd$ess_label))
  expect_true(is.na(pd$d))
  expect_equal(pd$model_label, "CTA")
})

test_that("E2-3: cta_plot_data() d matches cta_d_stat()", {
  skip_if(Sys.getenv("ODA_TEST_TIER", unset = Sys.getenv("ODACORE_TEST_TIER")) %in% c("", "cran"),
          "Skipping cta_plot_data d-stat check at cran tier")
  X <- data.frame(A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)))
  y <- c(rep(0L, 40), rep(1L, 20))
  tree <- cta_fit(X, y, mc_iter = 500L, mc_seed = 1L, loo = "off",
                  mindenom = 5L)
  if (isTRUE(tree$no_tree)) skip("no_tree")
  pd <- cta_plot_data(tree)
  expect_equal(pd$d, cta_d_stat(tree))
})

# ---------------------------------------------------------------------------
# E3: ort_plot_data() evidence fields
# ---------------------------------------------------------------------------

test_that("E3-1: ort_plot_data() returns evidence fields", {
  skip_if(Sys.getenv("ODA_TEST_TIER", unset = Sys.getenv("ODACORE_TEST_TIER")) %in% c("", "cran"),
          "Skipping ort_plot_data evidence fields at cran tier")
  X <- data.frame(
    A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
    B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
  )
  y <- c(rep(0L, 40), rep(1L, 20))
  lort <- lort_fit(X, y, mc_iter = 500L, mc_seed = 1L, loo = "off",
                   min_n = 5L)
  pd <- ort_plot_data(lort)
  expect_true("overall_ess" %in% names(pd))
  expect_true("ess_label"   %in% names(pd))
  expect_true("d"           %in% names(pd))
  expect_true("model_label" %in% names(pd))
  expect_true("training_n"  %in% names(pd))
  expect_equal(pd$model_label, "LORT")
  expect_true(is.numeric(pd$overall_ess))
})

test_that("E3-2: ort_plot_data() returns NA evidence for empty ort_nodes", {
  fake_ort <- structure(
    list(ort_nodes = list(), strata = NULL,
         overall_ess = NA_real_, has_weights = FALSE, n = 0L, C = 2L,
         no_tree = FALSE, nodes = list(), root_id = 1L,
         training_confusion = NULL, attr_names = character(0),
         miss_codes = NULL, priors_on = FALSE, alpha_split = 0.05,
         mindenom = 5L, prune_alpha = 0.05, max_depth = 8L, loo = "off",
         ort_settings = list(method = "lort", global_optimization = FALSE,
                             sda_anchored = FALSE)),
    class = c("cta_ort", "cta_tree")
  )
  pd <- ort_plot_data(fake_ort)
  expect_equal(pd$model_label, "LORT")
  expect_true(is.na(pd$overall_ess))
})

# ---------------------------------------------------------------------------
# E4: plot_cta_tree(show_metrics = TRUE) subtitle
# ---------------------------------------------------------------------------

test_that("E4-1: plot_cta_tree(show_metrics=TRUE) adds metrics to subtitle", {
  skip_if(!requireNamespace("ggplot2", quietly = TRUE),
          "ggplot2 not available")
  skip_if(Sys.getenv("ODA_TEST_TIER", unset = Sys.getenv("ODACORE_TEST_TIER")) %in% c("", "cran"),
          "Skipping plot_cta_tree show_metrics at cran tier")
  X <- data.frame(A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)))
  y <- c(rep(0L, 40), rep(1L, 20))
  tree <- cta_fit(X, y, mc_iter = 500L, mc_seed = 1L, loo = "off",
                  mindenom = 5L)
  if (isTRUE(tree$no_tree)) skip("no_tree")
  p <- plot_cta_tree(tree, show_metrics = TRUE, color_by = "none")
  expect_s3_class(p, "ggplot")
  # Subtitle should contain ESS
  gg_labels <- p$labels
  st <- gg_labels$subtitle
  expect_false(is.null(st))
  expect_true(grepl("ESS|WESS", st))
})

test_that("E4-2: plot_cta_tree(show_metrics=FALSE) leaves subtitle unchanged", {
  skip_if(!requireNamespace("ggplot2", quietly = TRUE),
          "ggplot2 not available")
  skip_if(Sys.getenv("ODA_TEST_TIER", unset = Sys.getenv("ODACORE_TEST_TIER")) %in% c("", "cran"),
          "Skipping plot_cta_tree at cran tier")
  X <- data.frame(A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)))
  y <- c(rep(0L, 40), rep(1L, 20))
  tree <- cta_fit(X, y, mc_iter = 500L, mc_seed = 1L, loo = "off",
                  mindenom = 5L)
  if (isTRUE(tree$no_tree)) skip("no_tree")
  p <- plot_cta_tree(tree, subtitle = "My note", show_metrics = FALSE,
                     color_by = "none")
  expect_equal(p$labels$subtitle, "My note")
})

test_that("E4-3: plot_cta_tree(show_metrics=TRUE) appends to existing subtitle", {
  skip_if(!requireNamespace("ggplot2", quietly = TRUE),
          "ggplot2 not available")
  skip_if(Sys.getenv("ODA_TEST_TIER", unset = Sys.getenv("ODACORE_TEST_TIER")) %in% c("", "cran"),
          "Skipping plot_cta_tree show_metrics append at cran tier")
  X <- data.frame(A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)))
  y <- c(rep(0L, 40), rep(1L, 20))
  tree <- cta_fit(X, y, mc_iter = 500L, mc_seed = 1L, loo = "off",
                  mindenom = 5L)
  if (isTRUE(tree$no_tree)) skip("no_tree")
  p <- plot_cta_tree(tree, subtitle = "My note", show_metrics = TRUE,
                     color_by = "none")
  st <- p$labels$subtitle
  expect_true(grepl("My note", st))
  expect_true(grepl("ESS|WESS", st))
})

# ---------------------------------------------------------------------------
# E5: plot_lort_tree(show_metrics = TRUE)
# ---------------------------------------------------------------------------

test_that("E5-1: plot_lort_tree(show_metrics=TRUE) adds metrics to subtitle", {
  skip_if(!requireNamespace("ggplot2", quietly = TRUE),
          "ggplot2 not available")
  skip_if(Sys.getenv("ODA_TEST_TIER", unset = Sys.getenv("ODACORE_TEST_TIER")) %in% c("", "cran"),
          "Skipping plot_lort_tree show_metrics at cran tier")
  X <- data.frame(
    A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
    B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
  )
  y <- c(rep(0L, 40), rep(1L, 20))
  lort <- lort_fit(X, y, mc_iter = 500L, mc_seed = 1L, loo = "off",
                   min_n = 5L)
  p <- plot_lort_tree(lort, show_metrics = TRUE, color_by = "none")
  expect_s3_class(p, "ggplot")
  st <- p$labels$subtitle
  # If LORT overall_ess is available, subtitle should contain ESS/WESS
  if (!is.na(lort$overall_ess)) {
    expect_true(grepl("ESS|WESS", st))
  }
})

# ---------------------------------------------------------------------------
# E6: No-fabrication / no-refitting invariants
# ---------------------------------------------------------------------------

test_that("E6-1: novo_boot_ci.oda_fit does not call oda_fit or cta_fit", {
  # Formal argument inspection: none of the fitting functions should appear
  # in the body of novo_boot_ci.oda_fit
  body_txt <- paste(deparse(body(novo_boot_ci.oda_fit)), collapse = " ")
  expect_false(grepl("oda_fit\\(|cta_fit\\(|lort_fit\\(", body_txt))
})

test_that("E6-2: novo_boot_ci.cta_tree does not call fitting functions", {
  body_txt <- paste(deparse(body(novo_boot_ci.cta_tree)), collapse = " ")
  expect_false(grepl("oda_fit\\(|cta_fit\\(|lort_fit\\(|oda_cta_fit\\(", body_txt))
})

test_that("E6-3: novo_boot_ci.cta_ort does not call fitting functions", {
  body_txt <- paste(deparse(body(novo_boot_ci.cta_ort)), collapse = " ")
  expect_false(grepl("oda_fit\\(|cta_fit\\(|lort_fit\\(|oda_cta_fit\\(", body_txt))
})

test_that("E6-4: cta_plot_data does not call fitting functions", {
  body_txt <- paste(deparse(body(cta_plot_data)), collapse = " ")
  expect_false(grepl("oda_fit\\(|cta_fit\\(|lort_fit\\(|oda_cta_fit\\(", body_txt))
})

test_that("E6-5: ort_plot_data does not call fitting functions", {
  body_txt <- paste(deparse(body(ort_plot_data)), collapse = " ")
  expect_false(grepl("oda_fit\\(|cta_fit\\(|lort_fit\\(|oda_cta_fit\\(", body_txt))
})

# ---------------------------------------------------------------------------
# E7: Reserved name checks
# ---------------------------------------------------------------------------

test_that("E7-1: sda_propensity_weights does not exist", {
  expect_false(existsFunction("sda_propensity_weights"))
})

test_that("E7-2: sort_propensity_weights does not exist", {
  expect_false(existsFunction("sort_propensity_weights"))
})

test_that("E7-3: gort_propensity_weights does not exist", {
  expect_false(existsFunction("gort_propensity_weights"))
})

test_that("E7-4: model_evidence generic does not exist", {
  expect_false(existsFunction("model_evidence"))
})

# ---------------------------------------------------------------------------
# E8: source_type / source_id / weighted provenance metadata
# ---------------------------------------------------------------------------

test_that("E8-1: default path carries source_type = 'matrix'", {
  m  <- matrix(c(146L, 40L, 36L, 33L), nrow = 2L, byrow = TRUE)
  ci <- novo_boot_ci(m, nboot = 50L, seed = 42L)
  expect_equal(ci$source_type, "matrix")
  expect_true(is.na(ci$source_id))
  expect_true(is.na(ci$weighted))
})

test_that("E8-2: oda_fit path carries source_type = 'oda_fit', weighted = FALSE", {
  fit <- oda_fit(x_bin, y_bin, attr_type = "ordered",
                 mcarlo = FALSE, loo = "off")
  ci <- novo_boot_ci(fit, nboot = 50L, seed = 42L)
  expect_equal(ci$source_type, "oda_fit")
  expect_true(is.na(ci$source_id))
  expect_false(ci$weighted)
})

test_that("E8-3: cta_tree full-tree path carries source_type = 'cta_tree'", {
  skip_if(Sys.getenv("ODA_TEST_TIER", unset = Sys.getenv("ODACORE_TEST_TIER")) %in% c("", "cran"),
          "Skipping cta_tree source_type at cran tier")
  X <- data.frame(
    A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
    B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
  )
  y <- c(rep(0L, 40), rep(1L, 20))
  tree <- cta_fit(X, y, mc_iter = 500L, mc_seed = 1L, loo = "off",
                  mindenom = 5L)
  if (isTRUE(tree$no_tree)) skip("no_tree")
  ci <- novo_boot_ci(tree, nboot = 50L, seed = 42L)
  expect_equal(ci$source_type, "cta_tree")
  expect_true(is.na(ci$source_id))
  expect_false(ci$weighted)
})

test_that("E8-4: cta_tree node_id path uses leaf class_counts and carries metadata", {
  skip_if(Sys.getenv("ODA_TEST_TIER", unset = Sys.getenv("ODACORE_TEST_TIER")) %in% c("", "cran"),
          "Skipping cta_tree node_id at cran tier")
  X <- data.frame(
    A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
    B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
  )
  y <- c(rep(0L, 40), rep(1L, 20))
  tree <- cta_fit(X, y, mc_iter = 500L, mc_seed = 1L, loo = "off",
                  mindenom = 5L)
  expect_false(isTRUE(tree$no_tree), label = "mc_seed=1 fixture must produce a tree")
  # Find a terminal node id - nodes is unnamed; iterate by position
  leaf_positions <- which(vapply(tree$nodes,
    function(nd) !is.null(nd) && isTRUE(nd$leaf), logical(1L)))
  expect_true(length(leaf_positions) > 0L, label = "mc_seed=1 fixture must have leaf nodes")
  nid <- tree$nodes[[leaf_positions[[1L]]]]$node_id
  ci <- novo_boot_ci(tree, node_id = nid, nboot = 50L, seed = 42L)
  expect_s3_class(ci, "novo_boot_ci")
  expect_equal(ci$source_type, "cta_tree_node")
  expect_equal(ci$source_id, nid)
  expect_false(ci$weighted)
  # Confusion column sums should be 0 on the non-predicted side
  nd       <- tree$nodes[[leaf_positions[[1L]]]]
  mc       <- nd$majority_class  # 0 or 1
  pred_col <- mc + 1L
  other_col <- 3L - pred_col
  expect_equal(sum(ci$confusion[, other_col]), 0L)
})

test_that("E8-5: cta_tree node_id errors on non-leaf node", {
  skip_if(Sys.getenv("ODA_TEST_TIER", unset = Sys.getenv("ODACORE_TEST_TIER")) %in% c("", "cran"),
          "Skipping at cran tier")
  X <- data.frame(
    A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
    B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
  )
  y <- c(rep(0L, 40), rep(1L, 20))
  tree <- cta_fit(X, y, mc_iter = 500L, mc_seed = 1L, loo = "off",
                  mindenom = 5L)
  expect_false(isTRUE(tree$no_tree), label = "mc_seed=1 fixture must produce a tree")
  # Root is the canonical non-leaf node - use root_id directly
  nid <- tree$root_id
  expect_false(isTRUE(tree$nodes[[nid]]$leaf), label = "root must be a split node")
  expect_error(novo_boot_ci(tree, node_id = nid, nboot = 50L),
               regexp = "not a terminal")
})

test_that("E8-6: cta_ort full-LORT path carries source_type = 'cta_ort'", {
  strata_df <- data.frame(
    stratum_id     = c(1L, 2L),
    node_id        = c(2L, 3L),
    terminal_class = c(0L, 1L),
    prop_class1    = c(0.1, 0.8),
    n              = c(20L, 20L),
    stringsAsFactors = FALSE
  )
  strata_df$class_counts <- list(c("0" = 18L, "1" = 2L),
                                  c("0" = 4L,  "1" = 16L))
  fake_ort <- structure(
    list(no_tree = FALSE, strata = strata_df, ort_nodes = list(),
         n_strata = 2L, overall_ess = NA_real_, has_weights = FALSE,
         n = 40L, C = 2L,
         ort_settings = list(method = "lort", global_optimization = FALSE,
                             sda_anchored = FALSE)),
    class = c("cta_ort", "cta_tree")
  )
  ci <- novo_boot_ci(fake_ort, nboot = 50L, seed = 1L)
  expect_equal(ci$source_type, "cta_ort")
  expect_true(is.na(ci$source_id))
  expect_false(ci$weighted)
})

test_that("E8-7: cta_ort stratum_id path uses single-stratum counts and carries metadata", {
  strata_df <- data.frame(
    stratum_id     = c(1L, 2L),
    node_id        = c(2L, 3L),
    terminal_class = c(0L, 1L),
    prop_class1    = c(0.1, 0.8),
    n              = c(20L, 20L),
    stringsAsFactors = FALSE
  )
  strata_df$class_counts <- list(c("0" = 18L, "1" = 2L),
                                  c("0" = 4L,  "1" = 16L))
  fake_ort <- structure(
    list(no_tree = FALSE, strata = strata_df, ort_nodes = list(),
         n_strata = 2L, overall_ess = NA_real_, has_weights = FALSE,
         n = 40L, C = 2L,
         ort_settings = list(method = "lort", global_optimization = FALSE,
                             sda_anchored = FALSE)),
    class = c("cta_ort", "cta_tree")
  )
  # Stratum 2: terminal_class=1, class_counts=c("0"=4,"1"=16)
  # Expected confusion: FP=4 (act=0,pred=1), TP=16 (act=1,pred=1), TN=0, FN=0
  ci <- novo_boot_ci(fake_ort, stratum_id = 2L, nboot = 50L, seed = 1L)
  expect_s3_class(ci, "novo_boot_ci")
  expect_equal(ci$source_type, "cta_ort_stratum")
  expect_equal(ci$source_id, 2L)
  expect_false(ci$weighted)
  expected <- matrix(c(0L, 4L, 0L, 16L), nrow = 2L, byrow = TRUE)
  expect_equal(ci$confusion, expected)
  expect_equal(ci$n, 20L)
})

test_that("E8-8: cta_ort stratum_id errors on missing stratum", {
  strata_df <- data.frame(
    stratum_id     = c(1L, 2L),
    terminal_class = c(0L, 1L),
    n              = c(20L, 20L),
    stringsAsFactors = FALSE
  )
  strata_df$class_counts <- list(c("0" = 18L, "1" = 2L),
                                  c("0" = 4L,  "1" = 16L))
  fake_ort <- structure(
    list(no_tree = FALSE, strata = strata_df, n_strata = 2L,
         overall_ess = NA_real_, has_weights = FALSE, n = 40L, C = 2L),
    class = c("cta_ort", "cta_tree")
  )
  expect_error(novo_boot_ci(fake_ort, stratum_id = 99L, nboot = 50L),
               regexp = "not found in strata")
})

# ---------------------------------------------------------------------------
# E9: as_confusion_matrix() - bridge from cta_confusion_table to 2x2 matrix
# ---------------------------------------------------------------------------

test_that("E9-1: as_confusion_matrix converts 2x2 tidy df to matrix", {
  df <- data.frame(
    actual    = c(0L, 0L, 1L, 1L),
    predicted = c(0L, 1L, 0L, 1L),
    n         = c(146L, 40L, 36L, 33L)
  )
  m <- as_confusion_matrix(df)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(2L, 2L))
  expect_equal(m[1L, 1L], 146L)   # TN
  expect_equal(m[1L, 2L], 40L)    # FP
  expect_equal(m[2L, 1L], 36L)    # FN
  expect_equal(m[2L, 2L], 33L)    # TP
})

test_that("E9-2: as_confusion_matrix round-trips cta_confusion_table", {
  skip_if(Sys.getenv("ODA_TEST_TIER", unset = Sys.getenv("ODACORE_TEST_TIER")) %in% c("", "cran"),
          "Skipping cta_confusion_table round-trip at cran tier")
  X <- data.frame(
    A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
    B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
  )
  y <- c(rep(0L, 40), rep(1L, 20))
  tree <- cta_fit(X, y, mc_iter = 500L, mc_seed = 1L, loo = "off",
                  mindenom = 5L)
  if (isTRUE(tree$no_tree)) skip("no_tree")
  ct <- cta_confusion_table(tree)
  m  <- as_confusion_matrix(ct)
  expect_equal(m, tree$training_confusion)
})

test_that("E9-3: as_confusion_matrix result passes novo_boot_ci.default", {
  df <- data.frame(
    actual    = c(0L, 0L, 1L, 1L),
    predicted = c(0L, 1L, 0L, 1L),
    n         = c(146L, 40L, 36L, 33L)
  )
  m  <- as_confusion_matrix(df)
  ci <- novo_boot_ci(m, nboot = 50L, seed = 1L)
  expect_s3_class(ci, "novo_boot_ci")
  expect_equal(ci$n, 255L)
})

test_that("E9-4: as_confusion_matrix errors on missing columns", {
  df <- data.frame(x = 1L, y = 2L)
  expect_error(as_confusion_matrix(df), regexp = "columns 'actual', 'predicted', 'n'")
})
