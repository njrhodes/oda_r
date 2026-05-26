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
#   A = 0 (n=20): all y = 0  → no discrimination; C=1 → no_tree at level 1
#   A = 1, B = 0 (n=20): y = 0
#   A = 1, B = 1 (n=20): y = 1
#
# Expected structure:
#   Root ORT node 1: MDSA finds split on A
#   Right child ORT node 2 (A > 0.5, n=40): MDSA finds split on B
#     Right grandchild ORT node 3 (A>0.5 AND B>0.5, n=20): no_tree (pure y=1)
#     Left  grandchild ORT node 4 (A>0.5 AND B<=0.5, n=20): no_tree (pure y=0)
#   Left  child ORT node 5 (A <= 0.5, n=20): no_tree (pure y=0)
#
# 3 terminal strata (nodes 3, 4, 5).
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
  ort <- do.call(cta_fit, syn_ort_args)
  expect_s3_class(ort, "cta_ort")
  expect_s3_class(ort, "cta_tree")
  expect_true(isTRUE(ort$recursive))
})

# ---------------------------------------------------------------------------
# Test 4: two-level synthetic tree — 3 terminal strata
# ---------------------------------------------------------------------------
test_that("cta_fit recursive: two-level synthetic tree has 3 strata", {
  ort <- do.call(cta_fit, syn_ort_args)
  expect_equal(ort$n_strata, 3L)
  expect_equal(nrow(ort$strata), 3L)
  # Strata sorted ascending by prop_class1
  expect_true(all(diff(ort$strata$prop_class1) >= 0))
})

# ---------------------------------------------------------------------------
# Test 5: strata n sums to total N
# ---------------------------------------------------------------------------
test_that("cta_fit recursive: strata n values sum to N", {
  ort <- do.call(cta_fit, syn_ort_args)
  expect_equal(sum(ort$strata$n), nrow(syn_X))
})

# ---------------------------------------------------------------------------
# Test 6: strata_check_passed is TRUE
# ---------------------------------------------------------------------------
test_that("cta_fit recursive: strata_check_passed is TRUE", {
  ort <- do.call(cta_fit, syn_ort_args)
  expect_true(isTRUE(ort$strata_check_passed))
})

# ---------------------------------------------------------------------------
# Test 7: predict type = "class" returns integer vector of length N
# ---------------------------------------------------------------------------
test_that("cta_fit recursive: predict type='class' has correct length and type", {
  ort  <- do.call(cta_fit, syn_ort_args)
  pred <- predict(ort, syn_X, type = "class")
  expect_true(is.integer(pred))
  expect_equal(length(pred), nrow(syn_X))
})

# ---------------------------------------------------------------------------
# Test 8: predict type = "all" returns data.frame with correct columns
# ---------------------------------------------------------------------------
test_that("cta_fit recursive: predict type='all' returns correct data.frame", {
  ort <- do.call(cta_fit, syn_ort_args)
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
  ort <- do.call(cta_fit, syn_ort_args)
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
  ort <- do.call(cta_fit, syn_ort_args)
  expect_silent(capture.output(print(ort)))
  sm  <- summary(ort)
  expect_s3_class(sm, "cta_ort_summary")
  expect_silent(capture.output(print(sm)))
})

# ---------------------------------------------------------------------------
# Test 12: ort_plot_data returns valid structure
# ---------------------------------------------------------------------------
test_that("ort_plot_data: returns list with nodes, edges, strata", {
  ort <- do.call(cta_fit, syn_ort_args)
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
  ort <- do.call(cta_fit, syn_ort_args)
  pdf(nullfile())
  on.exit(dev.off())
  expect_silent(plot(ort))
  expect_silent(plot(ort, target_class = 1L))
})
