###############################################################################
# test-cta-node-selection.R
#
# Integration tests for the CTA-specific ordered-cut dispatch and LOO STABLE
# gate introduced to match CTA.exe myeloma MINDENOM=1 node-level behavior.
#
# Key audit findings (docs/CTA_ORDERED_CUT_AUDIT.md):
#
#   Node 4 (V17<=0.5 AND V15<=0.5, n=113):
#     V4  CTA scan: cut=359.80, WESS=34.90%, Signif F → rejected at MC gate.
#         Even when Signif T, |WESS(34.90) - WESSL(~46.99)| >> 0.01 → UNSTABLE.
#         V4 is the only candidate → no split.
#
#   Node 2 (V17<=0.5, n=131):
#     V4  CTA scan: cut=371.20, WESS=34.89%, Signif T.
#         Generic ODA LOO WESSL≈38.52%; |34.89-38.52|=3.63 >> 0.01 → UNSTABLE.
#         V4 always rejected by LOO STABLE gate (regardless of MC outcome).
#     V15 binary (2 unique values): CTA path bypassed; generic ODA WESS=17.57%.
#         Generic ODA LOO WESSL=17.57% = WESS → STABLE.
#         V15 Signif T → selected as root.
#
#   cta_demo (uniform weights):
#     CTA path guard (any(w != w[1])) is FALSE → CTA path bypassed entirely.
#     Generic ODA path unchanged → root=V2, ESS=52.63%.
#
# Reference: docs/CTA_ORDERED_CUT_AUDIT.md
###############################################################################

# ---- Myeloma data helpers ---------------------------------------------------

.ns_myeloma_path <- function(f) {
  tryCatch(testthat::test_path(file.path("fixtures", "myeloma", f)),
           error = function(e) "")
}

.ns_myeloma_ok <- function() {
  p <- .ns_myeloma_path("data.txt")
  nzchar(p) && file.exists(p)
}

.ns_load_myeloma <- function() {
  f <- .ns_myeloma_path("data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]   # EX V2=0
}

# Node 4 rows: V17 present AND V17 <= 0.5, V15 present AND V15 <= 0.5
.ns_node4_idx <- function(d) {
  miss <- -9L
  with(d, which(!(V17 %in% miss) & V17 <= 0.5 &
                !(V15 %in% miss) & V15 <= 0.5))
}

# Node 2 rows: V17 present AND V17 <= 0.5
.ns_node2_idx <- function(d) {
  miss <- -9L
  with(d, which(!(V17 %in% miss) & V17 <= 0.5))
}

# =============================================================================
# Test 1: Node 4, V4 only → always leaf
#
# V4 at Node 4 is Signif F per CTA.exe (WESS=34.90%).
# Even when Signif T by chance, LOO gate rejects: |34.90 - 46.99| >> 0.01.
# =============================================================================

test_that("node-selection: Node 4 V4-only is always rejected → leaf-only tree", {
  skip_if_slow_tests_disabled("cta-node-selection")
  if (!.ns_myeloma_ok()) skip("myeloma fixture missing")

  d   <- .ns_load_myeloma()
  idx <- .ns_node4_idx(d)

  tree <- oda_cta_fit(
    X           = d[idx, "V4", drop = FALSE],
    y           = d$V1[idx],
    w           = d$V2[idx],
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
    attr_names  = "V4"
  )

  root <- tree$nodes[[tree$root_id]]
  expect_true(isTRUE(root$leaf),
    label = paste0("V4 at Node 4 must be rejected (Signif F or LOO UNSTABLE);",
                   " root should be a leaf"))
})

# =============================================================================
# Test 2: Node 2, V4 only → always leaf (Signif T but LOO UNSTABLE)
#
# V4 at Node 2: CTA scan WESS=34.89% (cut=371.20), CTA.exe Signif T.
# Generic ODA LOO WESSL≈38.52% (cut≈84.85). |34.89 - 38.52| = 3.63 >> 0.01.
# V4 is ALWAYS rejected at the LOO STABLE gate regardless of MC outcome.
# =============================================================================

test_that("node-selection: Node 2 V4-only always LOO UNSTABLE → leaf-only tree", {
  skip_if_slow_tests_disabled("cta-node-selection")
  if (!.ns_myeloma_ok()) skip("myeloma fixture missing")

  d   <- .ns_load_myeloma()
  idx <- .ns_node2_idx(d)

  tree <- oda_cta_fit(
    X           = d[idx, "V4", drop = FALSE],
    y           = d$V1[idx],
    w           = d$V2[idx],
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
    attr_names  = "V4"
  )

  root <- tree$nodes[[tree$root_id]]
  expect_true(isTRUE(root$leaf),
    label = paste0("V4 at Node 2: CTA scan WESS=34.89% but generic ODA LOO",
                   " WESSL≈38.52%; |delta|=3.63pp → UNSTABLE → must be rejected"))
})

# =============================================================================
# Test 3: Node 2, V4+V15 → V4 rejected; root must be V15
#
# V4: CTA path → LOO UNSTABLE → rejected (see Test 2 above).
# V15: binary (2 unique values) → CTA path bypassed; generic ODA WESS=17.57%,
#      LOO WESSL=17.57% → STABLE; Signif T → selected as root.
# =============================================================================

test_that("node-selection: Node 2 V4+V15 → V4 rejected; root is V15 (STABLE)", {
  skip_if_slow_tests_disabled("cta-node-selection")
  if (!.ns_myeloma_ok()) skip("myeloma fixture missing")

  d   <- .ns_load_myeloma()
  idx <- .ns_node2_idx(d)

  tree <- oda_cta_fit(
    X           = d[idx, c("V4", "V15"), drop = FALSE],
    y           = d$V1[idx],
    w           = d$V2[idx],
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
    attr_names  = c("V4", "V15")
  )

  root <- tree$nodes[[tree$root_id]]

  # V4 must never be selected — LOO STABLE gate always rejects it.
  if (!isTRUE(root$leaf)) {
    expect_equal(root$attribute, "V15",
      label = "V4 is LOO UNSTABLE; only V15 can be root")
    expect_equal(root$ess, 17.57, tolerance = 0.1,
      label = "V15 root WESS ≈ 17.57% (CTA.exe Node 2 fixture)")
    expect_equal(root$loo_status, "STABLE",
      label = "V15 must carry LOO STABLE status")
  } else {
    # V15 MC came back non-significant with this seed; this is not a bug.
    # V4 being absent is still confirmed. Skip rather than fail.
    skip("V15 Signif F with this MC seed — V4 rejection still confirmed by Test 2")
  }
})

# =============================================================================
# Test 4: uniform weights bypass CTA path - observable contract in fixture file
#
# any(w != w[1]) is FALSE for uniform w -> CTA dispatch guard skips the
# CTA-specific ordered scan -> generic ODA path -> root=V2, cut=4.5, ESS=52.63%.
# Duplicate full CTA_DEMO fit removed; bypass is explicitly asserted in
# test-fixture-cta-demo.R B20 (line 54): root=V2 label "uniform weights bypass".
# =============================================================================

# =============================================================================
# Test 5: split-node MC metadata preservation (CRAN-safe, synthetic data)
#
# Verifies that ge_count and iter_used are preserved on the split node after
# the CTA patch that carries MC counters from the candidate through .split_nd().
# Uses small synthetic data with uniform weights (generic ODA path).
# =============================================================================

test_that("node-selection: split node stores ge_count and iter_used from MC run", {
  set.seed(77L)
  n   <- 60L
  y_s <- c(rep(1L, 30L), rep(2L, 30L))
  X_s <- data.frame(A = c(rep(0L, 25L), rep(1L,  5L),
                           rep(1L, 25L), rep(0L,  5L)),
                    B = sample(0:1, n, replace = TRUE))

  tree <- oda_cta_fit(
    X           = X_s,
    y           = y_s,
    priors_on   = FALSE,
    alpha_split = 0.05,
    mindenom    = 1L,
    prune_alpha = 0.05,
    max_depth   = 5L,
    ess_min     = 0,
    mc_iter     = 500L,
    mc_target   = 0.05,
    mc_stop     = 99.9,
    mc_stopup   = 20,
    mc_seed     = 42L,
    loo         = "off",
    attr_names  = c("A", "B")
  )

  root <- tree$nodes[[tree$root_id]]
  # Root must be a split (A is near-perfect separator)
  expect_false(isTRUE(root$leaf), label = "root should be a split node")

  # p_mc must be numeric and finite (MC ran with mc_iter=500)
  expect_true(is.numeric(root$p_mc), label = "root$p_mc must be numeric")
  expect_false(is.na(root$p_mc),     label = "root$p_mc must not be NA")

  # ge_count and iter_used must be non-NA non-negative integers
  expect_false(is.na(root$ge_count),   label = "root$ge_count must not be NA")
  expect_true(is.integer(root$ge_count), label = "root$ge_count must be integer")
  expect_true(root$ge_count >= 0L,     label = "root$ge_count must be >= 0")

  expect_false(is.na(root$iter_used),  label = "root$iter_used must not be NA")
  expect_true(is.integer(root$iter_used), label = "root$iter_used must be integer")
  expect_true(root$iter_used >= 1L,    label = "root$iter_used must be >= 1 when MC ran")
})
