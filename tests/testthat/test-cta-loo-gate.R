###############################################################################
# test-cta-loo-gate.R
#
# Regression test: CTA LOO gate correctly rejects an ordered-cut attribute
# with higher training ESS but LOO instability in favour of a lower-ESS
# LOO-stable attribute.
#
# Synthetic dataset (n=43, binary class):
#   Trap   (V2): ordered, ESS~49.21% at cut<=7.5, LOO UNSTABLE
#                (tie boundary at Trap=7; dropping a class-2 obs there
#                 shifts the winning cut from 7.5 to 6.5, changing
#                 the held-out prediction)
#   Stable (V3): binary,  ESS~42.86% at cut<=0.5, LOO STABLE
#
# CTA.exe canonical output (ENUMERATE; LOO STABLE; MC ITER 5000):
#   V3  Node 1  Tree 3  ESS 42.86%  Best ESS 42.86%
#   INCLUDE V3>0.5 (i.e. cut_value=0.5, Stable=1 → class 1)
#   Selected None (no further splits)
#
# R result must agree exactly.
###############################################################################

.make_loo_gate_data <- function() {
  add <- function(cls, trap, stable, n)
    data.frame(Class = rep(cls, n), Trap = rep(trap, n), Stable = rep(stable, n))

  rbind(
    # Trap <= 6.5 region: c1=1, c2=22
    add(1L,  6, 0, 1), add(2L, 1, 0, 4), add(2L, 2, 0, 4),
    add(2L,  3, 0, 4), add(2L, 4, 0, 4), add(2L, 5, 0, 3),
    add(2L,  6, 0, 3),
    # Trap=7: tie-boundary bin that makes Trap LOO unstable (c1=1, c2=6)
    add(1L,  7, 0, 1), add(2L, 7, 0, 6),
    # Trap=8: c1=1, c2=0, Stable=1
    add(1L,  8, 1, 1),
    # Trap>8: c1=4, c2=8
    add(1L,  9, 1, 2), add(2L,  9, 0, 4),
    add(1L, 10, 0, 2), add(2L, 10, 0, 4)
  )
}

.loo_gate_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      d <- .make_loo_gate_data()
      fit <<- cta_fit(
        X           = d[, c("Trap", "Stable")],
        y           = d$Class,
        priors_on   = TRUE,
        alpha_split = 0.05,
        mindenom    = 1L,
        prune_alpha = 0.05,
        max_depth   = 3L,
        ess_min     = 0,
        mc_iter     = 5000L,
        mc_target   = 0.05,
        mc_stop     = 99.9,
        mc_stopup   = 99.9,
        mc_seed     = NULL,
        loo         = "stable",
        attr_names  = c("Trap", "Stable")
      )
    }
    fit
  }
})

# ---- root selection ---------------------------------------------------------

test_that("loo-gate: root is Stable, not LOO-unstable Trap", {
  skip_if_not_smoke("cta-loo-gate")
  tree <- .loo_gate_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_false(isTRUE(root$leaf))
  expect_equal(root$attribute, "Stable")
})

test_that("loo-gate: root cut value is 0.5 (matches CTA.exe V3>0.5)", {
  skip_if_not_smoke("cta-loo-gate")
  tree <- .loo_gate_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$rule$cut_value, 0.5, tolerance = 1e-9)
})

test_that("loo-gate: root ESS ~42.86% (matches CTA.exe)", {
  skip_if_not_smoke("cta-loo-gate")
  tree <- .loo_gate_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(round(root$ess, 2), 42.86)
})

test_that("loo-gate: root LOO status is STABLE", {
  skip_if_not_smoke("cta-loo-gate")
  tree <- .loo_gate_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$loo_status, "STABLE")
})

# ---- tree structure (CTA.exe: "Selected None" = no further splits) ----------

test_that("loo-gate: only 1 split node — no further splits below root", {
  skip_if_not_smoke("cta-loo-gate")
  tree <- .loo_gate_fit()
  tbl  <- cta_node_table(tree)
  expect_equal(sum(!tbl$leaf), 1L)
})

test_that("loo-gate: root n_obs = 43", {
  skip_if_not_smoke("cta-loo-gate")
  tree <- .loo_gate_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$n_obs, 43L)
})

# ---- LOO semantics: mode-aware loo_status ------------------------------------
#
# Fixes audited in Slice C:
#   Bug 3: .loo_info() returned "STABLE" for all allowed LOO regardless of mode.
#          With the fix, loo="stable" -> "STABLE", loo=numeric/"pvalue" -> "PVALUE".
#
# These tests are CRAN-safe (no skip_if_not_smoke).  They use n=20 synthetic
# data (x=1:20, y=10×0 + 10×1) confirmed by probe:
#   loo=0.99  → no_tree=FALSE, root_loo_status="PVALUE"  (LOO p=1e-4)
#   loo=1e-9  → no_tree=TRUE   (deterministic gate rejection; any p >= 1e-9)
#
# C-1 (loo="stable") is intentionally absent: synthetic integer data is always
# LOO-UNSTABLE (boundary shifts on every fold → ESS_LOO ≠ ESS_training).
# STABLE coverage is in smoke-gated C-5 (myeloma LOO gate fixture).

.loo_sem_cta_xy <- function() {
  list(x = seq_len(20L), y = c(rep(0L, 10L), rep(1L, 10L)))
}

test_that("C-2: cta_fit loo=numeric (permissive) produces loo_status='PVALUE' on root node", {
  # n=20 perfectly separable; LOO p=1e-4 < 0.99 threshold → gate passes.
  # Root must be a non-leaf split node with loo_status="PVALUE" (not "STABLE").
  d <- .loo_sem_cta_xy()
  tree <- cta_fit(
    data.frame(V1 = d$x), d$y,
    mindenom = 1L, loo = 0.99, mc_iter = 1000L, mc_seed = 1L,
    prune_alpha = 1.0
  )
  expect_false(isTRUE(tree$no_tree),
               label = "n=20 loo=0.99: tree must be built (probe confirmed p_loo=1e-4)")
  root <- tree$nodes[[tree$root_id]]
  expect_false(isTRUE(root$leaf))
  expect_equal(root$loo_status, "PVALUE",
               label = "root loo_status must be 'PVALUE' for numeric loo=0.99")
})

test_that("C-3: cta_fit loo='pvalue' string produces loo_status='PVALUE' on root node", {
  # loo="pvalue" uses default threshold 0.05; LOO p=1e-4 < 0.05 → gate passes.
  # Root loo_status must be "PVALUE", same as the numeric path.
  d <- .loo_sem_cta_xy()
  tree <- cta_fit(
    data.frame(V1 = d$x), d$y,
    mindenom = 1L, loo = "pvalue", mc_iter = 1000L, mc_seed = 1L,
    prune_alpha = 1.0
  )
  expect_false(isTRUE(tree$no_tree),
               label = "n=20 loo='pvalue': tree must be built (LOO p < 0.05)")
  root <- tree$nodes[[tree$root_id]]
  expect_false(isTRUE(root$leaf))
  expect_equal(root$loo_status, "PVALUE",
               label = "root loo_status must be 'PVALUE' for loo='pvalue'")
})

test_that("C-4a: CTA weighted ordered path with loo=0.99 produces PVALUE status on root node", {
  # Non-uniform weights trigger the CTA-specific weighted ordered path.
  # With loo=0.99 and n=20, the LOO gate passes (p_loo=1e-4 < 0.99).
  # Root must exist and report loo_status="PVALUE" (not "STABLE" or "OFF").
  y  <- c(rep(0L, 10L), rep(1L, 10L))
  w  <- c(rep(2, 10), rep(1, 10))    # non-uniform → CTA-specific path
  X  <- data.frame(V1 = 1:20)
  tree <- cta_fit(X, y, w = w,
                  loo = 0.99, mc_iter = 1000L, mc_seed = 1L,
                  prune_alpha = 1.0)
  expect_false(isTRUE(tree$no_tree),
               label = "n=20 weighted loo=0.99: tree must be built (probe confirmed)")
  root <- tree$nodes[[tree$root_id]]
  expect_false(isTRUE(root$leaf))
  expect_equal(root$loo_status, "PVALUE",
               label = "weighted CTA root must show 'PVALUE' for loo=0.99")
})

test_that("C-4b (deterministic): CTA weighted ordered with loo=1e-9 is rejected — no_tree", {
  # loo=1e-9 is impossible to pass: any LOO p-value is >= 1e-9.
  # The gate must reject all weighted ordered candidates deterministically.
  # Probe confirmed: n=20 weighted loo=1e-9 → no_tree=TRUE (seed-independent).
  y  <- c(rep(0L, 10L), rep(1L, 10L))
  w  <- c(rep(2, 10), rep(1, 10))
  X  <- data.frame(V1 = 1:20)
  tree <- cta_fit(X, y, w = w,
                  loo = 1e-9, mc_iter = 1000L, mc_seed = 1L,
                  prune_alpha = 1.0)
  expect_true(isTRUE(tree$no_tree),
              label = "loo=1e-9 must reject all candidates → no_tree (deterministic)")
})

test_that("C-5: existing LOO STABLE canon unaffected — loo='stable' still selects Stable over LOO-unstable Trap", {
  skip_if_not_smoke("cta-loo-gate")
  # Regression lock: the myeloma-derived LOO gate fixture must still work
  # exactly as before.  The .loo_gate_fit() dataset uses loo="stable".
  tree <- .loo_gate_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_false(isTRUE(root$leaf))
  expect_equal(root$attribute,   "Stable")
  expect_equal(root$loo_status,  "STABLE")
  expect_equal(round(root$ess, 2), 42.86)
})
