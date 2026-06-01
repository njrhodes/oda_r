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
#   INCLUDE V3>0.5 (i.e. cut_value=0.5, Stable=1 -> class 1)
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

test_that("loo-gate: only 1 split node -- no further splits below root", {
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
# data (x=1:20, y=10x0 + 10x1) confirmed by probe:
#   loo=0.99  -> no_tree=FALSE, root_loo_status="PVALUE"  (LOO p=1e-4)
#   loo=1e-9  -> no_tree=TRUE   (deterministic gate rejection; any p >= 1e-9)
#
# C-1 (loo="stable") is intentionally absent: synthetic integer data is always
# LOO-UNSTABLE (boundary shifts on every fold -> ESS_LOO != ESS_training).
# STABLE coverage is in smoke-gated C-5 (myeloma LOO gate fixture).

.loo_sem_cta_xy <- function() {
  list(x = seq_len(20L), y = c(rep(0L, 10L), rep(1L, 10L)))
}

test_that("C-2: cta_fit loo=numeric (permissive) produces loo_status='PVALUE' on root node", {
  # n=20 perfectly separable; LOO p=1e-4 < 0.99 threshold -> gate passes.
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
  # loo="pvalue" uses default threshold 0.05; LOO p=1e-4 < 0.05 -> gate passes.
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
  w  <- c(rep(2, 10), rep(1, 10))    # non-uniform -> CTA-specific path
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

test_that("C-4b (deterministic): CTA weighted ordered with loo=1e-9 is rejected -- no_tree", {
  # loo=1e-9 is impossible to pass: any LOO p-value is >= 1e-9.
  # The gate must reject all weighted ordered candidates deterministically.
  # Probe confirmed: n=20 weighted loo=1e-9 -> no_tree=TRUE (seed-independent).
  y  <- c(rep(0L, 10L), rep(1L, 10L))
  w  <- c(rep(2, 10), rep(1, 10))
  X  <- data.frame(V1 = 1:20)
  tree <- cta_fit(X, y, w = w,
                  loo = 1e-9, mc_iter = 1000L, mc_seed = 1L,
                  prune_alpha = 1.0)
  expect_true(isTRUE(tree$no_tree),
              label = "loo=1e-9 must reject all candidates -> no_tree (deterministic)")
})

test_that("C-5: existing LOO STABLE canon unaffected -- loo='stable' still selects Stable over LOO-unstable Trap", {
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

# ---- C-6: MC gate ordering -- LOO called only for MC-significant candidates --
#
# Canon invariant: EXE evaluates MC significance first; LOO only runs for
# Signif T candidates.  Running LOO unconditionally for all candidates
# (including Signif F) is non-canon and wastes n-fold LOO work.
#
# Synthetic dataset (n=20, binary class, equal sizes):
#   V_signal: 1:20  -- perfect separator (ESS=50% for C=2 equal class).
#             With mc_iter=10, mc_seed=1: ge_count=0 deterministically
#             (no permutation of y can sort 1:20 by accident in 10 tries;
#              P(ge_count>=1) < 0.01%).
#   V_noise:  rep(c(0L,1L), 10L) -- equal class distribution in both classes.
#             Cut <= 0.5 gives PAC=50%, ESS=0%.  ALL permutations of y achieve
#             ESS >= 0%, so ge_count=mc_iter=10 and p_mc=1.0 deterministically.
#
# With loo=0.99 (permissive numeric):
#   - V_signal: MC passes (p=0 < 0.05) -> LOO called -> appears in loo_log
#   - V_noise:  MC fails (p=1.0 >= 0.05) -> LOO NOT called -> absent from loo_log
#
# Regression: old code called oda_fit(..., loo=loo_arg, eval_order="loo_then_mc"),
# causing LOO to run inside unioda_core for all candidates unconditionally.
# That path did NOT write to diag_env$loo_log, so assertion (1) below would
# fail -- diag_env$loo_log would be empty even though LOO ran internally.

test_that("C-6: LOO called only for MC-significant generic_oda candidates", {
  # CRAN-safe: mc_iter=10, n=20, deterministic.
  y_c6  <- c(rep(0L, 10L), rep(1L, 10L))
  X_c6  <- data.frame(
    V_signal = 1:20,
    V_noise  = rep(c(0L, 1L), 10L)
  )
  de <- new.env(parent = emptyenv())

  tree <- cta_fit(
    X           = X_c6,
    y           = y_c6,
    mindenom    = 1L,
    priors_on   = FALSE,
    alpha_split = 0.05,
    ess_min     = 0,
    loo         = 0.99,
    mc_iter     = 10L,
    mc_seed     = 1L,
    mc_stop     = 99.9,
    mc_stopup   = 99.9,
    prune_alpha = 1.0,
    diag_env    = de
  )

  # Separate generic_oda entries from CTA-specific entries.
  gen_mc  <- Filter(function(x) identical(x$path, "generic_oda"), de$mc_log)
  gen_loo <- Filter(function(x) identical(x$path, "generic_oda"), de$loo_log)

  # (1) LOO was actually called for at least one generic_oda candidate.
  #     If regression re-introduces old code (LOO inside oda_fit, invisible to
  #     diag_env), this assertion fails -- loo_log would be empty.
  expect_gt(length(gen_loo), 0L,
            label = "LOO must be logged for at least one MC-significant candidate")

  # (2) Every attr in loo_log passed MC (signif_current=TRUE in mc_log).
  loo_attrs <- vapply(gen_loo, `[[`, character(1L), "attr_name")
  mc_pass_attrs <- vapply(
    Filter(function(x) isTRUE(x$signif_current), gen_mc),
    `[[`, character(1L), "attr_name"
  )
  expect_true(all(loo_attrs %in% mc_pass_attrs),
              label = "all loo_log entries must correspond to MC-significant candidates")

  # (3) V_noise failed MC -> must NOT appear in loo_log.
  noise_mc <- Filter(function(x) x$attr_name == "V_noise", gen_mc)
  if (length(noise_mc) > 0L) {
    expect_false(isTRUE(noise_mc[[1L]]$signif_current),
                 label = "V_noise must be Signif F (ESS=0%: all permutations match)")
    expect_false("V_noise" %in% loo_attrs,
                 label = "V_noise must NOT appear in loo_log (MC failed)")
  }

  # (4) Tree must be built (V_signal is the only admissible root).
  expect_false(isTRUE(tree$no_tree),
               label = "tree must be built -- V_signal is a valid root")
})

# ---- C-7: CTA family LOO gate -- family enumeration inherits MC-before-LOO --
#
# cta_descendant_family() passes ... to oda_cta_fit(), which uses .full_fit_one()
# with the fixed MC -> LOO gate ordering.  Each family member fit appends to the
# same diag_env, so mc_log and loo_log accumulate across all MINDENOM steps.
#
# Same V_signal/V_noise dataset as C-6 (n=20, deterministic).
# Family terminates quickly: V_signal stump until MINDENOM > 10 -> no_tree.

test_that("CTA family LOO gate: family enumeration inherits MC-before-LOO", {
  # CRAN-safe: mc_iter=10, n=20, deterministic.
  y_c7 <- c(rep(0L, 10L), rep(1L, 10L))
  X_c7 <- data.frame(
    V_signal = 1:20,
    V_noise  = rep(c(0L, 1L), 10L)
  )
  de <- new.env(parent = emptyenv())

  fam <- cta_descendant_family(
    X_c7, y_c7,
    alpha_split = 0.05,
    ess_min     = 0,
    loo         = 0.99,
    mc_iter     = 10L,
    mc_seed     = 1L,
    mc_stop     = 99.9,
    mc_stopup   = 99.9,
    prune_alpha = 1.0,
    priors_on   = FALSE,
    diag_env    = de
  )

  # Family must have produced at least one member.
  expect_gt(length(fam$members), 0L,
            label = "family must have at least one member")

  gen_mc  <- Filter(function(x) identical(x$path, "generic_oda"), de$mc_log)
  gen_loo <- Filter(function(x) identical(x$path, "generic_oda"), de$loo_log)

  # (1) At least one MC-failed candidate across all family members.
  expect_true(
    any(vapply(gen_mc, function(z) isFALSE(z$signif_current), logical(1L))),
    label = "at least one Signif F candidate must exist -- otherwise test is vacuous"
  )

  # (2) Every LOO log entry corresponds to an MC-significant candidate.
  if (length(gen_loo) > 0L) {
    loo_attrs <- vapply(gen_loo, `[[`, character(1L), "attr_name")
    mc_pass_attrs <- vapply(
      Filter(function(x) isTRUE(x$signif_current), gen_mc),
      `[[`, character(1L), "attr_name"
    )
    expect_true(all(loo_attrs %in% mc_pass_attrs),
                label = "all family loo_log entries must be MC-significant")
  }

  # (3) V_noise failed MC in every family member -> never in loo_log.
  noise_mc <- Filter(function(x) x$attr_name == "V_noise", gen_mc)
  if (length(noise_mc) > 0L) {
    all_noise_signif_f <- all(vapply(noise_mc, function(z) isFALSE(z$signif_current), logical(1L)))
    expect_true(all_noise_signif_f,
                label = "V_noise must be Signif F in every family member")
    loo_attrs_v <- if (length(gen_loo) > 0L)
      vapply(gen_loo, `[[`, character(1L), "attr_name") else character(0L)
    expect_false("V_noise" %in% loo_attrs_v,
                 label = "V_noise must NOT appear in any family member's loo_log")
  }
})

# ---- C-8: LORT LOO gate -- local CTA recursion inherits MC-before-LOO -------
#
# lort_fit() -> cta_fit(recursive=TRUE) -> .cta_ort_fit() ->
#   .cta_ort_fit_internal() -> cta_descendant_family() -> oda_cta_fit() ->
#   .full_fit_one() [FIXED].
#
# diag_env is now threaded through the full LORT chain so we can observe mc_log
# and loo_log accumulated across all ORT node family scans.
#
# Dataset (n=30): three groups of 10, two predictor structure.
#   V_signal: perfectly separating ordered predictor.
#   V_noise:  balanced binary predictor (ESS=0% -> Signif F deterministically).
# min_n=5 ensures recursion can proceed into sub-strata of n>=5.

test_that("LORT LOO gate: local CTA recursion inherits MC-before-LOO", {
  # CRAN-safe: mc_iter=10, n=30, deterministic.
  set.seed(NULL)  # do not inherit any global seed
  y_c8 <- c(rep(0L, 15L), rep(1L, 15L))
  X_c8 <- data.frame(
    V_signal = 1:30,
    V_noise  = rep(c(0L, 1L), 15L)
  )
  de <- new.env(parent = emptyenv())

  ort <- cta_fit(
    X                = X_c8,
    y                = y_c8,
    recursive        = TRUE,
    min_n            = 5L,
    max_depth        = 3L,
    max_nodes        = 15L,
    family_max_steps = 5L,
    alpha_split      = 0.05,
    prune_alpha      = 1.0,
    loo              = 0.99,
    mc_iter          = 10L,
    mc_seed          = 1L,
    mc_stop          = 99.9,
    mc_stopup        = 99.9,
    diag_env         = de
  )

  # LORT object must be valid.
  expect_true(inherits(ort, "cta_ort"),
              label = "LORT fit must return a cta_ort object")

  gen_mc  <- Filter(function(x) identical(x$path, "generic_oda"), de$mc_log)
  gen_loo <- Filter(function(x) identical(x$path, "generic_oda"), de$loo_log)

  # (1) At least one MC-failed candidate must have been evaluated (V_noise).
  expect_true(
    any(vapply(gen_mc, function(z) isFALSE(z$signif_current), logical(1L))),
    label = "at least one Signif F candidate must exist across all ORT node fits"
  )

  # (2) Every loo_log entry corresponds to an MC-significant candidate.
  if (length(gen_loo) > 0L) {
    loo_attrs <- vapply(gen_loo, `[[`, character(1L), "attr_name")
    mc_pass_attrs <- vapply(
      Filter(function(x) isTRUE(x$signif_current), gen_mc),
      `[[`, character(1L), "attr_name"
    )
    expect_true(all(loo_attrs %in% mc_pass_attrs),
                label = "all LORT loo_log entries must be MC-significant")
  }

  # (3) V_noise must never appear in loo_log across any ORT node fit.
  noise_mc <- Filter(function(x) x$attr_name == "V_noise", gen_mc)
  if (length(noise_mc) > 0L) {
    all_noise_f <- all(vapply(noise_mc, function(z) isFALSE(z$signif_current), logical(1L)))
    expect_true(all_noise_f,
                label = "V_noise must be Signif F in every ORT-node family fit")
    loo_attrs_v <- if (length(gen_loo) > 0L)
      vapply(gen_loo, `[[`, character(1L), "attr_name") else character(0L)
    expect_false("V_noise" %in% loo_attrs_v,
                 label = "V_noise must NOT appear in any ORT-node loo_log")
  }
})
