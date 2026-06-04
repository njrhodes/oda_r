###############################################################################
# test-cta-loo-gate.R
#
# Regression tests: CTA LOO gate correctly rejects an ordered-cut attribute
# with higher training ESS but LOO instability in favour of a lower-ESS
# LOO-stable attribute.
#
# Canon dataset (n=43, binary class):
#   Trap   (V2): ordered, ESS~49.21% at cut<=7.5, LOO UNSTABLE
#   Stable (V3): binary,  ESS~42.86% at cut<=0.5, LOO STABLE
#
# CTA.exe canonical output: V3, cut=0.5, ESS=42.86%, STABLE, no further splits.
###############################################################################

# ---- Helpers ----------------------------------------------------------------

.make_loo_gate_data <- function() {
  add <- function(cls, trap, stable, n)
    data.frame(Class = rep(cls, n), Trap = rep(trap, n), Stable = rep(stable, n))
  rbind(
    add(1L,  6, 0, 1), add(2L, 1, 0, 4), add(2L, 2, 0, 4),
    add(2L,  3, 0, 4), add(2L, 4, 0, 4), add(2L, 5, 0, 3),
    add(2L,  6, 0, 3),
    add(1L,  7, 0, 1), add(2L, 7, 0, 6),
    add(1L,  8, 1, 1),
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
        mc_seed     = 42L,
        loo         = "stable",
        attr_names  = c("Trap", "Stable")
      )
    }
    fit
  }
})

.loo_sem_cta_xy <- function() {
  list(x = seq_len(20L), y = c(rep(0L, 10L), rep(1L, 10L)))
}

# =============================================================================
# Smoke: LOO STABLE gate canon (n=43 synthetic fixture)
# =============================================================================

test_that("loo-gate: root is Stable/cut=0.5/ESS~42.86%/STABLE; no further splits", {
  testthat::skip_on_cran()
  tree <- .loo_gate_fit()
  root <- tree$nodes[[tree$root_id]]
  # Root selection
  expect_false(isTRUE(root$leaf))
  expect_equal(root$attribute,    "Stable")
  expect_equal(root$rule$cut_value, 0.5, tolerance = 1e-9)
  expect_equal(round(root$ess, 2), 42.86)
  expect_equal(root$loo_status,   "STABLE")
  expect_equal(root$n_obs,        43L)
  # "Selected None" = no further splits below root
  tbl <- cta_node_table(tree)
  expect_equal(sum(!tbl$leaf), 1L)
  # C-5 regression lock: redundant check; confirms fix is stable across runs
  expect_equal(root$attribute,  "Stable")
  expect_equal(root$loo_status, "STABLE")
})

# =============================================================================
# CRAN-safe: loo_status semantics
# =============================================================================

test_that("C-2/C-3: numeric and 'pvalue' loo both produce loo_status='PVALUE' on root", {
  d <- .loo_sem_cta_xy()
  # loo=0.99 (numeric permissive): p=1e-4 < 0.99 -> gate passes, status PVALUE
  t1 <- cta_fit(data.frame(V1 = d$x), d$y,
                mindenom = 1L, loo = 0.99, mc_iter = 1000L, mc_seed = 1L,
                prune_alpha = 1.0)
  expect_false(isTRUE(t1$no_tree))
  expect_equal(t1$nodes[[t1$root_id]]$loo_status, "PVALUE")
  # loo="pvalue" (string): same semantics
  t2 <- cta_fit(data.frame(V1 = d$x), d$y,
                mindenom = 1L, loo = "pvalue", mc_iter = 1000L, mc_seed = 1L,
                prune_alpha = 1.0)
  expect_false(isTRUE(t2$no_tree))
  expect_equal(t2$nodes[[t2$root_id]]$loo_status, "PVALUE")
})

test_that("C-4: weighted CTA ordered path - loo=0.99 passes (PVALUE); loo=1e-9 rejects (no_tree)", {
  y <- c(rep(0L, 10L), rep(1L, 10L))
  w <- c(rep(2, 10), rep(1, 10))    # non-uniform -> CTA-specific weighted path
  X <- data.frame(V1 = 1:20)
  # C-4a: loo=0.99 passes
  ta <- cta_fit(X, y, w = w, loo = 0.99, mc_iter = 1000L, mc_seed = 1L, prune_alpha = 1.0)
  expect_false(isTRUE(ta$no_tree))
  expect_equal(ta$nodes[[ta$root_id]]$loo_status, "PVALUE")
  # C-4b: loo=1e-9 deterministically rejected
  tb <- cta_fit(X, y, w = w, loo = 1e-9, mc_iter = 1000L, mc_seed = 1L, prune_alpha = 1.0)
  expect_true(isTRUE(tb$no_tree))
})

# =============================================================================
# CRAN-safe: MC-gate ordering (LOO called only for MC-significant candidates)
# =============================================================================

test_that("C-6: LOO called only for MC-significant generic_oda candidates", {
  # CRAN-safe: mc_iter=10, n=20, deterministic.
  y_c6 <- c(rep(0L, 10L), rep(1L, 10L))
  X_c6 <- data.frame(
    V_signal = 1:20,
    V_noise  = rep(c(0L, 1L), 10L)
  )
  de <- new.env(parent = emptyenv())

  tree <- cta_fit(
    X = X_c6, y = y_c6,
    mindenom = 1L, priors_on = FALSE, alpha_split = 0.05, ess_min = 0,
    loo = 0.99, mc_iter = 10L, mc_seed = 1L, mc_stop = 99.9, mc_stopup = 99.9,
    prune_alpha = 1.0, diag_env = de
  )

  gen_mc  <- Filter(function(x) identical(x$path, "generic_oda"), de$mc_log)
  gen_loo <- Filter(function(x) identical(x$path, "generic_oda"), de$loo_log)

  # LOO was actually called for at least one candidate
  expect_gt(length(gen_loo), 0L,
            label = "LOO must be logged for at least one MC-significant candidate")
  # Every loo_log entry corresponds to an MC-significant candidate
  loo_attrs    <- vapply(gen_loo, `[[`, character(1L), "attr_name")
  mc_pass      <- vapply(Filter(function(x) isTRUE(x$signif_current), gen_mc),
                         `[[`, character(1L), "attr_name")
  expect_true(all(loo_attrs %in% mc_pass),
              label = "all loo_log entries must be MC-significant")
  # V_noise failed MC -> NOT in loo_log
  noise_mc <- Filter(function(x) x$attr_name == "V_noise", gen_mc)
  if (length(noise_mc) > 0L) {
    expect_false(isTRUE(noise_mc[[1L]]$signif_current))
    expect_false("V_noise" %in% loo_attrs)
  }
  expect_false(isTRUE(tree$no_tree))
})

# =============================================================================
# CRAN-safe: CTA family LOO gate inherits MC-before-LOO
# =============================================================================

test_that("C-7: CTA family LOO gate - family enumeration inherits MC-before-LOO", {
  y_c7 <- c(rep(0L, 10L), rep(1L, 10L))
  X_c7 <- data.frame(V_signal = 1:20, V_noise = rep(c(0L, 1L), 10L))
  de <- new.env(parent = emptyenv())

  fam <- cta_descendant_family(
    X_c7, y_c7,
    alpha_split = 0.05, ess_min = 0, loo = 0.99,
    mc_iter = 10L, mc_seed = 1L, mc_stop = 99.9, mc_stopup = 99.9,
    prune_alpha = 1.0, priors_on = FALSE, diag_env = de
  )

  expect_gt(length(fam$members), 0L)

  gen_mc  <- Filter(function(x) identical(x$path, "generic_oda"), de$mc_log)
  gen_loo <- Filter(function(x) identical(x$path, "generic_oda"), de$loo_log)

  expect_true(any(vapply(gen_mc, function(z) isFALSE(z$signif_current), logical(1L))),
              label = "at least one Signif F candidate must exist")

  if (length(gen_loo) > 0L) {
    loo_attrs    <- vapply(gen_loo, `[[`, character(1L), "attr_name")
    mc_pass      <- vapply(Filter(function(x) isTRUE(x$signif_current), gen_mc),
                           `[[`, character(1L), "attr_name")
    expect_true(all(loo_attrs %in% mc_pass))
  }

  noise_mc <- Filter(function(x) x$attr_name == "V_noise", gen_mc)
  if (length(noise_mc) > 0L) {
    expect_true(all(vapply(noise_mc, function(z) isFALSE(z$signif_current), logical(1L))))
    loo_v <- if (length(gen_loo) > 0L)
      vapply(gen_loo, `[[`, character(1L), "attr_name") else character(0L)
    expect_false("V_noise" %in% loo_v)
  }
})

# =============================================================================
# CRAN-safe: LORT LOO gate inherits MC-before-LOO
# =============================================================================

test_that("C-8: LORT LOO gate - local CTA recursion inherits MC-before-LOO", {
  set.seed(NULL)
  y_c8 <- c(rep(0L, 15L), rep(1L, 15L))
  X_c8 <- data.frame(V_signal = 1:30, V_noise = rep(c(0L, 1L), 15L))
  de <- new.env(parent = emptyenv())

  ort <- cta_fit(
    X = X_c8, y = y_c8,
    recursive = TRUE, min_n = 5L, max_depth = 3L, max_nodes = 15L,
    family_max_steps = 5L, alpha_split = 0.05, prune_alpha = 1.0,
    loo = 0.99, mc_iter = 10L, mc_seed = 1L, mc_stop = 99.9, mc_stopup = 99.9,
    diag_env = de
  )

  expect_true(inherits(ort, "cta_ort"))

  gen_mc  <- Filter(function(x) identical(x$path, "generic_oda"), de$mc_log)
  gen_loo <- Filter(function(x) identical(x$path, "generic_oda"), de$loo_log)

  expect_true(any(vapply(gen_mc, function(z) isFALSE(z$signif_current), logical(1L))),
              label = "at least one Signif F candidate must exist")

  if (length(gen_loo) > 0L) {
    loo_attrs <- vapply(gen_loo, `[[`, character(1L), "attr_name")
    mc_pass   <- vapply(Filter(function(x) isTRUE(x$signif_current), gen_mc),
                        `[[`, character(1L), "attr_name")
    expect_true(all(loo_attrs %in% mc_pass))
  }

  noise_mc <- Filter(function(x) x$attr_name == "V_noise", gen_mc)
  if (length(noise_mc) > 0L) {
    expect_true(all(vapply(noise_mc, function(z) isFALSE(z$signif_current), logical(1L))))
    loo_v <- if (length(gen_loo) > 0L)
      vapply(gen_loo, `[[`, character(1L), "attr_name") else character(0L)
    expect_false("V_noise" %in% loo_v)
  }
})

# =============================================================================
# CRAN-safe: binary predictor LOO wiring through CTA generic ODA path
# =============================================================================

test_that("C-9: binary predictor LOO routes through oda_loo_binary_map_counts() via CTA generic path", {
  # B1-structure binary data: n00=8, n01=1, n10=1, n11=8 (ESS ~77.78%).
  # V_noise: alternating 0/1 balanced equally across classes (ESS = 0, Signif F).
  # Uniform weights -> CTA-specific ordered path bypassed; generic ODA path used.
  # Binary attribute (<=2 unique values) -> oda_fit() produces binary_map rule ->
  # oda_loo_for_rule() dispatches to oda_loo_binary_map_counts().
  # loo=0.99: permissive p-value gate; V_binary passes, loo_status="PVALUE".
  x_bin <- c(rep(0L, 9L), rep(1L, 9L))
  y_c9  <- c(rep(0L, 8L), 1L, 0L, rep(1L, 8L))   # n00=8, n01=1, n10=1, n11=8
  x_nz  <- rep(c(0L, 1L, 0L), 6L)                 # balanced noise: 0 ESS

  de <- new.env(parent = emptyenv())

  cta_fit(
    X           = data.frame(V_binary = x_bin, V_noise = x_nz),
    y           = y_c9,
    mindenom    = 1L, priors_on = TRUE, alpha_split = 0.05, ess_min = 0,
    loo         = 0.99,
    mc_iter     = 1000L, mc_seed = 1L, mc_stop = 99.9,
    prune_alpha = 1.0,
    diag_env    = de
  )

  # Step 1: LOO must have run for V_binary (it must pass MC gate)
  bin_loo <- Filter(function(z) z$attr_name == "V_binary", de$loo_log)
  expect_gt(length(bin_loo), 0L,
    label = "C-9: loo_log must have an entry for V_binary")
  bl <- bin_loo[[1L]]
  expect_identical(bl$loo_mode, "n_fold_or_algebraic",
    label = "C-9: binary LOO mode must be n_fold_or_algebraic (not failed)")
  expect_false(is.na(bl$ess_loo),
    label = "C-9: binary ess_loo must not be NA")

  # Step 2: oracle -- direct oda_loo_for_rule() on the same x/y
  fit_ref <- oda_fit(x = x_bin, y = y_c9, priors_on = TRUE,
                     mcarlo = FALSE, loo = "off")
  expect_true(isTRUE(fit_ref$ok),   label = "C-9: oracle oda_fit must succeed")
  expect_equal(fit_ref$rule$type, "binary_map",
    label = "C-9: binary attribute must produce binary_map rule")

  loo_ref <- oda:::oda_loo_for_rule(
    x         = x_bin,
    y         = y_c9,
    rule      = fit_ref$rule,
    attr_type = "binary",
    priors_on = TRUE
  )
  expect_true(isTRUE(loo_ref$allowed), label = "C-9: oracle loo_for_rule must be allowed")

  # Step 3: CTA internal ess_loo must equal oracle ess_loo.
  # Proves oda_loo_binary_map_counts() is exercised through the CTA generic path,
  # not only through the direct oda_loo_for_rule() proof tests (B1-B4).
  expect_equal(bl$ess_loo, loo_ref$ess_loo, tolerance = 1e-6,
    label = paste0(
      "C-9: CTA binary ess_loo (", round(bl$ess_loo, 4), "%) ",
      "must match oracle oda_loo_for_rule ess_loo (", round(loo_ref$ess_loo, 4), "%)"
    )
  )
})

# =============================================================================
# CRAN-safe: weighted binary predictor LOO wiring through CTA generic ODA path
# =============================================================================

test_that("C-10: weighted binary predictor LOO routes through oda_loo_binary_map_counts() via CTA generic path", {
  # B3-structure weighted data: same x/y as B3 (n00=8,n01=1,n10=1,n11=8),
  # w_c10: x=0 obs weight=2, x=1 obs weight=1 (WESS ~74.12%).
  # Non-uniform weights with binary predictor (<=2 unique x values):
  # CTA-specific ordered path requires >2 unique x values; binary falls through
  # to generic ODA path -> oda_loo_binary_map_counts() handles weighted case.
  # V_noise: balanced noise, expected Signif F.
  x_bin  <- c(rep(0L, 9L), rep(1L, 9L))
  y_c10  <- c(rep(0L, 8L), 1L, 0L, rep(1L, 8L))   # n00=8, n01=1, n10=1, n11=8
  w_c10  <- c(rep(2, 9L), rep(1, 9L))               # non-uniform: x=0 obs weight=2
  x_nz   <- rep(c(0L, 1L, 0L), 6L)                 # balanced noise: 0 ESS

  de <- new.env(parent = emptyenv())

  cta_fit(
    X           = data.frame(V_binary = x_bin, V_noise = x_nz),
    y           = y_c10,
    w           = w_c10,
    mindenom    = 1L, priors_on = TRUE, alpha_split = 0.05, ess_min = 0,
    loo         = 0.99,
    mc_iter     = 1000L, mc_seed = 1L, mc_stop = 99.9,
    prune_alpha = 1.0,
    diag_env    = de
  )

  # Step 1: LOO must have run for V_binary
  bin_loo <- Filter(function(z) z$attr_name == "V_binary", de$loo_log)
  expect_gt(length(bin_loo), 0L,
    label = "C-10: loo_log must have an entry for V_binary (weighted)")
  bl <- bin_loo[[1L]]
  expect_identical(bl$loo_mode, "n_fold_or_algebraic",
    label = "C-10: weighted binary LOO mode must be n_fold_or_algebraic (not failed)")
  expect_false(is.na(bl$ess_loo),
    label = "C-10: weighted binary ess_loo must not be NA")

  # Step 2: oracle -- direct oda_loo_for_rule() with same x/y/w
  fit_ref <- oda_fit(x = x_bin, y = y_c10, w = w_c10, priors_on = TRUE,
                     mcarlo = FALSE, loo = "off")
  expect_true(isTRUE(fit_ref$ok),   label = "C-10: oracle oda_fit must succeed")
  expect_equal(fit_ref$rule$type, "binary_map",
    label = "C-10: binary attribute must produce binary_map rule")

  loo_ref <- oda:::oda_loo_for_rule(
    x         = x_bin,
    y         = y_c10,
    w         = w_c10,
    rule      = fit_ref$rule,
    attr_type = "binary",
    priors_on = TRUE
  )
  expect_true(isTRUE(loo_ref$allowed), label = "C-10: oracle loo_for_rule must be allowed")

  # Step 3: CTA internal ess_loo must equal oracle ess_loo.
  # Proves oda_loo_binary_map_counts() handles case weights correctly through
  # the CTA generic path (not only through direct B3/B4 proof tests).
  expect_equal(bl$ess_loo, loo_ref$ess_loo, tolerance = 1e-6,
    label = paste0(
      "C-10: weighted CTA binary ess_loo (", round(bl$ess_loo, 4), "%) ",
      "must match oracle oda_loo_for_rule ess_loo (", round(loo_ref$ess_loo, 4), "%)"
    )
  )
})

# =============================================================================
# CRAN-safe: LORT recursive=TRUE binary LOO forwarding
# =============================================================================

test_that("C-11: LORT recursive=TRUE forwards binary LOO to oda_loo_binary_map_counts()", {
  # B1-structure binary data: n00=8, n01=1, n10=1, n11=8 (ESS ~77.78%).
  # V_noise: balanced, expected Signif F.
  # cta_fit(recursive=TRUE) -> .cta_ort_fit() -> .cta_ort_fit_internal() ->
  # cta_descendant_family() -> oda_cta_fit() -> .full_fit_one() (generic ODA path
  # for binary predictor) -> oda_loo_for_rule() -> oda_loo_binary_map_counts().
  # diag_env is forwarded intact through the full LORT chain.
  # Root-level entry (n_obs == 18) is the relevant one; after V_binary split each
  # child has only one unique V_binary value so V_binary is not re-evaluated.
  x_bin   <- c(rep(0L, 9L), rep(1L, 9L))
  y_c11   <- c(rep(0L, 8L), 1L, 0L, rep(1L, 8L))   # n00=8, n01=1, n10=1, n11=8
  x_nz    <- rep(c(0L, 1L, 0L), 6L)                 # balanced noise
  n_c11   <- length(y_c11)

  de <- new.env(parent = emptyenv())

  cta_fit(
    X                = data.frame(V_binary = x_bin, V_noise = x_nz),
    y                = y_c11,
    recursive        = TRUE,
    min_n            = 5L, max_depth = 3L, max_nodes = 15L,
    family_max_steps = 5L,
    priors_on = TRUE, alpha_split = 0.05, ess_min = 0,
    loo              = 0.99,
    mc_iter          = 1000L, mc_seed = 1L, mc_stop = 99.9,
    prune_alpha      = 1.0,
    diag_env         = de
  )

  # Root-level binary LOO entry (n_obs == n_c11 confirms root node evaluation)
  bin_loo <- Filter(
    function(z) z$attr_name == "V_binary" && isTRUE(z$n_obs == n_c11),
    de$loo_log
  )
  expect_gt(length(bin_loo), 0L,
    label = "C-11: LORT loo_log must contain a root-level V_binary entry")
  bl <- bin_loo[[1L]]
  expect_identical(bl$loo_mode, "n_fold_or_algebraic",
    label = "C-11: LORT binary LOO mode must be n_fold_or_algebraic (not failed)")
  expect_false(is.na(bl$ess_loo),
    label = "C-11: LORT binary ess_loo must not be NA")

  # Oracle -- same x/y, no weights
  fit_ref <- oda_fit(x = x_bin, y = y_c11, priors_on = TRUE,
                     mcarlo = FALSE, loo = "off")
  expect_true(isTRUE(fit_ref$ok), label = "C-11: oracle oda_fit must succeed")
  expect_equal(fit_ref$rule$type, "binary_map",
    label = "C-11: binary attribute must produce binary_map rule")

  loo_ref <- oda:::oda_loo_for_rule(
    x         = x_bin,
    y         = y_c11,
    rule      = fit_ref$rule,
    attr_type = "binary",
    priors_on = TRUE
  )
  expect_true(isTRUE(loo_ref$allowed), label = "C-11: oracle loo_for_rule must be allowed")

  # CTA-internal ess_loo (via LORT -> cta_descendant_family -> oda_cta_fit chain)
  # must equal oracle ess_loo, proving the algebraic path is not bypassed by LORT.
  expect_equal(bl$ess_loo, loo_ref$ess_loo, tolerance = 1e-6,
    label = paste0(
      "C-11: LORT binary ess_loo (", round(bl$ess_loo, 4), "%) ",
      "must match oracle oda_loo_for_rule ess_loo (", round(loo_ref$ess_loo, 4), "%)"
    )
  )
})
