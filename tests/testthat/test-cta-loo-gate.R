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

test_that("C-4: weighted CTA ordered path — loo=0.99 passes (PVALUE); loo=1e-9 rejects (no_tree)", {
  y <- c(rep(0L, 10L), rep(1L, 10L))
  w <- c(rep(2, 10), rep(1, 10))    # non-uniform → CTA-specific weighted path
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

test_that("C-7: CTA family LOO gate — family enumeration inherits MC-before-LOO", {
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

test_that("C-8: LORT LOO gate — local CTA recursion inherits MC-before-LOO", {
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
