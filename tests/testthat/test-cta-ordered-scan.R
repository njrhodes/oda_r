###############################################################################
# test-cta-ordered-scan.R
#
# Unit tests for .cta_ordered_scan() / .cta_fast_scan_perm().
#
# Rule (docs/CTA_ORDERED_CUT_AUDIT.md): maximise count-priors weighted
# accuracy sum(w x correct)/n_class over admissible cuts; MINDENOM gates the
# selected cut; reported ess is the usual WESS.  Anchors: every CTA.exe V4
# node report in the myeloma MINDENOM=1 output (11 nodes, Trees 1-3).
###############################################################################

# ---- Myeloma data helpers ---------------------------------------------------

.myeloma_cta_dat <- local({
  dat <- NULL
  function() {
    if (is.null(dat)) {
      fp <- tryCatch(
        testthat::test_path("fixtures", "myeloma", "data.txt"),
        error = function(e) NULL
      )
      if (is.null(fp) || !file.exists(fp)) return(NULL)
      raw <- read.table(fp, header = FALSE, col.names = paste0("V", 1:19))
      dat <<- raw[raw$V2 != 0, ]   # EX V2=0
    }
    dat
  }
})

.node4_rows <- function(dat) {
  miss <- -9L
  with(dat, which(!(V17 %in% miss) & V17 <= 0.5 &
                  !(V15 %in% miss) & V15 <= 0.5))
}

.node2_rows <- function(dat) {
  miss <- -9L
  with(dat, which(!(V17 %in% miss) & V17 <= 0.5))
}

.recode01 <- function(y) {
  labs <- sort(unique(as.integer(y)))
  ifelse(as.integer(y) == labs[1L], 0L, 1L)
}

# =============================================================================
# .cta_ordered_scan() - myeloma-specific audit evidence
# =============================================================================

test_that(".cta_ordered_scan: myeloma Node 4 V4 returns cut 359.80 ESS ~34.90%", {
  dat <- .myeloma_cta_dat()
  if (is.null(dat)) skip("myeloma fixture missing")

  idx <- .node4_rows(dat)
  x   <- dat$V4[idx]
  y01 <- .recode01(dat$V1[idx])
  w   <- dat$V2[idx]

  res <- oda:::.cta_ordered_scan(x, y01, w,
                                     priors_on  = TRUE,
                                     miss_codes = -9,
                                     mindenom   = 1L)

  expect_false(is.null(res),              label = "result is not NULL")
  expect_equal(res$rule$type,      "ordered_cut")
  expect_equal(res$rule$direction, "0->1")
  expect_equal(res$rule$cut_value, 359.80, tolerance = 0.01,
               label = "cut value matches CTA.exe myeloma MINDENOM=1 fixture")
  expect_equal(res$ess,            34.90,  tolerance = 0.01,
               label = "WESS matches CTA.exe myeloma MINDENOM=1 fixture")
  expect_gt(res$sens_wa, 0.5)
})

test_that(".cta_ordered_scan: myeloma Node 2 V4 returns cut 371.20 ESS ~34.89%", {
  dat <- .myeloma_cta_dat()
  if (is.null(dat)) skip("myeloma fixture missing")

  idx <- .node2_rows(dat)
  x   <- dat$V4[idx]
  y01 <- .recode01(dat$V1[idx])
  w   <- dat$V2[idx]

  res <- oda:::.cta_ordered_scan(x, y01, w,
                                     priors_on  = TRUE,
                                     miss_codes = -9,
                                     mindenom   = 1L)

  expect_false(is.null(res),              label = "result is not NULL")
  expect_equal(res$rule$cut_value, 371.20, tolerance = 0.01,
               label = "cut value - CTA.exe myeloma Node 2 V4 audit")
  expect_equal(res$ess,            34.89,  tolerance = 0.01,
               label = "WESS matches CTA.exe myeloma fixture")
  expect_gt(res$sens_wa, 0.5, label = "sens_wa > 0.5")
})

# =============================================================================
# Generic oda_univariate_core() now uses count priors for weighted data
# =============================================================================
# Before 0.1.5 the generic weighted path maximised WESS (cut 85.70, WESS 46.99).
# It now uses the same count-priors objective as CTA.exe / MegaODA, so it
# agrees with the CTA scan and the CTA.exe Node 4 V4 audit (cut 359.80, 34.90%).

test_that("generic oda_univariate_core: myeloma Node 4 V4 matches CTA.exe cut 359.80", {
  dat <- .myeloma_cta_dat()
  if (is.null(dat)) skip("myeloma fixture missing")

  idx <- .node4_rows(dat)
  x   <- dat$V4[idx]
  y01 <- .recode01(dat$V1[idx])
  w   <- dat$V2[idx]

  fit <- oda_univariate_core(
    x = x, y = y01, w = w,
    attr_type  = "ordered",
    priors_on  = TRUE,
    miss_codes = -9,
    mindenom   = 1L,
    mcarlo     = FALSE,
    loo        = "off"
  )

  expect_true(isTRUE(fit$ok), label = "generic ODA returns a result")
  expect_equal(fit$rule$cut_value, 359.80, tolerance = 0.01,
               label = "generic ODA cut = CTA.exe Node 4 V4 (359.80)")
  expect_equal(fit$ess, 34.90, tolerance = 0.05,
               label = "generic ODA WESS = CTA.exe Node 4 V4 (34.90%)")
  scan <- oda:::.cta_ordered_scan(x, y01, w, TRUE, -9, 1L)
  expect_equal(fit$rule$cut_value, scan$rule$cut_value, tolerance = 1e-9,
               label = "generic ODA and CTA scan choose the same cut")
  expect_identical(fit$rule$direction, scan$rule$direction)
})

# =============================================================================
# Edge cases for .cta_ordered_scan
# =============================================================================

test_that(".cta_ordered_scan: returns NULL for pure node", {
  x <- 1:6
  y <- rep(0L, 6)
  w <- rep(1, 6)
  expect_null(oda:::.cta_ordered_scan(x, y, w, TRUE, NULL, 1L))
})

test_that(".cta_ordered_scan: returns NULL when all obs missing", {
  x <- rep(-9, 6)
  y <- c(0L, 0L, 0L, 1L, 1L, 1L)
  w <- rep(1, 6)
  expect_null(oda:::.cta_ordered_scan(x, y, w, TRUE, -9, 1L))
})

test_that(".cta_ordered_scan: does not error under varying mindenom", {
  set.seed(1)
  x <- 1:8; y <- c(0L,0L,0L,0L,1L,1L,1L,1L); w <- rep(1, 8)
  r1 <- oda:::.cta_ordered_scan(x, y, w, TRUE, NULL, 1L)
  r5 <- oda:::.cta_ordered_scan(x, y, w, TRUE, NULL, 5L)
  expect_false(is.null(r1), label = "mindenom=1 finds a cut")
  expect_true(is.null(r5) || is.list(r5), label = "mindenom=5 returns NULL or list")
})

# =============================================================================
# CTA.exe parity: every weighted V4 node report, myeloma MINDENOM=1
# =============================================================================

test_that(".cta_ordered_scan: reproduces all 11 CTA.exe myeloma V4 node WESS values", {
  dat <- .myeloma_cta_dat()
  if (is.null(dat)) skip("myeloma fixture missing")
  le <- function(v) !(v %in% -9) & v <= 0.5
  gt <- function(v) !(v %in% -9) & v >  0.5
  all_rows <- rep(TRUE, nrow(dat))
  cases <- list(
    list("T1 N1", all_rows,                         6.93),
    list("T1 N2", le(dat$V17),                     34.89),
    list("T1 N4", le(dat$V17) & le(dat$V15),       34.90),
    list("T1 N5", le(dat$V17) & gt(dat$V15),       33.41),
    list("T1 N3", gt(dat$V17),                      3.42),
    list("T2 N2", le(dat$V14),                     11.40),
    list("T2 N4", le(dat$V14) & le(dat$V15),        1.70),
    list("T2 N5", le(dat$V14) & gt(dat$V15),       28.09),
    list("T2 N3", gt(dat$V14),                     20.22),
    list("T3 N2", le(dat$V15),                      8.34),
    list("T3 N3", gt(dat$V15),                     28.09)
  )
  for (cs in cases) {
    idx <- which(cs[[2]])
    res <- oda:::.cta_ordered_scan(dat$V4[idx], .recode01(dat$V1[idx]), dat$V2[idx],
                                   priors_on = TRUE, miss_codes = -9, mindenom = 1L)
    expect_false(is.null(res), label = cs[[1]])
    expect_equal(round(res$ess, 2), cs[[3]], tolerance = 0.006, label = cs[[1]])
  }
})

test_that(".cta_ordered_scan: MINDENOM gates the selected cut, it does not move it", {
  dat <- .myeloma_cta_dat()
  if (is.null(dat)) skip("myeloma fixture missing")
  idx <- .node2_rows(dat)                       # selected cut 371.20: nL = 102, nR = 29
  x <- dat$V4[idx]; y <- .recode01(dat$V1[idx]); w <- dat$V2[idx]
  r1  <- oda:::.cta_ordered_scan(x, y, w, TRUE, -9, 1L)
  r29 <- oda:::.cta_ordered_scan(x, y, w, TRUE, -9, 29L)
  r30 <- oda:::.cta_ordered_scan(x, y, w, TRUE, -9, 30L)
  expect_equal(r29$rule$cut_value, r1$rule$cut_value)
  expect_equal(r29$ess, r1$ess)
  expect_null(r30)
})

# =============================================================================
# Invariants (CRAN-safe, synthetic)
# =============================================================================

.scan_synth <- function(seed, n = 60L) {
  set.seed(seed)
  list(x = round(stats::rnorm(n), 1), y = stats::rbinom(n, 1L, 0.4),
       w = stats::runif(n, 0.2, 5))
}

test_that(".cta_ordered_scan: weight scale invariance and ESS within [-100, 100]", {
  for (s in 1:20) {
    d  <- .scan_synth(s)
    r  <- oda:::.cta_ordered_scan(d$x, d$y, d$w,      TRUE, NULL, 1L)
    rs <- oda:::.cta_ordered_scan(d$x, d$y, d$w * 10, TRUE, NULL, 1L)
    if (is.null(r)) { expect_null(rs); next }
    expect_equal(rs$rule$cut_value, r$rule$cut_value)
    expect_equal(rs$rule$direction, r$rule$direction)
    expect_equal(rs$ess, r$ess, tolerance = 1e-10)
    expect_true(r$ess >= -100 && r$ess <= 100)
    rf <- oda:::.cta_ordered_scan(d$x, d$y, d$w * 10, FALSE, NULL, 1L)
    if (!is.null(rf)) expect_true(rf$ess >= -100 && rf$ess <= 100)
  }
})

test_that(".cta_ordered_scan: uniform weights give the maximum admissible ESS", {
  for (s in 1:20) {
    d <- .scan_synth(s); w1 <- rep(1, length(d$x))
    r <- oda:::.cta_ordered_scan(d$x, d$y, w1, TRUE, NULL, 1L)
    u <- sort(unique(d$x)); best <- -Inf
    for (j in seq_len(length(u) - 1L)) {
      L <- d$x <= (u[j] + u[j + 1L]) / 2
      e <- (mean(!L[d$y == 1]) + mean(L[d$y == 0]) - 1) * 100
      if (any(L & d$y == 0) && any(!L & d$y == 1)) best <- max(best,  e)
      if (any(L & d$y == 1) && any(!L & d$y == 0)) best <- max(best, -e)
    }
    expect_equal(r$ess, best, tolerance = 1e-10)
  }
})

test_that(".cta_fast_scan_perm: identical to .cta_ordered_scan with no gate", {
  for (s in 1:20) {
    d  <- .scan_synth(s)
    pc <- oda:::.cta_mc_precomp_ord(d$x, d$w, NULL, 1L)
    for (pri in c(TRUE, FALSE)) {
      expect_identical(oda:::.cta_fast_scan_perm(d$y, pc, pri),
                       oda:::.cta_ordered_scan(d$x, d$y, d$w, pri, NULL, 1L))
    }
  }
})

test_that(".cta_fast_scan_perm: MINDENOM gate is not applied to permutations", {
  d  <- .scan_synth(3)
  pc <- oda:::.cta_mc_precomp_ord(d$x, d$w, NULL, 1000L)
  expect_false(is.null(oda:::.cta_fast_scan_perm(d$y, pc, TRUE)))
  expect_null(oda:::.cta_ordered_scan(d$x, d$y, d$w, TRUE, NULL, 1000L))
})

# =============================================================================
# Monte Carlo: attribute permutation with the count-priors objective statistic
# =============================================================================

.mc_ord <- function(x, y, w, seed, iter = 400L, priors_on = TRUE, miss = NULL)
  oda:::.cta_mc_ordered(x, y, w, obs_ess = NA_real_, priors_on = priors_on,
                        miss_codes = miss, mindenom = 1L, mc_iter = iter,
                        mc_target = 0.05, mc_stop = NA, mc_stopup = NA, mc_seed = seed)

test_that(".cta_scan_core: obj is the maximum admissible count-priors objective", {
  for (s in 1:10) {
    d <- .scan_synth(s)
    r <- oda:::.cta_ordered_scan(d$x, d$y, d$w, TRUE, NULL, 1L)
    n0 <- sum(d$y == 0); n1 <- sum(d$y == 1); u <- sort(unique(d$x)); best <- -Inf
    for (j in seq_len(length(u) - 1L)) {
      L <- d$x <= (u[j] + u[j + 1L]) / 2
      if (any(L & d$y == 0) && any(!L & d$y == 1))
        best <- max(best, sum(d$w[L & d$y == 0]) / n0 + sum(d$w[!L & d$y == 1]) / n1)
      if (any(L & d$y == 1) && any(!L & d$y == 0))
        best <- max(best, sum(d$w[L & d$y == 1]) / n1 + sum(d$w[!L & d$y == 0]) / n0)
    }
    expect_equal(r$obj, best, tolerance = 1e-12)
  }
})

test_that(".cta_mc_ordered: p-value is invariant to weight scale", {
  for (s in 1:5) {
    d  <- .scan_synth(s)
    p1 <- .mc_ord(d$x, d$y, d$w,      seed = 11L)
    p2 <- .mc_ord(d$x, d$y, d$w * 37, seed = 11L)
    expect_identical(p1$ge_count, p2$ge_count)
    expect_identical(p1$p_mc, p2$p_mc)
  }
})

test_that(".cta_mc_ordered: missing attribute rows are excluded from permutation", {
  d  <- .scan_synth(4)
  x9 <- c(d$x, rep(-9, 15)); y9 <- c(d$y, rep(1L, 15)); w9 <- c(d$w, rep(50, 15))
  expect_identical(.mc_ord(x9, y9, w9, seed = 5L, miss = -9)$p_mc,
                   .mc_ord(d$x, d$y, d$w, seed = 5L)$p_mc)
})

test_that(".cta_mc_ordered: reproduces CTA.exe Signif at the myeloma root (V15 T, V16 F, V17 T)", {
  skip_if_slow_tests_disabled("cta-mc-permutation")
  dat <- .myeloma_cta_dat()
  if (is.null(dat)) skip("myeloma fixture missing")
  p <- function(v) .mc_ord(dat[[v]], .recode01(dat$V1), dat$V2, seed = 2026L,
                           iter = 5000L, miss = -9)$p_mc
  expect_lt(p("V15"), 0.05)
  expect_gte(p("V16"), 0.05)   # WESS 13.57%: CTA.exe Signif F
  expect_lt(p("V17"), 0.05)    # WESS 16.51%: CTA.exe Signif T
})
