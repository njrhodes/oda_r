###############################################################################
# test-cta-ordered-scan.R
#
# Audit-evidence tests for .cta_ordered_scan().
#
# Scope: unit tests ONLY for the scan function itself.
#   The rule (rightmost sens_wa > 0.5) is empirically confirmed only for
#   myeloma V4 at Node 4 and Node 2 (direction 0->1, WEIGHT V2 active).
#   It has NOT been verified against other fixtures (e.g. cta_demo).
#   Do NOT integrate .cta_ordered_scan() into .fast_screen() or
#   .full_fit_one() until confirmed against all canonical CTA fixtures.
#
# Reference: docs/CTA_ORDERED_CUT_AUDIT.md
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
  expect_gt(res$sens_wa, 0.5,
            label = "sens_wa > 0.5: rightmost eligible cut confirmed")
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
# Generic oda_univariate_core() is unchanged
# =============================================================================

test_that("generic oda_univariate_core: myeloma Node 4 V4 still finds cut ~85.70", {
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
  expect_equal(fit$rule$cut_value, 85.70, tolerance = 0.01,
               label = "generic ODA cut unchanged at 85.70")
  expect_equal(fit$ess, 46.99, tolerance = 0.05,
               label = "generic ODA WESS unchanged near 46.99%")
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

test_that("weighted CTA scans are scale-invariant when priors are off", {
  x <- seq_len(8)
  y <- c(rep(0L, 4), rep(1L, 4))
  w <- c(0.3, 0.7, 1.2, 1.8, 0.2, 0.8, 1.4, 1.6)

  reference <- oda:::.cta_ordered_scan(
    x, y, w, priors_on = FALSE, miss_codes = NULL, mindenom = 1L
  )
  reference_mc <- oda:::.cta_mc_ordered(
    x = x, y_coded = y, w = w, obs_ess = reference$ess,
    priors_on = FALSE, miss_codes = NULL, mindenom = 1L,
    mc_iter = 200L, mc_target = NA_real_, mc_stop = NA_real_,
    mc_stopup = NA_real_, mc_seed = 42L
  )
  expect_gte(reference$ess, -100)
  expect_lte(reference$ess, 100)

  for (scale in c(0.1, 10)) {
    scaled <- oda:::.cta_ordered_scan(
      x, y, w * scale, priors_on = FALSE, miss_codes = NULL, mindenom = 1L
    )
    expect_equal(scaled$rule, reference$rule)
    expect_equal(scaled$ess, reference$ess)
    expect_gte(scaled$ess, -100)
    expect_lte(scaled$ess, 100)

    mc <- oda:::.cta_mc_ordered(
      x = x, y_coded = y, w = w * scale, obs_ess = scaled$ess,
      priors_on = FALSE, miss_codes = NULL, mindenom = 1L,
      mc_iter = 200L, mc_target = NA_real_, mc_stop = NA_real_,
      mc_stopup = NA_real_, mc_seed = 42L
    )
    expect_equal(mc$p_mc, reference_mc$p_mc)
  }
})

test_that("uniform case weights retain the unweighted CTA ESS when priors are off", {
  x <- seq_len(8)
  y <- c(rep(0L, 4), rep(1L, 4))

  uniform <- oda:::.cta_ordered_scan(
    x, y, rep(1, length(y)), priors_on = FALSE,
    miss_codes = NULL, mindenom = 1L
  )
  unweighted <- oda:::.cta_ordered_scan(
    x, y, rep(1, length(y)), priors_on = TRUE,
    miss_codes = NULL, mindenom = 1L
  )

  expect_equal(uniform$rule, unweighted$rule)
  expect_equal(uniform$ess, unweighted$ess)
  expect_gte(uniform$ess, -100)
  expect_lte(uniform$ess, 100)
})
