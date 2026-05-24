###############################################################################
# test-cta-ordered-scan.R
#
# Audit-evidence tests for .cta_ordered_scan().
#
# Scope: unit tests ONLY for the scan function itself.
#   The rule (rightmost sens_wa > 0.5) is empirically confirmed only for
#   myeloma V4 at Node 4 and Node 2 (direction 0→1, WEIGHT V2 active).
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
# .cta_ordered_scan() — myeloma-specific audit evidence
# =============================================================================

test_that(".cta_ordered_scan: myeloma Node 4 V4 returns cut 359.80 ESS ~34.90%", {
  skip_if_not_smoke("cta-ordered-scan")
  dat <- .myeloma_cta_dat()
  if (is.null(dat)) skip("myeloma fixture missing")

  idx <- .node4_rows(dat)
  x   <- dat$V4[idx]
  y01 <- .recode01(dat$V1[idx])
  w   <- dat$V2[idx]

  res <- odacore:::.cta_ordered_scan(x, y01, w,
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
  skip_if_not_smoke("cta-ordered-scan")
  dat <- .myeloma_cta_dat()
  if (is.null(dat)) skip("myeloma fixture missing")

  idx <- .node2_rows(dat)
  x   <- dat$V4[idx]
  y01 <- .recode01(dat$V1[idx])
  w   <- dat$V2[idx]

  res <- odacore:::.cta_ordered_scan(x, y01, w,
                                     priors_on  = TRUE,
                                     miss_codes = -9,
                                     mindenom   = 1L)

  expect_false(is.null(res),              label = "result is not NULL")
  expect_equal(res$rule$cut_value, 371.20, tolerance = 0.01,
               label = "cut value — CTA.exe myeloma Node 2 V4 audit")
  expect_equal(res$ess,            34.89,  tolerance = 0.01,
               label = "WESS matches CTA.exe myeloma fixture")
  expect_gt(res$sens_wa, 0.5, label = "sens_wa > 0.5")
})

# =============================================================================
# Generic oda_univariate_core() is unchanged
# =============================================================================

test_that("generic oda_univariate_core: myeloma Node 4 V4 still finds cut ~85.70", {
  skip_if_not_smoke("cta-ordered-scan")
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
  skip_if_not_smoke("cta-ordered-scan")
  x <- 1:6
  y <- rep(0L, 6)
  w <- rep(1, 6)
  expect_null(odacore:::.cta_ordered_scan(x, y, w, TRUE, NULL, 1L))
})

test_that(".cta_ordered_scan: returns NULL when all obs missing", {
  skip_if_not_smoke("cta-ordered-scan")
  x <- rep(-9, 6)
  y <- c(0L, 0L, 0L, 1L, 1L, 1L)
  w <- rep(1, 6)
  expect_null(odacore:::.cta_ordered_scan(x, y, w, TRUE, -9, 1L))
})

test_that(".cta_ordered_scan: does not error under varying mindenom", {
  skip_if_not_smoke("cta-ordered-scan")
  set.seed(1)
  x <- 1:8; y <- c(0L,0L,0L,0L,1L,1L,1L,1L); w <- rep(1, 8)
  r1 <- odacore:::.cta_ordered_scan(x, y, w, TRUE, NULL, 1L)
  r5 <- odacore:::.cta_ordered_scan(x, y, w, TRUE, NULL, 5L)
  expect_false(is.null(r1), label = "mindenom=1 finds a cut")
  expect_true(is.null(r5) || is.list(r5), label = "mindenom=5 returns NULL or list")
})
