###############################################################################
# test-oda-fit.R
# Tests for the unified oda_fit() entry point.
###############################################################################

test_that("oda_fit: routes to binary engine for C=2", {
  x <- c(1,2,3,4,5,6,7,8)
  y <- c(0L,0L,0L,0L,1L,1L,1L,1L)

  fit <- oda_fit(x=x, y=y, mcarlo=FALSE, loo="off")
  expect_true(fit$ok)
  expect_equal(fit$engine, "binary")
  expect_true(!is.null(fit$rule$cut_value) || !is.null(fit$rule$left_levels))
})

test_that("oda_fit: routes to multiclass engine for C=3", {
  x <- c(1,2,3,4,5,6,7,8,9)
  y <- c(1L,1L,1L,2L,2L,2L,3L,3L,3L)

  fit <- oda_fit(x=x, y=y, mcarlo=FALSE, loo="off")
  expect_true(fit$ok)
  expect_equal(fit$engine, "multiclass")
  expect_true(!is.null(fit$rule$cut_values))
})

test_that("oda_fit: missing_code alias works", {
  x <- c(1,2,3,4,5,-99)
  y <- c(0L,0L,0L,1L,1L,1L)

  fit_no_miss  <- oda_fit(x=x, y=y, mcarlo=FALSE)
  fit_with_miss <- oda_fit(x=x, y=y, missing_code=-99, mcarlo=FALSE)

  expect_true(fit_no_miss$ok)
  expect_true(fit_with_miss$ok)
  # With -99 excluded, 5 obs; cut should be in [1,5] not near -99
  expect_true(fit_with_miss$n_eff == 5L)
  expect_true(fit_with_miss$rule$cut_value > -10)
})

test_that("oda_fit: miss_codes vector also works", {
  x <- c(1,2,3,4,5,-99,-88)
  y <- c(0L,0L,0L,1L,1L,1L,0L)

  fit <- oda_fit(x=x, y=y, miss_codes=c(-99,-88), mcarlo=FALSE)
  expect_true(fit$ok)
  expect_equal(fit$n_eff, 5L)
})

test_that("oda_fit: pure node returns ok=FALSE", {
  x <- c(1,2,3,4)
  y <- c(0L,0L,0L,0L)

  fit <- oda_fit(x=x, y=y, mcarlo=FALSE)
  expect_false(fit$ok)
})

test_that("oda_fit: all_missing returns ok=FALSE", {
  x <- rep(NA_real_, 5)
  y <- c(0L,1L,0L,1L,0L)

  fit <- oda_fit(x=x, y=y, mcarlo=FALSE)
  expect_false(fit$ok)
})

test_that("oda_fit: binary and multiclass ESS are both finite positive for separated data", {
  # Binary
  x2 <- c(1,2,3,7,8,9); y2 <- c(0L,0L,0L,1L,1L,1L)
  f2 <- oda_fit(x=x2, y=y2, mcarlo=FALSE)
  expect_true(f2$ok)
  expect_true(is.finite(f2$ess) && f2$ess > 0)

  # Multiclass
  x3 <- c(1,2,3,4,5,6,7,8,9); y3 <- c(1L,1L,1L,2L,2L,2L,3L,3L,3L)
  f3 <- oda_fit(x=x3, y=y3, mcarlo=FALSE)
  expect_true(f3$ok)
  # multiclass exposes $ess (ess_pac retained as compat alias)
  expect_true(is.finite(f3$ess) && f3$ess > 0)
})

test_that("oda_fit multiclass: $ess and $ess_pac are both present and equal (Patch 1 transition)", {
  # Patch 1: $ess is the normalized public name; $ess_pac is the compat alias.
  # Both must exist on multiclass results and hold the same computed value
  # until $ess_pac is retired in a future patch.
  x <- c(1,2,3,4,5,6,7,8,9)
  y <- c(1L,1L,1L,2L,2L,2L,3L,3L,3L)
  fit <- oda_fit(x, y, mcarlo = FALSE)
  expect_true(fit$ok)
  expect_true(is.finite(fit$ess)     && fit$ess     > 0)
  expect_true(is.finite(fit$ess_pac) && fit$ess_pac > 0)
  expect_equal(fit$ess, fit$ess_pac)
})
