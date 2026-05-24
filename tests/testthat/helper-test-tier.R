###############################################################################
# tests/testthat/helper-test-tier.R
# Test tiering for slow canon / integration tests.
#
# Tier levels (ordered low → high):
#   cran  — CRAN-safe unit/contract tests only (default when unset or on CRAN)
#   fast  — same practical behavior as cran for slow fixtures; fast dev loop
#   smoke — CTA/MDSA/reporting/graphics production validation (myeloma, LOO
#           STABLE, translation stack, myeloma graphics)
#   full  — all MegaODA/CTA.exe fixture parity, mc_iter >= 10000, release gate
#
# Environment variable: ODACORE_TEST_TIER
#   unset / "cran" → cran tier (CRAN-safe; default)
#   "fast"         → fast dev loop (skips same slow tests as cran)
#   "smoke"        → smoke + cran; skips full-only tests
#   "full"         → runs all tests
#
# Usage commands:
#   # CRAN-safe / default check (no env var):
#   Rscript --vanilla -e "devtools::check(vignettes=FALSE)"
#
#   # Fast developer loop:
#   ODACORE_TEST_TIER=fast Rscript --vanilla -e "devtools::test(reporter='progress')"
#
#   # CTA/MDSA/reporting/graphics production smoke:
#   ODACORE_TEST_TIER=smoke Rscript --vanilla -e "devtools::test(reporter='progress')"
#
#   # Full canon / release gate:
#   ODACORE_TEST_TIER=full Rscript --vanilla -e "devtools::test(reporter='progress')"
#
# Guards:
#   skip_if_not_smoke(reason) — skip unless tier >= smoke; always skip on CRAN
#   skip_if_not_full(reason)  — skip unless tier == full; always skip on CRAN
#   skip_if_slow_tests_disabled(reason) — backward-compatible alias for
#                                         skip_if_not_smoke(); all existing
#                                         slow guards now skip by default/cran
###############################################################################

.tier_level <- function() {
  # On real CRAN (NOT_CRAN is unset), always cran-safe regardless of env var.
  if (!identical(Sys.getenv("NOT_CRAN"), "true")) return(1L)
  tier <- Sys.getenv("ODACORE_TEST_TIER", unset = "cran")
  switch(tier,
    "cran"  = 1L,
    "fast"  = 2L,
    "smoke" = 3L,
    "full"  = 4L,
    1L  # unknown value → cran-safe
  )
}

odacore_test_tier <- function() {
  if (!identical(Sys.getenv("NOT_CRAN"), "true")) return("cran")
  tier <- Sys.getenv("ODACORE_TEST_TIER", unset = "cran")
  if (tier %in% c("cran", "fast", "smoke", "full")) tier else "cran"
}

is_test_tier <- function(...) {
  odacore_test_tier() %in% c(...)
}

skip_unless_test_tier <- function(..., reason = NULL) {
  tiers <- c(...)
  if (!is_test_tier(tiers)) {
    msg <- if (!is.null(reason))
      paste0("ODACORE_TEST_TIER: ", reason, " skipped")
    else
      paste0("requires ODACORE_TEST_TIER in (", paste(tiers, collapse = "/"), ")")
    testthat::skip(msg)
  }
}

skip_if_not_smoke <- function(reason = NULL) {
  testthat::skip_on_cran()
  if (.tier_level() < 3L) {
    msg <- if (!is.null(reason))
      paste0("ODACORE_TEST_TIER<smoke: ", reason, " skipped")
    else
      "ODACORE_TEST_TIER<smoke: smoke test skipped"
    testthat::skip(msg)
  }
}

skip_if_not_full <- function(reason = NULL) {
  testthat::skip_on_cran()
  if (.tier_level() < 4L) {
    msg <- if (!is.null(reason))
      paste0("ODACORE_TEST_TIER<full: ", reason, " skipped")
    else
      "ODACORE_TEST_TIER<full: full-canon test skipped"
    testthat::skip(msg)
  }
}

# Backward-compatible alias: maps to skip_if_not_smoke().
# All existing skip_if_slow_tests_disabled() calls now require smoke or full tier.
skip_if_slow_tests_disabled <- function(reason = NULL) {
  skip_if_not_smoke(reason)
}
