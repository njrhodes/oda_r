###############################################################################
# tests/testthat/helper-test-tier.R
# Test tiering for slow canon / integration tests. See issue #9.
#
# Usage:
#   ODACORE_TEST_TIER=fast  -> slow tests are skipped with an explicit message
#   ODACORE_TEST_TIER=full  -> all tests run (same as unset)
#   (unset)                 -> all tests run (default-safe)
#
# Call skip_if_slow_tests_disabled() at the top of any slow test_that block.
###############################################################################

skip_if_slow_tests_disabled <- function(reason = "slow canon test") {
  tier <- Sys.getenv("ODACORE_TEST_TIER", unset = "full")
  if (identical(tier, "fast")) {
    testthat::skip(paste0("ODACORE_TEST_TIER=fast: ", reason, " skipped"))
  }
}
