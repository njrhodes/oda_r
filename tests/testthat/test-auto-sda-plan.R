###############################################################################
# test-auto-sda-plan.R - SDA-3 auto_sda_plan() tests
#
# Tier: CRAN-safe. No fitting; only planning/validation. Synthetic data only.
###############################################################################

# ---------------------------------------------------------------------------
# Synthetic test frame
#
# 40 rows, binary outcome (y), 8 candidate columns covering all edge cases:
#   A: numeric ordered predictor - clean
#   B: integer binary predictor  - clean
#   C: character predictor       - clean (categorical)
#   D: factor predictor          - clean (2 levels -> binary)
#   E: constant (zero variance)  - should be excluded
#   F: all-NA                    - should be excluded
#   G: list column               - should be excluded (invalid type)
#   H: exact duplicate of A      - collinear with A; should be flagged
#   id_col: ID column            - declared via role_map
# ---------------------------------------------------------------------------

set.seed(1)
.asp_n <- 40L
.asp_df <- data.frame(
  y      = c(rep(1L, 20L), rep(2L, 20L)),
  A      = c(rnorm(20, mean = 0), rnorm(20, mean = 2)),
  B      = c(rep(0L, 20L), rep(1L, 20L)),
  C      = rep(c("low", "high"), each = 20L),
  D      = factor(rep(c("no", "yes"), each = 20L)),
  E      = rep(5L, .asp_n),                  # constant
  F_col  = rep(NA_real_, .asp_n),            # all-missing
  id_col = seq_len(.asp_n),                  # ID
  stringsAsFactors = FALSE
)
# Add list column and exact duplicate after data.frame construction
.asp_df$G_list <- as.list(seq_len(.asp_n))   # list column - invalid type
.asp_df$H      <- .asp_df$A                  # exact duplicate of A

# ---------------------------------------------------------------------------
# Test 1: candidates=NULL excludes outcome; all other valid cols are candidates
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: candidates=NULL excludes outcome", {
  plan <- auto_sda_plan(.asp_df, outcome = "y")
  expect_false("y" %in% plan$candidate_names)
  expect_true("is_outcome" %in% plan$exclusion_reasons$reason)
})

# ---------------------------------------------------------------------------
# Test 2: explicit candidates list is respected
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: explicit candidates are respected", {
  plan <- auto_sda_plan(.asp_df, outcome = "y", candidates = c("A", "B"))
  # Only A and B (after any further exclusions) should be in candidate_names
  expect_true(all(plan$candidate_names %in% c("A", "B")))
  # outcome still excluded even if accidentally listed
  plan2 <- auto_sda_plan(.asp_df, outcome = "y", candidates = c("A", "B", "y"))
  expect_false("y" %in% plan2$candidate_names)
  expect_true(length(plan2$warnings) > 0L)  # warning about outcome in candidates
})

# ---------------------------------------------------------------------------
# Test 3: exclude removes variables with reason "force_excluded"
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: exclude removes with force_excluded reason", {
  plan <- auto_sda_plan(.asp_df, outcome = "y", candidates = c("A", "B", "C"),
                        exclude = "C")
  expect_false("C" %in% plan$candidate_names)
  reason_row <- plan$exclusion_reasons[plan$exclusion_reasons$name == "C", ]
  expect_equal(reason_row$reason, "force_excluded")
})

# ---------------------------------------------------------------------------
# Test 4: constant columns are excluded with reason "zero_variance"
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: constant columns excluded as zero_variance", {
  plan <- auto_sda_plan(.asp_df, outcome = "y")
  expect_false("E" %in% plan$candidate_names)
  reason_row <- plan$exclusion_reasons[plan$exclusion_reasons$name == "E", ]
  expect_equal(reason_row$reason, "zero_variance")
})

# ---------------------------------------------------------------------------
# Test 5: all-missing columns are excluded with reason "all_missing"
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: all-missing columns excluded as all_missing", {
  plan <- auto_sda_plan(.asp_df, outcome = "y")
  expect_false("F_col" %in% plan$candidate_names)
  reason_row <- plan$exclusion_reasons[plan$exclusion_reasons$name == "F_col", ]
  expect_equal(reason_row$reason, "all_missing")
})

# ---------------------------------------------------------------------------
# Test 6: list columns are excluded with reason "invalid_type"
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: list columns excluded as invalid_type", {
  plan <- auto_sda_plan(.asp_df, outcome = "y")
  expect_false("G_list" %in% plan$candidate_names)
  reason_row <- plan$exclusion_reasons[plan$exclusion_reasons$name == "G_list", ]
  expect_equal(reason_row$reason, "invalid_type")
})

# ---------------------------------------------------------------------------
# Test 7: exact duplicate columns are flagged and excluded as "collinear"
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: exact duplicate column excluded as collinear", {
  plan <- auto_sda_plan(.asp_df, outcome = "y")
  # H is an exact duplicate of A; one should survive, the other excluded
  a_in <- "A" %in% plan$candidate_names
  h_in <- "H" %in% plan$candidate_names
  # At least one should be excluded
  expect_false(a_in && h_in)  # they should not both be in the candidate set
  # The dropped one should have reason "collinear"
  collinear_rows <- plan$exclusion_reasons[
    plan$exclusion_reasons$reason == "collinear", ]
  expect_true(nrow(collinear_rows) >= 1L)
  # Warning should mention duplicate detection
  expect_true(any(grepl("duplicate", plan$warnings)))
})

# ---------------------------------------------------------------------------
# Test 8: role_map/time_map/stage_map are preserved on the returned object
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: role_map, time_map, stage_map preserved", {
  rm <- list(id_col = "id", B = "assignment_mechanism")
  tm <- c(A = 1, B = 2, y = 3)
  sm <- c(A = 1L, B = 2L)
  plan <- auto_sda_plan(.asp_df, outcome = "y",
                        candidates = c("A", "B"),
                        role_map = rm, time_map = tm, stage_map = sm)
  expect_identical(plan$role_map, rm)
  expect_identical(plan$time_map, tm)
  expect_identical(plan$stage_map, sm)
})

# ---------------------------------------------------------------------------
# Test 9: role_map "id" excludes the column; unknown map names go to warnings
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: role_map id excluded; unknown names -> warnings", {
  plan <- auto_sda_plan(.asp_df, outcome = "y",
                        role_map = list(id_col = "id",
                                        does_not_exist = "leakage"))
  # id_col excluded by role_id
  expect_false("id_col" %in% plan$candidate_names)
  reason_row <- plan$exclusion_reasons[
    plan$exclusion_reasons$name == "id_col", ]
  expect_equal(reason_row$reason, "role_id")
  # unknown name -> warning
  expect_true(any(grepl("does_not_exist", plan$warnings)))
})

test_that("auto_sda_plan: time_map unknown names -> warnings field", {
  plan <- auto_sda_plan(.asp_df, outcome = "y",
                        candidates = c("A", "B"),
                        time_map = c(A = 1, ghost_var = 99))
  expect_true(any(grepl("ghost_var", plan$warnings)))
})

test_that("auto_sda_plan: stage_map unknown names -> warnings field", {
  plan <- auto_sda_plan(.asp_df, outcome = "y",
                        candidates = c("A", "B"),
                        stage_map = c(A = 1L, phantom = 2L))
  expect_true(any(grepl("phantom", plan$warnings)))
})

# ---------------------------------------------------------------------------
# Test 10: dry_run = FALSE errors clearly
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: dry_run=FALSE errors with informative message", {
  expect_error(
    auto_sda_plan(.asp_df, outcome = "y", dry_run = FALSE),
    regexp = "dry_run = TRUE only"
  )
})

# ---------------------------------------------------------------------------
# Test 11: proposed_call contains candidate-names-compatible settings
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: proposed_call contains mode and settings", {
  plan <- auto_sda_plan(.asp_df, outcome = "y",
                        candidates = c("A", "B"),
                        mode = "unioda_max_ess")
  pc <- plan$proposed_call
  expect_true(is.list(pc))
  expect_equal(pc$mode, "unioda_max_ess")
  expect_true(!is.null(pc$mc_iter))
  expect_true(!is.null(pc$mc_seed))
  expect_true(!is.null(pc$alpha))
  # proposed_call should NOT contain X or y (data not stored there)
  expect_false("X" %in% names(pc))
  expect_false("y" %in% names(pc))
})

# ---------------------------------------------------------------------------
# Test 12: returned object has correct class
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: object has class c('auto_sda_plan','odacore_plan')", {
  plan <- auto_sda_plan(.asp_df, outcome = "y", candidates = c("A", "B"))
  expect_true(inherits(plan, "auto_sda_plan"))
  expect_true(inherits(plan, "odacore_plan"))
  expect_equal(class(plan), c("auto_sda_plan", "odacore_plan"))
})

# ---------------------------------------------------------------------------
# Test 13: print method smoke test
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: print runs without error", {
  plan <- auto_sda_plan(.asp_df, outcome = "y", candidates = c("A", "B", "C"))
  out <- capture.output(print(plan))
  expect_true(length(out) > 0L)
  expect_true(any(grepl("auto_sda_plan", out)))
})

# ---------------------------------------------------------------------------
# Test 14: attr_types inference
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: attr_types inferred correctly for clean candidates", {
  plan <- auto_sda_plan(.asp_df, outcome = "y",
                        candidates = c("A", "B", "C", "D"))
  at <- plan$attr_types
  expect_equal(at[["A"]], "ordered")    # numeric, >2 unique values
  expect_equal(at[["B"]], "binary")     # integer, 2 unique values
  expect_equal(at[["C"]], "binary")     # character, 2 unique values
  expect_equal(at[["D"]], "binary")     # factor, 2 levels
})

test_that("auto_sda_plan: declared attr_types override inference", {
  plan <- auto_sda_plan(.asp_df, outcome = "y",
                        candidates = c("A", "B"),
                        attr_types = c(A = "categorical"))
  expect_equal(plan$attr_types[["A"]], "categorical")
  expect_equal(plan$attr_types[["B"]], "binary")  # inferred
})

# ---------------------------------------------------------------------------
# Test 15: time_map leakage flag (warn; do not auto-exclude)
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: post-outcome time_map generates warning, not exclusion", {
  plan <- auto_sda_plan(.asp_df, outcome = "y",
                        candidates = c("A", "B"),
                        time_map = c(y = 1, A = 0, B = 2))  # B after outcome
  expect_true("B" %in% plan$candidate_names)   # NOT excluded
  expect_true(any(grepl("B", plan$warnings) & grepl("leakage", plan$warnings)))
})

# ---------------------------------------------------------------------------
# Test 16: outcome not in data errors clearly
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: outcome not in data errors clearly", {
  expect_error(
    auto_sda_plan(.asp_df, outcome = "not_a_column"),
    regexp = "not found in data"
  )
})

# ---------------------------------------------------------------------------
# Test 17: unknown candidates error (not warning)
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: unknown candidates error", {
  expect_error(
    auto_sda_plan(.asp_df, outcome = "y",
                  candidates = c("A", "does_not_exist")),
    regexp = "candidates not found in data"
  )
})

# ---------------------------------------------------------------------------
# Test 18: class_counts field correct
# ---------------------------------------------------------------------------
test_that("auto_sda_plan: class_counts field matches outcome distribution", {
  plan <- auto_sda_plan(.asp_df, outcome = "y", candidates = c("A", "B"))
  expect_equal(sum(plan$class_counts), .asp_n)
  expect_equal(as.integer(plan$class_counts), c(20L, 20L))
})
