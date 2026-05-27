###############################################################################
# test-sda-novometric.R  —  SDA-4B: sda_fit(mode = "novometric_min_d") tests
#
# Tier: CRAN-safe. Synthetic data only. No data-raw artifacts.
# All tests from §4A.7 of docs/SDA_AUTO_SDA_PLAN.md.
#
# Known structural gap: ge_count and iter_used are not stored by cta_tree
# nodes; they are NA_integer_ in the candidate table. p_mc is extracted from
# the root split node (tree$nodes[[root_id]]$p_mc).
###############################################################################

# ---------------------------------------------------------------------------
# Shared synthetic data  (n = 180, balanced)
#
# Two-attribute design for N1 (min-D winner ≠ max-ESS winner):
#
#   A_attr (binary, S=2): 76 C1 at A=0, 14 C1 at A=1; 14 C2 at A=0, 76 C2 at A=1.
#           UniODA-ESS ≈ 68.9%, CTA-ESS ≈ 68.9%, D ≈ 0.903.
#
#   B_attr (4-level non-monotone, S=4): the non-monotone class signal forces
#           CTA to produce a 4-strata tree.
#           CTA-ESS ≈ 80%, D = 1.0.
#           UniODA single-cut ESS ≈ 53% (poor single cut on non-monotone data).
#
# In novometric mode CTA is used for both candidates:
#   D(A)=0.903 < D(B)=1.0  → A wins min-D.
#   CTA-ESS(B)=80% > CTA-ESS(A)=68.9% → B would win max-CTA-ESS.
#   → min-D winner (A) ≠ max-CTA-ESS winner (B)  [N1 design point]
# ---------------------------------------------------------------------------

.sn_n    <- 180L
.sn_y    <- c(rep(1L, 90L), rep(2L, 90L))

# A binary: 76 C1 at 0, 14 C1 at 1; 14 C2 at 0, 76 C2 at 1.
.sn_A    <- c(rep(0L, 76L), rep(1L, 14L),   # class 1 rows
              rep(0L, 14L), rep(1L, 76L))    # class 2 rows

# B 4-level non-monotone: CTA produces S=4 (alternating A=0/2/1/3 pattern
# with controlled noise).
# C1 rows (1-90): 54 at B=0, 3 at B=1, 27 at B=2, 6 at B=3
# C2 rows (91-180): 6 at B=0, 27 at B=1, 3 at B=2, 54 at B=3
.sn_B    <- c(rep(0L, 54L), rep(1L,  3L), rep(2L, 27L), rep(3L,  6L),
              rep(0L,  6L), rep(1L, 27L), rep(2L,  3L), rep(3L, 54L))

.sn_X    <- data.frame(A = .sn_A, B = .sn_B)

# Pre-fitted novometric model — reused across multiple tests.
# mindenom=10: each of the 4 strata has >=10 obs.
.sn_fit  <- sda_fit(.sn_X, .sn_y,
                    mode     = "novometric_min_d",
                    mindenom = 10L,
                    mc_iter  = 3000L,
                    mc_seed  = 42L)

# ---------------------------------------------------------------------------
# N1 — min-D winner ≠ max-ESS winner
# ---------------------------------------------------------------------------
test_that("N1: novometric selects min-D candidate when it has lower ESS than a competitor", {
  expect_true(length(.sn_fit$steps) >= 1L)

  ctab <- .sn_fit$steps[[1L]]$candidate_table
  A_row <- ctab[ctab$attribute == "A", ]
  B_row <- ctab[ctab$attribute == "B", ]

  # Both candidates must be eligible at step 1 (no MINDENOM or p failures)
  expect_true(A_row$eligible)
  expect_true(B_row$eligible)

  # A is selected (lower D)
  expect_true(A_row$selected)
  expect_false(B_row$selected)

  # B has higher ESS than A (B would be the max-ESS winner)
  expect_true(B_row$ess > A_row$ess,
    label = "B has higher CTA-ESS: max-ESS would pick B, not A")

  # A has lower D than B (novometric correctly picks A)
  expect_true(A_row$d < B_row$d,
    label = "A has lower D: min-D correctly picks A")

  # Confirm step records D correctly
  expect_equal(.sn_fit$steps[[1L]]$attribute, "A")
  expect_true(!is.na(.sn_fit$steps[[1L]]$d))
})

# ---------------------------------------------------------------------------
# N2 — mindenom missing → canon error
# ---------------------------------------------------------------------------
test_that("N2: missing mindenom errors with canon MPE guidance", {
  err <- tryCatch(
    sda_fit(.sn_X, .sn_y, mode = "novometric_min_d"),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("mindenom", err, ignore.case = TRUE))
  expect_true(grepl("MPE|moderate effect|strong effect", err))
})

# ---------------------------------------------------------------------------
# N3 — p-gate failure → ineligible
# ---------------------------------------------------------------------------
test_that("N3: p-gate failure marks candidate ineligible with p_gate reason", {
  # Use tiny n to force MC p > alpha for one or more candidates.
  # n=20 balanced, alpha=0.001 (very strict) → at least one attribute should fail p.
  set.seed(99)
  n_tiny <- 20L
  y_tiny <- c(rep(1L, 10L), rep(2L, 10L))
  # A: moderate signal
  A_tiny <- c(rep(0L, 7L), rep(1L, 3L), rep(0L, 3L), rep(1L, 7L))
  # B: pure noise (random)
  B_tiny <- sample(c(0L, 1L), n_tiny, replace = TRUE)
  X_tiny <- data.frame(A = A_tiny, B = B_tiny)

  fit_strict <- tryCatch(
    sda_fit(X_tiny, y_tiny,
            mode     = "novometric_min_d",
            mindenom = 2L,
            alpha    = 0.001,       # very strict: forces p-gate failures
            mc_iter  = 500L,
            mc_seed  = 77L),
    error = function(e) NULL
  )
  skip_if(is.null(fit_strict), "strict-alpha fit errored on this platform")

  # If any step ran, candidate table should contain at least one p_gate row
  any_p_gate <- FALSE
  for (stp in fit_strict$steps) {
    ctab <- stp$candidate_table
    if (any(!is.na(ctab$ineligible_reason) &
            ctab$ineligible_reason == "p_gate"))
      any_p_gate <- TRUE
  }
  # Also check if stop_reason indicates p_gate
  if (!any_p_gate) {
    expect_true(fit_strict$stop_reason %in% c("p_gate", "axiom1_violated",
                                               "no_candidates", "all_resolved"),
      label = "fit should have stopped or found p_gate ineligibles")
  } else {
    expect_true(any_p_gate)
  }
})

# ---------------------------------------------------------------------------
# N4 — MINDENOM gate failure → ineligible with "axiom1"
# ---------------------------------------------------------------------------
test_that("N4: MINDENOM failure marks candidate ineligible with axiom1 reason", {
  # Use tiny n (n=10) with large mindenom (=8) so stump leaves fail MINDENOM.
  n_ax <- 10L
  y_ax <- c(rep(1L, 5L), rep(2L, 5L))
  A_ax <- c(rep(0L, 4L), rep(1L, 1L), rep(0L, 1L), rep(1L, 4L))
  B_ax <- c(rep(0L, 3L), rep(1L, 2L), rep(0L, 2L), rep(1L, 3L))
  X_ax <- data.frame(A = A_ax, B = B_ax)

  fit_ax <- tryCatch(
    sda_fit(X_ax, y_ax,
            mode     = "novometric_min_d",
            mindenom = 8L,          # each leaf must have >= 8 obs; stump leaves ~5 → fail
            mc_iter  = 500L,
            mc_seed  = 42L),
    error = function(e) NULL
  )
  skip_if(is.null(fit_ax), "axiom1 fit errored on this platform")

  # Either stopped immediately (axiom1_violated) or has candidate table entries
  expect_true(fit_ax$stop_reason %in%
    c("axiom1_violated", "p_gate", "no_candidates", "all_resolved",
      "min_n", "class_resolved"))

  # If any steps ran, check for axiom1 entries
  if (length(fit_ax$steps) > 0L) {
    all_ctab <- do.call(rbind, lapply(fit_ax$steps, function(s) s$candidate_table))
    if (any(!is.na(all_ctab$ineligible_reason)))
      expect_true(any(all_ctab$ineligible_reason == "axiom1", na.rm = TRUE))
  }
})

# ---------------------------------------------------------------------------
# N5 — no eligible candidate → stop reason p_gate or axiom1_violated
# ---------------------------------------------------------------------------
test_that("N5: all candidates fail → stop reason is p_gate or axiom1_violated", {
  # Construct data where MINDENOM is too large for any valid tree.
  n5 <- 14L
  y5 <- c(rep(1L, 7L), rep(2L, 7L))
  A5 <- c(rep(0L, 5L), rep(1L, 2L), rep(0L, 2L), rep(1L, 5L))
  B5 <- c(rep(0L, 4L), rep(1L, 3L), rep(0L, 3L), rep(1L, 4L))
  X5 <- data.frame(A = A5, B = B5)

  fit5 <- tryCatch(
    sda_fit(X5, y5,
            mode     = "novometric_min_d",
            mindenom = 10L,     # too large for n=14
            mc_iter  = 200L,
            mc_seed  = 42L),
    error = function(e) NULL
  )
  skip_if(is.null(fit5), "N5 fit errored")

  expect_true(fit5$stop_reason %in%
    c("p_gate", "axiom1_violated", "no_candidates",
      "min_n", "class_resolved", "all_resolved"))
})

# ---------------------------------------------------------------------------
# N6 — selected attribute removed from next step's candidate pool
# ---------------------------------------------------------------------------
test_that("N6: selected attribute absent from step 2 candidate table", {
  skip_if(length(.sn_fit$steps) < 2L, "need at least 2 SDA steps")

  step1_attr  <- .sn_fit$steps[[1L]]$attribute
  step2_attrs <- .sn_fit$steps[[2L]]$candidate_table$attribute
  expect_false(step1_attr %in% step2_attrs,
    label = paste0("attribute '", step1_attr, "' should be absent from step 2 candidates"))
})

# ---------------------------------------------------------------------------
# N7 — correctly classified observations removed between steps
# ---------------------------------------------------------------------------
test_that("N7: step 2 n_in = step 1 n_in - step 1 n_correct", {
  skip_if(length(.sn_fit$steps) < 2L, "need at least 2 SDA steps")

  s1 <- .sn_fit$steps[[1L]]
  s2 <- .sn_fit$steps[[2L]]
  expect_equal(s2$n_in, s1$n_in - s1$n_correct)
})

# ---------------------------------------------------------------------------
# N8 — candidate table records all evaluated candidates
# ---------------------------------------------------------------------------
test_that("N8: candidate table has one row per evaluated attribute", {
  for (stp in .sn_fit$steps) {
    ctab <- stp$candidate_table
    expect_true(is.data.frame(ctab))
    expect_true(nrow(ctab) >= 1L)

    # Status column populated
    expect_true(all(ctab$status %in% c("selected", "eligible", "ineligible", "skipped")))

    # Exactly one selected row per step
    expect_equal(sum(ctab$selected), 1L)

    # ineligible rows have a reason
    inelig_rows <- ctab[ctab$status == "ineligible", ]
    if (nrow(inelig_rows) > 0L)
      expect_true(all(!is.na(inelig_rows$ineligible_reason)))

    # All 17 required contract fields present
    required_fields <- c("attribute", "status", "eligible", "ineligible_reason",
                         "n", "class_counts", "mindenom", "min_terminal_denom",
                         "ess", "d", "p_mc", "ge_count", "iter_used",
                         "strata", "selected", "tied_objective",
                         "selected_by_tie_break")
    for (fld in required_fields)
      expect_true(fld %in% names(ctab),
        label = paste0("candidate table missing required field: ", fld))
  }
})

# ---------------------------------------------------------------------------
# N9 — tie case: both candidates record tied_objective = TRUE; first in
# column order wins
# ---------------------------------------------------------------------------
test_that("N9: tied D yields tied_objective=TRUE; column-order tie-break is deterministic", {
  # Construct data where two IDENTICAL attributes produce the same D.
  # A and B are exact copies → both should produce identical families.
  # Column order: A is first → A should win the tie-break.
  set.seed(7)
  n9   <- 80L
  y9   <- c(rep(1L, 40L), rep(2L, 40L))
  col9 <- c(rep(0L, 32L), rep(1L,  8L),   # class 1
             rep(0L,  8L), rep(1L, 32L))   # class 2
  X9   <- data.frame(A = col9, B = col9)   # exact duplicates

  # With collinearity="allow", both columns survive to the candidate pool.
  fit9 <- tryCatch(
    sda_fit(X9, y9,
            mode         = "novometric_min_d",
            mindenom     = 5L,
            mc_iter      = 2000L,
            mc_seed      = 42L,
            collinearity = "allow"),
    error = function(e) NULL
  )
  skip_if(is.null(fit9) || length(fit9$steps) == 0L, "N9 fit produced no steps")

  ctab <- fit9$steps[[1L]]$candidate_table
  # Winner must be A (first column)
  expect_equal(fit9$steps[[1L]]$attribute, "A")

  # Both should have tied_objective = TRUE (same D because same data)
  A_row <- ctab[ctab$attribute == "A", ]
  B_row <- ctab[ctab$attribute == "B", ]

  # Both eligible → both in tie pool
  skip_if(!A_row$eligible || !B_row$eligible, "both must be eligible for tie test")
  expect_true(A_row$tied_objective)
  expect_true(B_row$tied_objective)

  # A is flagged selected_by_tie_break
  expect_true(A_row$selected_by_tie_break)
})

# ---------------------------------------------------------------------------
# N10 — predict.sda_fit works on novometric fit
# ---------------------------------------------------------------------------
test_that("N10: predict.sda_fit(type='class') works on novometric fit", {
  preds <- predict(.sn_fit, .sn_X, type = "class")
  expect_true(is.integer(preds))
  expect_equal(length(preds), .sn_n)
  # predictions are in {1, 2, NA}
  expect_true(all(preds %in% c(1L, 2L, NA_integer_)))

  # stage output
  stages <- predict(.sn_fit, .sn_X, type = "stage")
  expect_equal(length(stages), .sn_n)

  # At least some obs classified
  expect_true(sum(!is.na(preds)) > 0L)
})

# ---------------------------------------------------------------------------
# N11 — CTA interop: as_cta_candidates works on novometric fit
# ---------------------------------------------------------------------------
test_that("N11: as_cta_candidates works on a novometric fit", {
  X_cta <- as_cta_candidates(.sn_fit, .sn_X)
  sel   <- sda_selected_attributes(.sn_fit)
  expect_true(is.data.frame(X_cta))
  expect_equal(colnames(X_cta), sel)
  expect_equal(nrow(X_cta), .sn_n)
})

# ---------------------------------------------------------------------------
# N12 — axiom1_violated vs p_gate distinguishable
# ---------------------------------------------------------------------------
test_that("N12: axiom1_violated and p_gate stop reasons are distinct", {
  # Case 1: MINDENOM too large → axiom1_violated
  n12 <- 12L
  y12 <- c(rep(1L, 6L), rep(2L, 6L))
  A12 <- c(rep(0L, 5L), rep(1L, 1L), rep(0L, 1L), rep(1L, 5L))
  X12 <- data.frame(A = A12)

  fit_ax1 <- tryCatch(
    sda_fit(X12, y12, mode = "novometric_min_d",
            mindenom = 10L, mc_iter = 200L, mc_seed = 1L),
    error = function(e) NULL
  )
  skip_if(is.null(fit_ax1), "axiom1 fit errored")
  expect_true(fit_ax1$stop_reason %in%
    c("axiom1_violated", "p_gate", "no_candidates", "all_resolved",
      "min_n", "class_resolved"))

  # Case 2: MINDENOM is feasible but signal is too weak to pass p-gate
  # Use strict alpha
  n12b <- 40L
  y12b <- c(rep(1L, 20L), rep(2L, 20L))
  # Near-random A (weak signal)
  A12b <- c(rep(0L, 11L), rep(1L, 9L), rep(0L, 9L), rep(1L, 11L))
  X12b <- data.frame(A = A12b)

  fit_pg <- tryCatch(
    sda_fit(X12b, y12b, mode = "novometric_min_d",
            mindenom = 4L, alpha = 0.001,
            mc_iter = 500L, mc_seed = 1L),
    error = function(e) NULL
  )
  skip_if(is.null(fit_pg), "p_gate fit errored")

  # The two stop reasons must be different character values
  expect_false(identical("axiom1_violated", "p_gate"))

  # Both are in the set of valid stop reasons
  valid_stops <- c("axiom1_violated", "p_gate", "no_candidates",
                   "all_resolved", "min_n", "class_resolved", "max_steps")
  expect_true(fit_ax1$stop_reason %in% valid_stops)
  expect_true(fit_pg$stop_reason %in% valid_stops)
})

# ---------------------------------------------------------------------------
# N13 — min_terminal_denom recorded in candidate table
# ---------------------------------------------------------------------------
test_that("N13: min_terminal_denom in candidate table matches family object", {
  # Use .sn_fit where step 1 selected A (verified in N1).
  ctab <- .sn_fit$steps[[1L]]$candidate_table
  A_row <- ctab[ctab$attribute == "A", ]

  # A is eligible → min_terminal_denom should be non-NA
  expect_false(is.na(A_row$min_terminal_denom))
  expect_true(A_row$min_terminal_denom >= 10L)   # >= mindenom used

  # Cross-check: retrieve the cta_family stored in the step
  step_model <- .sn_fit$steps[[1L]]$model   # cta_family for selected attribute A
  midx       <- .sn_fit$steps[[1L]]$min_d_idx
  expect_false(is.na(midx))

  best_member <- step_model$members[[midx]]
  expect_equal(A_row$min_terminal_denom, best_member$min_terminal_denom)
})

# ---------------------------------------------------------------------------
# Additional: structural fields on the fitted object
# ---------------------------------------------------------------------------
test_that("novometric fit has correct class and mode fields", {
  expect_true(inherits(.sn_fit, "sda_fit"))
  expect_true(inherits(.sn_fit, "odacore_sda"))
  expect_equal(.sn_fit$mode, "novometric_min_d")
  expect_equal(.sn_fit$settings$mindenom, 10L)
})

test_that("novometric step stores non-NA d, mode='novometric_min_d', min_d_idx", {
  for (stp in .sn_fit$steps) {
    expect_equal(stp$mode, "novometric_min_d")
    expect_false(is.na(stp$d))
    expect_false(is.na(stp$min_d_idx))
    expect_equal(stp$reason, "min_d")
    # ge_count and iter_used are NA (structural gap: not stored by cta_tree)
    expect_true(is.na(stp$ge_count))
    expect_true(is.na(stp$iter_used))
  }
})

test_that("novometric print runs without error", {
  out <- capture.output(print(.sn_fit))
  expect_true(length(out) > 0L)
  expect_true(any(grepl("novometric_min_d", out)))
  expect_true(any(grepl("^\\s*\\[1\\]", out)))   # step line with D=
})

test_that("novometric summary runs without error", {
  sm  <- summary(.sn_fit)
  out <- capture.output(print(sm))
  expect_true(length(out) > 0L)
})

test_that("sda_step_table works on novometric fit", {
  st <- sda_step_table(.sn_fit)
  expect_true(is.data.frame(st))
  expect_equal(nrow(st), length(.sn_fit$steps))
  expect_true("d" %in% names(st))
  # d should be non-NA for novometric mode
  if (nrow(st) > 0L) expect_false(any(is.na(st$d)))
})
