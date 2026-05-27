###############################################################################
# test-sda-interop.R — SDA-2 CTA interop accessor tests
#
# Tier: CRAN-safe. All fits use small mc_iter on synthetic data.
###############################################################################

# Shared synthetic fit (reused across tests)
.sda_iop_X <- data.frame(
  A = c(rep(0L, 40), rep(1L, 10), rep(1L, 40), rep(0L, 10)),
  B = c(rep(0L, 30), rep(1L, 20), rep(1L, 30), rep(0L, 20)),
  C = c(rep(0L, 50), rep(1L, 50))
)
.sda_iop_y <- c(rep(1L, 50), rep(2L, 50))

.sda_iop_fit <- sda_fit(
  .sda_iop_X, .sda_iop_y,
  mode    = "unioda_max_ess",
  mc_iter = 500L,
  mc_seed = 42L,
  loo     = "off"
)

# ---------------------------------------------------------------------------
# sda_selected_attributes
# ---------------------------------------------------------------------------
test_that("sda_selected_attributes: returns character vector in step order", {
  sel <- sda_selected_attributes(.sda_iop_fit)
  expect_true(is.character(sel))
  expect_equal(length(sel), length(.sda_iop_fit$steps))
  # First selected attribute must match step 1
  if (length(.sda_iop_fit$steps) > 0L)
    expect_equal(sel[1], .sda_iop_fit$steps[[1]]$attribute)
})

test_that("sda_selected_attributes: errors on non-sda_fit", {
  expect_error(sda_selected_attributes(list()), regexp = "sda_fit")
})

# ---------------------------------------------------------------------------
# sda_step_table
# ---------------------------------------------------------------------------
test_that("sda_step_table: returns data.frame with required columns", {
  st <- sda_step_table(.sda_iop_fit)
  expect_true(is.data.frame(st))
  required <- c("step_id", "attribute", "n_in", "n_correct", "n_incorrect",
                "ess", "d", "p_mc", "mindenom")
  for (col in required)
    expect_true(col %in% names(st), label = paste0("sda_step_table missing: ", col))
})

test_that("sda_step_table: one row per step; step_id sequential", {
  st <- sda_step_table(.sda_iop_fit)
  expect_equal(nrow(st), length(.sda_iop_fit$steps))
  expect_equal(st$step_id, seq_len(nrow(st)))
})

test_that("sda_step_table: d is NA for unioda_max_ess mode", {
  st <- sda_step_table(.sda_iop_fit)
  expect_true(all(is.na(st$d)))
})

test_that("sda_step_table: empty data.frame for zero-step fit", {
  # Construct a minimal zero-step fit by using a pure-node dataset
  X_pure <- data.frame(A = rep(0L, 10))
  y_pure  <- rep(1L, 10)
  expect_error(
    sda_fit(X_pure, y_pure, mode = "unioda_max_ess", mc_iter = 100L),
    regexp = "2 distinct"  # validation error: y must have 2 classes
  )
  # Directly test with an empty steps list
  fake_fit <- structure(list(steps = list()), class = c("sda_fit", "odacore_sda"))
  st <- sda_step_table(fake_fit)
  expect_equal(nrow(st), 0L)
  expect_true(is.data.frame(st))
})

# ---------------------------------------------------------------------------
# sda_candidate_table
# ---------------------------------------------------------------------------
test_that("sda_candidate_table: step=NULL returns named list", {
  tabs <- sda_candidate_table(.sda_iop_fit)
  expect_true(is.list(tabs))
  expect_equal(length(tabs), length(.sda_iop_fit$steps))
  expect_true(all(grepl("^step_", names(tabs))))
})

test_that("sda_candidate_table: single step returns data.frame with step_id col", {
  if (length(.sda_iop_fit$steps) == 0L) skip("no steps")
  ctab <- sda_candidate_table(.sda_iop_fit, step = 1L)
  expect_true(is.data.frame(ctab))
  expect_true("step_id" %in% names(ctab))
  expect_true(all(ctab$step_id == 1L))
})

test_that("sda_candidate_table: out-of-range step errors", {
  n_steps <- length(.sda_iop_fit$steps)
  expect_error(
    sda_candidate_table(.sda_iop_fit, step = n_steps + 1L),
    regexp = "integer in 1"
  )
})

test_that("sda_candidate_table: each table has one selected=TRUE row per step", {
  tabs <- sda_candidate_table(.sda_iop_fit)
  for (tab in tabs)
    expect_equal(sum(tab$selected), 1L)
})

# ---------------------------------------------------------------------------
# as_cta_candidates
# ---------------------------------------------------------------------------
test_that("as_cta_candidates: returns X subset to selected attributes", {
  X_cta <- as_cta_candidates(.sda_iop_fit, .sda_iop_X)
  sel   <- sda_selected_attributes(.sda_iop_fit)
  expect_true(is.data.frame(X_cta))
  expect_equal(ncol(X_cta), length(sel))
  expect_equal(colnames(X_cta), sel)
  expect_equal(nrow(X_cta), nrow(.sda_iop_X))
})

test_that("as_cta_candidates: wide X (extra columns) is accepted", {
  wide_X <- cbind(.sda_iop_X, extra = rnorm(100))
  X_cta  <- as_cta_candidates(.sda_iop_fit, wide_X)
  sel    <- sda_selected_attributes(.sda_iop_fit)
  expect_equal(colnames(X_cta), sel)
})

test_that("as_cta_candidates: errors when required column missing", {
  sel      <- sda_selected_attributes(.sda_iop_fit)
  missing  <- sel[1]
  X_bad    <- .sda_iop_X[, setdiff(colnames(.sda_iop_X), missing), drop = FALSE]
  expect_error(as_cta_candidates(.sda_iop_fit, X_bad),
               regexp = "missing SDA-selected")
})

# ---------------------------------------------------------------------------
# sda_to_cta_data
# ---------------------------------------------------------------------------
test_that("sda_to_cta_data: returns list with X_cta and y_cta", {
  out <- sda_to_cta_data(.sda_iop_fit, .sda_iop_X, .sda_iop_y)
  expect_true(is.list(out))
  expect_true(all(c("X_cta", "y_cta") %in% names(out)))
})

test_that("sda_to_cta_data: X_cta has selected columns; y_cta is full length", {
  out <- sda_to_cta_data(.sda_iop_fit, .sda_iop_X, .sda_iop_y)
  sel <- sda_selected_attributes(.sda_iop_fit)
  expect_equal(colnames(out$X_cta), sel)
  expect_equal(length(out$y_cta), 100L)
  expect_true(is.integer(out$y_cta))
})

test_that("sda_to_cta_data: nrow(X) != length(y) errors", {
  expect_error(
    sda_to_cta_data(.sda_iop_fit, .sda_iop_X, .sda_iop_y[1:50]),
    regexp = "nrow\\(X\\)"
  )
})

# ---------------------------------------------------------------------------
# End-to-end: SDA -> cta_fit pipeline
# ---------------------------------------------------------------------------
test_that("sda -> cta_fit end-to-end: runs without error", {
  out     <- sda_to_cta_data(.sda_iop_fit, .sda_iop_X, .sda_iop_y)
  # cta_fit on SDA-constrained candidate frame (non-recursive, fast settings)
  cta_result <- cta_fit(
    out$X_cta, out$y_cta,
    mindenom = 1L,
    mc_iter  = 200L,
    mc_seed  = 42L,
    loo      = "off"
  )
  expect_true(inherits(cta_result, "cta_tree"))
})

test_that("sda -> cta_fit: selected candidates are the only CTA columns", {
  sel    <- sda_selected_attributes(.sda_iop_fit)
  X_cta  <- as_cta_candidates(.sda_iop_fit, .sda_iop_X)
  expect_equal(colnames(X_cta), sel)
  # CTA receives exactly the SDA-constrained frame — no extra columns
  expect_equal(ncol(X_cta), length(sel))
})
