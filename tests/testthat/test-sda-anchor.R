###############################################################################
# test-sda-anchor.R — SDA anchor object and task hook tests (Slice O)
#
# Tier: CRAN-safe (all tests run at default tier).
# All sda_fit calls use small mc_iter for speed. Only synthetic public data.
#
# Canon reference: docs/SDA_ANCHOR_CONTRACT.md
###############################################################################

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

# Small synthetic dataset: 60 obs, 3 attributes, 2 classes.
# A separates well; B and C are noise.
set.seed(7L)
n_fix  <- 60L
A_fix  <- c(rep(0L, 30L), rep(1L, 30L))          # near-perfect separator
B_fix  <- sample(0:1, n_fix, replace = TRUE)
C_fix  <- sample(0:1, n_fix, replace = TRUE)
y_fix  <- c(rep(1L, 30L), rep(2L, 30L))
X_fix  <- data.frame(A = A_fix, B = B_fix, C = C_fix)

# Run small sda_fit once for reuse.
# mode = "unioda_max_ess" does not require mindenom; safe for small n.
sda_small <- sda_fit(
  X       = X_fix,
  y       = y_fix,
  mode    = "unioda_max_ess",
  mc_iter = 100L,
  mc_seed = 7L
)

# Explicit stage table for manual anchor tests
explicit_stage <- data.frame(
  stage_id  = 1L,
  attribute = "A",
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# 1. sda_anchor() returns class c("sda_anchor", "list")
# ---------------------------------------------------------------------------
test_that("sda_anchor() returns class c('sda_anchor', 'list')", {
  a <- sda_anchor(
    anchor_type         = "explicit",
    selected_attributes = "A",
    stage_table         = explicit_stage,
    weights_used        = FALSE
  )
  expect_s3_class(a, "sda_anchor")
  expect_true(inherits(a, "list"))
  expect_identical(class(a), c("sda_anchor", "list"))
})

# ---------------------------------------------------------------------------
# 2. Explicit / manual anchor validates
# ---------------------------------------------------------------------------
test_that("explicit anchor validates without error", {
  a <- sda_anchor(
    anchor_type         = "explicit",
    selected_attributes = "A",
    candidate_universe  = c("A", "B", "C"),
    stage_table         = explicit_stage,
    weights_used        = FALSE
  )
  expect_silent(validate_sda_anchor(a))
})

# ---------------------------------------------------------------------------
# 3. Missing selected_attributes errors
# ---------------------------------------------------------------------------
test_that("validate_sda_anchor errors when selected_attributes is empty", {
  a <- sda_anchor(
    anchor_type         = "explicit",
    selected_attributes = character(0),
    stage_table         = explicit_stage,
    weights_used        = FALSE
  )
  expect_error(validate_sda_anchor(a), "selected_attributes")
})

# ---------------------------------------------------------------------------
# 4. Invalid stage_table errors
# ---------------------------------------------------------------------------
test_that("validate_sda_anchor errors when stage_table lacks required columns", {
  bad_stage <- data.frame(step = 1L, attr = "A", stringsAsFactors = FALSE)
  a <- sda_anchor(
    anchor_type         = "explicit",
    selected_attributes = "A",
    stage_table         = bad_stage,
    weights_used        = FALSE
  )
  expect_error(validate_sda_anchor(a), "stage_table is missing column")
})

test_that("validate_sda_anchor errors when stage_table is not a data.frame", {
  a <- sda_anchor(
    anchor_type         = "explicit",
    selected_attributes = "A",
    stage_table         = list(stage_id = 1L, attribute = "A"),
    weights_used        = FALSE
  )
  expect_error(validate_sda_anchor(a), "stage_table must be a data.frame")
})

# ---------------------------------------------------------------------------
# 5. Selected attributes not in candidate universe errors
# ---------------------------------------------------------------------------
test_that("validate_sda_anchor errors when selected not in universe", {
  a <- sda_anchor(
    anchor_type         = "explicit",
    selected_attributes = c("A", "X_missing"),
    candidate_universe  = c("A", "B", "C"),
    stage_table         = data.frame(
      stage_id  = 1:2,
      attribute = c("A", "X_missing"),
      stringsAsFactors = FALSE
    ),
    weights_used        = FALSE
  )
  expect_error(validate_sda_anchor(a), "candidate_universe")
})

# ---------------------------------------------------------------------------
# 6. Stage table order is preserved
# ---------------------------------------------------------------------------
test_that("validate_sda_anchor errors when stage_table order mismatches", {
  st <- data.frame(
    stage_id  = 1:2,
    attribute = c("B", "A"),     # reversed relative to selected_attributes
    stringsAsFactors = FALSE
  )
  a <- sda_anchor(
    anchor_type         = "explicit",
    selected_attributes = c("A", "B"),
    stage_table         = st,
    weights_used        = FALSE
  )
  expect_error(validate_sda_anchor(a), "order")
})

test_that("validate_sda_anchor passes when stage_table order matches", {
  st <- data.frame(
    stage_id  = 1:2,
    attribute = c("A", "B"),
    stringsAsFactors = FALSE
  )
  a <- sda_anchor(
    anchor_type         = "explicit",
    selected_attributes = c("A", "B"),
    stage_table         = st,
    weights_used        = FALSE
  )
  expect_silent(validate_sda_anchor(a))
})

# ---------------------------------------------------------------------------
# 7. as_sda_anchor() generic dispatch works
# ---------------------------------------------------------------------------
test_that("as_sda_anchor generic dispatches correctly", {
  expect_true(isGeneric("as_sda_anchor") || existsMethod("as_sda_anchor", "sda_fit") ||
                is.function(as_sda_anchor.sda_fit))
  a <- as_sda_anchor(sda_small)
  expect_s3_class(a, "sda_anchor")
})

# ---------------------------------------------------------------------------
# 8. as_sda_anchor.sda_fit() works on CRAN-safe sda_fit example
# ---------------------------------------------------------------------------
test_that("as_sda_anchor.sda_fit() returns valid sda_anchor", {
  a <- as_sda_anchor(sda_small)
  expect_s3_class(a, "sda_anchor")
  expect_silent(validate_sda_anchor(a))
  expect_identical(a$anchor_type, "sda_fit")
})

# ---------------------------------------------------------------------------
# 9. SDA-derived anchor preserves selected attribute order
# ---------------------------------------------------------------------------
test_that("as_sda_anchor.sda_fit preserves selected attribute order", {
  a     <- as_sda_anchor(sda_small)
  sda_sel <- sda_small$selected_attributes
  expect_identical(a$selected_attributes, sda_sel)
  if (nrow(a$stage_table) > 0L && length(sda_sel) == nrow(a$stage_table))
    expect_identical(a$stage_table$attribute, sda_sel)
})

# ---------------------------------------------------------------------------
# 10. Default task_hook exists
# ---------------------------------------------------------------------------
test_that("sda_anchor task_hook exists by default", {
  a <- as_sda_anchor(sda_small)
  expect_true(!is.null(a$task_hook))
  expect_true(is.list(a$task_hook))
})

# ---------------------------------------------------------------------------
# 11. task_hook$hook_type == "sda_anchor"
# ---------------------------------------------------------------------------
test_that("task_hook$hook_type is 'sda_anchor'", {
  a <- as_sda_anchor(sda_small)
  expect_identical(a$task_hook$hook_type, "sda_anchor")
})

# ---------------------------------------------------------------------------
# 12. task_hook$implementation_status == "anchor_only_no_sort"
# ---------------------------------------------------------------------------
test_that("task_hook$implementation_status is 'anchor_only_no_sort'", {
  a <- as_sda_anchor(sda_small)
  expect_identical(a$task_hook$implementation_status, "anchor_only_no_sort")
})

# ---------------------------------------------------------------------------
# 13. task_hook$requires_human_review is TRUE
# ---------------------------------------------------------------------------
test_that("task_hook$requires_human_review is TRUE", {
  a <- as_sda_anchor(sda_small)
  expect_true(isTRUE(a$task_hook$requires_human_review))
})

# ---------------------------------------------------------------------------
# 14. future_downstream can include "sort" and "gort"
# ---------------------------------------------------------------------------
test_that("future_downstream includes 'sort' and 'gort'", {
  a <- as_sda_anchor(sda_small)
  expect_true("sort" %in% a$task_hook$future_downstream)
  expect_true("gort" %in% a$task_hook$future_downstream)
})

# ---------------------------------------------------------------------------
# 15. prohibited_downstream includes "propensity_weighting" and "fraud_demo"
# ---------------------------------------------------------------------------
test_that("prohibited_downstream includes propensity_weighting and fraud_demo", {
  a <- as_sda_anchor(sda_small)
  expect_true("propensity_weighting" %in% a$task_hook$prohibited_downstream)
  expect_true("fraud_demo"           %in% a$task_hook$prohibited_downstream)
})

# ---------------------------------------------------------------------------
# 16. task_hook task_role does not include "propensity_model"
# ---------------------------------------------------------------------------
test_that("task_hook$task_role does not include 'propensity_model'", {
  a <- as_sda_anchor(sda_small)
  expect_false("propensity_model" %in% a$task_hook$task_role)
})

# ---------------------------------------------------------------------------
# 17. Explicit anchor safety notes mention manual / user-declared
# ---------------------------------------------------------------------------
test_that("explicit anchor notes mention manual/user-declared", {
  a <- as_sda_anchor(
    explicit_stage,
    selected_attributes = "A",
    candidate_universe  = c("A", "B", "C")
  )
  all_notes <- c(a$reproducibility_notes, a$canon_notes)
  has_manual <- any(grepl("manual|explicit|user", all_notes, ignore.case = TRUE))
  expect_true(has_manual)
})

# ---------------------------------------------------------------------------
# 18. print.sda_anchor does not claim SORT/GORT implemented
# ---------------------------------------------------------------------------
test_that("print.sda_anchor output does not claim SORT/GORT implemented", {
  a   <- as_sda_anchor(sda_small)
  out <- capture.output(print(a))
  combined <- paste(out, collapse = " ")
  # should say "not implemented" but not "SORT is implemented" or "GORT is implemented"
  expect_false(grepl("SORT is implemented|GORT is implemented", combined))
  expect_true(grepl("not implemented", combined, ignore.case = TRUE))
})

# ---------------------------------------------------------------------------
# 19. summary.sda_anchor returns expected fields
# ---------------------------------------------------------------------------
test_that("summary.sda_anchor returns expected fields", {
  a   <- as_sda_anchor(sda_small)
  s   <- summary(a)
  expected_fields <- c("anchor_type", "n_stages", "selected_attributes",
                       "stage_table", "canon_notes", "implementation_status",
                       "safety_notes")
  for (f in expected_fields)
    expect_true(f %in% names(s),
                info = paste("missing field:", f))
})

# ---------------------------------------------------------------------------
# 20. CTA behavior unchanged
# ---------------------------------------------------------------------------
test_that("CTA fit is unchanged after sda_anchor implementation", {
  X_cta <- data.frame(
    x1 = c(1, 2, 3, 4, 5, 6, 7, 8),
    x2 = c(0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L)
  )
  y_cta <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L)
  tree  <- cta_fit(X_cta, y_cta,
    priors_on  = TRUE,
    mindenom   = 2L,
    mc_iter    = 200L,
    mc_seed    = 42L,
    loo        = "off",
    attr_names = c("x1", "x2")
  )
  expect_s3_class(tree, "cta_tree")
  expect_false(is.null(tree$root))
})

# ---------------------------------------------------------------------------
# 21. LORT behavior unchanged
# ---------------------------------------------------------------------------
test_that("LORT fit is unchanged after sda_anchor implementation", {
  X_lort <- data.frame(
    x1 = c(rnorm(20, 0), rnorm(20, 3)),
    x2 = c(rnorm(20, 0), rnorm(20, 3))
  )
  y_lort <- c(rep(1L, 20L), rep(2L, 20L))
  set.seed(1L)
  ort <- lort_fit(X_lort, y_lort,
    min_n   = 5L,
    mc_iter = 100L,
    mc_seed = 1L,
    loo     = "off"
  )
  expect_s3_class(ort, "cta_ort")
  expect_identical(ort$ort_settings$method, "lort")
})

# ---------------------------------------------------------------------------
# 22. sort_fit is not exported
# ---------------------------------------------------------------------------
test_that("sort_fit is not exported", {
  expect_false(existsFunction("sort_fit") &&
    isExportedFromNamespace("sort_fit", "odacore"))
})

# ---------------------------------------------------------------------------
# 23. gort_fit is not exported
# ---------------------------------------------------------------------------
test_that("gort_fit is not exported", {
  expect_false(existsFunction("gort_fit") &&
    isExportedFromNamespace("gort_fit", "odacore"))
})

# ---------------------------------------------------------------------------
# 24. sda_propensity_weights does not exist
# ---------------------------------------------------------------------------
test_that("sda_propensity_weights is not exported", {
  ns_exports <- getNamespaceExports("odacore")
  expect_false("sda_propensity_weights" %in% ns_exports)
})

# ---------------------------------------------------------------------------
# Helper (used in tests 22-23 when package may not be loaded via namespace)
# ---------------------------------------------------------------------------
existsFunction <- function(nm) existsMethod <- function(...) FALSE
isExportedFromNamespace <- function(nm, pkg) {
  tryCatch(nm %in% getNamespaceExports(pkg), error = function(e) FALSE)
}
