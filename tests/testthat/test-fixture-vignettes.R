###############################################################################
# test-fixture-vignettes.R
#
# Training gold tests for the four legacy MegaODA vignette examples.
#
# Fixtures live under:
#   tests/testthat/fixtures/vignettes/data/example-{1..4}/
# Each folder contains data*.csv (observation-level, no header) and
# MODEL1.OUT (MegaODA gold output, re-run 2026-05-15).
#
# Coverage: deterministic training parity only.
#   - rule type and cut / level grouping
#   - confusion matrix (raw integer counts, actual × predicted)
#   - ESS
#
# Deliberately excluded from this commit:
#   - MC p-value assertions: low-mc_iter CRAN calls produce noisy estimates;
#     engine agreement with MegaODA verified at mc_iter=25000 (within MC error).
#   - LOO assertions: deferred to follow-on commit once binary LOO ok-gate
#     behavior is verified.
#
# Example summaries:
#   Ex1 (n=407) — binary class (vote Con/Pro), binary ordered attr (party affil)
#   Ex2 (n=720) — binary class (motivation Ind/Com), 4-cat nominal attr (adjustment)
#   Ex3 (n=325) — 4-class (biological protein type), 4-cat nominal attr (amino acid type)
#   Ex4 (n=67)  — binary class (treatment arm), ordered attr (migraine attacks 0–7)
###############################################################################

# ---- helpers ----------------------------------------------------------------

# Build an actual-rows × predicted-cols integer matrix from a binary fit's
# confusion list.  fit$confusion returns a list of TP/TN/FP/FN scalars.
# Orientation (byrow = TRUE):
#   row 0 (actual 0): [ TN, FP ]   (pred 0 = TN, pred 1 = FP)
#   row 1 (actual 1): [ FN, TP ]   (pred 0 = FN, pred 1 = TP)
.binary_confusion_matrix <- function(conf) {
  matrix(
    c(conf$TN, conf$FP,
      conf$FN, conf$TP),
    nrow = 2L,
    byrow = TRUE
  )
}

.read_vignette_csv <- function(n, colnames) {
  read.csv(
    test_path(sprintf("fixtures/vignettes/data/example-%d/data%d.csv", n, n)),
    header = FALSE, col.names = colnames
  )
}

# ---- Example 1: binary class, binary attribute (ordered) --------------------
#
# Refugee Act 1980, U.S. House (n=407).
# Class v1: vote (0=Con, 1=Pro).  Attribute v2: party (0=Rep, 1=Dem).
# MegaODA: Categorical=OFF (ordered scan).
# Gold rule: V2 <= 0.5 → V1=0;  V2 > 0.5 → V1=1.
# Confusion: [[118, 78], [34, 177]]   ESS: 44.09%

test_that("vignette Ex1: rule type and cut value", {
  df  <- .read_vignette_csv(1L, c("v1", "v2"))
  fit <- oda_fit(x = df$v2, y = df$v1,
                 attr_type = "ordered",
                 mcarlo    = FALSE,
                 loo       = "off")

  expect_true(fit$ok)
  expect_equal(fit$rule$type,      "ordered_cut")
  expect_equal(fit$rule$cut_value, 0.5)
})

test_that("vignette Ex1: confusion matrix matches MegaODA gold", {
  df  <- .read_vignette_csv(1L, c("v1", "v2"))
  fit <- oda_fit(x = df$v2, y = df$v1,
                 attr_type = "ordered",
                 mcarlo    = FALSE,
                 loo       = "off")

  gold <- matrix(c(118L, 78L,
                    34L, 177L), nrow = 2L, byrow = TRUE)
  expect_equal(.binary_confusion_matrix(fit$confusion), gold,
               label = "Ex1 confusion (actual × predicted)")
})

test_that("vignette Ex1: ESS matches MegaODA gold (44.09)", {
  df  <- .read_vignette_csv(1L, c("v1", "v2"))
  fit <- oda_fit(x = df$v2, y = df$v1,
                 attr_type = "ordered",
                 mcarlo    = FALSE,
                 loo       = "off")

  expect_equal(fit$ess, 44.09, tolerance = 0.05, label = "Ex1 ESS")
})

# ---- Example 2: binary class, multicategorical attribute --------------------
#
# Gully erosion adjustment study (n=720).
# Class v2: motivation (0=Individual, 1=Community).
# Attribute v1: adjustment type (1=Ridges, 2=Shifting, 3=Relocation, 4=Intensified).
# MegaODA: Categorical=ON, Degen=ON.
# Note: Degen=ON in MegaODA is binary-display only; odacore's degen param is
#   multiclass-only and is not passed to oda_univariate_core for C=2.
# Gold rule: V1∈{1,2}→V2=1;  V1∈{3,4}→V2=0.
# Confusion: [[217, 150], [10, 343]]   ESS: 56.30%

test_that("vignette Ex2: rule type and level grouping", {
  df  <- .read_vignette_csv(2L, c("v1", "v2"))
  fit <- oda_fit(x = df$v1, y = df$v2,
                 attr_type = "categorical",
                 mcarlo    = FALSE,
                 loo       = "off")

  expect_true(fit$ok)
  expect_equal(fit$rule$type, "nominal_cut")
  # {3,4} → Individual (class 0 = left side)
  # {1,2} → Community  (class 1 = right side)
  expect_equal(sort(as.integer(fit$rule$left_levels)),  c(3L, 4L),
               label = "Ex2 left levels → Individual (0)")
  expect_equal(sort(as.integer(fit$rule$right_levels)), c(1L, 2L),
               label = "Ex2 right levels → Community (1)")
})

test_that("vignette Ex2: confusion matrix matches MegaODA gold", {
  df  <- .read_vignette_csv(2L, c("v1", "v2"))
  fit <- oda_fit(x = df$v1, y = df$v2,
                 attr_type = "categorical",
                 mcarlo    = FALSE,
                 loo       = "off")

  gold <- matrix(c(217L, 150L,
                    10L, 343L), nrow = 2L, byrow = TRUE)
  expect_equal(.binary_confusion_matrix(fit$confusion), gold,
               label = "Ex2 confusion (actual × predicted)")
})

test_that("vignette Ex2: ESS matches MegaODA gold (56.30)", {
  df  <- .read_vignette_csv(2L, c("v1", "v2"))
  fit <- oda_fit(x = df$v1, y = df$v2,
                 attr_type = "categorical",
                 mcarlo    = FALSE,
                 loo       = "off")

  expect_equal(fit$ess, 56.30, tolerance = 0.05, label = "Ex2 ESS")
})

# ---- Example 3: 4-class variable, multicategorical attribute ----------------
#
# Protein classification — biological type vs. amino acid type (n=325).
# Class v1: biological type (1–4).  Attribute v2: amino acid type (1–4).
# MegaODA: Categorical=ON (CAT v2), 4 classes.
# Gold rule: V2=k → V1=k for k∈{1,2,3,4} (identity mapping).
# Confusion 4×4 (actual rows × predicted cols, class labels {1,2,3,4}):
#   [[98,16,5,3],[13,50,2,8],[6,4,23,12],[7,19,14,45]]
# ESS: 50.96%
# LOO: not asserted — per Rmd, LOO is not conducted for multicategorical class.

test_that("vignette Ex3: rule type is multiclass_nominal", {
  df  <- .read_vignette_csv(3L, c("v1", "v2"))
  fit <- oda_fit(x = df$v2, y = df$v1,
                 attr_type = "categorical",
                 mcarlo    = FALSE,
                 loo       = "off")

  expect_true(fit$ok)
  expect_equal(fit$rule$type, "multiclass_nominal")
})

test_that("vignette Ex3: confusion matrix matches MegaODA gold", {
  df  <- .read_vignette_csv(3L, c("v1", "v2"))
  fit <- oda_fit(x = df$v2, y = df$v1,
                 attr_type = "categorical",
                 mcarlo    = FALSE,
                 loo       = "off")

  gold <- matrix(c(
    98L, 16L,  5L,  3L,
    13L, 50L,  2L,  8L,
     6L,  4L, 23L, 12L,
     7L, 19L, 14L, 45L
  ), nrow = 4L, byrow = TRUE)
  expect_equal(unname(fit$confusion), gold,
               label = "Ex3 confusion (actual × predicted)")
})

test_that("vignette Ex3: ESS matches MegaODA gold (50.96)", {
  df  <- .read_vignette_csv(3L, c("v1", "v2"))
  fit <- oda_fit(x = df$v2, y = df$v1,
                 attr_type = "categorical",
                 mcarlo    = FALSE,
                 loo       = "off")

  expect_equal(fit$ess, 50.96, tolerance = 0.05, label = "Ex3 ESS")
})

# ---- Example 4: binary class, ordered attribute (multicategory) -------------
#
# Migraine attacks clinical trial (n=67).
# Class v2: treatment arm (0=Treatment1, 1=Treatment2).
# Attribute v1: number of migraine attacks (ordered integer, 0–7).
# MegaODA: Categorical=OFF (ordered scan).
# Gold rule: V1 <= 0.5 → V2=0;  V1 > 0.5 → V2=1.
# Confusion: [[13, 20], [5, 29]]   ESS: 24.69%
# Note: MegaODA MC p ≈ 0.0859 (NS). MC p-bucket assertion deferred — see #5.

test_that("vignette Ex4: rule type and cut value", {
  df  <- .read_vignette_csv(4L, c("v1", "v2"))
  fit <- oda_fit(x = df$v1, y = df$v2,
                 attr_type = "ordered",
                 mcarlo    = FALSE,
                 loo       = "off")

  expect_true(fit$ok)
  expect_equal(fit$rule$type,      "ordered_cut")
  expect_equal(fit$rule$cut_value, 0.5)
})

test_that("vignette Ex4: confusion matrix matches MegaODA gold", {
  df  <- .read_vignette_csv(4L, c("v1", "v2"))
  fit <- oda_fit(x = df$v1, y = df$v2,
                 attr_type = "ordered",
                 mcarlo    = FALSE,
                 loo       = "off")

  gold <- matrix(c(13L, 20L,
                    5L, 29L), nrow = 2L, byrow = TRUE)
  expect_equal(.binary_confusion_matrix(fit$confusion), gold,
               label = "Ex4 confusion (actual × predicted)")
})

test_that("vignette Ex4: ESS matches MegaODA gold (24.69)", {
  df  <- .read_vignette_csv(4L, c("v1", "v2"))
  fit <- oda_fit(x = df$v1, y = df$v2,
                 attr_type = "ordered",
                 mcarlo    = FALSE,
                 loo       = "off")

  expect_equal(fit$ess, 24.69, tolerance = 0.05, label = "Ex4 ESS")
})
