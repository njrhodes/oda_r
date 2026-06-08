###############################################################################
# test-fixture-vignettes.R - training gold for four legacy MegaODA vignettes
#
# Fixtures: tests/testthat/fixtures/vignettes/data/example-{1..4}/
# Coverage: rule type/cut, confusion matrix (raw integer, actual x predicted), ESS,
#           LOO ESS and LOO confusion parity.
# All CRAN-safe (mcarlo = FALSE).
#
# Ex1 (n=407): binary class (vote Con/Pro), binary ordered attr (party affil)
# Ex2 (n=720): binary class (motivation Ind/Com), 4-cat nominal attr (adjustment)
#              Gold: nondirectional categorical ODA (no DIRECTION command in EXE)
# Ex3 (n=325): 4-class (protein type), 4-cat nominal attr (amino acid type)
#              Gold: nondirectional categorical ODA (no DIRECTION command in EXE)
# Ex4 (n= 67): binary class (treatment arm), ordered attr (migraine attacks 0-7)
###############################################################################

.binary_confusion_matrix <- function(conf) {
  matrix(
    c(conf$TN, conf$FP,
      conf$FN, conf$TP),
    nrow = 2L, byrow = TRUE
  )
}

.read_vignette_csv <- function(n, colnames) {
  read.csv(
    test_path(sprintf("fixtures/vignettes/data/example-%d/data%d.csv", n, n)),
    header = FALSE, col.names = colnames
  )
}

# =============================================================================
# Ex1 - binary class, binary ordered attribute (MegaODA: Categorical=OFF)
# Gold rule: V2 <= 0.5 -> V1=0;  V2 > 0.5 -> V1=1
# Confusion: [[118,78],[34,177]]  ESS: 44.09%
# =============================================================================

test_that("vignette Ex1: ordered_cut at 0.5, confusion [[118,78],[34,177]], ESS=44.09", {
  df  <- .read_vignette_csv(1L, c("v1", "v2"))
  fit <- oda_fit(x = df$v2, y = df$v1,
                 attr_type = "ordered", mcarlo = FALSE, loo = "off")
  expect_true(fit$ok)
  expect_equal(fit$rule$type,      "ordered_cut")
  expect_equal(fit$rule$cut_value,  0.5)
  gold <- matrix(c(118L, 78L, 34L, 177L), nrow = 2L, byrow = TRUE)
  expect_equal(.binary_confusion_matrix(fit$confusion), gold)
  expect_equal(fit$ess, 44.09, tolerance = 0.05)
})

# =============================================================================
# Ex2 - binary class, 4-cat nominal attribute (MegaODA: Categorical=ON)
# Gold rule: V1 in {1,2}->V2=1 (Community);  V1 in {3,4}->V2=0 (Individual)
# Confusion: [[217,150],[10,343]]  ESS: 56.30%
# =============================================================================

test_that("vignette Ex2: nominal_cut {3,4}->0/{1,2}->1, confusion [[217,150],[10,343]], ESS=56.30", {
  df  <- .read_vignette_csv(2L, c("v1", "v2"))
  fit <- oda_fit(x = df$v1, y = df$v2,
                 attr_type = "categorical", mcarlo = FALSE, loo = "off")
  expect_true(fit$ok)
  expect_equal(fit$rule$type, "nominal_cut")
  expect_equal(sort(as.integer(fit$rule$left_levels)),  c(3L, 4L),
               label = "left levels -> Individual (0)")
  expect_equal(sort(as.integer(fit$rule$right_levels)), c(1L, 2L),
               label = "right levels -> Community (1)")
  gold <- matrix(c(217L, 150L, 10L, 343L), nrow = 2L, byrow = TRUE)
  expect_equal(.binary_confusion_matrix(fit$confusion), gold)
  expect_equal(fit$ess, 56.30, tolerance = 0.05)
})

test_that("vignette Ex2 LOO: ess_loo=56.30 equals training ESS (stable nondirectional rule)", {
  # Gold: MegaODA LOO same as training (optimal mapping stable across folds).
  # No direction_map: nondirectional categorical search per fold.
  df  <- .read_vignette_csv(2L, c("v1", "v2"))
  fit <- oda_fit(x = df$v1, y = df$v2,
                 attr_type = "categorical", mcarlo = FALSE, loo = "on")
  expect_true(isTRUE(fit$loo$allowed), label = "Ex2 LOO allowed")
  expect_equal(fit$loo$ess_loo, 56.30, tolerance = 0.05,
               label = "Ex2 LOO ESS equals training ESS 56.30")
  expect_lt(fit$loo$p_value, 0.001, label = "Ex2 LOO p < 0.001")
})

# =============================================================================
# Ex3 - 4-class, 4-cat nominal attribute (MegaODA: Categorical=ON, CAT v2)
# Gold rule: V2=k -> V1=k (identity mapping)
# Confusion 4x4 (actual x predicted): [[98,16,5,3],[13,50,2,8],[6,4,23,12],[7,19,14,45]]
# ESS: 50.96%
# =============================================================================

test_that("vignette Ex3: multiclass_nominal, confusion 4x4 gold, ESS=50.96", {
  df  <- .read_vignette_csv(3L, c("v1", "v2"))
  fit <- oda_fit(x = df$v2, y = df$v1,
                 attr_type = "categorical", mcarlo = FALSE, loo = "off")
  expect_true(fit$ok)
  expect_equal(fit$rule$type, "multiclass_nominal")
  gold <- matrix(c(
    98L, 16L,  5L,  3L,
    13L, 50L,  2L,  8L,
     6L,  4L, 23L, 12L,
     7L, 19L, 14L, 45L
  ), nrow = 4L, byrow = TRUE)
  expect_equal(unname(fit$confusion), gold)
  expect_equal(fit$ess, 50.96, tolerance = 0.05)
})

test_that("vignette Ex3 LOO: ess_loo=50.96 equals training ESS; LOO confusion equals training", {
  # Gold: MegaODA LOO same as training (identity mapping optimal and stable across folds).
  # No direction command: nondirectional categorical search per fold.
  df  <- .read_vignette_csv(3L, c("v1", "v2"))
  fit <- oda_fit(x = df$v2, y = df$v1,
                 attr_type = "categorical", mcarlo = FALSE, loo = "on")
  expect_true(isTRUE(fit$loo$allowed), label = "Ex3 LOO allowed")
  s <- summary(fit)
  expect_equal(s$loo$ess_loo, 50.96, tolerance = 0.05,
               label = "Ex3 LOO ESS equals training ESS 50.96")
  gold <- matrix(c(
    98L, 16L,  5L,  3L,
    13L, 50L,  2L,  8L,
     6L,  4L, 23L, 12L,
     7L, 19L, 14L, 45L
  ), nrow = 4L, byrow = TRUE)
  expect_equal(unname(s$loo$confusion), gold,
               label = "Ex3 LOO confusion equals training confusion (identity rule stable)")
})

# =============================================================================
# Ex4 - binary class, ordered attribute (multicategory, values 0-7)
# Gold rule: V1 <= 0.5 -> V2=0;  V1 > 0.5 -> V2=1
# Confusion: [[13,20],[5,29]]  ESS: 24.69%
# =============================================================================

# =============================================================================
# Ex2 LOO directional: binary class, 4-cat nominal, fixed direction_map
# Gold: same ESS/confusion as nondirectional (optimal mapping is the fixed map).
# Directional LOO applies the fixed training rule per fold; binary Fisher p
# is one-tailed per MPE p.34 hold-out canon.
# =============================================================================

test_that("vignette Ex2 directional LOO: fixed direction_map, ESS=56.30, Fisher p one-tailed", {
  df  <- .read_vignette_csv(2L, c("v1", "v2"))
  # v1 = adjustment (attribute); v2 = motivation (class)
  # direction_map a priori: {1,2}->Community(1); {3,4}->Individual(0)
  fit <- oda_fit(x = df$v1, y = df$v2,
                 attr_type     = "categorical",
                 direction_map = c("1" = 1L, "2" = 1L, "3" = 0L, "4" = 0L),
                 mcarlo = FALSE, loo = "on")
  expect_true(fit$ok)
  expect_true(isTRUE(fit$loo$allowed), label = "Ex2 directional LOO allowed")
  expect_equal(fit$ess, 56.30, tolerance = 0.05,
               label = "Ex2 directional training ESS=56.30")
  expect_equal(fit$loo$ess_loo, 56.30, tolerance = 0.05,
               label = "Ex2 directional LOO ESS=56.30")
  # Binary LOO confusion is stored as PAC fractions (hold-out contract).
  # Gold actual x predicted: [[TN=217, FP=150], [FN=10, TP=343]].
  # Recover raw counts: TN_frac * n_class0 and TP_frac * n_class1.
  n0 <- sum(df$v2 == 0L)   # 367 Individual
  n1 <- sum(df$v2 == 1L)   # 353 Community
  expect_equal(round(fit$loo$confusion$TN * n0), 217L,
               label = "Ex2 directional LOO TN=217")
  expect_equal(round(fit$loo$confusion$TP * n1), 343L,
               label = "Ex2 directional LOO TP=343")
  # Binary Fisher LOO p: one-tailed per MPE p.34 hold-out canon.
  expect_true(!is.null(fit$loo$p_value) && is.finite(fit$loo$p_value),
              label = "Ex2 directional LOO p finite")
  expect_lt(fit$loo$p_value, 0.001,     label = "Ex2 directional LOO p < 0.001")
  expect_equal(fit$loo$alternative, "greater",
               label = "Ex2 directional LOO alternative=greater (one-tailed)")
  # summary() exposes p_method and the one-tailed label.
  s <- summary(fit)
  expect_match(s$loo$p_method, "one-tailed",
               label = "Ex2 directional summary LOO p_method contains one-tailed")
  # print() must display p(LOO).  MC p and LOO p are separate calculations.
  expect_output(print(fit), "p\\(LOO\\)",
                label = "Ex2 directional print shows p(LOO)")
})

# =============================================================================
# Ex3 LOO directional: 4-class, 4-cat nominal, direction="ascending"
# Gold: same ESS/confusion as nondirectional (identity map is globally optimal).
# LOO folds search nondirectionally; no C x C Fisher p reported.
# =============================================================================

test_that("vignette Ex3 directional LOO: direction=ascending, ESS=50.96, no LOO Fisher p", {
  df  <- .read_vignette_csv(3L, c("v1", "v2"))
  # v2 = amino_acid_type (attribute); v1 = biological_type (class)
  fit <- oda_fit(x = df$v2, y = df$v1,
                 attr_type = "categorical",
                 direction  = "ascending",
                 mcarlo = FALSE, loo = "on")
  expect_true(fit$ok)
  expect_equal(fit$ess, 50.96, tolerance = 0.05,
               label = "Ex3 directional training ESS=50.96")
  s <- summary(fit)
  expect_equal(s$loo$ess_loo, 50.96, tolerance = 0.05,
               label = "Ex3 directional LOO ESS=50.96")
  # LOO confusion equals training confusion (identity map stable across folds).
  gold <- matrix(c(
    98L, 16L,  5L,  3L,
    13L, 50L,  2L,  8L,
     6L,  4L, 23L, 12L,
     7L, 19L, 14L, 45L
  ), nrow = 4L, byrow = TRUE)
  expect_equal(unname(s$loo$confusion), gold,
               label = "Ex3 directional LOO confusion equals gold")
  # Multicategorical LOO Fisher p: undefined; must be absent/NA.
  expect_true(is.na(s$loo$p_value) || is.null(s$loo$p_value),
              label = "Ex3 directional LOO p_value NA for multicategorical")
  expect_equal(s$loo$p_method, "none",
               label = "Ex3 directional LOO p_method=none")
  # print() must explicitly state p(LOO) is not reported for multicategorical.
  expect_output(print(fit), "not reported for multicategorical ODA",
                label = "Ex3 directional print says p(LOO) not reported")
})

# =============================================================================
# Ex4 - binary class, ordered attribute (multicategory, values 0-7)
# Gold rule: V1 <= 0.5 -> V2=0;  V1 > 0.5 -> V2=1
# Confusion: [[13,20],[5,29]]  ESS: 24.69%
# =============================================================================

test_that("vignette Ex4: ordered_cut at 0.5, confusion [[13,20],[5,29]], ESS=24.69", {
  df  <- .read_vignette_csv(4L, c("v1", "v2"))
  fit <- oda_fit(x = df$v1, y = df$v2,
                 attr_type = "ordered", mcarlo = FALSE, loo = "off")
  expect_true(fit$ok)
  expect_equal(fit$rule$type,      "ordered_cut")
  expect_equal(fit$rule$cut_value,  0.5)
  gold <- matrix(c(13L, 20L, 5L, 29L), nrow = 2L, byrow = TRUE)
  expect_equal(.binary_confusion_matrix(fit$confusion), gold)
  expect_equal(fit$ess, 24.69, tolerance = 0.05)
})
