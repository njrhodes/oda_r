###############################################################################
# test-fixture-vignettes.R - training gold for four legacy MegaODA vignettes
#
# Fixtures: tests/testthat/fixtures/vignettes/data/example-{1..4}/
# Coverage: rule type/cut, confusion matrix (raw integer, actual x predicted), ESS.
# All CRAN-safe (mcarlo = FALSE).
#
# Ex1 (n=407): binary class (vote Con/Pro), binary ordered attr (party affil)
# Ex2 (n=720): binary class (motivation Ind/Com), 4-cat nominal attr (adjustment)
# Ex3 (n=325): 4-class (protein type), 4-cat nominal attr (amino acid type)
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
