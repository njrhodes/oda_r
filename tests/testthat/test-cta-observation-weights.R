###############################################################################
# test-cta-observation-weights.R — cta_observation_weights()
#
# Per-observation propensity weight assignment.
# Fast tier: synthetic trees only.
# Slow tier: myeloma fixture regression anchors.
###############################################################################

# ---- Shared synthetic fixtures -----------------------------------------------

# 8-observation binary dataset; stump with cut at x=4.5
#   node2 (x<=4.5): obs 1–4, class0=4, class1=0   => endpoint 1
#   node3 (x>4.5):  obs 5–8, class0=0, class1=4   => endpoint 2
.ow_stump_fit <- function() {
  X <- data.frame(x = 1:8)
  y <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  suppressMessages(
    oda_cta_fit(X, y, mindenom = 2L, mc_iter = 300L, mc_seed = 1L, loo = "off")
  )
}

.ow_no_tree_fit <- function() {
  X <- data.frame(x = 1:8)
  y <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  suppressMessages(
    oda_cta_fit(X, y, mindenom = 999L, mc_iter = 100L, mc_seed = 1L,
                loo = "off")
  )
}

# ---- Row count invariant -----------------------------------------------------

test_that("output has exactly nrow(newdata) rows for stump", {
  tree <- .ow_stump_fit()
  X    <- data.frame(x = 1:8)
  y    <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  ow   <- cta_observation_weights(tree, X, y)
  expect_equal(nrow(ow), 8L)
})

test_that("output has exactly nrow(newdata) rows for no-tree fit", {
  tree <- .ow_no_tree_fit()
  X    <- data.frame(x = 1:8)
  y    <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  ow   <- cta_observation_weights(tree, X, y)
  expect_equal(nrow(ow), 8L)
})

# ---- Column presence ---------------------------------------------------------

test_that("output has all 11 required columns", {
  tree <- .ow_stump_fit()
  X    <- data.frame(x = 1:8)
  y    <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  ow   <- cta_observation_weights(tree, X, y)
  expected_cols <- c("row_id", "actual_class", "endpoint_node_id",
                     "endpoint_id", "target_class", "propensity_weight",
                     "adjusted_propensity_weight", "undefined_empirical",
                     "perfectly_predicted_endpoint", "adjusted", "assigned")
  expect_true(all(expected_cols %in% names(ow)))
})

# ---- row_id column -----------------------------------------------------------

test_that("row_id is sequential 1:n", {
  tree <- .ow_stump_fit()
  X    <- data.frame(x = 1:8)
  y    <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  ow   <- cta_observation_weights(tree, X, y)
  expect_identical(ow$row_id, 1:8)
})

# ---- actual_class column -----------------------------------------------------

test_that("actual_class is character coercion of y", {
  tree <- .ow_stump_fit()
  X    <- data.frame(x = 1:8)
  y    <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  ow   <- cta_observation_weights(tree, X, y)
  expect_equal(ow$actual_class, as.character(y))
})

# ---- assigned column — no-tree fit -------------------------------------------

test_that("assigned is FALSE for all rows in no-tree fit", {
  tree <- .ow_no_tree_fit()
  X    <- data.frame(x = 1:8)
  y    <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  ow   <- cta_observation_weights(tree, X, y)
  expect_true(all(!ow$assigned))
})

test_that("weight columns are NA for no-tree fit", {
  tree <- .ow_no_tree_fit()
  X    <- data.frame(x = 1:8)
  y    <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  ow   <- cta_observation_weights(tree, X, y)
  expect_true(all(is.na(ow$propensity_weight)))
  expect_true(all(is.na(ow$adjusted_propensity_weight)))
})

# ---- assigned column — stump (perfectly predicted endpoints) -----------------

test_that("all 8 stump observations are assigned when endpoints are perfect", {
  tree <- .ow_stump_fit()
  X    <- data.frame(x = 1:8)
  y    <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  ow   <- cta_observation_weights(tree, X, y, adjusted = FALSE)
  # endpoint 1: all class 0 → undefined_empirical for class1; class0 assigned
  # endpoint 2: all class 1 → undefined_empirical for class0; class1 assigned
  # Each obs has a matching (endpoint, actual_class) cell
  expect_true(all(ow$assigned))
})

# ---- NA y label → assigned = FALSE ------------------------------------------

test_that("NA y gives assigned=FALSE and NA weight columns", {
  tree <- .ow_stump_fit()
  X    <- data.frame(x = 1:8)
  y    <- c(NA_integer_, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  ow   <- cta_observation_weights(tree, X, y)
  expect_false(ow$assigned[1L])
  expect_true(is.na(ow$propensity_weight[1L]))
  # Other rows still assigned
  expect_true(all(ow$assigned[-1L]))
})

# ---- length(y) mismatch error ------------------------------------------------

test_that("length(y) != nrow(newdata) is an error", {
  tree <- .ow_stump_fit()
  X    <- data.frame(x = 1:8)
  expect_error(
    cta_observation_weights(tree, X, y = 1:5),
    regexp = "length\\(y\\)"
  )
})

# ---- non-data.frame newdata coercion -----------------------------------------

test_that("matrix newdata is coerced to data.frame", {
  tree <- .ow_stump_fit()
  X    <- matrix(1:8, ncol = 1L, dimnames = list(NULL, "x"))
  y    <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  expect_no_error(cta_observation_weights(tree, X, y))
})

# ---- missing_action propagation ----------------------------------------------

test_that("missing_action='majority' produces same nrow as 'na'", {
  tree <- .ow_stump_fit()
  X    <- data.frame(x = c(NA, 1, 2, 3, 5, 6, 7, 8))
  y    <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  ow_na  <- cta_observation_weights(tree, X, y, missing_action = "na")
  ow_maj <- cta_observation_weights(tree, X, y, missing_action = "majority")
  expect_equal(nrow(ow_na),  8L)
  expect_equal(nrow(ow_maj), 8L)
  # Row 1: NA under "na", routed under "majority"
  expect_true(is.na(ow_na$endpoint_id[1L]))
  expect_false(is.na(ow_maj$endpoint_id[1L]))
})

# ---- target_class annotation -------------------------------------------------

test_that("target_class is annotation-only: all obs remain assigned", {
  # cta_propensity_weights target_class is annotation — does not filter rows.
  # The pw table always has all (endpoint, class) combinations.
  tree <- .ow_stump_fit()
  X    <- data.frame(x = 1:8)
  y    <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  ow1 <- cta_observation_weights(tree, X, y, target_class = 1L,
                                  adjusted = FALSE)
  # All 8 obs classified and matched regardless of target_class
  expect_true(all(ow1$assigned))
  # target_class column is integer and reflects the annotation
  expect_type(ow1$target_class, "integer")
  expect_true(all(ow1$target_class[ow1$assigned] == 1L))
})

# ---- unmatched classified obs warning ----------------------------------------

test_that("unmatched classified obs triggers warning", {
  tree <- .ow_stump_fit()
  X    <- data.frame(x = 1:8)
  # Class "2" does not appear in the fitted tree
  y    <- c(2L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  expect_warning(
    cta_observation_weights(tree, X, y),
    regexp = "could not be matched"
  )
})

# ---- not a cta_tree error ----------------------------------------------------

test_that("non-cta_tree tree argument is an error", {
  expect_error(
    cta_observation_weights(list(), data.frame(x = 1:4), y = 1:4),
    regexp = "cta_tree"
  )
})

# ==============================================================================
# Slow-tier myeloma regression anchors
# ==============================================================================

# data.txt has no header; columns are named V1..V19 after loading.
# V1 = class, V2 = weight (case weight; EX V2=0 excludes zero-weight rows).
# Attributes: V4 V9 V11 V12 V14 V15 V16 V17 V18 V19.
.ow_load_myeloma <- function() {
  fixture_dir <- testthat::test_path("fixtures", "myeloma")
  raw <- read.table(file.path(fixture_dir, "data.txt"), header = FALSE)
  names(raw) <- paste0("V", seq_len(ncol(raw)))   # V1..V19
  raw <- raw[raw$V2 != 0, ]                        # EX V2=0
  raw
}

.ow_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.ow_myeloma_fit <- function(mindenom) {
  d <- .ow_load_myeloma()
  suppressMessages(
    oda_cta_fit(
      X           = d[, .ow_myeloma_attrs],
      y           = as.integer(d$V1),
      w           = d$V2,
      mindenom    = mindenom,
      mc_iter     = 25000L,
      mc_seed     = 42L,
      miss_codes  = -9,
      loo         = "stable",
      attr_names  = .ow_myeloma_attrs
    )
  )
}

test_that("myeloma MINDENOM=1 observation weights — row count and assignment", {
  skip_if_slow_tests_disabled("fixture-myeloma-obs-weights")

  d    <- .ow_load_myeloma()
  tree <- .ow_myeloma_fit(1L)
  X    <- d[, .ow_myeloma_attrs]
  y    <- as.integer(d$V1)
  ow   <- cta_observation_weights(tree, X, y, adjusted = TRUE)

  # nrow invariant
  expect_equal(nrow(ow), nrow(X))
  # row_id is sequential
  expect_identical(ow$row_id, seq_len(nrow(X)))
  # assigned must be logical
  expect_type(ow$assigned, "logical")
  # Obs with NA endpoint_node_id are not assigned
  na_ep <- is.na(ow$endpoint_node_id)
  expect_true(all(!ow$assigned[na_ep]))
  # Assigned obs with defined empirical cell have non-NA propensity_weight
  assigned_rows <- ow[ow$assigned, ]
  non_undef <- assigned_rows[!isTRUE(assigned_rows$undefined_empirical), ]
  expect_true(all(!is.na(non_undef$propensity_weight)))
})

test_that("myeloma MINDENOM=1 assigned count matches classified n=255", {
  skip_if_slow_tests_disabled("fixture-myeloma-obs-weights")

  d    <- .ow_load_myeloma()
  tree <- .ow_myeloma_fit(1L)
  X    <- d[, .ow_myeloma_attrs]
  y    <- as.integer(d$V1)
  ow   <- cta_observation_weights(tree, X, y, adjusted = TRUE)

  # MINDENOM=1 root=V14; obs with V14=-9 (missing) are unroutable.
  # Canon: 255 obs classified.
  expect_equal(sum(ow$assigned), 255L)
})

test_that("myeloma MINDENOM=30 assigned count matches classified n=186", {
  skip_if_slow_tests_disabled("fixture-myeloma-obs-weights")

  d    <- .ow_load_myeloma()
  tree <- .ow_myeloma_fit(30L)
  X    <- d[, .ow_myeloma_attrs]
  y    <- as.integer(d$V1)
  ow   <- cta_observation_weights(tree, X, y, adjusted = TRUE)

  # MINDENOM=30 root=V17; obs with V17=-9 are unroutable.
  # Canon: 186 obs classified.
  expect_equal(sum(ow$assigned), 186L)
})

test_that("myeloma MINDENOM=56 all unassigned (no-tree)", {
  skip_if_slow_tests_disabled("fixture-myeloma-obs-weights")

  d    <- .ow_load_myeloma()
  tree <- .ow_myeloma_fit(56L)
  X    <- d[, .ow_myeloma_attrs]
  y    <- as.integer(d$V1)
  ow   <- cta_observation_weights(tree, X, y, adjusted = TRUE)

  expect_true(all(!ow$assigned))
  expect_true(all(is.na(ow$endpoint_id)))
})
