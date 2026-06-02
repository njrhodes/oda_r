###############################################################################
# test-cta-observation-weights.R — cta_observation_weights()
#
# Per-observation propensity weight assignment.
###############################################################################

# ---- Helpers ----------------------------------------------------------------

.ow_stump_fit <- function() {
  suppressMessages(
    oda_cta_fit(data.frame(x = 1:8),
                c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
                mindenom = 2L, mc_iter = 300L, mc_seed = 1L, loo = "off")
  )
}

.ow_no_tree_fit <- function() {
  suppressMessages(
    oda_cta_fit(data.frame(x = 1:8),
                c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
                mindenom = 999L, mc_iter = 100L, mc_seed = 1L, loo = "off")
  )
}

# =============================================================================
# Contract tests
# =============================================================================

test_that("ow: schema + no-tree — nrow, columns, no-tree all unassigned/NA", {
  tree <- .ow_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  X <- data.frame(x = 1:8)
  y <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  ow <- cta_observation_weights(tree, X, y)
  # Row count
  expect_equal(nrow(ow), 8L)
  # Required columns
  expected_cols <- c("row_id", "actual_class", "endpoint_node_id",
                     "endpoint_id", "target_class", "propensity_weight",
                     "adjusted_propensity_weight", "undefined_empirical",
                     "perfectly_predicted_endpoint", "adjusted", "assigned")
  expect_true(all(expected_cols %in% names(ow)))
  # row_id is sequential 1:n
  expect_identical(ow$row_id, 1:8)
  # actual_class is character coercion of y
  expect_equal(ow$actual_class, as.character(y))
  # No-tree: all unassigned, weight columns NA
  nt <- .ow_no_tree_fit()
  ow_nt <- cta_observation_weights(nt, X, y)
  expect_equal(nrow(ow_nt), 8L)
  expect_true(all(!ow_nt$assigned))
  expect_true(all(is.na(ow_nt$propensity_weight)))
  expect_true(all(is.na(ow_nt$adjusted_propensity_weight)))
})

test_that("ow: stump contract — assigned, missing_action, NA y, target_class annotation", {
  tree <- .ow_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  X <- data.frame(x = 1:8)
  y <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  # All 8 obs assigned when using adjusted=FALSE on perfectly separated data
  ow <- cta_observation_weights(tree, X, y, adjusted = FALSE)
  expect_true(all(ow$assigned))
  # missing_action='na' → row 1 unassigned; 'majority' → row 1 assigned
  X_miss <- data.frame(x = c(NA, 1, 2, 3, 5, 6, 7, 8))
  y_miss  <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  ow_na  <- cta_observation_weights(tree, X_miss, y_miss, missing_action = "na")
  ow_maj <- cta_observation_weights(tree, X_miss, y_miss, missing_action = "majority")
  expect_true(is.na(ow_na$endpoint_id[1L]))
  expect_false(is.na(ow_maj$endpoint_id[1L]))
  expect_equal(nrow(ow_na), 8L); expect_equal(nrow(ow_maj), 8L)
  # NA y → assigned=FALSE and NA weight; other rows still assigned
  y_na <- c(NA_integer_, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  ow_yn <- cta_observation_weights(tree, X, y_na)
  expect_false(ow_yn$assigned[1L])
  expect_true(is.na(ow_yn$propensity_weight[1L]))
  expect_true(all(ow_yn$assigned[-1L]))
  # target_class annotation: all obs remain assigned; column is integer
  ow1 <- cta_observation_weights(tree, X, y, target_class = 1L, adjusted = FALSE)
  expect_true(all(ow1$assigned))
  expect_type(ow1$target_class, "integer")
  expect_true(all(ow1$target_class[ow1$assigned] == 1L))
})

test_that("ow: errors — y length mismatch, non-cta_tree, unmatched obs warning", {
  tree <- .ow_stump_fit()
  skip_if(isTRUE(tree$no_tree), "mc sampling missed")
  X <- data.frame(x = 1:8)
  expect_error(cta_observation_weights(tree, X, y = 1:5), regexp = "length\\(y\\)")
  expect_error(cta_observation_weights(list(), data.frame(x=1:4), y=1:4), regexp = "cta_tree")
  # Class "2" not in fitted tree → warning about unmatched obs
  y_bad <- c(2L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  expect_warning(cta_observation_weights(tree, X, y_bad), regexp = "could not be matched")
  # matrix newdata coerced to data.frame without error
  X_mat <- matrix(1:8, ncol = 1L, dimnames = list(NULL, "x"))
  expect_no_error(cta_observation_weights(tree, X_mat, y = rep(0L, 8L)))
})

# =============================================================================
# Smoke: myeloma canon
# =============================================================================

.ow_load_myeloma <- function() {
  raw <- read.table(testthat::test_path("fixtures", "myeloma", "data.txt"), header = FALSE)
  names(raw) <- paste0("V", seq_len(ncol(raw)))
  raw[raw$V2 != 0, ]
}
.ow_myeloma_attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.ow_myeloma_fit <- local({
  fits <- list()
  function(mindenom) {
    key <- as.character(mindenom)
    if (is.null(fits[[key]])) {
      d <- .ow_load_myeloma()
      fits[[key]] <<- suppressMessages(oda_cta_fit(
        X = d[, .ow_myeloma_attrs], y = as.integer(d$V1), w = d$V2,
        mindenom = mindenom, mc_iter = 25000L, mc_seed = 42L,
        miss_codes = -9, loo = "stable", attr_names = .ow_myeloma_attrs
      ))
    }
    fits[[key]]
  }
})

test_that("ow: myeloma — assigned counts and no-tree for MINDENOM=56/1/30", {
  skip_if_slow_tests_disabled("fixture-myeloma-obs-weights")
  d <- .ow_load_myeloma()
  X <- d[, .ow_myeloma_attrs]
  y <- as.integer(d$V1)

  # MINDENOM=56: no-tree → all unassigned, all NA endpoint_id
  ow56 <- cta_observation_weights(.ow_myeloma_fit(56L), X, y, adjusted = TRUE)
  expect_true(all(!ow56$assigned))
  expect_true(all(is.na(ow56$endpoint_id)))

  # MINDENOM=1: classified = 255; nrow = nrow(X); assigned obs have non-NA weights
  ow1 <- cta_observation_weights(.ow_myeloma_fit(1L), X, y, adjusted = TRUE)
  expect_equal(nrow(ow1), nrow(X))
  expect_identical(ow1$row_id, seq_len(nrow(X)))
  expect_equal(sum(ow1$assigned), 255L)
  assigned1 <- ow1[ow1$assigned, ]
  non_undef1 <- assigned1[!isTRUE(assigned1$undefined_empirical), ]
  expect_true(all(!is.na(non_undef1$propensity_weight)))

  # MINDENOM=30: classified = 186
  ow30 <- cta_observation_weights(.ow_myeloma_fit(30L), X, y, adjusted = TRUE)
  expect_equal(sum(ow30$assigned), 186L)
})
