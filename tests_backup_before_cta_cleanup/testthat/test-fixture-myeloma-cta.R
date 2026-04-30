###############################################################################
# test-fixture-myeloma-cta.R
#
# Gold regression for oda_cta_fit() against CTA.exe output (MODEL1.TXT).
#
# Source: fixtures/myeloma/data.txt  (256 obs, space-delimited, no header)
#         fixtures/myeloma/cta.pgm   (exact command file used)
#
# MegaODA/CTA command summary:
#   VARS V1 TO V19; CLASS V1; ATTRIBUTE V4 V9 V11 V12 V14 TO V19;
#   EX V2=0; MISSING ALL (-9); WEIGHT V2;
#   MC ITER 5000 CUTOFF 0.05 STOP 99.9; MINDENOM 1; PRUNE 0.05;
#   ENUMERATE; LOO STABLE;
#
# Canonical result: ENUMERATED tree (our oda_cta_fit always enumerates).
#
#   NODE  ATTR  OBS   WESS     LOO     CUT
#      1  V14   255  14.06%  STABLE   <=0.5-->0; >0.5-->1
#      2  V15   211  20.71%  STABLE   <=0.5-->0; >0.5-->1
#
# Training confusion (raw counts, rows=actual, cols=predicted):
#   [[146, 40],
#    [ 36, 33]]
# Weighted ESS (using V2 as case weight) = 27.69%
#
# p-values are NOT asserted: CTA.exe uses system-time MC seeds.
# Structure, cut values, ESS, confusion, and LOO status are deterministic.
###############################################################################

# ---- Fixture helpers --------------------------------------------------------

.myeloma_fixture_path <- function(f) {
  tryCatch(testthat::test_path(file.path("fixtures","myeloma", f)),
           error = function(e) "")
}

.myeloma_fixtures_ok <- function() {
  p <- .myeloma_fixture_path("data.txt")
  nzchar(p) && file.exists(p)
}

.load_myeloma <- function() {
  f <- .myeloma_fixture_path("data.txt")
  d <- read.table(f, header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))   # V1..V19
  d <- d[d$V2 != 0, ]                          # EX V2=0
  d
}

# Attributes used (MegaODA ATTRIBUTE command):
.myeloma_attr_names <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

# Fit the enumerated CTA tree once; cache across all tests in this file.
.myeloma_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      if (!.myeloma_fixtures_ok()) return(NULL)
      d  <- .load_myeloma()
      X  <- d[, .myeloma_attr_names]
      y  <- d$V1
      w  <- d$V2
      fit <<- oda_cta_fit(
        X           = X,
        y           = y,
        w           = w,
        priors_on   = TRUE,
        miss_codes  = -9,
        alpha_split = 0.05,
        mindenom    = 1L,
        prune_alpha = 0.05,
        max_depth   = 20L,
        ess_min     = 0,
        mc_iter     = 5000L,
        mc_target   = 0.05,
        mc_stop     = 99.9,
        mc_stopup   = 20,
        mc_seed     = NULL,
        loo         = "stable",
        attr_names  = .myeloma_attr_names
      )
    }
    fit
  }
})

# =============================================================================
# Structure tests
# =============================================================================

test_that("myeloma gold: root attribute is V14", {
  if (!.myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .myeloma_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_false(isTRUE(root$leaf))
  expect_equal(root$attribute, "V14")
})

test_that("myeloma gold: root cut value is 0.5", {
  if (!.myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .myeloma_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$rule$cut_value, 0.5)
})

test_that("myeloma gold: root n_obs = 255", {
  if (!.myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .myeloma_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$n_obs, 255L)
})

test_that("myeloma gold: root WESS within 0.1% of 14.06%", {
  if (!.myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .myeloma_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$ess, 14.06, tolerance = 0.1, label = "root WESS")
})

test_that("myeloma gold: root LOO status is STABLE", {
  if (!.myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .myeloma_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$loo_status, "STABLE")
})

test_that("myeloma gold: depth-2 split is V15 with cut 0.5", {
  if (!.myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .myeloma_fit()
  tbl  <- cta_node_table(tree)
  d2   <- tbl[tbl$depth == 2L & !tbl$leaf, ]
  expect_equal(nrow(d2), 1L, label = "one depth-2 split node")
  expect_equal(d2$attribute[[1L]], "V15")
  # Cut value from node object
  nid <- d2$node_id[[1L]]
  nd  <- tree$nodes[[nid]]
  expect_equal(nd$rule$cut_value, 0.5)
})

test_that("myeloma gold: node 2 n_obs = 211", {
  if (!.myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .myeloma_fit()
  tbl  <- cta_node_table(tree)
  d2   <- tbl[tbl$depth == 2L & !tbl$leaf, ]
  expect_equal(d2$n_obs[[1L]], 211L)
})

test_that("myeloma gold: node 2 WESS within 0.1% of 20.71%", {
  if (!.myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .myeloma_fit()
  tbl  <- cta_node_table(tree)
  d2   <- tbl[tbl$depth == 2L & !tbl$leaf, ]
  nid  <- d2$node_id[[1L]]
  nd   <- tree$nodes[[nid]]
  expect_equal(nd$ess, 20.71, tolerance = 0.1, label = "node-2 WESS")
})

test_that("myeloma gold: node 2 LOO status is STABLE", {
  if (!.myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .myeloma_fit()
  tbl  <- cta_node_table(tree)
  d2   <- tbl[tbl$depth == 2L & !tbl$leaf, ]
  nid  <- d2$node_id[[1L]]
  nd   <- tree$nodes[[nid]]
  expect_equal(nd$loo_status, "STABLE")
})

test_that("myeloma gold: all split nodes have LOO = STABLE", {
  if (!.myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .myeloma_fit()
  tbl  <- cta_node_table(tree)
  split_loo <- tbl$loo_status[!tbl$leaf]
  expect_true(all(split_loo == "STABLE"),
              info = paste("non-STABLE:", paste(split_loo[split_loo != "STABLE"], collapse = ",")))
})

# =============================================================================
# Prediction and confusion tests
# =============================================================================

test_that("myeloma gold: full-tree training confusion matches [[146,40],[36,33]]", {
  if (!.myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  d    <- .load_myeloma()
  tree <- .myeloma_fit()
  X    <- d[, .myeloma_attr_names]
  y    <- as.integer(d$V1)
  preds <- as.integer(predict(tree, X))

  conf <- matrix(0L, 2L, 2L)   # rows = actual {0,1}+1; cols = predicted
  for (i in seq_along(y)) {
    a <- y[i] + 1L; p <- preds[i] + 1L
    if (a %in% 1:2 && p %in% 1:2) conf[a, p] <- conf[a, p] + 1L
  }

  expect_equal(conf[1L, 1L], 146L, label = "TN (actual 0, pred 0)")
  expect_equal(conf[1L, 2L],  40L, label = "FP (actual 0, pred 1)")
  expect_equal(conf[2L, 1L],  36L, label = "FN (actual 1, pred 0)")
  expect_equal(conf[2L, 2L],  33L, label = "TP (actual 1, pred 1)")
})

test_that("myeloma gold: weighted ESS within 0.5% of 27.69%", {
  if (!.myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  d     <- .load_myeloma()
  tree  <- .myeloma_fit()
  X     <- d[, .myeloma_attr_names]
  y     <- as.integer(d$V1)
  w     <- d$V2
  preds <- as.integer(predict(tree, X))

  # Weighted PAC per class (case-weight-adjusted sensitivity/specificity)
  m0    <- (y == 0L)
  m1    <- (y == 1L)
  wpac0 <- sum(w[m0 & preds == 0L]) / sum(w[m0])
  wpac1 <- sum(w[m1 & preds == 1L]) / sum(w[m1])
  wess  <- 2 * ((wpac0 + wpac1) / 2 - 0.5) * 100

  expect_equal(wess, 27.69, tolerance = 0.5, label = "weighted ESS")
})
