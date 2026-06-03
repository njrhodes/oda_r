###############################################################################
# test-fixture-myeloma-cta.R — gold regression vs CTA.exe (MODEL1/30/56.TXT)
#
# Source: fixtures/myeloma/data.txt (256 obs, space-delimited, no header)
# Command: EX V2=0; MISSING ALL (-9); WEIGHT V2; MC ITER 5000 STOP 99.9;
#   MINDENOM <1|30|56>; PRUNE 0.05; ENUMERATE; LOO STABLE;
#
# MINDENOM=1 canon (MODEL1.TXT Enumerated):
#   Node 1: V14 cut=0.5 n=255 WESS=14.06% STABLE
#   Node 2: V15 cut=0.5 n=211 WESS=20.71% STABLE
#   Confusion: [[146,40],[36,33]]  WEIGHTED ESS=27.69%
#
# MINDENOM=30 canon (MODEL30.TXT): V17 stump, n=186, confusion [[101,34],[30,21]]
#   OVERALL ESS=15.99%, WEIGHTED ESS=16.51%
#
# MINDENOM=56 canon (MODEL56.TXT): no admissible root → leaf only.
###############################################################################

# ---- fixture helpers ---------------------------------------------------------

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
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]
}

.myeloma_attr_names <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")

.myeloma_fit1 <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      if (!.myeloma_fixtures_ok()) return(NULL)
      d  <- .load_myeloma()
      fit <<- oda_cta_fit(
        X           = d[, .myeloma_attr_names],
        y           = d$V1,
        w           = d$V2,
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
        mc_seed     = NULL,
        loo         = "stable",
        attr_names  = .myeloma_attr_names
      )
    }
    fit
  }
})

.myeloma_fit30 <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      if (!.myeloma_fixtures_ok()) return(NULL)
      d  <- .load_myeloma()
      fit <<- oda_cta_fit(
        X           = d[, .myeloma_attr_names],
        y           = d$V1,
        w           = d$V2,
        priors_on   = TRUE,
        miss_codes  = -9,
        alpha_split = 0.05,
        mindenom    = 30L,
        prune_alpha = 0.05,
        max_depth   = 20L,
        ess_min     = 0,
        mc_iter     = 5000L,
        mc_target   = 0.05,
        mc_stop     = 99.9,
        mc_seed     = NULL,
        loo         = "stable",
        attr_names  = .myeloma_attr_names
      )
    }
    fit
  }
})

# =============================================================================
# 1. MINDENOM=1 — tree structure
# =============================================================================

test_that("myeloma MINDENOM=1: root V14/cut0.5/n255/WESS≈14.06%/STABLE; node2 V15/n211/WESS≈20.71%/STABLE; all STABLE", {
  skip_if_slow_tests_disabled("fixture-myeloma-cta")
  if (!.myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  tree <- .myeloma_fit1()
  if (is.null(tree)) skip("myeloma fit returned NULL")

  root <- tree$nodes[[tree$root_id]]
  expect_false(isTRUE(root$leaf))
  expect_equal(root$attribute,      "V14")
  expect_equal(root$rule$cut_value,  0.5)
  expect_equal(root$n_obs,          255L)
  expect_equal(root$ess, 14.06, tolerance = 0.1, label = "root WESS")
  expect_equal(root$loo_status,    "STABLE")

  tbl <- cta_node_table(tree)
  d2  <- tbl[tbl$depth == 2L & !tbl$leaf, ]
  expect_equal(nrow(d2), 1L)
  nid <- d2$node_id[[1L]]
  nd2 <- tree$nodes[[nid]]
  expect_equal(nd2$attribute,        "V15")
  expect_equal(nd2$rule$cut_value,    0.5)
  expect_equal(nd2$n_obs,           211L)
  expect_equal(nd2$ess, 20.71, tolerance = 0.1, label = "node-2 WESS")
  expect_equal(nd2$loo_status,      "STABLE")

  split_loo <- tbl$loo_status[!tbl$leaf]
  expect_true(all(split_loo == "STABLE"),
              info = paste("non-STABLE:", paste(split_loo[split_loo != "STABLE"], collapse = ",")))
})

# =============================================================================
# 2. MINDENOM=1 — training performance
# =============================================================================

test_that("myeloma MINDENOM=1: confusion [[146,40],[36,33]] and WEIGHTED ESS ≈ 27.69%", {
  skip_if_slow_tests_disabled("fixture-myeloma-cta")
  if (!.myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  d    <- .load_myeloma()
  tree <- .myeloma_fit1()
  if (is.null(tree)) skip("myeloma fit returned NULL")

  X     <- d[, .myeloma_attr_names]
  y     <- as.integer(d$V1)
  w     <- as.numeric(d$V2)
  preds <- as.integer(predict(tree, X))

  conf <- matrix(0L, 2L, 2L)
  for (i in seq_along(y)) {
    a <- y[i] + 1L; p <- preds[i] + 1L
    if (a %in% 1:2 && p %in% 1:2) conf[a, p] <- conf[a, p] + 1L
  }
  expect_equal(conf[1L,1L], 146L, label = "TN")
  expect_equal(conf[1L,2L],  40L, label = "FP")
  expect_equal(conf[2L,1L],  36L, label = "FN")
  expect_equal(conf[2L,2L],  33L, label = "TP")

  m0    <- (y == 0L); m1 <- (y == 1L)
  wpac0 <- sum(w[m0 & preds == 0L]) / sum(w[m0])
  wpac1 <- sum(w[m1 & preds == 1L]) / sum(w[m1])
  wess  <- 2 * ((wpac0 + wpac1) / 2 - 0.5) * 100
  expect_equal(wess, 27.69, tolerance = 0.5, label = "weighted ESS")

  # cta_confusion_table() API: column names, row count, exact values
  # (absorbed from test-cta-confusion-table.R smoke block)
  ct1 <- cta_confusion_table(tree)
  expect_equal(names(ct1), c("actual", "predicted", "n"))
  expect_equal(sum(ct1$n), 255L)
  expect_equal(ct1$n[ct1$actual == 0L & ct1$predicted == 0L], 146L)
  expect_equal(ct1$n[ct1$actual == 0L & ct1$predicted == 1L],  40L)
  expect_equal(ct1$n[ct1$actual == 1L & ct1$predicted == 0L],  36L)
  expect_equal(ct1$n[ct1$actual == 1L & ct1$predicted == 1L],  33L)
})

# =============================================================================
# 3. MINDENOM=30 — V17 stump, n=186, confusion, OVERALL and WEIGHTED ESS
# =============================================================================

test_that("myeloma MINDENOM=30: V17 stump, n=186, confusion [[101,34],[30,21]], OVERALL≈15.99%, WEIGHTED≈16.51%", {
  # Path-local scoring: V17 stump (n=186, WESS=16.51%) wins ENUMERATE over
  # V14 stump (n=255, WESS=14.06%). 69 V17-missing obs excluded by NA prediction.
  skip_if_slow_tests_disabled("fixture-myeloma-cta")
  if (!.myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  d    <- .load_myeloma()
  tree <- .myeloma_fit30()
  if (is.null(tree)) skip("myeloma fit returned NULL")

  root <- tree$nodes[[tree$root_id]]
  tbl  <- cta_node_table(tree)
  expect_false(isTRUE(root$leaf))
  expect_equal(root$attribute, "V17",
    label = "MINDENOM=30 root must be V17 (path-local ENUMERATE wins over V14)")
  expect_equal(sum(!tbl$leaf), 1L, label = "stump: exactly one split node")

  X     <- d[, .myeloma_attr_names]
  y     <- as.integer(d$V1)
  w     <- as.numeric(d$V2)
  preds <- predict(tree, X, missing_action = "na")
  ok    <- !is.na(preds)

  expect_equal(sum(ok), 186L, label = "69 V17-missing obs excluded")

  conf <- matrix(0L, 2L, 2L)
  for (i in which(ok)) {
    a <- y[i] + 1L; p <- preds[i] + 1L
    if (a %in% 1:2 && p %in% 1:2) conf[a, p] <- conf[a, p] + 1L
  }
  expect_equal(conf[1L,1L], 101L, label = "TN")
  expect_equal(conf[1L,2L],  34L, label = "FP")
  expect_equal(conf[2L,1L],  30L, label = "FN")
  expect_equal(conf[2L,2L],  21L, label = "TP")

  n0   <- sum(y[ok] == 0L); n1 <- sum(y[ok] == 1L)
  pac0 <- sum(y[ok] == 0L & preds[ok] == 0L) / n0
  pac1 <- sum(y[ok] == 1L & preds[ok] == 1L) / n1
  ess  <- (pac0 + pac1 - 1) * 100
  expect_equal(ess, 15.99, tolerance = 0.1, label = "OVERALL ESS")

  m0    <- ok & y == 0L; m1 <- ok & y == 1L
  wpac0 <- sum(w[m0 & preds == 0L]) / sum(w[m0])
  wpac1 <- sum(w[m1 & preds == 1L]) / sum(w[m1])
  wess  <- (wpac0 + wpac1 - 1) * 100
  expect_equal(wess, 16.51, tolerance = 0.1, label = "WEIGHTED ESS")

  # cta_confusion_table() API: column names, row count, exact values
  # (absorbed from test-cta-confusion-table.R smoke block)
  ct30 <- cta_confusion_table(tree)
  expect_equal(names(ct30), c("actual", "predicted", "n"))
  expect_equal(sum(ct30$n), 186L)
  expect_equal(cta_strata(tree), 2L)
  expect_equal(ct30$n[ct30$actual == 0L & ct30$predicted == 0L], 101L)
  expect_equal(ct30$n[ct30$actual == 0L & ct30$predicted == 1L],  34L)
  expect_equal(ct30$n[ct30$actual == 1L & ct30$predicted == 0L],  30L)
  expect_equal(ct30$n[ct30$actual == 1L & ct30$predicted == 1L],  21L)
})

# =============================================================================
# 4. MINDENOM=56 — no admissible root → no tree
# =============================================================================

test_that("myeloma MINDENOM=56: all candidates fail MINDENOM gate → no tree", {
  skip_if_slow_tests_disabled("fixture-myeloma-cta")
  if (!.myeloma_fixtures_ok()) skip("myeloma fixture files missing")
  d   <- .load_myeloma()
  fit <- oda_cta_fit(
    X           = d[, .myeloma_attr_names],
    y           = d$V1,
    w           = d$V2,
    priors_on   = TRUE,
    miss_codes  = -9,
    alpha_split = 0.05,
    mindenom    = 56L,
    prune_alpha = 0.05,
    max_depth   = 20L,
    ess_min     = 0,
    mc_iter     = 5000L,
    mc_target   = 0.05,
    mc_stop     = 99.9,
    mc_seed     = NULL,
    loo         = "stable",
    attr_names  = .myeloma_attr_names
  )
  tbl <- cta_node_table(fit)
  expect_equal(sum(!tbl$leaf), 0L, label = "no split nodes when MINDENOM=56")

  # cta_confusion_table() API: no-tree -> zero rows with correct column names
  # (absorbed from test-cta-confusion-table.R smoke block)
  ct56 <- cta_confusion_table(fit)
  expect_equal(nrow(ct56), 0L)
  expect_equal(names(ct56), c("actual", "predicted", "n"))
})
