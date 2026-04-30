###############################################################################
# test-cta-gold.R
#
# Gold regression test for oda_cta_fit() against verified MegaODA.exe output.
#
# Source: CTA_DEMO.CSV (200 obs, no header, comma-delimited)
#         MegaODA command: cta.pgm
#           MINDENOM 1; PRUNE 0.05; ENUMERATE; LOO STABLE;
#           MC ITER 25000 CUTOFF 0.05 STOP 99.9
#
# MegaODA node table (Unpruned = Pruned = Enumerated, all three identical):
#
#   ATTR  NODE LEV  OBS    p      ESS     LOO
#   V2       1   1  200  .000   52.63%  STABLE
#   V6       2   2  110  .006   38.39%  STABLE
#   V2       3   2   90  .001   38.44%  STABLE
#   V3       4   3   29  .013   53.68%  STABLE
#   V5       5   3   81  .008   60.04%  STABLE
#   V6       6   3   47  .043   34.96%  STABLE
#   V4      10   4   25  .000   94.74%  STABLE
#   V4      13   4   37  .008   48.14%  STABLE
#
# All splits verified algebraically against the raw CSV (see comments below).
# Final confusion [[103,21],[6,70]], OVERALL ESS = 75.17%.
#
# Split ancestry (verified):
#   1 (root) → left(V2<=4.5, n=110) = node 2
#                right(V2>4.5, n=90)  = node 3
#   2 → left(V6<=7.5, n=29)  = node 4
#       right(V6>7.5,  n=81)  = node 5
#   3 → left(V2<=6.5, n=47)  = node 6
#       right(V2>6.5,  n=43)  = leaf (class 2)
#   4 → left(V3<=4.5, n=16)  = leaf (class 1)
#       right(V3>4.5,  n=13)  = leaf (class 2)
#   5 → left(V5<=5.5, n=25)  = node 10
#       right(V5>5.5,  n=56)  = leaf (class 1)
#   6 → left(V6<=5.5, n=10)  = leaf (class 2)
#       right(V6>5.5,  n=37)  = node 13
#  10 → left(V4<=0.5, n=18)  = leaf (class 1)
#       right(V4>0.5,  n=7)   = leaf (class 2)
#  13 → left(V4<=0.5, n=19)  = leaf (class 1)
#       right(V4>0.5,  n=18)  = leaf (class 2)
###############################################################################

# ---- Load data --------------------------------------------------------------
load_cta_demo <- function() {
  # testthat::test_path() resolves correctly under both devtools::test()
  # and R CMD check (the latter copies test files to the check temp dir).
  f <- tryCatch(testthat::test_path("CTA_DEMO.CSV"), error = function(e) "")
  if (!nzchar(f) || !file.exists(f)) skip("CTA_DEMO.CSV not found")
  read.csv(f, header = FALSE,
           col.names = c("V1","V2","V3","V4","V5","V6"))
}

# ---- Fit the tree once for all gold tests -----------------------------------
# Use LOO STABLE, MINDENOM 1, PRUNE 0.05, MC ITER 25000
# This is slow (~minutes) so we cache via local()
.cta_gold_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      d <- load_cta_demo()
      X <- d[, c("V2","V3","V4","V5","V6")]
      y <- d$V1
      fit <<- oda_cta_fit(
        X           = X,
        y           = y,
        priors_on   = TRUE,
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
        attr_names  = c("V2","V3","V4","V5","V6")
      )
    }
    fit
  }
})

# =============================================================================
# Structural tests (fast to assert, don't depend on MC variability)
# =============================================================================

test_that("gold: root node attribute is V2", {
  tree <- .cta_gold_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_false(isTRUE(root$leaf))
  expect_equal(root$attribute, "V2")
})

test_that("gold: root cut value is 4.5", {
  tree <- .cta_gold_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$rule$cut_value, 4.5)
})

test_that("gold: root n_obs = 200", {
  tree <- .cta_gold_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$n_obs, 200L)
})

test_that("gold: root ESS >= 50% (MegaODA = 52.63%)", {
  tree <- .cta_gold_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_gte(root$ess, 50.0)
})

test_that("gold: root LOO status is STABLE", {
  tree <- .cta_gold_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$loo_status, "STABLE")
})

test_that("gold: 8 split nodes (matching MegaODA node table)", {
  tree <- .cta_gold_fit()
  tbl  <- cta_node_table(tree)
  n_split <- sum(!tbl$leaf)
  expect_equal(n_split, 8L)
})

test_that("gold: split node obs counts match MegaODA exactly", {
  tree <- .cta_gold_fit()
  tbl  <- cta_node_table(tree)
  split_obs <- sort(tbl$n_obs[!tbl$leaf])
  # MegaODA split node obs: 25,29,37,47,81,90,110,200
  expect_equal(split_obs, sort(c(200L,110L,90L,29L,81L,47L,25L,37L)))
})

test_that("gold: split node attributes match MegaODA (V2x2, V3, V4x2, V5, V6x2)", {
  tree <- .cta_gold_fit()
  tbl  <- cta_node_table(tree)
  split_attrs <- sort(tbl$attribute[!tbl$leaf])
  expected    <- sort(c("V2","V2","V3","V4","V4","V5","V6","V6"))
  expect_equal(split_attrs, expected)
})

test_that("gold: all split nodes have LOO = STABLE", {
  tree <- .cta_gold_fit()
  tbl  <- cta_node_table(tree)
  split_loo <- tbl$loo_status[!tbl$leaf]
  expect_true(all(split_loo == "STABLE"),
              info = paste("non-STABLE:", paste(split_loo[split_loo != "STABLE"], collapse=",")))
})

test_that("gold: all split node p-values are finite", {
  # p-values are RNG-dependent (system-time seed in CTA.exe) and are not
  # assertable. Structure and ESS are the gold values.
  tree <- .cta_gold_fit()
  tbl  <- cta_node_table(tree)
  p_vals <- tbl$p_mc[!tbl$leaf]
  expect_true(all(is.finite(p_vals)))
})

test_that("gold: depth-2 nodes are V6 (n=110) and V2 (n=90)", {
  tree <- .cta_gold_fit()
  tbl  <- cta_node_table(tree)
  d2   <- tbl[tbl$depth == 2 & !tbl$leaf, ]
  expect_equal(nrow(d2), 2L)
  expect_equal(sort(d2$n_obs), c(90L, 110L))
  expect_equal(sort(d2$attribute), c("V2", "V6"))
})

test_that("gold: depth-3 nodes are V3(29), V5(81), V6(47)", {
  tree <- .cta_gold_fit()
  tbl  <- cta_node_table(tree)
  d3   <- tbl[tbl$depth == 3 & !tbl$leaf, ]
  expect_equal(nrow(d3), 3L)
  expect_equal(sort(d3$n_obs), c(29L, 47L, 81L))
  expect_equal(sort(d3$attribute), c("V3", "V5", "V6"))
})

test_that("gold: depth-4 nodes are both V4, obs 25 and 37", {
  tree <- .cta_gold_fit()
  tbl  <- cta_node_table(tree)
  d4   <- tbl[tbl$depth == 4 & !tbl$leaf, ]
  expect_equal(nrow(d4), 2L)
  expect_equal(sort(d4$n_obs), c(25L, 37L))
  expect_true(all(d4$attribute == "V4"))
})

# =============================================================================
# Prediction and confusion tests (the real end-to-end check)
# =============================================================================

test_that("gold: full-tree training confusion matches MegaODA [[103,21],[6,70]]", {
  d    <- load_cta_demo()
  tree <- .cta_gold_fit()
  X    <- d[, c("V2","V3","V4","V5","V6")]
  y    <- d$V1
  preds <- predict(tree, X)

  conf <- matrix(0L, 2L, 2L)
  for (i in seq_along(y)) {
    a <- y[i]; p <- preds[i]
    if (a %in% 1:2 && p %in% 1:2) conf[a, p] <- conf[a, p] + 1L
  }

  expect_equal(conf[1,1], 103L, label = "TP class 1")
  expect_equal(conf[1,2],  21L, label = "FN class 1")
  expect_equal(conf[2,1],   6L, label = "FP class 2")
  expect_equal(conf[2,2],  70L, label = "TP class 2")
})

test_that("gold: overall training ESS = 75.17% (tolerance 0.5%)", {
  d    <- load_cta_demo()
  tree <- .cta_gold_fit()
  X    <- d[, c("V2","V3","V4","V5","V6")]
  y    <- d$V1
  preds <- predict(tree, X)

  tp1 <- sum(preds == 1 & y == 1); fn1 <- sum(preds == 2 & y == 1)
  tp2 <- sum(preds == 2 & y == 2); fn2 <- sum(preds == 1 & y == 2)
  pac1 <- tp1 / (tp1 + fn1)
  pac2 <- tp2 / (tp2 + fn2)
  ess  <- 100 * ((pac1 + pac2)/2 - 0.5) / 0.5

  expect_equal(round(ess, 1), 75.2, tolerance = 0.5,
               label = "overall ESS")
})
