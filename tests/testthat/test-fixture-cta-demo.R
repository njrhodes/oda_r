###############################################################################
# test-fixture-cta-demo.R
#
# Gold regression test for oda_cta_fit() ENUMERATE against
# tests/testthat/fixtures/CTA_DEMO_output.txt (MegaODA.exe canonical output).
#
# Canon command:
#   OPEN CTA_DEMO.CSV;
#   VARS V1 TO V6; CLASS V1; ATTRIBUTE V2 V3 V4 V5 V6;
#   MC ITER 25000 CUTOFF 0.05 STOP 99.9;
#   MINDENOM 1; PRUNE 0.05; ENUMERATE; LOO STABLE;
#
# ENUMERATE output from fixture (Tree 1, the winner with Best ESS 75.17%):
#   Node 1: V2  cut=4.5  n=200  ESS=52.63%  STABLE
#   Node 2: V6  cut=7.5  n=110  ESS=38.39%  STABLE
#   Node 3: V2  cut=6.5  n=90   ESS=38.44%  STABLE
#   Full-tree confusion: class1={103 correct, 21 wrong}
#                        class2={6 wrong, 70 correct}
#   Overall ESS = 75.17%
###############################################################################

# ---- helpers ----------------------------------------------------------------

load_cta_demo <- function() {
  f <- tryCatch(testthat::test_path("fixtures/cta_demo/CTA_DEMO.CSV"), error = function(e) "")
  if (!nzchar(f) || !file.exists(f)) skip("CTA_DEMO.CSV not found")
  read.csv(f, header = FALSE,
           col.names = c("V1", "V2", "V3", "V4", "V5", "V6"))
}

.cta_demo_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      d <- load_cta_demo()
      X <- d[, c("V2", "V3", "V4", "V5", "V6")]
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
        mc_stopup   = NA_real_,
        mc_seed     = NULL,
        loo         = "stable",
        attr_names  = c("V2", "V3", "V4", "V5", "V6")
      )
    }
    fit
  }
})

# ---- top-three-node structure -----------------------------------------------

test_that("cta-demo gold: enumerated root is V2", {
  tree <- .cta_demo_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_false(isTRUE(root$leaf))
  expect_equal(root$attribute, "V2")
})

test_that("cta-demo gold: root cut value is 4.5", {
  tree <- .cta_demo_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$rule$cut_value, 4.5)
})

test_that("cta-demo gold: root n_obs = 200", {
  tree <- .cta_demo_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$n_obs, 200L)
})

test_that("cta-demo gold: root LOO status is STABLE", {
  tree <- .cta_demo_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$loo_status, "STABLE")
})

test_that("cta-demo gold: depth-2 nodes are V6 (n=110) and V2 (n=90), both STABLE", {
  tree <- .cta_demo_fit()
  tbl  <- cta_node_table(tree)
  d2   <- tbl[tbl$depth == 2 & !tbl$leaf, ]
  expect_equal(nrow(d2), 2L,
               label = "two split nodes at depth 2")
  expect_equal(sort(d2$n_obs), c(90L, 110L),
               label = "depth-2 obs counts")
  expect_equal(sort(d2$attribute), c("V2", "V6"),
               label = "depth-2 attributes")
  expect_true(all(d2$loo_status == "STABLE"),
              info = paste("non-STABLE:", paste(d2$loo_status[d2$loo_status != "STABLE"], collapse = ",")))
})

# ---- full-tree performance --------------------------------------------------

test_that("cta-demo gold: full-tree training confusion matches [[103,21],[6,70]]", {
  d    <- load_cta_demo()
  tree <- .cta_demo_fit()
  X    <- d[, c("V2", "V3", "V4", "V5", "V6")]
  y    <- d$V1
  preds <- predict(tree, X)

  conf <- matrix(0L, 2L, 2L)
  for (i in seq_along(y)) {
    a <- y[i]; p <- preds[i]
    if (a %in% 1:2 && p %in% 1:2) conf[a, p] <- conf[a, p] + 1L
  }

  expect_equal(conf[1, 1], 103L, label = "TP class 1")
  expect_equal(conf[1, 2],  21L, label = "FN class 1")
  expect_equal(conf[2, 1],   6L, label = "FP class 2")
  expect_equal(conf[2, 2],  70L, label = "TP class 2")
})

test_that("cta-demo gold: overall training ESS = 75.17% (tolerance 0.5%)", {
  d    <- load_cta_demo()
  tree <- .cta_demo_fit()
  X    <- d[, c("V2", "V3", "V4", "V5", "V6")]
  y    <- d$V1
  preds <- predict(tree, X)

  tp1 <- sum(preds == 1 & y == 1); fn1 <- sum(preds == 2 & y == 1)
  tp2 <- sum(preds == 2 & y == 2); fn2 <- sum(preds == 1 & y == 2)
  pac1 <- tp1 / (tp1 + fn1)
  pac2 <- tp2 / (tp2 + fn2)
  ess  <- 100 * ((pac1 + pac2) / 2 - 0.5) / 0.5

  expect_equal(round(ess, 1), 75.2, tolerance = 0.5,
               label = "overall ESS")
})
