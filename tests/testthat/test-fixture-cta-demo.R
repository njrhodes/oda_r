###############################################################################
# test-fixture-cta-demo.R — gold regression vs CTA_DEMO_output.txt (MegaODA)
#
# Canon: OPEN CTA_DEMO.CSV; VARS V1-V6; CLASS V1; ATTRIBUTE V2-V6;
#   MC ITER 25000 CUTOFF 0.05 STOP 99.9; MINDENOM 1; PRUNE 0.05;
#   ENUMERATE; LOO STABLE;
#
# Winning tree (ENUMERATE Tree 1, Best ESS 75.17%):
#   Node 1: V2 cut=4.5 n=200 ESS=52.63% STABLE
#   Node 2: V6 cut=7.5 n=110 ESS=38.39% STABLE
#   Node 3: V2 cut=6.5 n= 90 ESS=38.44% STABLE
#   Confusion: [[103,21],[6,70]]   Overall ESS = 75.17%
###############################################################################

# ---- fixture helpers ---------------------------------------------------------

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
      fit <<- oda_cta_fit(
        X           = d[, c("V2", "V3", "V4", "V5", "V6")],
        y           = d$V1,
        priors_on   = TRUE,
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
        attr_names  = c("V2", "V3", "V4", "V5", "V6")
      )
    }
    fit
  }
})

# =============================================================================
# 1. Root node contract
# =============================================================================

test_that("cta-demo gold: root attr=V2, cut=4.5, n=200, LOO=STABLE, ESS>=50% (uniform-weight bypass confirmed)", {
  # No WEIGHT command -> w defaults to all-1 -> any(w != w[1]) is FALSE ->
  # CTA ordered-scan path is bypassed; generic ODA path used throughout.
  # Observable bypass contract: root=V2, cut=4.5, ESS=52.63% matches CTA.exe.
  skip_if_not_smoke("fixture-cta-demo")
  tree <- .cta_demo_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_false(isTRUE(root$leaf))
  expect_equal(root$attribute,      "V2",
    label = "root=V2 confirms uniform-weight CTA path bypass")
  expect_equal(root$rule$cut_value,  4.5)
  expect_equal(root$n_obs,          200L)
  expect_equal(root$loo_status,    "STABLE")
  expect_gte(root$ess,              50.0)
})

# =============================================================================
# 2. Tree structure
# =============================================================================

test_that("cta-demo gold: 8 splits with EXE-canonical attrs, obs, depth breakdowns, all STABLE", {
  skip_if_not_smoke("fixture-cta-demo")
  tree   <- .cta_demo_fit()
  tbl    <- cta_node_table(tree)
  splits <- tbl[!tbl$leaf, ]

  expect_equal(nrow(splits), 8L)
  expect_equal(sort(splits$attribute),
               sort(c("V2", "V2", "V3", "V4", "V4", "V5", "V6", "V6")))
  expect_equal(sort(splits$n_obs),
               sort(c(200L, 110L, 90L, 81L, 47L, 29L, 37L, 25L)))
  expect_true(all(splits$loo_status == "STABLE"),
              info = paste("non-STABLE:", paste(splits$loo_status[splits$loo_status != "STABLE"], collapse = ",")))
  expect_true(all(is.finite(splits$p_mc)))

  d2 <- tbl[tbl$depth == 2L & !tbl$leaf, ]
  expect_equal(nrow(d2), 2L)
  expect_equal(sort(d2$attribute), c("V2", "V6"))
  expect_equal(sort(d2$n_obs),     c(90L, 110L))
  expect_true(all(d2$loo_status == "STABLE"))

  d3 <- tbl[tbl$depth == 3L & !tbl$leaf, ]
  expect_equal(nrow(d3), 3L)
  expect_equal(sort(d3$attribute), c("V3", "V5", "V6"))
  expect_equal(sort(d3$n_obs),     c(29L, 47L, 81L))

  d4 <- tbl[tbl$depth == 4L & !tbl$leaf, ]
  expect_equal(nrow(d4), 2L)
  expect_true(all(d4$attribute == "V4"))
  expect_equal(sort(d4$n_obs), c(25L, 37L))
})

# =============================================================================
# 3. Training performance
# =============================================================================

test_that("cta-demo gold: confusion [[103,21],[6,70]] and overall ESS ≈ 75.17%", {
  skip_if_not_smoke("fixture-cta-demo")
  d     <- load_cta_demo()
  tree  <- .cta_demo_fit()
  preds <- predict(tree, d[, c("V2", "V3", "V4", "V5", "V6")])
  y     <- d$V1

  conf <- matrix(0L, 2L, 2L)
  for (i in seq_along(y)) {
    a <- y[i]; p <- preds[i]
    if (a %in% 1:2 && p %in% 1:2) conf[a, p] <- conf[a, p] + 1L
  }
  expect_equal(conf[1, 1], 103L, label = "TP class 1")
  expect_equal(conf[1, 2],  21L, label = "FN class 1")
  expect_equal(conf[2, 1],   6L, label = "FP class 2")
  expect_equal(conf[2, 2],  70L, label = "TP class 2")

  pac1 <- conf[1,1] / (conf[1,1] + conf[1,2])
  pac2 <- conf[2,2] / (conf[2,2] + conf[2,1])
  ess  <- 100 * ((pac1 + pac2) / 2 - 0.5) / 0.5
  expect_equal(round(ess, 1), 75.2, tolerance = 0.5, label = "overall ESS")
})
