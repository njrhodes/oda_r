###############################################################################
# test-loo-ties.R - ordered-cut LOO equals explicit per-fold refits, including
# folds whose optimum is tied (MAXSENS / SAMPLEREP) and directional rules.
###############################################################################

.loo_refit_ess <- function(x, y, direction) {
  n <- length(y); p <- integer(n)
  for (i in seq_len(n)) {
    f <- oda_fit(x = x[-i], y = y[-i], attr_type = "ordered", priors_on = TRUE,
                 mcarlo = FALSE, loo = "off", direction = direction)
    p[i] <- oda:::oda_rule_predict(x[i], f$rule)
  }
  (mean(p[y == 1] == 1) + mean(p[y == 0] == 0) - 1) * 100
}

test_that("ordered-cut LOO matches per-fold refits on tied folds", {
  for (s in 1:25) {
    set.seed(s)
    n <- 22L
    x <- sample(40:80, n, replace = TRUE)          # coarse integers -> frequent ties
    y <- rbinom(n, 1L, 0.5); if (length(unique(y)) < 2L) next
    for (dir in c("both", "greater", "less")) {
      f <- oda_fit(x = x, y = y, attr_type = "ordered", priors_on = TRUE,
                   mcarlo = FALSE, loo = "off", direction = dir)
      if (!isTRUE(f$ok)) next
      lr <- oda:::oda_loo_for_rule(x = x, y = y, rule = f$rule, attr_type = "ordered",
                                   priors_on = TRUE, direction = dir)
      expect_equal(lr$ess_loo, .loo_refit_ess(x, y, dir), tolerance = 1e-9,
                   label = paste0("seed ", s, " direction ", dir))
    }
  }
})

test_that("tied fold without a refit function makes the shortcut return NULL", {
  x <- c(1, 2, 3, 4, 5, 6); y <- c(0L, 1L, 0L, 1L, 0L, 1L)
  rule <- list(type = "ordered_cut", cut_value = 1.5, direction = "0->1")
  expect_null(oda:::oda_loo_ordered_cut_counts(x, y, rep(1, 6), TRUE, rule))
})

# CTA node LOO refits each fold with the CTA scan (the optimizer that chose the
# node rule).  Anchors from CTA.exe outputs.
test_that("CTA-scan LOO: cta_demo MODEL8 node 6 V6 is STABLE (CTA.exe)", {
  f <- tryCatch(testthat::test_path("fixtures/cta_demo/CTA_DEMO.CSV"), error = function(e) "")
  if (!nzchar(f) || !file.exists(f)) skip("CTA_DEMO.CSV not found")
  d <- read.csv(f, header = FALSE, col.names = c("V1", "V2", "V3", "V4", "V5", "V6"))
  d <- d[d$V2 > 4.5 & d$V2 <= 6.5, ]                 # node 6, n = 47, 5 tied LOO folds
  y <- as.integer(d$V1 == 2L); w <- rep(1, nrow(d))
  s <- oda:::.cta_ordered_scan(d$V6, y, w, TRUE, NULL, 1L)
  l <- oda:::.cta_loo_scan(d$V6, y, w, TRUE, s$rule)
  expect_equal(s$rule$cut_value, 5.5)
  expect_equal(l$ess_loo, s$ess, tolerance = 0.01)
})

test_that("CTA-scan LOO: myeloma Tree 1 node 2 V4 UNSTABLE, V15 STABLE (CTA.exe selects V15)", {
  fp <- tryCatch(testthat::test_path("fixtures", "myeloma", "data.txt"), error = function(e) "")
  if (!nzchar(fp) || !file.exists(fp)) skip("myeloma fixture missing")
  d <- read.table(fp, header = FALSE, col.names = paste0("V", 1:19))
  d <- d[d$V2 != 0 & d$V17 != -9 & d$V17 <= 0.5, ]   # node 2
  for (v in c("V4", "V15")) {
    k <- d[[v]] != -9; x <- d[[v]][k]; y <- as.integer(d$V1[k]); w <- d$V2[k]
    s <- oda:::.cta_ordered_scan(x, y, w, TRUE, NULL, 1L)
    l <- oda:::.cta_loo_scan(x, y, w, TRUE, s$rule)
    if (v == "V4")  expect_gt(abs(l$ess_loo - s$ess), 0.01)        # 34.89 vs 27.29
    if (v == "V15") expect_equal(l$ess_loo, s$ess, tolerance = 0.01) # 17.57
  }
})
