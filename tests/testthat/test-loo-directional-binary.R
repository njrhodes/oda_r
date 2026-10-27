###############################################################################
# test-loo-directional-binary.R - LOO honours a declared `direction`.
#
# Canonical LOO (MPE Ch. 2, pp. 32-33): hold out each observation, refit ODA
# on the rest, classify the held-out observation. With direction "greater" or
# "less" every refit is restricted to one rule direction; for a binary
# attribute that leaves one admissible rule, so ess_loo == ess.
# Ordered-attribute equivalent: test-loo-ties.R.
###############################################################################

.loo_oracle_ess_binary <- function(x, y, direction) {
  n <- length(y); p <- integer(n)
  for (i in seq_len(n)) {
    f <- oda_fit(x = x[-i], y = y[-i], attr_type = "binary", priors_on = TRUE,
                 mcarlo = FALSE, loo = "off", direction = direction)
    stopifnot(isTRUE(f$ok))
    p[i] <- oda:::oda_rule_predict(x[i], f$rule)
  }
  (mean(p[y == 1] == 1) + mean(p[y == 0] == 0) - 1) * 100
}

# x = 1 goes with class 0 (5/6) and x = 0 with class 1 (5/6): the hypothesis
# "x = 1 predicts class 1" (direction = "greater") is wrong for these data.
# Training ESS = (1/6 + 1/6 - 1) * 100 = -66.67; "less" gives +66.67.
.xb <- c(1, 1, 1, 1, 1, 0,  0, 0, 0, 0, 0, 1)
.yb <- c(0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 1)

test_that("directional binary LOO keeps the declared direction", {
  for (dir in c("greater", "less")) {
    f <- oda_fit(x = .xb, y = .yb, attr_type = "binary", priors_on = TRUE,
                 mcarlo = FALSE, loo = "off", direction = dir)
    expect_equal(f$ess, if (dir == "greater") -200 / 3 else 200 / 3,
                 tolerance = 1e-9, label = paste(dir, "training ESS"))
    lr <- oda:::oda_loo_for_rule(x = .xb, y = .yb, rule = f$rule,
                                 attr_type = "binary", priors_on = TRUE,
                                 direction = dir)
    expect_true(isTRUE(lr$allowed), label = paste(dir, "LOO allowed"))
    expect_equal(lr$ess_loo, .loo_oracle_ess_binary(.xb, .yb, dir),
                 tolerance = 1e-9, label = paste(dir, "LOO == delete-one refits"))
    expect_equal(lr$ess_loo, f$ess, tolerance = 1e-9,
                 label = paste(dir, "LOO ESS == training ESS"))
  }
})

test_that("directional binary fit with loo = 'stable' is not spuriously unstable", {
  f <- oda_fit(x = .xb, y = .yb, attr_type = "binary", priors_on = TRUE,
               mcarlo = FALSE, loo = "stable", direction = "greater")
  expect_false(identical(f$reason, "loo_not_stable"))
  expect_equal(f$loo$ess_loo, f$ess, tolerance = 1e-9)
})

test_that("weighted directional binary LOO keeps the declared direction", {
  # priors off: weighted 2x2 shortcut path (was affected); priors on with
  # non-uniform weights: count-priors per-fold refits (already correct).
  wt <- c(2, 1, 1, 3, 1, 1,  1, 2, 1, 1, 3, 1)
  for (pri in c(FALSE, TRUE)) {
    for (dir in c("greater", "less")) {
      f <- oda_fit(x = .xb, y = .yb, w = wt, attr_type = "binary", priors_on = pri,
                   mcarlo = FALSE, loo = "off", direction = dir)
      lr <- oda:::oda_loo_for_rule(x = .xb, y = .yb, rule = f$rule, w = wt,
                                   attr_type = "binary", priors_on = pri,
                                   direction = dir)
      lab <- paste0("priors_on=", pri, " ", dir)
      expect_true(isTRUE(lr$allowed), label = paste(lab, "LOO allowed"))
      expect_equal(lr$ess_loo, f$ess, tolerance = 1e-9,
                   label = paste(lab, "weighted LOO ESS == training ESS"))
    }
  }
})

test_that("count-table LOO shortcuts refuse directional rules", {
  fb <- oda_fit(x = .xb, y = .yb, attr_type = "binary", priors_on = TRUE,
                mcarlo = FALSE, loo = "off", direction = "greater")
  w1 <- rep(1, length(.yb))
  expect_null(oda:::oda_loo_binary_map_counts(.xb, .yb, w1, TRUE, fb$rule,
                                              direction = "greater"))
  expect_false(is.null(oda:::oda_loo_binary_map_counts(.xb, .yb, w1, TRUE, fb$rule)))

  xo <- c(1, 2, 3, 4, 5, 6, 7, 8); yo <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  fo <- oda_fit(x = xo, y = yo, attr_type = "ordered", priors_on = TRUE,
                mcarlo = FALSE, loo = "off", direction = "greater")
  expect_null(oda:::oda_loo_ordered_cut_counts(xo, yo, rep(1, 8), TRUE, fo$rule,
                                               direction = "greater"))
  expect_false(is.null(oda:::oda_loo_ordered_cut_counts(xo, yo, rep(1, 8), TRUE, fo$rule)))
})

test_that("oda_loo_for_rule stops when the rule contradicts the declared direction", {
  f <- oda_fit(x = .xb, y = .yb, attr_type = "binary", priors_on = TRUE,
               mcarlo = FALSE, loo = "off", direction = "less")      # rule "1->0"
  expect_error(
    oda:::oda_loo_for_rule(x = .xb, y = .yb, rule = f$rule, attr_type = "binary",
                           priors_on = TRUE, direction = "greater"),
    "contradicts direction")
})

test_that("binary LOO matches delete-one refits (random data, all directions)", {
  for (s in 1:10) {
    set.seed(s)
    x <- rbinom(16L, 1L, 0.5); y <- rbinom(16L, 1L, 0.5)
    if (min(table(factor(x, 0:1))) < 2L || min(table(factor(y, 0:1))) < 2L) next
    for (dir in c("both", "greater", "less")) {
      f <- oda_fit(x = x, y = y, attr_type = "binary", priors_on = TRUE,
                   mcarlo = FALSE, loo = "off", direction = dir)
      if (!isTRUE(f$ok)) next
      lr <- oda:::oda_loo_for_rule(x = x, y = y, rule = f$rule, attr_type = "binary",
                                   priors_on = TRUE, direction = dir)
      expect_equal(lr$ess_loo, .loo_oracle_ess_binary(x, y, dir), tolerance = 1e-9,
                   label = paste0("seed ", s, " direction ", dir))
    }
  }
})
