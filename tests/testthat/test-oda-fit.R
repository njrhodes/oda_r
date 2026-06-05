###############################################################################
# test-oda-fit.R
# Tests for the unified oda_fit() entry point.
###############################################################################

test_that("oda_fit: routes to binary engine for C=2", {
  x <- c(1,2,3,4,5,6,7,8)
  y <- c(0L,0L,0L,0L,1L,1L,1L,1L)

  fit <- oda_fit(x=x, y=y, mcarlo=FALSE, loo="off")
  expect_true(fit$ok)
  expect_equal(fit$engine, "binary")
  expect_true(!is.null(fit$rule$cut_value) || !is.null(fit$rule$left_levels))
})

test_that("oda_fit: routes to multiclass engine for C=3", {
  x <- c(1,2,3,4,5,6,7,8,9)
  y <- c(1L,1L,1L,2L,2L,2L,3L,3L,3L)

  fit <- oda_fit(x=x, y=y, mcarlo=FALSE, loo="off")
  expect_true(fit$ok)
  expect_equal(fit$engine, "multiclass")
  expect_true(!is.null(fit$rule$cut_values))
})

test_that("oda_fit: missing_code alias works", {
  x <- c(1,2,3,4,5,-99)
  y <- c(0L,0L,0L,1L,1L,1L)

  fit_no_miss  <- oda_fit(x=x, y=y, mcarlo=FALSE)
  fit_with_miss <- oda_fit(x=x, y=y, missing_code=-99, mcarlo=FALSE)

  expect_true(fit_no_miss$ok)
  expect_true(fit_with_miss$ok)
  # With -99 excluded, 5 obs; cut should be in [1,5] not near -99
  expect_true(fit_with_miss$n_eff == 5L)
  expect_true(fit_with_miss$rule$cut_value > -10)
})

test_that("oda_fit: miss_codes vector also works", {
  x <- c(1,2,3,4,5,-99,-88)
  y <- c(0L,0L,0L,1L,1L,1L,0L)

  fit <- oda_fit(x=x, y=y, miss_codes=c(-99,-88), mcarlo=FALSE)
  expect_true(fit$ok)
  expect_equal(fit$n_eff, 5L)
})

test_that("oda_fit: pure node returns ok=FALSE", {
  x <- c(1,2,3,4)
  y <- c(0L,0L,0L,0L)

  fit <- oda_fit(x=x, y=y, mcarlo=FALSE)
  expect_false(fit$ok)
})

test_that("oda_fit: all_missing returns ok=FALSE", {
  x <- rep(NA_real_, 5)
  y <- c(0L,1L,0L,1L,0L)

  fit <- oda_fit(x=x, y=y, mcarlo=FALSE)
  expect_false(fit$ok)
})

test_that("oda_fit: binary and multiclass ESS are both finite positive for separated data", {
  # Binary
  x2 <- c(1,2,3,7,8,9); y2 <- c(0L,0L,0L,1L,1L,1L)
  f2 <- oda_fit(x=x2, y=y2, mcarlo=FALSE)
  expect_true(f2$ok)
  expect_true(is.finite(f2$ess) && f2$ess > 0)

  # Multiclass
  x3 <- c(1,2,3,4,5,6,7,8,9); y3 <- c(1L,1L,1L,2L,2L,2L,3L,3L,3L)
  f3 <- oda_fit(x=x3, y=y3, mcarlo=FALSE)
  expect_true(f3$ok)
  # multiclass exposes $ess (ess_pac retained as compat alias)
  expect_true(is.finite(f3$ess) && f3$ess > 0)
})

test_that("oda_fit multiclass: $ess and $ess_pac are both present and equal (Patch 1 transition)", {
  # Patch 1: $ess is the normalized public name; $ess_pac is the compat alias.
  # Both must exist on multiclass results and hold the same computed value
  # until $ess_pac is retired in a future patch.
  x <- c(1,2,3,4,5,6,7,8,9)
  y <- c(1L,1L,1L,2L,2L,2L,3L,3L,3L)
  fit <- oda_fit(x, y, mcarlo = FALSE)
  expect_true(fit$ok)
  expect_true(is.finite(fit$ess)     && fit$ess     > 0)
  expect_true(is.finite(fit$ess_pac) && fit$ess_pac > 0)
  expect_equal(fit$ess, fit$ess_pac)
})

# ---- Phase 6C: direction="ascending"/"descending" and direction_map -----------

test_that("direction='ascending' on binary gives error (use 'greater'/'less')", {
  x <- c(1,2,3,7,8,9); y <- c(0L,0L,0L,1L,1L,1L)
  expect_error(
    oda_fit(x=x, y=y, direction="ascending", mcarlo=FALSE),
    regexp = "ascending.*multiclass|greater.*less",
    ignore.case = TRUE
  )
})

test_that("direction='descending' on binary gives error", {
  x <- c(1,2,3,7,8,9); y <- c(0L,0L,0L,1L,1L,1L)
  expect_error(
    oda_fit(x=x, y=y, direction="descending", mcarlo=FALSE),
    regexp = "descending.*multiclass|greater.*less",
    ignore.case = TRUE
  )
})

test_that("direction='greater'/'less' on multiclass gives warning and is ignored", {
  x <- c(1,2,3,4,5,6,7,8,9); y <- c(1L,1L,1L,2L,2L,2L,3L,3L,3L)
  # greater on multiclass warns and falls back to nondirectional
  expect_warning(oda_fit(x=x, y=y, direction="greater", mcarlo=FALSE),
                 regexp = "greater.*binary|Chapter 2", ignore.case = TRUE)
  expect_warning(oda_fit(x=x, y=y, direction="less",    mcarlo=FALSE),
                 regexp = "less.*binary|Chapter 2",    ignore.case = TRUE)
})

test_that("multiclass ordered ascending matches nondirectional when identity rule is optimal", {
  # Separated 3-class ordered data: ascending rule IS the optimal rule.
  x <- c(1,2,3,4,5,6,7,8,9); y <- c(1L,1L,1L,2L,2L,2L,3L,3L,3L)
  fit_nd  <- oda_fit(x=x, y=y, mcarlo=FALSE, direction="both")
  fit_asc <- oda_fit(x=x, y=y, mcarlo=FALSE, direction="ascending")
  expect_true(fit_nd$ok && fit_asc$ok)
  expect_equal(fit_asc$ess, fit_nd$ess, tolerance = 1e-6)
  expect_equal(fit_asc$rule$seg_classes, 1:3)
})

test_that("multiclass ordered descending finds reversed assignment", {
  # Reversed 3-class data: class 3 at low end, class 1 at high end.
  # descending should find seg_classes = c(3,2,1).
  x <- c(1,2,3,4,5,6,7,8,9); y <- c(3L,3L,3L,2L,2L,2L,1L,1L,1L)
  fit_desc <- oda_fit(x=x, y=y, mcarlo=FALSE, direction="descending")
  fit_nd   <- oda_fit(x=x, y=y, mcarlo=FALSE)
  expect_true(fit_desc$ok)
  expect_equal(fit_desc$rule$seg_classes, c(3L,2L,1L))
  # ESS should match nondirectional when descending is the global optimum
  expect_equal(fit_desc$ess, fit_nd$ess, tolerance = 1e-6)
})

test_that("multiclass categorical ascending with L==C auto-creates identity map", {
  # 3 categories, 3 classes: ascending -> identity direction_map
  x <- c(1L,1L,1L, 2L,2L,2L, 3L,3L,3L)
  y <- c(1L,1L,2L, 2L,2L,3L, 3L,3L,1L)
  fit <- oda_fit(x=x, y=y, attr_type="categorical", mcarlo=FALSE, direction="ascending")
  expect_true(fit$ok)
  expect_equal(fit$rule$type, "multiclass_nominal")
  # Identity map: level 1->class1, level 2->class2, level 3->class3
  expect_equal(as.integer(fit$rule$level_class), c(1L, 2L, 3L))
})

test_that("multiclass categorical ascending with L!=C and no direction_map gives error", {
  # 4 attribute levels, 3 classes: L != C without direction_map -> error
  x <- c(1L,1L, 2L,2L, 3L,3L, 4L,4L)
  y <- c(1L,1L, 2L,2L, 3L,3L, 1L,2L)
  expect_error(
    oda_fit(x=x, y=y, attr_type="categorical", direction="ascending", mcarlo=FALSE),
    regexp = "direction_map|L.*C|direction.*ascending",
    ignore.case = TRUE
  )
})

test_that("binary categorical direction_map: fixed-partition evaluated correctly", {
  # 4-level categorical, binary class: use direction_map to specify partition
  x <- c(1L,1L,2L,2L,3L,3L,4L,4L)
  y <- c(1L,1L,1L,1L,0L,0L,0L,0L)
  # direction_map: levels 1,2 -> class 1; levels 3,4 -> class 0
  fit <- oda_fit(x=x, y=y, attr_type="categorical", mcarlo=FALSE,
                 direction_map = c("1"=1L, "2"=1L, "3"=0L, "4"=0L))
  expect_true(fit$ok)
  expect_equal(fit$rule$type, "nominal_cut")
  # ESS = 100% (perfect separation)
  expect_equal(fit$ess, 100, tolerance = 1e-6)
})

test_that("binary direction_map with mismatched levels gives ok=FALSE", {
  x <- c(1L,2L,3L,4L); y <- c(0L,0L,1L,1L)
  fit <- oda_fit(x=x, y=y, attr_type="categorical", mcarlo=FALSE,
                 direction_map = c("1"=0L, "2"=0L, "99"=1L))  # "99" not a level
  expect_false(fit$ok)
  expect_match(fit$reason %||% "direction_map_levels_mismatch",
               "direction_map_levels_mismatch|ok")
})

test_that("multiclass categorical direction_map: fixed-partition overrides search", {
  # 4 levels, 4 classes: supply a non-identity direction_map, verify it's used
  x <- c(1L,1L, 2L,2L, 3L,3L, 4L,4L)
  y <- c(1L,2L, 2L,3L, 3L,4L, 4L,1L)
  # Reverse map: level 1->class4, level 2->class3, level 3->class2, level 4->class1
  dmap <- c("1"=4L, "2"=3L, "3"=2L, "4"=1L)
  fit <- oda_fit(x=x, y=y, attr_type="categorical", mcarlo=FALSE,
                 direction_map=dmap)
  expect_true(fit$ok)
  expect_equal(as.integer(fit$rule$level_class), c(4L, 3L, 2L, 1L))
})

# ---- LOO semantics: numeric loo must not become "off" -----------------------
#
# Bug fixed: oda_fit(..., loo = numeric) was silently converted to "off" by
# switch(as.character(0.05), ...) hitting the default branch.  After the fix,
# numeric loo maps to loo = "pvalue" with loo_alpha = the numeric value.
#
# Threshold semantics:
#   loo = "pvalue"  -> pvalue gate; default threshold 0.05; pass requires p < 0.05
#   loo = numeric   -> pvalue gate; threshold = supplied value; pass requires p < loo
#   loo = "stable"  -> STABLE gate; |WESSL - WESS| <= 0.01 pp
#   loo = "off"     -> no LOO; fit$loo will be NULL
#
# Fixture: n=20, x=1:20, y=10x0 + 10x1. Probe confirms:
#   loo=0.99 -> ok=TRUE, p_loo=1e-4 (passes)
#   loo="pvalue" -> ok=TRUE, p_loo=1e-4 (passes default 0.05)
#   loo=1e-9  -> ok=FALSE, reason=loo_p_not_significant (gate rejects deterministically)
#   loo="stable" -> ok=FALSE (LOO genuinely unstable on integer-separable data - not a bug)
#   loo="off" -> fit$loo = NULL (LOO never runs)

.loo_sem_xy <- function() {
  list(x = seq_len(20L), y = c(rep(0L, 10L), rep(1L, 10L)))
}

test_that("oda_fit: loo='off' leaves fit$loo NULL (no LOO computation)", {
  d   <- .loo_sem_xy()
  fit <- oda_fit(d$x, d$y, loo = "off", mcarlo = TRUE, mc_iter = 1000L, mc_seed = 1L)
  expect_true(isTRUE(fit$ok), label = "loo='off' fit must succeed")
  expect_null(fit$loo,        label = "loo='off' must leave fit$loo NULL")
})

test_that("oda_fit: numeric loo=0.99 runs LOO and produces non-NULL fit$loo", {
  # loo=0.99 passes the gate (p_loo=1e-4 < 0.99) - proves LOO ran, not silently 'off'.
  d   <- .loo_sem_xy()
  fit <- oda_fit(d$x, d$y, loo = 0.99, mcarlo = TRUE, mc_iter = 1000L, mc_seed = 1L)
  expect_true(isTRUE(fit$ok),
              label = "loo=0.99 must succeed on n=20 fixture (p_loo=1e-4 << 0.99)")
  expect_false(is.null(fit$loo),
               label = "fit$loo must be non-NULL when numeric loo ran (not 'off')")
})

test_that("oda_fit: loo=1e-9 rejects via LOO p-value gate with correct reason", {
  # p_loo=1e-4 for this fixture; 1e-4 >= 1e-9 is TRUE -> gate rejects.
  # ok=FALSE with reason="loo_p_not_significant" proves the gate ran and rejected.
  # Rule and confusion are present even when gated (ok=FALSE is not "no model").
  d   <- .loo_sem_xy()
  fit <- oda_fit(d$x, d$y, loo = 1e-9, mcarlo = TRUE, mc_iter = 1000L, mc_seed = 1L)
  expect_false(isTRUE(fit$ok),
               label = "loo=1e-9 must be rejected by LOO gate (p_loo=1e-4 >= 1e-9)")
  expect_equal(fit$reason, "loo_p_not_significant",
               label = "rejection reason must be 'loo_p_not_significant'")
  expect_false(is.null(fit$rule),
               label = "rule must be present even when LOO gated")
  expect_false(is.null(fit$confusion),
               label = "confusion must be present even when LOO gated")
  expect_false(is.null(fit$loo),
               label = "loo result must be present when LOO gated")
})

test_that("oda_fit: loo=0.99 and loo='pvalue' produce same LOO ESS (same computation, different threshold)", {
  # Both pass on n=20 fixture (p_loo=1e-4 < min(0.99, 0.05)).
  # They must produce identical ess_loo because the underlying LOO computation
  # is the same regardless of threshold.
  d       <- .loo_sem_xy()
  fit_num <- oda_fit(d$x, d$y, loo = 0.99,     mcarlo = TRUE, mc_iter = 1000L, mc_seed = 1L)
  fit_pv  <- oda_fit(d$x, d$y, loo = "pvalue", mcarlo = TRUE, mc_iter = 1000L, mc_seed = 1L)
  expect_true(isTRUE(fit_num$ok), label = "loo=0.99 must succeed on n=20 fixture")
  expect_true(isTRUE(fit_pv$ok),  label = "loo='pvalue' must succeed on n=20 fixture")
  expect_equal(fit_num$loo$ess_loo, fit_pv$loo$ess_loo, tolerance = 1e-9,
               label = "loo=0.99 and loo='pvalue' must produce same LOO ESS")
})

# ---- Numeric loo validation: invalid inputs must error ----------------------
#
# oda_fit() validates numeric loo: must be length 1, non-NA, finite, in (0, 1).
# Rejection message contains "strictly in (0, 1)".

test_that("oda_fit: loo=NA_real_ errors", {
  x <- 1:8; y <- c(0L,0L,0L,0L,1L,1L,1L,1L)
  expect_error(oda_fit(x, y, loo = NA_real_), regexp = "strictly in .0, 1.")
})

test_that("oda_fit: loo=c(0.01, 0.05) (length > 1) errors", {
  x <- 1:8; y <- c(0L,0L,0L,0L,1L,1L,1L,1L)
  expect_error(oda_fit(x, y, loo = c(0.01, 0.05)), regexp = "strictly in .0, 1.")
})

test_that("oda_fit: loo=-0.1 (negative) errors", {
  x <- 1:8; y <- c(0L,0L,0L,0L,1L,1L,1L,1L)
  expect_error(oda_fit(x, y, loo = -0.1), regexp = "strictly in .0, 1.")
})

test_that("oda_fit: loo=1 (boundary, not strictly < 1) errors", {
  x <- 1:8; y <- c(0L,0L,0L,0L,1L,1L,1L,1L)
  expect_error(oda_fit(x, y, loo = 1), regexp = "strictly in .0, 1.")
})

test_that("oda_fit: loo=Inf errors", {
  x <- 1:8; y <- c(0L,0L,0L,0L,1L,1L,1L,1L)
  expect_error(oda_fit(x, y, loo = Inf), regexp = "strictly in .0, 1.")
})

# ---- summary(multiclass LOO) — per-class PAC visible in print output ----
#
# Requirement: summary.oda_fit must carry 'classes' so that
# print.oda_fit_summary can render the LOO CLASS/PAC table for multiclass
# and multicategorical fits.  Tests lock this behaviour for both ordered
# and categorical multiclass paths.

test_that("summary: iris ordered multiclass LOO includes classes and per-class PAC", {
  data(iris)
  fit <- oda_fit(x = iris$Petal.Length, y = as.integer(iris$Species),
                 attr_type = "ordered", loo = "on", mcarlo = FALSE)
  s   <- summary(fit)

  # classes carried into summary
  expect_equal(s$classes, 1:3, label = "iris summary carries classes 1:3")

  # LOO section present and allowed
  expect_true(isTRUE(s$loo$allowed), label = "iris LOO allowed")

  # per-class PAC present (proportion scale: 1, 0.88, 0.94).
  # pac_by_class is a named vector (class labels as names); ignore_attr = TRUE.
  pb <- s$loo$pac_by_class
  expect_equal(length(pb), 3L, label = "iris LOO pac_by_class has 3 elements")
  expect_equal(pb[[1]], 1.00, tolerance = 1e-6, ignore_attr = TRUE, label = "iris LOO PAC class 1 = 100%")
  expect_equal(pb[[2]], 0.88, tolerance = 1e-6, ignore_attr = TRUE, label = "iris LOO PAC class 2 = 88%")
  expect_equal(pb[[3]], 0.94, tolerance = 1e-6, ignore_attr = TRUE, label = "iris LOO PAC class 3 = 94%")

  # aggregate LOO metrics (percent scale)
  expect_equal(s$loo$mean_pac, 94.00, tolerance = 0.01, label = "iris LOO Mean PAC = 94%")
  expect_equal(s$loo$ess_loo,  91.00, tolerance = 0.01, label = "iris LOO ESS = 91%")

  # no LOO p-value for multicategorical
  expect_equal(s$loo$p_status, "not_computed",
               label = "iris LOO p_status not_computed (no canon C>2 Fisher p)")

  # print output contains per-class table and p(LOO) note
  out <- capture.output(print(s))
  expect_true(any(grepl("CLASS.*PAC", out)),        label = "iris summary print: CLASS/PAC header")
  expect_true(any(grepl("1.*100\\.0%", out)),        label = "iris summary print: class 1 = 100%")
  expect_true(any(grepl("2.*88\\.0%",  out)),        label = "iris summary print: class 2 = 88%")
  expect_true(any(grepl("3.*94\\.0%",  out)),        label = "iris summary print: class 3 = 94%")
  expect_true(any(grepl("LOO ESS.*91", out)),        label = "iris summary print: LOO ESS = 91%")
  expect_true(any(grepl("p\\(LOO\\).*not reported", out)),
              label = "iris summary print: p(LOO) not reported note")
})

test_that("summary: protein categorical multiclass LOO includes classes, rule, per-class PAC", {
  biological_type <- c(
    rep(1L, 98), rep(2L, 13), rep(3L,  6), rep(4L,  7),
    rep(1L, 16), rep(2L, 50), rep(3L,  4), rep(4L, 19),
    rep(1L,  5), rep(2L,  2), rep(3L, 23), rep(4L, 14),
    rep(1L,  3), rep(2L,  8), rep(3L, 12), rep(4L, 45)
  )
  amino_acid_type <- c(rep(1L, 124), rep(2L, 89), rep(3L, 44), rep(4L, 68))

  fit <- oda_fit(x = amino_acid_type, y = biological_type,
                 attr_type = "categorical", direction = "ascending",
                 mc_iter = 500L, mc_seed = 42L, loo = "on")
  s   <- summary(fit)

  # classes carried into summary
  expect_equal(s$classes, 1:4, label = "protein summary carries classes 1:4")

  # rule renders as identity mapping (not placeholder)
  expect_true(grepl("1 --> 1", s$rule_string), label = "protein rule_string contains identity map")
  expect_false(grepl("<nominal", s$rule_string), label = "protein rule_string not placeholder")

  # LOO section present and allowed
  expect_true(isTRUE(s$loo$allowed), label = "protein LOO allowed")

  # per-class PAC present (proportion scale).
  # pac_by_class is a named vector (class labels as names); ignore_attr = TRUE.
  pb <- s$loo$pac_by_class
  expect_equal(length(pb), 4L, label = "protein LOO pac_by_class has 4 elements")
  expect_equal(pb[[1]], 98/122, tolerance = 1e-4, ignore_attr = TRUE, label = "protein LOO PAC class 1")
  expect_equal(pb[[2]], 50/73,  tolerance = 1e-4, ignore_attr = TRUE, label = "protein LOO PAC class 2")
  expect_equal(pb[[3]], 23/45,  tolerance = 1e-4, ignore_attr = TRUE, label = "protein LOO PAC class 3")
  expect_equal(pb[[4]], 45/85,  tolerance = 1e-4, ignore_attr = TRUE, label = "protein LOO PAC class 4")

  # aggregate LOO metrics (percent scale)
  expect_equal(s$loo$mean_pac, 63.22, tolerance = 0.01, label = "protein LOO Mean PAC")
  expect_equal(s$loo$ess_loo,  50.96, tolerance = 0.01, label = "protein LOO ESS")

  # no LOO p-value
  expect_equal(s$loo$p_status, "not_computed",
               label = "protein LOO p_status not_computed")

  # print output contains per-class table and p(LOO) note
  out <- capture.output(print(s))
  expect_true(any(grepl("CLASS.*PAC", out)),              label = "protein summary print: CLASS/PAC header")
  expect_true(any(grepl("1.*80\\.3%", out)),              label = "protein summary print: class 1 = 80.3%")
  expect_true(any(grepl("2.*68\\.5%", out)),              label = "protein summary print: class 2 = 68.5%")
  expect_true(any(grepl("3.*51\\.1%", out)),              label = "protein summary print: class 3 = 51.1%")
  expect_true(any(grepl("4.*52\\.9%", out)),              label = "protein summary print: class 4 = 52.9%")
  expect_true(any(grepl("p\\(LOO\\).*not reported", out)),
              label = "protein summary print: p(LOO) not reported note")
})
