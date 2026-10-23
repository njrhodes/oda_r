###############################################################################
# test-unioda-weighted-megaoda.R - weighted ordered UniODA vs MegaODA.exe
#
# Fixture: fixtures/myeloma/data.txt (same file as the CTA fixtures).
# Programs (run MegaODA.exe on each in that folder):
#   megaoda_v4.pgm    -> UNIV4.OUT    cut only
#     OPEN data.txt; OUTPUT UNIV4.OUT; VARS v1 to v19; CLASS v1; ATTRIBUTE v4;
#     EX v2=0; MISSING ALL (-9); WEIGHT v2; GO;
#   megaoda_v4_mc.pgm -> UNIV4MC.OUT  cut, WESS, MC p, LOO
#     (as above plus MCARLO ITER 25000 TARGET 0.05 STOP 99.9; LOO;)
# Each test skips until its .OUT file is present.
#
# MegaODA (like CTA.exe) chooses weighted cuts with count priors: maximise
# sum(w x correct) / n_class over RAW class counts, not WESS. On this fixture
# count priors give the V4 cut 890.25 (WESS 6.93); maximising WESS would give
# 64.0 (WESS 16.16). The MC p comes from permuting (class, weight) pairs jointly
# with the count-priors objective as the statistic; LOO refits every fold.
###############################################################################

.univ4_path <- function(f) {
  tryCatch(testthat::test_path(file.path("fixtures", "myeloma", f)),
           error = function(e) "")
}

.univ4_out <- function(f) {
  p <- .univ4_path(f)
  if (!nzchar(p) || !file.exists(p)) return(NULL)
  readLines(p, warn = FALSE)
}

.univ4_num <- function(s) as.numeric(sub("%", "", s))

.univ4_data <- function() {
  d <- read.table(.univ4_path("data.txt"), header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d[d$V2 != 0, ]                                   # EX v2=0
}

.univ4_cut <- function(l) {
  r <- grep("^\\s*IF\\s+V4\\s*<=", l, value = TRUE)[1L]
  if (is.na(r)) return(NULL)
  as.numeric(sub("^\\s*IF\\s+V4\\s*<=\\s*([-0-9.eE+]+).*$", "\\1", r))
}

test_that("weighted UniODA: myeloma V4 cut matches MegaODA.exe (UNIV4.OUT)", {
  l <- .univ4_out("UNIV4.OUT")
  if (is.null(l)) skip("fixtures/myeloma/UNIV4.OUT missing: run MegaODA on megaoda_v4.pgm")
  exe_cut <- .univ4_cut(l)
  d <- .univ4_data()
  fit <- oda_fit(x = d$V4, y = d$V1, w = d$V2, priors_on = TRUE,
                 miss_codes = -9, mcarlo = FALSE, loo = "off")
  expect_true(isTRUE(fit$ok))
  expect_equal(fit$rule$cut_value, exe_cut, tolerance = 0.01,
               label = "oda_fit weighted V4 cut vs MegaODA.exe")
})

test_that("weighted UniODA: myeloma V4 WESS, MC p and LOO match MegaODA.exe (UNIV4MC.OUT)", {
  l <- .univ4_out("UNIV4MC.OUT")
  if (is.null(l)) skip("fixtures/myeloma/UNIV4MC.OUT missing: run MegaODA on megaoda_v4_mc.pgm")
  exe_cut  <- .univ4_cut(l)
  exe_pmc  <- .univ4_num(sub(".*:", "", grep("^\\s*Estimated p:", l, value = TRUE)[1L]))
  row      <- grep("^\\s*Effect Strength Wtd PAC\\s", l, value = TRUE)[1L]   # Train  LOO
  exe_wess <- .univ4_num(strsplit(trimws(sub("Effect Strength Wtd PAC", "", row)), "\\s+")[[1L]])
  exe_loop <- .univ4_num(sub(".*p\\s*=\\s*", "", grep("Fisher's exact test", l, value = TRUE)[1L]))

  d <- .univ4_data()
  fit <- oda_univariate_core(x = d$V4, y = d$V1, w = d$V2, attr_type = "ordered",
                             priors_on = TRUE, miss_codes = -9,
                             mcarlo = TRUE, mc_iter = 25000L, mc_target = 0.05,
                             mc_stop = 99.9, mc_seed = 20261023L,
                             loo = "pvalue", loo_alpha = 0.99)
  expect_true(isTRUE(fit$ok))
  expect_equal(fit$rule$cut_value, exe_cut, tolerance = 0.01, label = "cut")
  expect_equal(fit$ess, exe_wess[1L], tolerance = 0.01, label = "training WESS")
  # MC p is stochastic: 25000 iterations each side, SE of the difference ~0.004.
  expect_lt(abs(fit$p_mc - exe_pmc), 0.015, label = "MC p vs MegaODA")
  expect_equal(fit$loo$ess_loo, exe_wess[2L], tolerance = 0.01, label = "LOO WESS")
  expect_equal(fit$loo$p_value, exe_loop, tolerance = 1e-4, label = "LOO Fisher p")
})
