###############################################################################
# test-unioda-weighted-megaoda.R - weighted ordered UniODA vs MegaODA.exe
#
# Fixture: fixtures/myeloma/data.txt (same file as the CTA fixtures).
# Program: fixtures/myeloma/megaoda_v4.pgm
#   OPEN data.txt; OUTPUT UNIV4.OUT; VARS v1 to v19; CLASS v1; ATTRIBUTE v4;
#   EX v2=0; MISSING ALL (-9); WEIGHT v2; GO;
# Expected output: fixtures/myeloma/UNIV4.OUT (run MegaODA.exe on the program
# above in that folder). The test skips until UNIV4.OUT is present.
#
# Why: CTA.exe chooses weighted cuts with count priors (sum of correct case
# weights / raw class count). oda_apply_priors() divides by class weight totals
# instead, which maximises WESS. The two rules give different V4 cuts on this
# fixture (count priors: 890.25, WESS 6.93; max WESS: 64.0, WESS 16.16).
# This test records whichever MegaODA actually does.
###############################################################################

.univ4_path <- function(f) {
  tryCatch(testthat::test_path(file.path("fixtures", "myeloma", f)),
           error = function(e) "")
}

.univ4_exe_cut <- function() {
  f <- .univ4_path("UNIV4.OUT")
  if (!nzchar(f) || !file.exists(f)) return(NULL)
  l <- grep("^\\s*IF\\s+V4\\s*<=", readLines(f, warn = FALSE), value = TRUE)[1L]
  if (is.na(l)) return(NULL)
  as.numeric(sub("^\\s*IF\\s+V4\\s*<=\\s*([-0-9.eE+]+).*$", "\\1", l))
}

test_that("weighted UniODA: myeloma V4 cut matches MegaODA.exe (UNIV4.OUT)", {
  exe_cut <- .univ4_exe_cut()
  if (is.null(exe_cut)) skip("fixtures/myeloma/UNIV4.OUT missing: run MegaODA on megaoda_v4.pgm")
  d <- read.table(.univ4_path("data.txt"), header = FALSE)
  names(d) <- paste0("V", seq_len(ncol(d)))
  d <- d[d$V2 != 0, ]
  fit <- oda_fit(x = d$V4, y = d$V1, w = d$V2, priors_on = TRUE,
                 miss_codes = -9, mcarlo = FALSE, loo = "off")
  expect_true(isTRUE(fit$ok))
  expect_equal(fit$rule$cut_value, exe_cut, tolerance = 0.01,
               label = "oda_fit weighted V4 cut vs MegaODA.exe")
})
