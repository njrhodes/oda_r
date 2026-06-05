# data-raw/prepare_datasets.R
# Run from the package root to regenerate data/myeloma.rda and data/cta_demo.rda.
#
# Usage:
#   Rscript --vanilla data-raw/prepare_datasets.R
#
# These datasets mirror the fixture files in tests/testthat/fixtures/ exactly.
# Do not modify column names or row order without updating fixture expectations.

# ---- myeloma ----------------------------------------------------------------
myeloma <- read.table(
  "tests/testthat/fixtures/myeloma/data.txt",
  header = FALSE
)
colnames(myeloma) <- paste0("V", seq_len(ncol(myeloma)))
save(myeloma, file = "data/myeloma.rda", compress = "bzip2")
message("Saved data/myeloma.rda  [", nrow(myeloma), " x ", ncol(myeloma), "]")

# ---- cta_demo ---------------------------------------------------------------
cta_demo <- read.csv(
  "tests/testthat/fixtures/cta_demo/CTA_DEMO.CSV",
  header = FALSE
)
colnames(cta_demo) <- paste0("V", seq_len(ncol(cta_demo)))
save(cta_demo, file = "data/cta_demo.rda", compress = "bzip2")
message("Saved data/cta_demo.rda  [", nrow(cta_demo), " x ", ncol(cta_demo), "]")
