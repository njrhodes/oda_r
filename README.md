# odacore

Pure-R reimplementation of the MegaODA classification engine — binary-class
UniODA and multiclass MultiODA for ordered and categorical attributes.

## Installation (development)

```r
# install.packages("devtools")
devtools::install_github("njrhodes/ODA", subdir = ".", ref = "rcore")
```

Or locally:

```r
devtools::install("path/to/oda_rcore")
```

## Quick start

```r
library(odacore)

# 3-class ordered example (iris Petal.Length)
y <- as.integer(factor(iris$Species,
                        levels = c("setosa","versicolor","virginica")))
fit <- oda_multiclass_unioda_core(
  x          = iris$Petal.Length,
  y          = y,
  attr_type  = "ordered",
  priors_on  = TRUE,
  K_segments = 3L,
  mcarlo     = TRUE,
  loo        = "on"
)

fit$rule$cut_values   # c(2.45, 4.75)
fit$rule$seg_classes  # c(1, 2, 3)
fit$mean_pac          # ~95.3%
```

## Architecture

```
R/
├── utils.R          %||%, tick(), fmt helpers
├── unioda_core.R    Binary-class ODA (oda_univariate_core)
├── multioda_core.R  Multiclass ODA  (oda_multiclass_unioda_core)
└── harness_utils.R  Diagnostic harness (confusion_raw, harness_loo_refit_*)
```

## Tie-breaking spec (MegaODA-faithful)

```
PRIMARY   = MAXSENS   (overall PAC in priors-weighted objective space)
SECONDARY = SAMPLEREP (L1 distance: predicted vs observed class frequencies)
TERTIARY  = FIRST IDENTIFIED (enumeration order — tick() integer quantisation)
LOO       = refit-per-fold (no global rule reuse)
```

## Running tests

```r
devtools::test()           # all tests
devtools::test(filter = "iris")          # iris gold only
devtools::test(filter = "tie-breaking")  # SAMPLEREP isolation
```

## CI

GitHub Actions runs `R CMD check` on ubuntu/windows/macos × R-release/R-devel
on every push to `main` and `rcore` branches.
