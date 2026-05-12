# odacore

Pure-R reimplementation of the MegaODA / CTA classification engine.

- **UniODA** (`oda_univariate_core`) — binary-class ODA for ordered, categorical, and binary attributes
- **MultiODA** (`oda_multiclass_unioda_core`) — multiclass ODA for ordered and categorical attributes
- **`oda_fit`** — unified dispatcher: calls UniODA when C = 2, MultiODA when C ≥ 3
- **CTA** (`oda_cta_fit`) — Classification Tree Analysis with ENUMERATE, LOO STABLE, pruning, and fixture-driven CTA.exe parity tests

Core UniODA/MultiODA behavior and covered CTA fixtures are validated against MegaODA.exe / CTA.exe golden outputs.

## Installation (development)

```r
# install.packages("devtools")
devtools::install_github("njrhodes/oda_rcore")
```

Or locally:

```r
devtools::install("path/to/oda_rcore")
```

## Quick start

### Binary class (C = 2)

```r
library(odacore)

x <- c(1, 2, 3, 4, 5, 6, 7, 8)
y <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)

fit <- oda_fit(x, y, mcarlo = FALSE)
fit$rule$cut_value   # 4.5
fit$rule$direction   # "0->1"
fit$ess              # 100
fit$confusion        # list: TP, TN, FP, FN, sensitivity, specificity, mean_pac
```

### Multiclass (C ≥ 3)

```r
library(odacore)

y <- as.integer(factor(iris$Species,
                        levels = c("setosa", "versicolor", "virginica")))
fit <- oda_fit(
  x         = iris$Petal.Length,
  y         = y,
  attr_type = "ordered",
  priors_on = TRUE,
  mcarlo    = FALSE,
  loo       = "off"
)

fit$rule$cut_values   # c(2.45, 4.75)
fit$rule$seg_classes  # c(1, 2, 3)
fit$ess_pac           # ~93  (PAC-based ESS)
fit$mean_pac          # ~95.3
```

### Classification Tree Analysis

`oda_cta_fit()` supports binary-class CTA with ENUMERATE, LOO STABLE, PRUNE,
and MINDENOM endpoint constraints. Predictions use `predict(tree, newdata)`.

```r
library(odacore)

# Binary synthetic example
set.seed(1)
n  <- 60
X  <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
y  <- as.integer(X$x1 + X$x2 > 0)

tree <- oda_cta_fit(
  X           = X,
  y           = y,
  priors_on   = TRUE,
  alpha_split = 0.05,
  prune_alpha = 0.05,
  loo         = "stable",
  mc_iter     = 5000L
)

preds    <- predict(tree, X)                           # majority fallback for missing path attributes
preds_na <- predict(tree, X, missing_action = "na")    # path-local NA for missing path attributes
```

## Architecture

```
R/
├── utils.R          — %||%, tick(), fmt helpers (loaded first)
├── unioda_core.R    — Binary-class engine: oda_univariate_core()
├── multioda_core.R  — Multiclass engine: oda_multiclass_unioda_core()
├── oda_fit.R        — Unified dispatcher: oda_fit() routes on C
└── cta_core.R       — Classification tree: oda_cta_fit()
```

## Tie-breaking spec (MegaODA-faithful)

```
PRIMARY   = MAXSENS   (overall PAC in priors-weighted objective space)
SECONDARY = SAMPLEREP (L1 distance: predicted vs observed class frequencies)
TERTIARY  = FIRST IDENTIFIED (enumeration order — tick() integer quantisation)
```

LOO validity is a post-selection check, not a tie-break criterion.

## Running tests

```r
devtools::test()                          # all tests
devtools::test(filter = "iris")           # iris multiclass gold
devtools::test(filter = "fixture-cta")    # CTA_DEMO gold fixture
devtools::test(filter = "tie-breaking")   # SAMPLEREP isolation
devtools::check()                         # full R CMD check
```

## CI

GitHub Actions runs `R CMD check` on:

| OS | R |
|---|---|
| ubuntu-latest | release, devel |
| windows-latest | release |
| macos-latest | release |

Triggered on push/PR to `main`.
