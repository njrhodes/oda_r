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

### Classification Tree Analysis and MDSA reporting

`oda_cta_fit()` supports binary-class CTA with ENUMERATE, LOO STABLE, PRUNE,
and MINDENOM endpoint constraints. Predictions use `predict(tree, newdata)`.

```r
library(odacore)

# Fit myeloma CTA family (MINDENOM=1 → termination)
fam <- cta_descendant_family(
  X = df[, attrs], y = y, w = w,
  start_mindenom = 1L, max_steps = 20L,
  priors_on = TRUE, miss_codes = -9,
  alpha_split = 0.05, prune_alpha = 0.05,
  mc_iter = 5000L, loo = "stable", attr_names = attrs
)

# Family summary with D-statistic and minimum-D selection
cta_family_table(fam)
#   index mindenom     status no_tree strata min_terminal_denom next_mindenom
# 1     1        1 valid_tree   FALSE      3                 29            30
# 2     2       30      stump   FALSE      2                 55            56
# 3     3       56    no_tree    TRUE     NA                 NA            NA
#   overall_ess has_weights        d selected_min_d
# 1    27.68636        TRUE  7.83566           TRUE
# 2    16.51269        TRUE 10.11190          FALSE
# 3          NA        TRUE       NA          FALSE

# Selected member (minimum D) — MINDENOM=1
t1 <- fam$members[[fam$min_d_idx]]$tree

# Staging table (ordered by target-class propensity)
cta_staging_table(t1, target_class = 1L)[,
  c("stage", "path", "terminal_prediction", "target_proportion")]

# Endpoint × class propensity weights
cta_propensity_weights(t1, target_class = 1L)[,
  c("endpoint_id", "class", "class_n", "propensity_weight")]

# Assign individual observations to endpoints
cta_assign_endpoints(t1, newdata = df[1:5, attrs])

# Observation-level propensity weights
cta_observation_weights(t1, newdata = df[1:5, attrs], y = y[1:5])[,
  c("row_id", "actual_class", "endpoint_id", "propensity_weight", "assigned")]

# Novometric bootstrap CI on the selected tree's confusion
conf <- matrix(c(146, 40,
                  36, 33), nrow = 2, byrow = TRUE)  # actual x predicted
novo_boot_ci(conf, nboot = 5000L, seed = 42L)
```

See `docs/CTA_TRANSLATION_STACK.md` for the full reporting pipeline and
`docs/myeloma-cta-translation.md` for a complete walkthrough with computed output.

## Architecture

```
R/
├── utils.R          — %||%, tick(), fmt helpers (loaded first)
├── unioda_core.R    — Binary-class engine: oda_univariate_core()
├── multioda_core.R  — Multiclass engine: oda_multiclass_unioda_core()
├── oda_fit.R        — Unified dispatcher: oda_fit() routes on C
├── oda_s3.R         — ODA S3 methods: predict/print/summary
├── cta_core.R       — Classification tree: oda_cta_fit()
├── cta_s3.R         — CTA S3 + translation layer: staging/propensity/endpoints
└── cta_family.R     — MDSA family: cta_descendant_family(), cta_family_table()
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
