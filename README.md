# odacore

Pure-R reimplementation of the MegaODA / CTA classification engine.

- **`oda_fit()`** - ODA dispatcher: UniODA (C = 2) or MultiODA (C >= 3) for ordered, categorical, and binary attributes
- **`cta_fit()`** - Classification Tree Analysis: ENUMERATE, LOO STABLE, PRUNE, MINDENOM endpoint constraints
- **`lort_fit()`** - LORT (Locally Optimal Recursive Tree): greedy recursive CTA/MDSA; adjacent workflow using canon CTA/MDSA components
- **`sda_fit()`** - SDA (Structural Decomposition Analysis): canonical staged attribute identification; novometric `min_d` and legacy `unioda_max_ess` modes
- **`cta_descendant_family()`** - MDSA descendant family: D/parsimony evidence, family-level comparison
- **`novo_boot_ci()`** - Novometric bootstrap CI; S3 generic dispatching on ODA / CTA / LORT fit objects
- **Balance and propensity:** `oda_balance_table()`, `smd_balance_table()`, `cta_balance_table()`, `oda_propensity_weights()`, `lort_propensity_weights()`
- **Graphics v3 (ggplot2):** `plot_cta_tree()`, `plot_lort_tree()`, `plot_oda_balance()`, `plot_smd_balance()`, `plot_cta_balance()`

Core ODA/CTA behavior and covered fixtures are validated against MegaODA.exe / CTA.exe golden outputs.

## Installation (development)

```r
# install.packages("devtools")
devtools::install_github("njrhodes/oda_rcore")
```

Or locally:

```r
devtools::install("path/to/oda_rcore")
```

## Collaborator setup

Preferred dependency installation (validated on macOS arm64 / R 4.3.2):

```r
# Install pak if not already present
install.packages("pak")

# Install all package dependencies
pak::local_install_deps(dependencies = TRUE)
```

Then load and verify:

```r
devtools::load_all()
```

Smoke check (CTA, MDSA, reporting, graphics):

```bash
ODACORE_TEST_TIER=smoke Rscript --vanilla -e "devtools::test(reporter='progress')"
```

Package check:

```bash
Rscript --vanilla -e "devtools::check(vignettes=FALSE)"
```

Expected result: 0 errors / 0 warnings / at most 1 note.
The note about future file timestamps can appear on cloud-synced or
network-mounted folders (e.g., Dropbox). It is a filesystem timestamp
artifact, not a package failure.

`devtools::install_deps(dependencies = TRUE)` is an alternative if pak is
unavailable, but pak is preferred on macOS.

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

### Multiclass (C >= 3)

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

`cta_fit()` is the public entry point for binary-class CTA with ENUMERATE,
LOO STABLE, PRUNE, and MINDENOM endpoint constraints.
Predictions use `predict(tree, newdata)`.
(`oda_cta_fit()` is the internal engine name retained for backward compatibility.)

```r
library(odacore)

# Fit myeloma CTA family (MINDENOM=1 -> termination)
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

# Selected member (minimum D) - MINDENOM=1
t1 <- fam$members[[fam$min_d_idx]]$tree

# Staging table (ordered by target-class propensity)
cta_staging_table(t1, target_class = 1L)[,
  c("stage", "path", "terminal_prediction", "target_proportion")]

# Endpoint x class propensity weights
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

# Tree plot - structural
plot(t1)

# Target-class enriched plot (endpoint fill color encodes relative
# target-class proportion; does not imply clinical thresholds)
plot(t1, target_class = 1L, class_labels = c("0" = "Alive", "1" = "Deceased"))
```

See `docs/CTA_TRANSLATION_STACK.md` for the full reporting pipeline and
`docs/examples/myeloma-cta-translation.md` for a complete walkthrough with computed output.

### Graphics v3 (ggplot2)

Graphics v3 provides direct ggplot2 renderers for trees and balance diagnostics.
`ggplot2 >= 3.4.0` is in `Suggests`; a clear error is raised if it is absent.
Base `plot.*` methods are unchanged and do not require ggplot2.

#### Tree diagrams

```r
library(odacore)
if (requireNamespace("ggplot2", quietly = TRUE)) {

  X <- data.frame(
    A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
    B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
  )
  y <- c(rep(0L, 40), rep(1L, 20))

  # CTA tree
  tree <- cta_fit(X, y, mindenom = 5L, mc_iter = 200L, mc_seed = 42L, loo = "off")
  p <- plot_cta_tree(tree, color_by = "target_rate")
  print(p)

  # LORT tree
  lort <- lort_fit(X, y, mc_iter = 200L, mc_seed = 42L, loo = "off", min_n = 5L)
  p2 <- plot_lort_tree(lort, color_by = "prediction")
  print(p2)

  # All renderers return ggplot objects - save with ggsave
  ggplot2::ggsave("tree.png", p, width = 8, height = 5, dpi = 300)
}
```

#### Covariate balance plots

Plot functions are pure renderers - they do not fit models.

```r
if (requireNamespace("ggplot2", quietly = TRUE)) {

  group <- c(rep(0L, 30), rep(1L, 30))
  X_bl  <- data.frame(
    age   = c(rep(40L, 30), rep(60L, 30)),
    score = rnorm(60, 50, 10)
  )

  # ODA balance: ESS/WESS dot plot
  bt <- oda_balance_table(group, X_bl, mcarlo = TRUE, mc_iter = 500L)
  pd <- oda_balance_plot_data(bt)
  print(plot_oda_balance(pd))

  # SMD / Love plot: |SMD| with reference lines at 0.10 and 0.20
  smd <- smd_balance_table(group, X_bl)
  print(plot_smd_balance(smd, ref_020 = TRUE))
  print(plot_balance_love(smd))   # identical to plot_smd_balance()

  # CTA multivariate balance
  ct  <- cta_balance_table(group, X_bl, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
  cpd <- cta_balance_plot_data(ct)
  p   <- plot_cta_balance(cpd)
  print(p)
  # status = "no_tree" renders a message panel:
  # "No combination of covariates predicted group membership...
  #  This is favorable evidence of multivariable balance..."
  # This is the desired result for well-balanced study arms.
}
```

See `docs/GRAPHICS_V3.md` for the full function reference and interpretation
guide for CTA balance `no_tree` results.

## Architecture

```
R/
+-- utils.R             - %||%, tick(), fmt helpers, .validate_case_weights() (loaded first)
+-- unioda_core.R       - Binary-class engine: oda_univariate_core()
+-- multioda_core.R     - Multiclass engine: oda_multiclass_unioda_core()
+-- oda_fit.R           - Unified dispatcher: oda_fit() routes on C
+-- oda_s3.R            - ODA S3 methods: predict/print/summary
+-- cta_core.R          - CTA engine: cta_fit() / oda_cta_fit() (internal)
+-- cta_s3.R            - CTA S3 + translation layer: staging/propensity/endpoints
+-- cta_family.R        - MDSA family: cta_descendant_family(), cta_family_table()
+-- cta_ort.R           - LORT: lort_fit(), ort_plot_data(), predict.cta_ort()
+-- sda_core.R          - SDA engine: sda_fit(), sda_selected_attributes(), sda_step_table()
+-- sda_anchor.R        - SDA anchor: sda_anchor(), as_sda_anchor(), validate_sda_anchor()
+-- balance.R           - Balance diagnostics: oda/smd/cta balance tables + plot data
+-- production_tools.R  - Readiness checks: oda_readiness_check()
\-- graphics_v3.R       - ggplot2 renderers: plot_cta_tree(), plot_lort_tree(),
                          plot_oda_balance(), plot_smd_balance(),
                          plot_balance_love(), plot_cta_balance()
```

## Tie-breaking spec (MegaODA-faithful)

```
PRIMARY   = MAXSENS   (overall PAC in priors-weighted objective space)
SECONDARY = SAMPLEREP (L1 distance: predicted vs observed class frequencies)
TERTIARY  = FIRST IDENTIFIED (enumeration order - tick() integer quantisation)
```

LOO validity is a post-selection check, not a tie-break criterion.

## Running tests

Tests are gated by `ODACORE_TEST_TIER` (see `tests/testthat/helper-test-tier.R`).
The default (unset) runs CRAN-safe tests only.

```bash
# CRAN-safe / default package check (no env var required):
Rscript --vanilla -e "devtools::check(vignettes=FALSE)"

# Fast developer loop (skips all slow canon fixtures):
ODACORE_TEST_TIER=fast Rscript --vanilla -e "devtools::test(reporter='progress')"

# CTA/MDSA/reporting/graphics production gate (required before those commits):
ODACORE_TEST_TIER=smoke Rscript --vanilla -e "devtools::test(reporter='progress')"

# Full canon / release gate (required before release tags or canon parity work):
ODACORE_TEST_TIER=full Rscript --vanilla -e "devtools::test(reporter='progress')"
```

Targeted filter examples:

```r
devtools::test(filter = "iris")           # iris multiclass gold
devtools::test(filter = "fixture-cta")    # CTA_DEMO gold fixture (smoke tier)
devtools::test(filter = "tie-breaking")   # SAMPLEREP isolation
```

## CI

GitHub Actions runs `R CMD check` on:

| OS | R |
|---|---|
| ubuntu-latest | release, devel |
| windows-latest | release |
| macos-latest | release |

Triggered on push/PR to `main`.
