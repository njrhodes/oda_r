# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package is

`odacore` is a pure-R reimplementation of the **MegaODA** classification engine. It implements:

- **UniODA** (`oda_univariate_core`): univariate binary-class ODA for ordered, categorical, and binary attributes.
- **MultiODA** (`oda_multiclass_unioda_core`): univariate multiclass ODA for ordered and categorical attributes.
- **`oda_fit`**: unified dispatcher — calls UniODA when C=2, MultiODA when C≥3. This is the entry point CTA (Classification Tree Analysis) nodes should use.

All results are tested for exact parity against MegaODA.exe golden outputs.

## Commands

```r
devtools::test()                          # run all tests
devtools::test(filter = "iris")           # run a single test file (e.g., test-iris.R)
devtools::test(filter = "tie-breaking")   # run test-tie-breaking.R
devtools::check()                         # full R CMD check
devtools::document()                      # regenerate docs from roxygen (NAMESPACE is hand-maintained)
devtools::install(".")                    # install locally
```

## Architecture

Files are sourced alphabetically by R's package loader, so `utils.R` is always first.

```
R/
├── utils.R          — %||%, tick(), fmt helpers (loaded first)
├── unioda_core.R    — Binary-class engine: oda_univariate_core() and helpers
├── multioda_core.R  — Multiclass engine: oda_multiclass_unioda_core() and helpers
├── oda_fit.R        — Unified dispatcher: oda_fit() routes on C
└── harness_utils.R  — Diagnostic/parity harness (not production API; exported for tests)
```

### Key internal concepts

**`tick()`** (`utils.R`): integer quantization at 1e9 precision. Used for deterministic floating-point tie comparisons throughout the multiclass search. Must never return `NA_integer_`.

**Tie-breaking spec (MegaODA-faithful):**
```
PRIMARY   = MAXSENS   (overall PAC in priors-weighted objective space)
SECONDARY = SAMPLEREP (L1 distance: predicted vs. observed class frequencies — always uses raw counts, regardless of priors_on)
TERTIARY  = FIRST IDENTIFIED (enumeration order)
```
In multiclass ordered rules: SAMPLEREP operates **within** a cut position (across segment assignments), not **across** cut positions. Across cuts, first-enumerated cut wins.

**Priors weighting:** When `priors_on = TRUE`, each class's weights are normalized so every class contributes equally to the objective. For binary class, this is in `oda_apply_priors()`; for multiclass, in `oda_apply_priors_multiclass()`. ESS and confusion matrices are computed on both weighted (`confusion_wt`) and raw (`confusion`) counts — the returned `$confusion` field always holds raw integer counts.

**LOO (leave-one-out):** True refit-per-fold — each fold drops one observation and refits the full model on n−1 cases, then predicts the held-out case. There is no global rule reuse. Weighted categorical LOO is explicitly forbidden.

**Monte Carlo p-value:** Fisher randomization with Clopper-Pearson early stopping (`STOP`/`STOPUP` confidence bounds). ESS of permuted labels compared against observed ESS.

**Rule object structure** (consistent across UniODA and MultiODA):
- `$rule$type`: `"ordered_cut"` | `"binary_map"` | `"nominal_cut"` | `"multiclass_ordered"` | `"multiclass_nominal"`
- UniODA ordered: `$rule$cut_value`, `$rule$direction`
- MultiODA ordered: `$rule$cut_values` (vector), `$rule$seg_classes` (integer vector)

**`degen` flag (MultiODA):** When `FALSE` (default), all C classes must appear in the predicted labels. When `TRUE`, degenerate solutions are allowed and `priors_on` is forced `FALSE`.

**`boundary_mode`:** Controls how ordered cut values map to segment boundaries.
- `"right_closed"` — matches MegaODA.exe golden outputs used in tests.
- `"megaoda_halfopen"` — alternative convention.

### NAMESPACE note

`NAMESPACE` is hand-maintained. The harness utilities (`same_int_mat`, `pac_overall`, `loo_match_status`) are exported so test scripts can call them directly, but they are not part of the production inference API.

## Tests

Test files in `tests/testthat/`:
- `test-iris.R` — MegaODA.exe gold regression across all 4 iris attributes, K=3 ordered multiclass
- `test-unioda.R` — Binary-class ODA unit tests
- `test-tie-breaking.R` — SAMPLEREP isolation tests
- `test-synthetic-multiclass.R` — Synthetic multiclass scenarios
- `test-oda-fit.R` — Unified dispatcher tests
- `helper-odacore.R` — Loaded automatically before every test file; duplicates key harness functions so they work under both `devtools::test()` and `R CMD check`

Gold values in tests come from MegaODA.exe output. When a test checks confusion matrices, it uses `confusion_raw()` (raw integer counts, C×C matrix) and `pac_overall()` (percentage of correct classifications).

## Upcoming work (from `inst/dev/oda-faithful-spec.md`)

The next major feature is **CTA (Classification Tree Analysis)**, which chains ODA models as tree nodes. Before that branch:
1. ~~`oda_fit()` is the required dispatcher~~ — done.
2. ~~`miss_codes` alias needs to be consistent across both engines~~ — done (`missing_code` alias added to `oda_multiclass_unioda_core`).
3. ~~LOO parameters on `oda_multiclass_unioda_core` should be consolidated into an `loo_opts` list~~ — done.
4. ~~Several currently-exported internal functions should be moved off the export list~~ — already clean; none of the listed functions were in NAMESPACE.

Implementation policy:
First reproduce CTA.exe canonical behavior exactly.
Extensions beyond CTA.exe are allowed only behind explicit new options.

Canonical:
- ENUMERATE = evaluate valid combinations in top three nodes.
- Retain best tree by full-tree WESS.
- Do not enumerate deeper in canonical mode.

Future extension:
- deeper enumeration may be implemented as a separate non-canonical mode,
  e.g. enumerate_depth > 3 or search = "beam"/"global".
