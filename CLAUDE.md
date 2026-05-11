# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package is

`odacore` is a pure-R reimplementation of the **MegaODA** classification engine. It implements:

- **UniODA** (`oda_univariate_core`): univariate binary-class ODA for ordered, categorical, and binary attributes.
- **MultiODA** (`oda_multiclass_unioda_core`): univariate multiclass ODA for ordered and categorical attributes.
- **`oda_fit`**: unified dispatcher — calls UniODA when C=2, MultiODA when C≥3. This is the entry point CTA nodes should use.
- **`oda_cta_fit`** (`cta_core.R`): CTA engine — implemented and under fixture-parity refinement.

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
└── cta_core.R       — CTA engine: oda_cta_fit(), predict.cta_tree(), helpers
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

**Path-local missingness (CTA):** `predict.cta_tree()` supports `missing_action = c("na", "majority")`.
- `"na"` — canonical: observation is excluded (returns `NA_integer_`) when the split attribute is missing on its actual traversal path.
- `"majority"` — legacy: route missing obs to the split node's majority class.

### NAMESPACE note

`NAMESPACE` is hand-maintained. The harness utilities (`same_int_mat`, `pac_overall`, `loo_match_status`) are exported so test scripts can call them directly, but they are not part of the production inference API.

## Terminology — must match CTA.exe labels exactly

| WEIGHT command | Node-level | Full-tree | Best candidate |
|----------------|-----------|-----------|----------------|
| Not active     | ESS / ESSL | OVERALL ESS | Best ESS |
| Active         | WESS / WESSL | WEIGHTED ESS | Best WESS |

Never call unit-weight scores WESS. Never call weighted scores generic ESS when quoting CTA output.

## Canon fixtures

### CTA_DEMO (no WEIGHT command → ESS)

```
tests/testthat/fixtures/cta_demo/
  CTA_DEMO.CSV
  cta.pgm            — MINDENOM=1
  MODEL1.TXT
  CTA_DEMO_output.txt
  cta_8.pgm          — MINDENOM=8, mc_iter=25000 required for parity
  MODEL8.TXT
  CTA_DEMO__8_output.txt
```

CTA_DEMO MINDENOM=8 canon anchor:
- Exact mc_iter=25000 matches CTA.exe.
- Selected root: V2.
- Enumerated/pruned OVERALL ESS: 68.08%.

### Myeloma (WEIGHT V2 active → WESS)

```
tests/testthat/fixtures/myeloma/
  data.txt
  data.csv
  cta_1.pgm          — MINDENOM=1
  MODEL1.TXT
  myeloma_1_output.txt
  cta_30.pgm         — MINDENOM=30
  MODEL30.TXT
  myeloma_30_output.txt
  cta_56.pgm         — MINDENOM=56
  MODEL56.TXT
  myeloma_56_output.txt
```

Myeloma settings (all MINDENOM variants):
- `EX V2=0` — exclude zero-weight rows.
- `MISSING ALL (-9)` — miss code is -9.
- `WEIGHT V2` — case weights.
- Attributes: V4 V9 V11 V12 V14 V15 V16 V17 V18 V19.

Myeloma canon anchors:

**MINDENOM=1 — enumerated tree (MODEL1.TXT Enumerated section):**
- Root: V14. Child on V14≤0.5 side: V15.
- Classified n = 255 (no missing on V14 path).
- Confusion: [[146, 40], [36, 33]] (actual × predicted, classes 0 and 1).
- OVERALL ESS = 26.32%. WEIGHTED ESS = 27.69%.

**MINDENOM=1 — pruned tree (MODEL1.TXT Pruned section):**
- Root: V17. Child on V17≤0.5 side: V15.
- Classified n = 186 (69 obs have V17=−9, excluded by path-local).
- Confusion: [[92, 43], [21, 30]].
- OVERALL ESS = 26.97%. WEIGHTED ESS = 25.43%.

**MINDENOM=30 — stump (MODEL30.TXT all sections):**
- Root: V17 only (V15 right child would be n=18 < 30; all children fail).
- Classified n = 186.
- Confusion: [[101, 34], [30, 21]].
- OVERALL ESS = 15.99%. WEIGHTED ESS = 16.51%.

**MINDENOM=56 — no tree (MODEL56.TXT all sections):**
- V17 right child n=55 < 56; V14 right child n=44 < 56. All candidates fail.
- Result: leaf node only.

## Tests

Test files in `tests/testthat/`:
- `test-iris.R` — MegaODA.exe gold regression across all 4 iris attributes, K=3 ordered multiclass
- `test-unioda.R` — Binary-class ODA unit tests
- `test-tie-breaking.R` — SAMPLEREP isolation tests
- `test-synthetic-multiclass.R` — Synthetic multiclass scenarios
- `test-oda-fit.R` — Unified dispatcher tests
- `test-cta.R` — CTA path-local prediction deterministic tests
- `helper-odacore.R` — Loaded automatically before every test file

Gold values in tests come from MegaODA.exe output. When a test checks confusion matrices, it uses `confusion_raw()` (raw integer counts, C×C matrix) and `pac_overall()` (percentage of correct classifications).

## Canon references

- `docs/ODA_CANON.md` — canonical `oda_fit()` / UniODA / MultiODA behavior spec.
- `docs/CTA_CANON.md` — canonical `oda_cta_fit()` behavior spec, including ENUMERATE, pruning, LOO, and current gold fixture status.

## Current known-good state

- MINDENOM child-size enforcement: fixed and merged.
- Verbose CTA reporting: merged.
- Public CTA prediction: `predict.cta_tree(missing_action = c("na", "majority"))` implemented.
  - `"na"` is canonical path-local missingness.
  - `"majority"` is legacy compatibility.
- The old full-tree path-local ENUMERATE patch was reverted. **Do not reapply blindly.**
- **Active issue:** myeloma MINDENOM=30 — ENUMERATE scoring uses majority-fallback (`missing_action="majority"` equivalent), which routes 69 V17-missing obs to class 0 and deflates V17's WEIGHTED ESS below V14's, causing wrong root selection. CTA.exe excludes those obs entirely (path-local). Patch not yet applied; fixture arithmetic pre-validation required first.

## Rules for Claude

- Do not run CTA.exe unless explicitly instructed.
- Do not infer from stale fixtures; read the fixture matching the exact MINDENOM setting.
- Do not patch R code before reproducing exact fixture arithmetic with a temp script.
- Do not call unit-weight scores WESS.
- Do not call weighted WESS "ESS" when quoting CTA output.
- Use temp scripts for probes; delete them after reporting.
- Stop after reports when asked — do not proceed to patch without approval.
- Do not touch MINDENOM child-size enforcement.
- Do not touch verbose reporting.

## Implementation policy

Reproduce CTA.exe canonical behavior exactly. Extensions beyond CTA.exe are allowed only behind explicit new options.

Canonical ENUMERATE: evaluate each valid root candidate, grow the CTA.exe-compatible HO-CTA candidate tree below it, compute the full-tree score using CTA.exe's classified/scored universe, and retain best. Path-local missingness is required for prediction/scoring, but do not patch ENUMERATE scoring until CTA node-level fitting matches CTA.exe, especially ordered-cut eligibility.
