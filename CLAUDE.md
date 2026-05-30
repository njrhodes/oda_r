# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package is

`odacore` is a pure-R reimplementation of the **MegaODA** classification engine. It implements:

- **UniODA** (`oda_univariate_core`): univariate binary-class ODA for ordered, categorical, and binary attributes.
- **MultiODA** (`oda_multiclass_unioda_core`): univariate multiclass ODA for ordered and categorical attributes.
- **`oda_fit`**: unified dispatcher — calls UniODA when C=2, MultiODA when C≥3. This is the entry point CTA nodes should use.
- **`cta_fit`** (`oda_fit.R`): public CTA entry point — thin wrapper around the internal engine; preferred public name.
- **`oda_cta_fit`** (`cta_core.R`): CTA engine (internal name; retained for backward compatibility; prefer `cta_fit()` in user-facing code).

Covered fixtures are tested against MegaODA.exe / CTA.exe outputs. Extension behavior must be explicitly marked and must not disturb covered parity. Roadmap details are in `docs/ODACORE_VISION.md`; engine canon is in `docs/ODA_CANON.md`, `docs/CTA_CANON.md`, and `docs/CTA_ORDERED_CUT_AUDIT.md`.

## Commands

```r
devtools::document()   # regenerate docs from roxygen (NAMESPACE is hand-maintained)
devtools::install(".") # install locally
```

### Test tiers

Tests are gated by the `ODACORE_TEST_TIER` environment variable (`tests/testthat/helper-test-tier.R`).

| Tier | When to use |
|------|-------------|
| unset / `cran` | CRAN-safe unit/contract tests only — default when env var is unset |
| `fast` | Fast dev loop — same slow-skip behavior as cran |
| `smoke` | Production gate before any CTA/MDSA/reporting/graphics commit |
| `full` | Release gate before release tags or deep canon parity work |

```bash
# CRAN-safe / default check (no env var — runs only CRAN-safe tests):
Rscript --vanilla -e "devtools::check(vignettes=FALSE)"

# Fast developer loop (skip all slow canon fixtures):
ODACORE_TEST_TIER=fast Rscript --vanilla -e "devtools::test(reporter='progress')"

# CTA/MDSA/reporting/graphics production smoke (required before every commit):
ODACORE_TEST_TIER=smoke Rscript --vanilla -e "devtools::test(reporter='progress')"

# Full canon / release gate (required before release tags):
ODACORE_TEST_TIER=full Rscript --vanilla -e "devtools::test(reporter='progress')"
```

**Do not run `ODACORE_TEST_TIER=full` by default.** Use it only before release tags or when
explicit canon parity work requires it.

### Production-smoke gate

**Required before pushing any CTA, MDSA, reporting, or graphics changes.**
`ODACORE_TEST_TIER=fast` (or unset/cran) skips all myeloma and CTA canon fixtures and is
NOT sufficient for production validation.

**Important:** `devtools::test(filter=...)` matches **test file names**, not skip-reason strings.
`"cta-family"` and `"cta-myeloma-chain"` are skip-reason labels, not file names. The slow tests
they gate live inside `test-cta.R` (matched by `filter='cta'`).

```bash
# Targeted smoke (myeloma + confusion + CTA + graphics):
ODACORE_TEST_TIER=smoke Rscript --vanilla -e \
  "devtools::test(filter='fixture-myeloma-cta|cta-confusion-table|cta|cta-plot|cta-loo-gate|cta-ordered-scan', reporter='progress')"

# Full broad smoke (all files):
ODACORE_TEST_TIER=smoke Rscript --vanilla -e "devtools::test(reporter='progress')"
```

## Architecture

Files are sourced alphabetically by R's package loader, so `utils.R` is always first.

```
R/
├── utils.R          — %||%, tick(), fmt helpers, .validate_case_weights() (loaded first)
├── unioda_core.R    — Binary-class engine: oda_univariate_core() and helpers
├── multioda_core.R  — Multiclass engine: oda_multiclass_unioda_core() and helpers
├── oda_fit.R        — Unified dispatcher: oda_fit() routes on C
├── oda_s3.R         — ODA S3 methods: predict/print/summary for oda_fit
├── cta_core.R       — CTA engine: oda_cta_fit(), predict.cta_tree(), helpers
├── cta_s3.R         — CTA S3 + translation layer: summary/print/accessors/staging/
│                      propensity/assign_endpoints/observation_weights
└── cta_family.R     — MDSA family: cta_descendant_family(), cta_family_table(),
                       summary/print methods for cta_family
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

**`degen` flag (UniODA/MultiODA only):** When `FALSE` (default), all C classes must appear in the predicted labels — degenerate solutions (where only one class is ever predicted) are rejected. When `TRUE`, degenerate solutions are allowed and `priors_on` is forced `FALSE`. `degen = TRUE` is available only for `oda_fit()`, `oda_univariate_core()`, and `oda_multiclass_unioda_core()`. It does not exist for `cta_fit()` or `oda_cta_fit()`. CTA never produces degenerate trees: a candidate split that predicts only one class is ineligible at the node gate. A CTA or recursive CTA (LORT) result in which all endpoints predict the same class is not a valid tree — it is a model failure and must surface as `no_tree`.

**`boundary_mode`:** Controls how ordered cut values map to segment boundaries.
- `"right_closed"` — matches MegaODA.exe golden outputs used in tests.
- `"megaoda_halfopen"` — alternative convention.

**Path-local missingness (CTA):** `predict.cta_tree()` supports `missing_action = c("na", "majority")`.
- `"na"` — canonical: observation is excluded (returns `NA_integer_`) when the split attribute is missing on its actual traversal path.
- `"majority"` — legacy: route missing obs to the split node's majority class.

**CTA lean-fit invariant:** `oda_cta_fit()` stores only what is needed for on-demand reporting:
- Tree nodes (split rules, child pointers, `n_obs`)
- `training_confusion` — C×C raw integer confusion from the final selected tree
- Per-leaf `class_counts_raw` and `class_counts_weighted` vectors

Nothing else is stored at fit time: no training `X` or `y`, no row indices, no endpoint membership, no staging tables, no propensity tables, no observation-weight tables, no cached reporting artifacts of any kind.
**Do not add fit-time storage to `cta_tree` for reporting or translation convenience.** All translation/reporting artifacts are computed on explicit function call. This invariant must be maintained as the translation stack grows.

**Global case-weight guard:** `.validate_case_weights(w, n, arg = "w")` in `utils.R` is the single canonical weight validator. It rejects non-numeric, wrong length, NA/NaN, Inf/-Inf, zero, and negative weights. It is wired into `oda_fit()`, `oda_cta_fit()`, and `cta_descendant_family()`. Do not add redundant weight checks elsewhere.

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
- `docs/CTA_ORDERED_CUT_AUDIT.md` — audit evidence and canon for weighted ordered-cut selection and LOO STABLE law; MPE.pdf anchors and myeloma V4/V15 empirical evidence.
- `docs/CTA_TRANSLATION_STACK.md` — navigation map for the CTA reporting and translation pipeline (lean-fit principle, function map, pipeline overview).
- `docs/myeloma-cta-translation.md` — canonical myeloma CTA walkthrough with actual computed values covering MINDENOM=1/30/56.

## Current known-good state

- MINDENOM child-size enforcement: fixed and merged.
- Verbose CTA reporting: merged.
- Public CTA prediction: `predict.cta_tree(missing_action = c("na", "majority"))` implemented.
  - `"na"` is canonical path-local missingness.
  - `"majority"` is legacy compatibility.
- The old full-tree path-local ENUMERATE patch was reverted. **Do not reapply blindly.**
- **Global case-weight guard (production):** `.validate_case_weights()` in `utils.R` wired into `oda_fit()`, `oda_cta_fit()`, `cta_descendant_family()`. Rejects non-numeric, wrong-length, NA/NaN, Inf/-Inf, zero, and negative weights.
- **CTA translation stack — first production version (complete):**
  - `cta_endpoint_summary()` — endpoint structure and denominators
  - `cta_endpoint_counts()` — endpoint × class raw/weighted counts
  - `cta_staging_table()` — target-class count, proportion, odds, perfect-endpoint flags, adjusted values
  - `cta_propensity_weights()` — stabilized weights (raw counts); adjusted variant; `undefined_empirical` flag
  - `cta_assign_endpoints()` — row → endpoint traversal on demand
  - `cta_observation_weights()` — row → endpoint-level weight assignment on demand
  - `cta_confusion_table()` — actual × predicted confusion with PAC/ESS from stored `training_confusion`
  - `cta_family_table()` — MINDENOM family summary with D-statistic and min-D selection
  - All on-demand: lean-fit invariant maintained (see Architecture section above).
- **Weighted ordered predictor node-level fitting (production, commit 85459a4):**
  - `.cta_ordered_scan()`, `.cta_mc_ordered()`, `.cta_full_fit_ordered()` are implemented helpers in `R/cta_core.R`.
  - `.full_fit_one()` dispatches to the CTA path for non-uniform weights + binary class + >2 unique non-missing attribute values.
  - Binary attributes (≤2 unique values) and uniform-weight datasets always use generic ODA — CTA path bypassed.
  - The CTA path is coupled with the LOO STABLE gate:
    - Full-model WESS uses `.cta_ordered_scan()` (rightmost cut with class-1 right-branch priors-adjusted PAC > 0.5).
    - LOO STABLE uses true generic ODA per-fold refit (`oda_loo_for_rule()`).
    - Candidates are rejected if WESSL ≠ WESS (|delta| > 0.01 pp). Signif T alone is insufficient.
  - Regression locks (all passing):
    - myeloma Node 2 (V17≤0.5, n=131): V4 rejected as UNSTABLE (CTA WESS=34.89% vs generic LOO WESSL=38.52%, |Δ|=3.63 pp); V15 STABLE and selected.
    - myeloma Node 4 (V17≤0.5 AND V15≤0.5, n=113): V4 scores at cut=359.80/WESS=34.90% but is not selected (Signif F or LOO UNSTABLE).
    - cta_demo root remains V2, cut=4.5, ESS=52.63% (uniform weights → CTA path bypassed).
  - All CTA tests passed after integration: 100/100.
  - Full test suite passed after integration.
- **Root-only ENUMERATE stump phase (production, commit a2e2a9d):**
  - After the expanded ENUMERATE loop, a second loop evaluates each root candidate as a stump (root split + two leaves), scored path-locally.
  - Path-local scoring: `.apply_cand(A_cand, seq_len(n))$y_pred` already carries `NA_integer_` for obs whose root attribute is in `miss_codes`. WESS is computed over `ok = !is.na(stump_preds)` only — missing-root obs are excluded, not majority-routed.
  - Root-only candidates compete against expanded candidates; the overall best WESS wins.
  - Regression locks (all passing):
    - MINDENOM=1: expanded V14→V15 WESS=27.69% > all root-only stumps → V14→V15 tree selected. Classified n=255, confusion [[146,40],[36,33]].
    - MINDENOM=30: V17 stump (n_classified=186, WESS=16.51%) > V14 stump (n_classified=255, WESS=14.06%) → V17 stump selected. Confusion [[101,34],[30,21]], OVERALL ESS=15.99%.
    - MINDENOM=56: no admissible root candidate (all child sizes < 56) → leaf node only.
    - cta_demo: uniform weights → CTA path bypassed entirely → root=V2, cut=4.5, ESS=52.63%. Unaffected by stump phase.
  - **Do not globally apply path-local scoring to expanded ENUMERATE candidates.** Doing so causes V17's deeper MC-grown tree to score > 27.69% for MINDENOM=1, incorrectly displacing V14→V15. This was attempted twice and reverted.

## Recursive CTA method taxonomy — LORT / SORT / GORT

Full contract in `docs/LORT_SORT_GORT_TAXONOMY.md`. Summary rules:

- **LORT** = Locally Optimal Recursive Tree. Current `cta_fit(recursive = TRUE)`.
  Greedy local min-D at each node. No lookahead. No SDA anchor.
  **Adjacent workflow-layer composition** using canon CTA/MDSA components at each
  node; not MPE Chapter 12 canon. `lort_fit()` is public workflow API, not a canon
  claim. Legacy ORT naming is compatibility only.
  Object metadata: `ort_settings$method == "lort"`, `global_optimization == FALSE`, `sda_anchored == FALSE`.
- **SORT** = Sequentially Optimal Recursive Tree. Reserved. Requires SDA source object.
  Do not implement SORT inside `recursive = TRUE`. Reserved entry: `sort_fit()`.
- **GORT** = Globally Optimal Recursive Tree. Reserved. Future design only. Reserved entry: `gort_fit()`.

Agent rules for recursive CTA tasks:
1. Any task touching recursive CTA must explicitly say LORT, SORT, or GORT.
2. LORT task: do not add lookahead, SDA anchoring, or global search to `recursive = TRUE`.
3. SORT task: require an SDA result or equivalent anchor before starting. Do not begin without it.
4. GORT task: do not modify LORT or SORT behavior. Require approved design doc first.
5. Do not use "ORT" alone in new docs. Prefer LORT/SORT/GORT.
6. Do not change `ort_settings$method` from `"lort"` for current `recursive = TRUE` fits.
7. Do not export or stub SORT/GORT reserved names without explicit approval.
8. Do not rename S3 class `cta_ort` or break existing dispatch in this or future slices.

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
- Do not globally apply path-local scoring to expanded ENUMERATE candidates. Only root-only stump candidates are scored path-locally. Expanded candidates use majority-fallback `.predict_all()`. This split is canon (MODEL1.TXT Trees 2–4 vs Trees 5–7).
- Do not touch the weighted ordered scan / LOO STABLE gate (`cta_ordered_scan`, `.cta_mc_ordered`, `.cta_full_fit_ordered`, `.full_fit_one`) unless a regression in the node-selection tests proves it is involved.
- Do not add fit-time storage to `cta_tree` for reporting or translation convenience. The lean-fit invariant (tree nodes + `training_confusion` + per-leaf `class_counts_raw`/`class_counts_weighted` only) must be preserved. All translation/reporting artifacts are computed on explicit function call.
- Do not add redundant weight validation outside `.validate_case_weights()` in `utils.R`. It is already wired into all public fit entrypoints.
- Do not accept or label a CTA or recursive CTA (LORT) result as valid when all terminal endpoints predict the same class. That is a degenerate solution. CTA never produces degenerate trees; such a result means no admissible split existed and the correct output is `no_tree`. `degen = TRUE` does not exist for CTA or LORT — it is a UniODA/MultiODA-only option.

## Implementation policy

Reproduce CTA.exe canonical behavior exactly. Extensions beyond CTA.exe are allowed only behind explicit new options.

Canonical ENUMERATE: evaluate each valid root candidate, grow the CTA.exe-compatible HO-CTA candidate tree below it, compute the full-tree score using CTA.exe's classified/scored universe, and retain best. Path-local missingness is required for prediction/scoring.

**ENUMERATE architecture (production):** Two distinct phases, matching MODEL1.TXT:
1. **Expanded phase** (Trees 2–4): grow full HO-CTA below each root candidate; score with majority-fallback `.predict_all()`. Do not apply path-local scoring here.
2. **Root-only stump phase** (Trees 5–7): score each root candidate as a stump path-locally (missing-root obs excluded). Runs after the expanded phase; competes for overall best WESS.

Node-level weighted ordered fitting and LOO STABLE gate are implemented and passing (commit 85459a4). Root-only ENUMERATE stump phase is implemented and passing (commit a2e2a9d). All canon fixtures pass: 765/765.

**Remaining known risk — MC stochasticity at Node 4 (Tree 4):**
- Node 4 = V17≤0.5 ∧ V15≤0.5, n=113. CTA.exe reports V14 Signif F (WESS=19.19%).
- R's MC may accept V14 (p≈0.03) in some runs, causing R to grow V17→V15→V14 (3 splits) vs CTA.exe's V17→V15 (2 splits) for the ENUMERATE expanded candidate.
- The root-only stump phase resolves the MINDENOM=30 fixture without depending on Node 4 behavior.
- This risk only matters if a future fixture requires exact sub-tree parity for the MINDENOM=1 expanded candidate (Tree 4). No current test asserts Tree 4 internal structure.
