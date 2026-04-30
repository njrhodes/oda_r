# RCORE_CLAUDE.md — Claude Code Brief: odacore Pure R Package

Maintainer: Nathaniel J. Rhodes, PharmD, MSc
Repository: https://github.com/njrhodes/ODA
Branch: rcore/phase1-validate-unioda (off processx-modernization)
Last updated: 2026-04-05

> **Authority:** This document defers to CONTRACT.md → RCORE.md → CLAUDE.md
> in that order. When in conflict, the higher document governs.
> This document contains operational instructions specific to the odacore
> pure R package. It does not restate invariants already in CONTRACT.md.

---

## CRITICAL: Read this first

This document was written based on a package scaffold built in a design
session. The scaffold in `rcore/` may be out of date with the current
state of the uploaded GPT/Gemini source files.

Before writing any code, Claude Code must:

1. Read `CONTRACT.md` fully
2. Read `RCORE.md` fully
3. Read `CLAUDE.md` fully
4. **Diff the current `rcore/R/` files against what is described here**
5. If the current files have evolved, the current files take precedence
   over any code examples in this document

Do not assume this document reflects the current state of the repo.
Verify before acting.

---

## What odacore is

A pure R package implementing the ODA and CTA algorithms without the
Windows x86 exe dependency. It lives in `rcore/` within the main ODA
repository and is developed on `rcore/*` branches.

The exe (`MegaODA.exe`, `CTA.exe`) remains the ground truth. Every
function in this package must produce output that matches exe output
within the tolerances in `RCORE.md`. Divergence is a bug in odacore,
not an alternative interpretation.

---

## Package structure

```
rcore/
  DESCRIPTION
  NAMESPACE
  README.md
  odacore.Rproj
  .gitignore
  .github/workflows/R-CMD-check.yml
  R/
    utils.R           %||%, tick(), fmt helpers — defined exactly once
    unioda_core.R     Binary-class UniODA engine
    multioda_core.R   Multiclass UniODA engine
    harness_utils.R   Diagnostic harness utilities
    oda_switch.R      Binary/multiclass dispatch router
  tests/
    VALIDATION_STATUS.md    Per-attribute exe comparison table
    testthat/
      helper-odacore.R             Auto-loaded, dual-mode loader
      test-unioda.R                Binary-class unit tests
      test-tie-breaking.R          MAXSENS→SAMPLEREP→FIRST tests
      test-synthetic-multiclass.R  Multiclass consistency tests
      test-iris.R                  Iris gold regression (all 4 attrs)
      test-fixture-myeloma.R       Exe fixture tests (skip until fixtures)
      fixtures/
        README.md
        myeloma_oda_parsed.rds   (generated on Windows — not committed yet)
        myeloma_cta_parsed.rds   (generated on Windows — not committed yet)
```

---

## Source files: origin and known issues

The R/ files originated from GPT-generated code across two passes and
one Gemini refinement pass. The following issues were confirmed and
should have been fixed in the scaffold — verify they are actually fixed
before proceeding:

### Confirmed fixed in scaffold (verify still true)

| Issue | File | Fix |
|---|---|---|
| `source()` calls in package context | multioda_core.R, oda_switch.R | Removed — package loads via namespace |
| `%||%` defined 4 times | all files | Moved to utils.R, removed elsewhere |
| `loo="on"` inside per-fold loop | harness_utils.R | Changed to `loo="off"` — halves O(n²) to O(n) |
| `stats::setNames` false conflict | multioda_core.R | Reverted — `setNames` is base R, not stats |
| `local({}) for` loop in iris test | test-iris.R | Replaced with 12 explicit named test_that blocks |
| `where="package:..."` in test helper | helper-odacore.R | Replaced with asNamespace() + source fallback |

### Known remaining issues (not yet fixed — work items)

See the validation framework section below.

---

## What functions exist and what they do

### R/utils.R
- `` `%||%` `` — null coalesce, defined once
- `tick()` — progress bar for long loops
- `fmt2()`, `fmt6()` — numeric formatters for comparison
- `identical_fmt2()`, `identical_fmt6()` — formatted equality checks
- `p_bucket()` — bins p-value into reporting bucket (<.001, <.01, etc.)

### R/unioda_core.R — binary class only
- `oda_univariate_core()` — main entry point, binary class
- `oda_clean_xy()` — per-attribute missing value removal
- `oda_resolve_attr_type()` — auto-detect ordered/categorical/binary
- `oda_apply_priors()` — inverse-frequency class reweighting
- `oda_rule_side()`, `oda_rule_predict()` — apply a rule to data
- `oda_confusion_binary()` — TP/TN/FP/FN from predictions
- `oda_make_blocks_ordered()` — cumulative sweep structure
- `oda_mc_p_value()` — Monte Carlo with ITER/TARGET/STOP/STOPUP
- `oda_loo_for_rule()` — jackknife LOO for a fixed rule
- `oda_apply_primary_secondary()` — tie-breaking selector

### R/multioda_core.R — multiclass (C > 2)
- `oda_multiclass_unioda_core()` — main entry, C > 2
- `oda_enforce_weighting_policy()` — weighted categorical LOO guard
- `oda_apply_priors_multiclass()` — class-normalize weights
- `oda_confusion_multiclass()` — C×C confusion matrix
- `oda_best_ordered_multiclass_partition()` — exhaustive segment search
- `oda_rule_predict_multiclass()` — apply multiclass rule
- `oda_mc_p_value_multiclass()` — MC for multiclass
- `oda_loo_multiclass()`, `oda_loo_multiclass_ordered()` — multiclass LOO
- `tick()` — progress for enumeration loop

### R/harness_utils.R — diagnostic only, not for production
- `confusion_raw()`, `confusion_weighted()` — build confusion matrices
- `extract_confusion_matrix()` — safe extractor from various formats
- `extract_core_loo_conf_raw_int()` — pull integer LOO confusion from fit
- `pac_overall()`, `mean_pac_from_conf()`, `ess_pac_from_conf()` — metrics
- `same_int_mat()`, `same_num_mat()` — matrix equality with tolerance
- `harness_loo_refit_ordered_raw()` — independent LOO recomputation
- `rule_is_nondegenerate()` — confirms all C classes predicted
- `safe_fit_mode()` — robust wrapper for fold-level refits
- `check_fold()`, `fold_debug_audit()` — fold-level diagnostics

### R/oda_switch.R
- `oda_univariate_dispatch()` — routes to binary or multiclass core

---

## Known correctness risks (from RCORE.md)

These are the areas where the R core is most likely to diverge from exe
output. Every one must be validated before the function is trusted.

### Risk 1 — Categorical attribute level ordering (HIGHEST RISK)
**Location:** `unioda_core.R` categorical branch, `multioda_core.R`
**Issue:** The GPT implementation orders categorical levels by weighted
class-1 rate before sweeping. The exe may use a different convention.
Any ordering difference produces different cutpoints and different ESS.
**Test:** `test-iris.R` V2 (Sepal.Width is the only categorical-adjacent
attribute in the iris gold set). The myeloma fixture is the definitive test.
**Do not fix without fixture validation.**

### Risk 2 — MC p-value early stopping divergence
**Location:** `oda_mc_p_value()` in `unioda_core.R`
**Issue:** STOP/STOPUP thresholds use a binomial CI approximation that
may diverge from the exe at small iteration counts.
**Test:** `test-fixture-myeloma.R` p-value bucket comparison.
**Acceptable behavior:** Same bucket result with mc_iter=25000 and matching seed.

### Risk 3 — Weighted LOO weight distinction
**Location:** `oda_loo_for_rule()` in `unioda_core.R`
**Issue:** LOO refit must use case weights, not objective weights.
The myeloma fixture is weighted with LOO — this is the critical test.
**Test:** `test-fixture-myeloma.R` LOO confusion comparison.

### Risk 4 — CTA core is NOT implemented
**What exists:** `cta_core.R` in the uploaded files is a usage harness
that `source()`s `oda_cta_core.R` which was not uploaded and does not
exist in the package. The full CTA tree recursion — greedy top-down
growth, alpha pruning, exhaustive enumeration, LOO STABLE per node — is
not yet implemented.
**Do not attempt CTA implementation without explicit instruction.**
**Do not create placeholder CTA functions that silently return NULL.**
If asked to implement CTA, follow the build sequence in RCORE.md Phase 3.

---

## Iris gold values

These are confirmed from MegaODA.exe output. They are used in test-iris.R.
If you modify any function and the iris tests fail, it is a regression.

```
V1 (Sepal.Length): cuts=[5.45, 6.15], segs=[1,2,3]
  train: PAC=74.67%, conf=[[45,5,0],[6,28,16],[1,10,39]]
  LOO:   PAC=73.33%, conf=[[45,5,0],[6,28,16],[1,12,37]]

V2 (Sepal.Width):  cuts=[2.95, 3.35], segs=[2,3,1]
  train: PAC=59.33%, conf=[[31,2,17],[1,34,15],[5,21,24]]
  LOO:   PAC=59.33%, conf=[[31,2,17],[1,34,15],[5,21,24]]

V3 (Petal.Length): cuts=[2.45, 4.75], segs=[1,2,3]
  train: PAC=95.33%, conf=[[50,0,0],[0,44,6],[0,1,49]]
  LOO:   PAC=94.00%, conf=[[50,0,0],[0,44,6],[0,3,47]]

V4 (Petal.Width):  cuts=[0.80, 1.65], segs=[1,2,3]
  train: PAC=96.00%, conf=[[50,0,0],[0,48,2],[0,4,46]]
  LOO:   PAC=95.33%, conf=[[50,0,0],[0,48,2],[0,5,45]]
```

All iris tests run with: priors_on=TRUE, degen=FALSE, K=3,
boundary=right_closed, loo_grid_mode=refit, loo_priors_mode=fold.

---

## Validation framework

### Phase 1 (current phase): Validate UniODA against myeloma fixture

**What Phase 1 requires:**
1. Generate `tests/testthat/fixtures/myeloma_oda_parsed.rds` on Windows
2. For each attribute in MODEL1.OUT: compare ESS, cutpoint, confusion, LOO
3. Record results in `tests/VALIDATION_STATUS.md`
4. Fix any discrepancies — one branch per failing attribute

**The fixture generation script** (run on Windows with ODA package installed):
```r
# Run from the ODA repository root on Windows Intel
library(ODA)
ref_oda <- ODAparse(
  run          = 1,
  mod          = 1,
  assign_global = FALSE,
  base_path    = here::here("vignettes/myeloma")
)
dir.create("rcore/tests/testthat/fixtures", showWarnings=FALSE)
saveRDS(ref_oda, "rcore/tests/testthat/fixtures/myeloma_oda_parsed.rds")
message("Fixture written.")
```

**What happens when you run tests without fixtures:**
`test-fixture-myeloma.R` skips cleanly with message:
`"Myeloma exe fixture not found — run fixture generation script on Windows first"`

### Phase 2: Fix UniODA discrepancies

One branch per failing attribute. One commit per fix. The failing attribute
name goes in the branch name: `rcore/fix-categorical-ordering-V2`.

### Phase 3: CTA implementation

Do not start until Phase 2 passes. CTA build sequence is in RCORE.md.

---

## Claude Code session rules for this package

### Session scoping
Each session targets ONE of:
- A single correctness fix with test coverage
- A single new utility function with test
- Phase 1 fixture generation script
- VALIDATION_STATUS.md population (after fixtures exist)

### Required before any commit to rcore/

```r
setwd("rcore/")
devtools::load_all()        # no errors
devtools::test()            # all non-skipped tests pass
devtools::check_man()       # no documentation warnings
```

The fixture tests will skip on Linux/macOS — that is correct behavior,
not a failure. They only run when fixtures are present.

### What Claude Code must not do

- Must not add `source()` calls anywhere in R/
- Must not define `%||%` anywhere except utils.R
- Must not implement CTA without explicit instruction
- Must not create fake passing tests that don't actually call the functions
- Must not modify iris gold values without exe confirmation
- Must not claim a function is correct without fixture comparison
- Must not call `loo="on"` inside a per-fold refit loop (O(n²) bug)
- Must not combine weighted and unweighted confusion matrices
- Must not treat `loo_suppressed=TRUE` as equivalent to `loostab=0`

### Commit message convention (same as CLAUDE.md, rcore scope)

```
type(rcore/scope): description [RISK-N or PHASE-N]

Examples:
  fix(rcore/unioda): correct categorical level ordering [RISK-1]
  test(rcore/fixtures): add myeloma fixture comparison for V7-V14 [PHASE-1]
  feat(rcore/harness): add generate_fixtures.R script
  fix(rcore/harness): remove nested loo="on" in fold loop [PERF]
```

### Branch naming for rcore work

```
rcore/phase1-validate-unioda      (current integration branch)
  └── rcore/fix-{issue-desc}      (one per fix)
  └── rcore/feat-{feature-desc}   (new functions)
  └── rcore/test-{test-desc}      (test additions only)
```

PRs from `rcore/*` branches target `rcore/phase1-validate-unioda`.
The phase branch merges to `processx-modernization` when Phase 1 passes.

---

## The one test that matters most

When the myeloma fixture exists, this is the single most important test:

```r
# V7 (or whichever attribute has the largest weighted ESS in the fixture)
# must match exe within 0.01%

test_that("myeloma V7: ESS matches exe within 0.01%", {
  ref    <- readRDS(test_path("fixtures", "myeloma_oda_parsed.rds"))
  data   <- read.csv(here::here("vignettes/myeloma/ODA/1/inputs/data.csv"))
  result <- oda_univariate_core(
    x = data[["V7"]], y = data[[cls_col]], w = data[[wt_col]],
    attr_type = "ordered", priors_on = TRUE,
    mcarlo = FALSE, loo = "off"
  )
  ref_ess <- as.numeric(gsub("[^0-9.-]", "",
    ref$model[ref$model$Attribute == "V7", "Overall.ESS."]))
  expect_equal(result$ess, ref_ess, tolerance = 0.01)
})
```

If this single test passes for all 14 attributes in the myeloma fixture,
Phase 1 is complete. Everything else is scaffolding.
