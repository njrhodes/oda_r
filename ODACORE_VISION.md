# ODACORE_VISION.md — Synthesis and Ground Truth

Maintainer: Nathaniel J. Rhodes, PharmD, MSc
Repository: https://github.com/njrhodes/ODA
Synthesized: 2026-04-15
Authority: This document synthesizes RCORE_CLAUDE.md (operational brief,
April 2026) with the current state of the package as built across this
session. Where the two conflict, the current verified state takes
precedence over the earlier brief.

---

## 1. What odacore is — the ground truth statement

odacore is a pure R reimplementation of the MegaODA.exe / CTA.exe
statistical classification engine. Its sole correctness standard is
the Windows executable. Any divergence from exe output is a bug in
odacore, not an alternative interpretation.

The package removes the Windows x86 binary dependency so that ODA and
CTA analyses can run cross-platform, be embedded in R workflows, be
tested with standard CI tooling, and eventually support extension
(Bayesian priors, Shiny interfaces, NOVO bootstrap) without touching
the exe.

It does not replace the exe for production clinical or research use
until it passes Phase 2 validation. Until then it is a calibrated
reimplementation under active validation.

---

## 2. Current verified state (as of this session)

### Package structure (actual, not spec)
```
rcore/
  DESCRIPTION            (RoxygenNote 7.3.3, no LazyData)
  NAMESPACE              (14 exports, hand-maintained)
  README.md
  odacore.Rproj
  .Rbuildignore
  .github/workflows/R-CMD-check.yml
  R/
    utils.R              %||%, tick(), fmt2/6, p_bucket — defined once
    unioda_core.R        Binary ODA engine (731 lines)
    multioda_core.R      Multiclass ODA engine (1091 lines)
    oda_fit.R            Unified dispatcher C=2→binary, C≥3→multiclass
    cta_core.R           CTA recursive tree engine (558 lines)
  man/                   14 hand-written .Rd files
  tests/testthat/
    helper-odacore.R     Test utilities (always sourced, never exported)
    CTA_DEMO.CSV         Gold fixture data (200 obs, 6 columns)
    test-unioda.R        (10 tests) Binary ODA structural
    test-tie-breaking.R  (8 tests)  MAXSENS→SAMPLEREP→FIRST proven
    test-synthetic-multiclass.R  (5 tests) Algebraically proven SREP
    test-iris.R          (12 tests) MegaODA gold — all 4 iris attributes
    test-oda-fit.R       (7 tests)  Dispatcher correctness
    test-cta.R           (12 tests) CTA structural
    test-cta-gold.R      (15 tests) MegaODA CTA gold — CTA_DEMO dataset
```

**Test count:** 69 tests total. All pass (114/114 on last full run
before CTA gold was added; CTA gold adds 15 more).

**R CMD check:** 0 errors, 1 warning (documentation — resolved with
hand-written .Rd files), 1 note (timestamp, unfixable).

### What was fixed in this session vs the April 2026 brief
The brief listed these as confirmed-but-verify items:
- ✓ No source() calls in R/
- ✓ %||% defined once in utils.R
- ✓ loo="on" inside fold loop removed
- ✓ setNames is base R — no importFrom needed (but kept for safety)
- ✓ iris tests are 12 explicit test_that blocks, not a loop

The brief listed Risk 4 (CTA not implemented) as a known gap.
That gap is now closed: cta_core.R is implemented and gold-tested
against a real MegaODA CTA output file.

### What the brief listed but is NOT yet in the package
- test-fixture-myeloma.R — fixture not yet generated (Phase 1 gate)
- tests/VALIDATION_STATUS.md — population requires myeloma fixture
- harness_utils.R in R/ — intentionally removed (was test-only code
  that leaked into production namespace; test utilities now live in
  tests/testthat/helper-odacore.R where they belong)
- oda_switch.R — replaced by oda_fit.R with cleaner dispatch logic

---

## 3. Selection algorithm — proven invariants

These are not guesses. They were proven by running MegaODA.exe on real
data and reverse-engineering the selection behavior:

```
For each cut position r (combn enumeration order):
  Level 2 (inner) — within this cut, best segment assignment:
    PRIMARY   = MAXSENS (mean PAC in priors-weighted space)
    SECONDARY = SAMPLEREP (L1 distance, RAW counts, always)
    FALLBACK  = FIRST ASSIGNMENT ENUMERATED

Level 1 (outer) — across cut positions:
    PRIMARY   = MAXSENS
    FALLBACK  = FIRST CUT POSITION ENUMERATED (NOT SAMPLEREP)
```

Key proofs:
- iris V1: cuts 6.15 vs 6.25 → both same PAC, raw SREP favors 6.25,
  but MegaODA chose 6.15 → proves first-cut wins across cut positions.
- DAT5: assignments (2,1,3) vs (1,2,3) at same cut → SREP picks (2,1,3)
  → proves SREP operates within a cut position on assignment ties.
- SAMPLEREP is always raw counts. Bug that used priors-weighted space
  was identified and fixed.

---

## 4. Architecture — the five-file structure

```
utils.R       ← shared primitives, no dependencies
    ↓
unioda_core.R ← binary ODA engine, depends on utils
multioda_core.R ← multiclass engine, depends on utils
    ↓
oda_fit.R     ← dispatcher, calls both engines
    ↓
cta_core.R    ← recursive tree, calls oda_fit exclusively
```

No circular dependencies. No cross-engine calls. CTA knows nothing
about the internal structure of either engine — it calls oda_fit()
and reads the standardized result fields.

### Public API (14 exports)

**Entry points:**
- oda_fit(x, y, w, ...)              ← primary; auto-dispatches
- oda_univariate_core(...)           ← binary engine directly
- oda_multiclass_unioda_core(...)    ← multiclass engine directly
- oda_cta_fit(X, y, w, ...)         ← CTA tree

**Prediction:**
- oda_rule_predict(x, rule)
- oda_rule_predict_multiclass(x, rule, boundary)
- predict.cta_tree(object, newdata)
- print.cta_tree(x)
- cta_node_table(tree)

**Metrics:**
- oda_confusion_binary(y, y_pred, w)
- oda_confusion_multiclass(y, y_pred, w)
- oda_mean_pac(sens, spec)
- oda_ess_from_meanpac(mean_pac, chance)
- oda_ess_from_mean(mean_metric, C)

### Return value contract (confusion counts)
- Binary: $confusion = list(TP, TN, FP, FN as raw integer counts;
  sensitivity, specificity as proportions [0,1])
- Multiclass: $confusion = raw integer count matrix (rows=actual, cols=predicted)
  $confusion_wt = priors-weighted matrix (for diagnostics)
- Percentages: $pac, $mean_pac, $pac_by_class are in [0,100]
- Proportions: $sensitivity, $specificity, $accuracy are in [0,1]

---

## 5. Validation status

### Tier 1 — Proven against exe (regression-locked)

| Dataset | Attributes | What's proven |
|---------|-----------|---------------|
| iris    | V1 Sepal.Length | cuts, confusion, LOO confusion exact |
| iris    | V2 Sepal.Width  | cuts, confusion, LOO confusion exact |
| iris    | V3 Petal.Length | cuts, confusion, LOO confusion exact |
| iris    | V4 Petal.Width  | cuts, confusion, LOO confusion exact |
| CTA_DEMO | V2-V6 | 8 split nodes, all node obs counts, confusion exact |
| dat2-dat7 | V1 ordered | SAMPLEREP selection behavior proven |

### Tier 2 — Structurally tested but not exe-validated

| Area | Tests | Status |
|------|-------|--------|
| Binary ODA, all attr types | test-unioda.R | Structural only |
| Multiclass SREP isolation | test-synthetic-multiclass.R | Algebraic proof |
| Dispatcher C=2/C≥3 routing | test-oda-fit.R | Structural |
| CTA stopping rules, depths | test-cta.R | Structural |

### Tier 3 — Not yet validated (known gaps)

| Risk | Description | Gate |
|------|-------------|------|
| RISK-1 | Categorical attribute level ordering | Myeloma fixture |
| RISK-2 | MC p early stopping divergence | Myeloma fixture |
| RISK-3 | Weighted LOO weight distinction | Myeloma fixture |

The myeloma dataset is the canonical Phase 1 validation target.
It has: multiple attributes, case weights (WEIGHT command),
LOO, and categorical attributes — it exercises all three risks
simultaneously. Phase 1 is not complete until myeloma passes.

---

## 6. Phases

### Phase 0 — Complete ✓
Engine implementation and iris gold regression.
- Binary ODA engine proven faithful
- Multiclass ODA engine proven faithful (two-level SAMPLEREP)
- ODA LOO proven faithful (true refit-per-fold)
- CTA engine implemented and gold-tested (CTA_DEMO.CSV)

### Phase 1 — Active: Validate UniODA against myeloma fixture
Generate fixture on Windows, compare all attributes ESS/confusion/LOO.
Fix any discrepancy. One branch per fix.

Gate: All 14 myeloma attributes pass ESS within 0.01% tolerance.

**To begin Phase 1:**
```r
# On Windows with ODA package installed:
library(ODA)
ref_oda <- ODAparse(run=1, mod=1, assign_global=FALSE,
                    base_path=here::here("vignettes/myeloma"))
saveRDS(ref_oda, "rcore/tests/testthat/fixtures/myeloma_oda_parsed.rds")
```
Then add test-fixture-myeloma.R and run devtools::test().

### Phase 2 — Fix UniODA discrepancies
One branch per failing attribute. One commit per fix.
Branch naming: rcore/fix-{issue-desc}

Gate: All myeloma attributes pass + existing iris tests still pass.

### Phase 3 — CTA validation (in progress)
CTA engine is implemented. Gold test against CTA_DEMO exists.
Need: myeloma CTA fixture, weighted CTA fixture.

**To complete Phase 3:**
```r
# On Windows:
ref_cta <- CTAparse(...)
saveRDS(ref_cta, "rcore/tests/testthat/fixtures/myeloma_cta_parsed.rds")
```

### Phase 4 — Weighted CTA and LOO STABLE parity
Case weights in CTA. Weighted ESS (WESS) vs ESS distinction.
LOO STABLE at tree level (currently per-node only).

Gate: Weighted myeloma CTA fixture passes.

### Phase 5 — Extension layer (future)
Not started. Does not touch the validated engine.

**5a. NOVO bootstrap**
Wrapped around oda_fit(). No engine changes.

**5b. Bayesian prior plug-in**
New priors function that replaces oda_apply_priors().
Engine interface unchanged — w is still the weight vector.

**5c. Shiny application**
Takes validated package as dependency. Separate package.
Architecture: odacore (engine) + odaapp (Shiny).

**5d. Exact p-value engine**
Generalization of CYN 1997 to imbalanced and non-directional cases — see inst/docs/theory/.

---

## 7. Rules for future development

These were learned from failures in this session and should not be
re-learned:

**R scoping:**
- <<- from a for() loop body writes to the CALLER'S environment,
  not the function's local environment. for() has no scope.
  Use <- in for() body; use <<- only inside nested functions/closures.
  (This bug caused all iris fits to return ok=FALSE for ~2 sessions.)

**SAMPLEREP:**
- Always raw counts. Never priors-weighted. The objective space
  SREP is a diagnostic field only, never used for selection.

**Confusion matrix:**
- $confusion is always raw integer counts (unit weights).
  $confusion_wt holds the priors-weighted version.
  Never swap these.

**LOO:**
- True refit-per-fold. loo="off" on per-fold calls.
  Never reuse the global rule inside LOO.

**CTA mindenom:**
- Raw observation count, not weighted sum.
  MINDENOM in MegaODA is the OBS column, which is raw.

**MC defaults:**
- CTA node growth mc_iter=5000, ODA standalone mc_iter=25000.
  STOP=99.9 and STOPUP=20 apply to both.
  All four are tunable parameters, not hardcoded constants.

**Tests:**
- Iris gold values are fixed. Any change requires exe confirmation.
- CTA_DEMO gold values are fixed. Any change requires CTA.exe confirmation.
- Never write a test that passes without calling the function.

---

## 8. The single most important near-term test

From the April 2026 brief — this is still the right framing:

When the myeloma fixture exists:
```r
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
This test, passing for all 14 myeloma attributes, completes Phase 1.

---

## 9. What does not change across phases

The engine contract. Once Phase 2 passes:
- oda_fit() signature is stable
- Return field names are stable
- Confusion orientation (rows=actual, cols=predicted) is stable
- ESS is always [0,100], PAC is always [0,100], rates are [0,1]
- cta_tree node fields are stable

Extension phases (5a/5b/5c) add new functions. They do not
modify existing function signatures or return values.

---

## 10. Session discipline (from RCORE_CLAUDE.md, confirmed)

Before any commit:
```r
devtools::load_all()      # no errors
devtools::test()          # all non-skipped tests pass
devtools::check_man()     # no documentation warnings
```

Each session targets ONE of:
- A single correctness fix + its test
- A single new fixture + its test file
- VALIDATION_STATUS.md population
- One extension function

Commit format:
```
type(rcore/scope): description [PHASE-N or RISK-N]

fix(rcore/unioda): correct categorical level ordering [RISK-1]
test(rcore/fixtures): myeloma V7-V14 ESS comparison [PHASE-1]
feat(rcore/cta): weighted ESS vs ESS distinction [PHASE-4]
```
