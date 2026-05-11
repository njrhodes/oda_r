# ODACORE_VISION.md

Maintainer: Nathaniel J. Rhodes, PharmD, MSc
Repository: https://github.com/njrhodes/oda_rcore
Updated: 2026-05-11

---

## 1. What odacore is

odacore is a pure-R reimplementation of the MegaODA.exe / CTA.exe statistical
classification engine. Where gold executable coverage exists, the correctness
standard is the Windows executable: any divergence from exe output on a covered
fixture is a bug in odacore, not an alternative interpretation.

The package removes the Windows x86 binary dependency so that ODA and CTA
analyses can run cross-platform, be embedded in R workflows, tested with
standard CI tooling, and extended (novometric MDSA, interpretability artifacts)
without touching the exe.

Extensions beyond exe behavior are allowed only when clearly marked as
extensions and must not alter parity behavior on covered fixtures.

---

## 2. Current production state

### Architecture (five-file structure)

```
utils.R          ← shared primitives (%||%, tick(), fmt helpers)
    ↓
unioda_core.R    ← binary ODA engine (oda_univariate_core)
multioda_core.R  ← multiclass ODA engine (oda_multiclass_unioda_core)
    ↓
oda_fit.R        ← unified dispatcher: C=2 → binary, C≥3 → multiclass
    ↓
cta_core.R       ← CTA recursive tree (oda_cta_fit, predict.cta_tree, helpers)
```

No circular dependencies. CTA calls `oda_fit()` and reads standardized result
fields; it knows nothing about the internal structure of either ODA engine.

### Public API

**Entry points:**
- `oda_fit(x, y, w, ...)` — primary dispatcher
- `oda_univariate_core(...)` — binary ODA engine directly
- `oda_multiclass_unioda_core(...)` — multiclass ODA engine directly
- `oda_cta_fit(X, y, w, ...)` — CTA tree

**Prediction / inspection:**
- `predict.cta_tree(object, newdata, missing_action = c("majority", "na"))`
- `print.cta_tree(x)`
- `cta_node_table(tree)`
- `oda_rule_predict(x, rule)`
- `oda_rule_predict_multiclass(x, rule, boundary)`

**Metrics:**
- `oda_confusion_binary`, `oda_confusion_multiclass`
- `oda_mean_pac`, `oda_ess_from_meanpac`, `oda_ess_from_mean`

### Return value contract

- `$confusion` — always raw integer counts (rows = actual, cols = predicted)
- `$confusion_wt` — priors-weighted matrix (diagnostics only)
- `$pac`, `$mean_pac`, `$pac_by_class` — in [0, 100]
- `$sensitivity`, `$specificity`, `$accuracy` — in [0, 1]
- `$ess`, `$wess` — in [0, 100]; WESS only when WEIGHT command is active

### Covered fixture validation status

| Dataset    | MINDENOM | Type               | Status   | Key anchor                        |
|------------|----------|--------------------|----------|-----------------------------------|
| iris V1–V4 | —        | MultiODA           | ✓ green  | All 4 attributes, K=3 ordered     |
| CTA_DEMO   | 1        | CTA (no weights)   | ✓ green  | Root V2, cut 4.5, ESS 52.63%      |
| CTA_DEMO   | 8        | CTA (no weights)   | ✓ green  | mc_iter=25000, ESS 68.08%         |
| myeloma    | 1        | CTA (WEIGHT V2)    | ✓ green  | V14→V15, WESS 27.69%, n=255       |
| myeloma    | 30       | CTA (WEIGHT V2)    | ✓ green  | V17 stump, WESS 16.51%, n=186     |
| myeloma    | 56       | CTA (WEIGHT V2)    | ✓ green  | No tree (all child sizes < 56)    |

**Full test suite: 224/224 passing.**

**devtools::check(): 0 errors / 0 warnings / 1 note** (network timestamp
check — environment issue, not a package issue).

### Key CTA implementation details

**Weighted ordered scan + LOO STABLE gate** (commit 85459a4):
- `.cta_ordered_scan()` selects the rightmost cut where the class-1
  right-branch priors-adjusted PAC > 0.5.
- LOO STABLE requires WESSL = WESS (|delta| ≤ 0.01 pp). Signif T alone is
  insufficient.
- Binary attributes (≤ 2 unique values) and uniform-weight datasets bypass the
  CTA path and use generic ODA.

**Root-only ENUMERATE stump phase** (commit a2e2a9d):
- After the expanded ENUMERATE loop (candidate trees grown below each root),
  a second loop evaluates each root candidate as a stump scored path-locally.
- Observations missing the root attribute are excluded (NA), not majority-routed.
- Root-only stump candidates compete against expanded candidates; best WESS wins.
- This matches CTA.exe Trees 5–7 behavior in MODEL1.TXT.
- **Do not globally apply path-local scoring to expanded ENUMERATE candidates.**
  This was attempted twice and reverted; it incorrectly displaces correct
  expanded trees for MINDENOM=1.

**MINDENOM:**
- Raw child-node row-count admissibility. Unweighted, no priors adjustment.

**`missing_action` in `predict.cta_tree()`:**
- `"majority"` (default) — majority-fallback; used in expanded ENUMERATE scoring.
- `"na"` — canonical path-local missingness; returns NA_integer_ for missing obs.

---

## 3. Selection algorithm — proven invariants

Reverse-engineered from MegaODA.exe runs and locked by fixture tests.

```
For each cut position r (enumeration order):
  Within-cut (inner) — best segment assignment:
    PRIMARY   = MAXSENS  (mean PAC in priors-weighted space)
    SECONDARY = SAMPLEREP  (L1 distance, RAW counts always)
    FALLBACK  = FIRST ASSIGNMENT ENUMERATED

Across cut positions (outer):
    PRIMARY   = MAXSENS
    FALLBACK  = FIRST CUT POSITION ENUMERATED  (not SAMPLEREP)
```

Key proofs:
- iris V1: cuts 6.15 vs 6.25 — same PAC, raw SREP favors 6.25, MegaODA
  chose 6.15 → first-cut wins across cut positions.
- DAT5: assignments (2,1,3) vs (1,2,3) at same cut — SREP picks (2,1,3)
  → SREP operates within a cut position on assignment ties.
- SAMPLEREP is always raw counts. Using priors-weighted space is a bug.

---

## 4. Novometric axioms

### Axiom 1 — Sample size / exact distribution

For binary class × binary attribute, unit-weighted UniODA exact distributions
converge to Fisher's exact test. More general ODA exact p-values are obtained
through Fisher randomization / permutation or require further exact-distribution
theory. Additional theory is needed for sample-size calculation for the ODA
exact distribution.

### Axiom 2 — Structural Decomposition Analysis (SDA)

SDA is used in multiattribute applications to identify the attribute subset
producing the GO-CTA model. It is analogous to PCA but maximizes predictive
accuracy rather than explained variance. It selects attributes successively
over monotonically diminishing sample partitions:

In step 1, EO-CTA is applied and the attribute yielding the minimum D statistic
with exact p < 0.05 is selected. Correctly classified observations are then
removed and the selected attribute is omitted from later steps. The process
repeats on the remaining misclassified observations until all observations in
either class are correctly classified, p > 0.05, or too few observations remain
to satisfy Axiom 1.

D is the parsimony-normalized criterion:

```
D = 100 / ((ESS_or_WESS) / strata_length) - strata_length
```

Where `strata_length` counts terminal leaf endpoints only; internal and
intersection nodes are not counted. MPE defines D as [100 / (ESS / strata)] − strata,
where strata is the number of strata/endpoints. Use WESS in place of ESS only
as odacore weighted-extension language when case weights are declared.

### Axiom 3 — MDSA descendant family

MDSA operates on EO-CTA models configured on the SDA-selected attribute subset.
Starting from MINDENOM 1, the descendant family is traced by computing the
minimum terminal endpoint denominator of the current model and stepping to:

```
Next MINDENOM = current model's minimum terminal endpoint denominator + 1
```

Myeloma example:
- MINDENOM 1 model has minimum terminal endpoint denominator 29 → next MINDENOM 30.
- MINDENOM 30 stump has minimum terminal endpoint denominator 55 → next MINDENOM 56.
- MINDENOM 56 yields no feasible tree → terminate.

No-tree states have no terminal endpoint denominators and terminate the
descendant family. The resulting sequence {MINDENOM 1, 30, 56} is the descendant
family; the minD model is selected from the feasible members.

### Axiom 4 — Reliability

MPE canon: hold-out, LOO/jackknife, bootstrap, or test-retest validity analyses
estimate cross-generalizability of D, ESS, ESP, and related performance indices.

odacore planning policy: reliability status should be represented explicitly;
when a workflow declares reliability as an acceptance criterion, only reliable
complete models should compete. If multiple reliable models exist, theoretical
or empirical use considerations must be weighed.

---

## 5. Phase map

### Phase 0 — Completed: parity stabilization ✓

- UniODA / MultiODA parity on covered fixtures (iris, synthetic tie-breaking).
- Binary CTA fixture parity for CTA_DEMO (MINDENOM 1 and 8) and myeloma
  (MINDENOM 1, 30, 56).
- Weighted ordered scan + LOO STABLE gate.
- Root-only ENUMERATE stump phase.
- Package hygiene: 0 errors / 0 warnings on devtools::check().
- Dev-only theory assets removed from package tracking and build.
- Copilot review instructions established.

### Phase 1 — Completed: documentation and canon alignment ✓

- CLAUDE.md, data-raw/README.md, ODACORE_VISION.md aligned to current state.
- Copilot instructions cover all known parity invariants and scope limits.
- No behavior changes.

### Phase 2 — Next: novometric / MDSA functionalization

Functionalize MDSA across a descendant family of CTA trees.

Initial example: myeloma MINDENOM 1 → 30 → 56.
- Each model/stump/no-tree state records terminal endpoint denominators.
- Next MINDENOM = current model's minimum terminal endpoint denominator + 1.
- No-tree states terminate the descendant family.
- Family supports minD model selection across feasible members.

Requires safe no-tree state representation and a clean family container object.

### Phase 3 — Interpretability artifacts

- Tree plots with meaningful node labels.
- Confusion matrices (raw and priors-weighted).
- ESS/WESS/PAC metrics with MC/LOO/reliability annotations.
- Within-tree and across-family comparisons.
- Safe no-tree reporting.
- Clear distinction between UniODA/MultiODA single-rule outputs, CTA tree
  outputs, and MDSA family outputs.

### Phase 4 — Permutation / bootstrap / 95% CI performance

- Review current MC permutation and CI implementation.
- Compare with ODA repository R implementation.
- Define statistical target before optimizing.
- Preserve seed/reproducibility policy.
- Correctness tests before speed tests.

### Phase 5 — Multiclass CTA extension

- Most far-reaching extension; no gold executable benchmark exists.
- Must be explicitly documented as extension behavior, not MegaODA/CTA parity.
- Fixture strategy must be synthetic/property-based, not gold-exe parity.
- Must not alter binary CTA parity behavior.

---

## 6. Invariants that must not regress

These were learned from debugging and are locked by fixture tests.

**SAMPLEREP:** Always raw counts. Never priors-weighted. Using the objective
space SREP for selection is a bug.

**Confusion matrix:** `$confusion` is always raw integer counts. `$confusion_wt`
holds the priors-weighted version. Never swap.

**LOO:** True refit-per-fold. Never reuse the global rule inside a LOO fold.

**MINDENOM:** Raw observation count (OBS column in CTA.exe output), not weighted
sum. Do not use priors-adjusted or classified-only counts.

**MC defaults:** CTA node growth mc_iter=5000; standalone ODA mc_iter=25000.
STOP/STOPUP are tunable. These are parameters, not hardcoded constants.

**R scoping in loops:** `<<-` from a `for()` body writes to the caller's
environment. `for()` has no scope. Use `<-` in loop bodies; use `<<-` only
inside nested closures.

---

## 7. Reference index

| Document | Purpose |
|----------|---------|
| `CLAUDE.md` | Operational guardrails for Claude Code sessions |
| `README.md` | Public package overview and quick-start |
| `docs/ODA_CANON.md` | ODA engine canonical behavior spec |
| `docs/CTA_CANON.md` | CTA engine canonical behavior spec |
| `docs/CTA_ORDERED_CUT_AUDIT.md` | Weighted ordered scan / LOO STABLE audit evidence |
| `.github/copilot-instructions.md` | AI code review policy (10 rules) |
| `data-raw/README.md` | Fixture provenance (MegaODA.exe and CTA.exe run settings) |
