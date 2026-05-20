# Production Roadmap: CTA / MDSA / odacore

Date: 2026-05-19
Main branch: `eca9395 docs: regenerate Rd for new_cta_family and cta_descendant_family`

---

## Current State

### Closed issues

| Issue | Description |
|---|---|
| #3 | Canonical CTA top-three A×B×C ENUMERATE — commit `01ab7fd` |
| #4 | Canonical PRUNE (Sidak-Bonferroni + maximum-accuracy) — commit `7927cfe` |
| #5 | ODA cleanup (binary MC raw p, multiclass MC mean-PAC p, LOO Fisher direction) |
| #7 | Multiclass categorical MC p parity fixed |

### Open issues

| Issue | Description | Priority |
|---|---|---|
| #8 | CTA MC TARGET/STOP operational parity | **Next — prod blocker** |
| #6 | DIRECTIONAL support for categorical and ordered ODA | Deferred |

### What is solid (parity-tested, regression-locked)

- UniODA binary and multiclass core
- MC p-value: binary raw proportion, multiclass mean-PAC — both correct
- LOO Fisher direction behavior
- CTA weighted ordered scan + LOO STABLE gate
- Canonical Sidak-Bonferroni PRUNE
- Canonical top-three A×B×C ENUMERATE
- Endpoint denominators, strata, D statistic — read-only, no recomputation
- CTA descendant family: MDSA chain loop, min-D selection, termination handling
- Fixture coverage: CTA_DEMO (MINDENOM 1, 8), Myeloma (MINDENOM 1, 30, 56), iris K=3 — all passing
- Full test suite: 531 pass / 0 fail / 2 pre-existing DIRECTIONAL skips
- Check: 0 errors / 0 warnings / 1 pre-existing timestamp note

### What is newly landed and needs hardening

- `cta_descendant_family()`: passes myeloma chain but has no `print.cta_family()`, no `summary.cta_family()`, no `@examples`
- A×B×C ENUMERATE: canonical and correct, but adds an uncharacterized MC call cost multiplier
- PRUNE: correct but no diagnostic record of what was pruned

### What is still experimental or partial

- DIRECTIONAL: binary ordered partially wired; categorical/multiclass deferred (#6)
- LOO p-value contract: design notes only
- Bootstrap CI: not started
- Report table objects: none exist
- Print/summary for `cta_tree`, `cta_family`: minimal
- Graphics: not started

---

## Production Blockers (ranked)

### 1. #8 CTA MC TARGET/STOP operational parity — CRITICAL

MC STOP is implemented in both MC helpers (`oda_mc_p_value()` in `unioda_core.R` and `.cta_mc_ordered()` in `cta_core.R`). Parameters `mc_target`, `mc_stop`, `mc_stopup` are propagated correctly through the call chain. The unresolved question is empirical: with `STOP = 99.9` and `mc_iter = 5000`, how many CTA candidate MC calls in a full A×B×C ENUMERATE actually stop early?

`oda_mc_p_value()` has `min_check = 50L` / `check_every = 50L` — stopping is not even attempted for the first 49 iterations. For attributes with moderate or near-null ESS in branch partitions, the 99.9% confidence threshold may not be crossed until thousands of iterations.

With A×B×C ENUMERATE, `.all_cands()` runs MC+LOO for all candidate attributes at both left and right branches for every root. For myeloma (~10 attributes, ~10 root candidates), this is potentially O(10 roots × 18 branch-attribute calls × 5000 iterations) = ~900,000 MC iterations for candidate generation alone, before any `.ho_grow()` depth-3+ calls.

Every downstream feature — descendant family, bootstrap, report tables — inherits this cost if unresolved.

### 2. Test tiering — HIGH

Full suite runs 8–10 min. Slow components: `cta-node-selection` (2.5–3.5 min), `fixture-cta-demo` (2–3 min), `cta` (1–1.5 min). Without tiering, every small edit costs 10 minutes to validate. Adding bootstrap CI tests without tiering makes the suite unusably slow.

### 3. Public API/docs gaps — MEDIUM

`cta_descendant_family()` has an Rd file but no `@examples`. No `print.cta_family()` or `summary.cta_family()`. A user who calls the function gets an R list with no display contract.

### 4. Report table objects — MEDIUM (prerequisite for bootstrap and graphics)

Bootstrap CI needs stable output fields to attach CIs to. Graphics need stable endpoint structures to render. Designing bootstrap and graphics before report tables means retrofitting — the wrong order.

### 5. Bootstrap CI design — LOW-MEDIUM (deferred)

Open design questions: row resampling vs weighted resampling with `WEIGHT V2`, miss-code behavior during resamples, seed increment policy across family steps, performance under CTA family (each resample runs the full MINDENOM chain). Not blocked by #8, but impractical to implement until #8 is resolved and runtime is sane.

### 6. #6 DIRECTIONAL — LOW unless release explicitly promises DIRECTIONAL parity

Binary ordered direction is partially implemented. Categorical DIRECTIONAL (`DIRECTIONAL < 1 2 3 4` on TABLE input) requires hard design decisions and could consume weeks. Two fixture tests are cleanly labeled as skipped/deferred. Does not block any other current feature.

---

## Phase Roadmap

### Phase 0 — #8 MC TARGET/STOP audit and fix

**Goal:** Determine whether CTA MC STOP is operationally correct and whether runtime is canon cost or implementation gap.

This is a read/instrument/report task before any code change.

**Evidence to inspect:**
- `oda_mc_p_value()` stopping logic and `min_check`/`check_every` guards (`unioda_core.R`)
- `.cta_mc_ordered()` stopping logic (`cta_core.R`)
- `iter_used` distribution for myeloma MINDENOM=1 under full ENUMERATE
- Proportion of MC calls that reach all `mc_iter` iterations vs stop early

**Three possible outcomes:**

| Decision | Finding | Response |
|---|---|---|
| A | STOP is wrong or missing | Targeted fix; canon fixtures must still pass |
| B | STOP is correct but diagnostics absent | Add `iter_used`/`ge_count`/`stop_reason` to diagnostic output |
| C | STOP is correct; runtime is true canon cost | Mark slow tests; implement tiering; document expectation |

**Files touched:** `R/unioda_core.R` (diagnostic helper), `R/cta_core.R` (diagnostic + possible fix)

**Exit criteria:** All 10 questions in #8 answered. Decision documented. If fix applied, all canon fixtures pass and check is clean.

---

### Phase 0B — Test tiering (concurrent with or immediately after Phase 0)

**Goal:** Fast synthetic/unit test path (< 30s) separate from slow canon/integration tests.

**Approach:** Annotate slow tests with `skip_if(Sys.getenv("ODACORE_TESTS") == "fast")` or move canon/integration tests to a dedicated prefix (`test-integration-*.R`) with documented runner convention.

**Files touched:** `tests/testthat/test-fixture-cta-demo.R`, `test-fixture-myeloma-cta.R`, `test-cta.R` (family tests), `test-cta-node-selection.R`

**Exit criteria:** Fast path < 30s. Full path still runs everything. Convention documented.

---

### Phase 1 — CTA public API stabilization

**Goal:** `cta_descendant_family()` has a usable public surface.

**Deliverables:**
- `print.cta_family()`: compact display of family chain, status, ESS/WESS, D, min-D member highlighted
- `summary.cta_family()`: structured list; reads from stored fields only — no refitting
- `@examples` in `cta_descendant_family.Rd` that runs in < 5s
- Safe no-tree display in all methods

**Files touched:** `R/cta_family.R`, `man/cta_descendant_family.Rd`, `NAMESPACE`

**Exit criteria:** `print(family)` and `summary(family)` work; `@examples` passes in check; 0 new warnings.

---

### Phase 2 — Report table objects

**Goal:** Stable report table data frames that bootstrap and graphics will consume.

**Candidate functions:**
- `cta_endpoint_summary(tree)`: one row per leaf — path, prediction, n, class counts, target-class proportion, denominator
- `cta_family_table(family)`: formatted family comparison data frame
- `oda_rule_table(fit)`: attribute, rule type, cut, direction, ESS/WESS, p_mc, LOO status

**Constraint:** Every column must come from stored fit-time state. No silent recomputation. No LOO p-value columns until design is settled.

**Files touched:** New `R/cta_reports.R`, Rd files for new functions

**Exit criteria:** Each function returns a clean data frame; `print.data.frame` renders usably; no field requires recomputation.

---

### Phase 3 — Bootstrap CI

**Goal:** Bootstrap confidence intervals for ESS/WESS and PAC at ODA and CTA levels.

**Design decisions to settle before code:**
1. Row resampling vs weighted resampling for `WEIGHT V2` datasets
2. Miss-code behavior during bootstrap resamples
3. Seed increment policy across CTA family steps
4. Performance guardrail: B=50 for tests, B=1000 for production

**Prerequisite:** #8 resolved and runtime sane. Each bootstrap resample runs a full CTA fit; CTA family bootstraps run the full MINDENOM chain per resample.

**Files touched:** New `R/oda_bootstrap.R`, Rd files

**Exit criteria:** `oda_boot_ci(fit, B=1000, seed=1L)` reproducible; `cta_family_boot_ci(family, B=50)` returns per-member CIs; no-tree member returns NA CI.

---

### Phase 4 — Graphics / reporting polish

**Goal:** `plot.cta_tree()`, `plot.cta_family()` (D curve over MINDENOM), optional Mermaid text export.

**Prerequisite:** Report table objects from Phase 2. Graphics consume `cta_endpoint_summary()` output, not raw `tree$nodes`.

**Files touched:** New `R/cta_plot.R`, Rd files

**Exit criteria:** `plot(tree)` and `plot(family)` render safely on all fixture trees including no-tree; Mermaid text export is a separate optional function, not the rendering engine.

---

### Phase 5 — Deferred canon expansion (#6 DIRECTIONAL)

**Goal:** Categorical DIRECTIONAL (Bowker, political affiliation Chapter 4 examples) and ordered multiclass DIRECTIONAL. Promote Chapter 4 skip-stubs to passing tests.

**Risk:** Categorical `DIRECTIONAL < 1 2 3 4` on TABLE input is ordinal, not nominal categorical scanning — design decisions are genuinely hard and could consume weeks. Binary ordered direction is partially in place but the LOO Fisher `alternative` propagation has an open design question.

**Files touched:** `R/unioda_core.R`, `R/multioda_core.R`, `R/oda_fit.R`

**Exit criteria:** Both Chapter 4 skip-stubs pass; all existing canon fixtures unaffected; check clean.

---

## Production MVP Definition

```
- oda_fit() + predict.oda_fit() + summary.oda_fit()
- oda_cta_fit() + predict.cta_tree() + cta_node_table()
- cta_descendant_family() + print.cta_family() + summary.cta_family()
- cta_endpoint_summary() returning a clean data frame
- #8 resolved (MC STOP correct + runtime sane)
- Test tiering in place
- 0 errors / 0 warnings in devtools::check()
- No DIRECTIONAL promise
- No bootstrap promise
```

## Beta / Deferred (explicitly not MVP)

```
- Bootstrap CI
- Graphics / Mermaid export
- #6 DIRECTIONAL categorical/multiclass
- LOO p-value implementation (design not settled)
- SDA / auto_SDA
- Multiclass CTA
```

---

## Next Single Task

**#8 iter_used / TARGET / STOP diagnostic for Myeloma MINDENOM=1 under full A×B×C ENUMERATE.**

Instrument or add a non-exported diagnostic helper to capture per-candidate-MC-call:
- node / position
- attribute
- cutoff / target
- p_mc
- ge_count
- iter_used
- stop reason (early-STOP, early-STOPUP, or max-iter)
- final Signif T/F decision

Run under canon-like settings:

```r
mc_iter    = 5000L
mc_target  = 0.05
mc_stop    = 99.9
mc_stopup  = 20
loo        = "stable"
prune_alpha = 0.05
```

Report the distribution of `iter_used`. Determine whether most calls run all 5000 iterations or stop early. Make decision A/B/C from Phase 0. Implement fix only if warranted. Canon fixtures must still pass.

This is the only task that reduces production risk before any other work begins.
