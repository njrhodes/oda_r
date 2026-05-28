# DOCS_INDEX.md

Navigation index for odacore documentation. **This is the authoritative starting
point for understanding active specs, design history, and implementation status.**

Do not use `prompts/` files as active roadmap. They are historical archive
material. This file supersedes them.

---

## Active Canon Docs

These files define "correct" behavior. A divergence from a covered fixture or
from these specs is a **bug**, not an alternative interpretation.

| File | Covers | Status |
|------|--------|--------|
| [ODA_CANON.md](ODA_CANON.md) | `oda_fit()` / UniODA / MultiODA: degeneracy policy, objective, ordered/categorical/binary attributes, weights, missing values, MC, LOO, DIRECTION parameters (Chapter 2 binary and Chapter 4 multiclass/categorical) | Active canon |
| [CTA_CANON.md](CTA_CANON.md) | `cta_fit()` / `oda_cta_fit()`: degeneracy policy, node-level construction, HO-CTA, EO-CTA/ENUMERATE, pruning (Sidak-Bonferroni + max-accuracy), weights, missing values, LOO, MC, gold fixtures (Myeloma + CTA_DEMO) | Active canon |

**Degeneracy boundary (added 2026-05-27):** Post-pruning degeneracy gate is
now production. CTA never returns a tree where all terminals predict the same
class — such candidates are rejected in the ENUMERATE loop before entering
best-tree competition. `degen = TRUE` is not an option for CTA or LORT (recursive CTA).

---

## CTA / MDSA Translation and Reporting Docs

These files cover the on-demand reporting pipeline that runs on top of fitted
`cta_tree` objects. All translation is on-demand; the lean-fit invariant
(no training X/y stored at fit time) must be preserved.

| File | Covers | Status |
|------|--------|--------|
| [CTA_TRANSLATION_STACK.md](CTA_TRANSLATION_STACK.md) | Navigation map: lean-fit principle, function map for all reporting/translation functions, pipeline overview | Active reference — first production version complete |
| [myeloma-cta-translation.md](myeloma-cta-translation.md) | Worked walkthrough of the full reporting stack using myeloma MINDENOM=1/30/56; all values computed from actual package runs | Completed reference / worked example |
| [CTA_ORDERED_CUT_AUDIT.md](CTA_ORDERED_CUT_AUDIT.md) | Audit evidence for weighted ordered-cut selection rule and LOO STABLE gate; MPE.pdf anchors and myeloma V4/V15 empirical evidence | Active audit / implementation reference |

---

## LORT / SORT / GORT Method Taxonomy

| File | Covers | Status |
|------|--------|--------|
| [LORT_SORT_GORT_TAXONOMY.md](LORT_SORT_GORT_TAXONOMY.md) | **Agent handoff contract:** LORT/SORT/GORT definitions, method relationships, confirmed LORT object metadata, reserved names, what agents may/may not modify, future SORT/GORT task templates, canon boundary | **Active — read before any recursive CTA task** |

Key distinctions:
- **LORT** (implemented): greedy local min-D recursion; `method = "lort"`; `global_optimization = FALSE`; `sda_anchored = FALSE`.
- **SORT** (reserved): SDA-anchored sequential recursive CTA/MDSA; requires SDA source object; replaces generic "staged CTA" as a specific method name.
- **GORT** (future design only): global recursive search over all configurations; not implemented; reserved namespace.

---

## Recursive CTA Design Docs

| File | Covers | Status |
|------|--------|--------|
| [CTA_ORT_DESIGN.md](CTA_ORT_DESIGN.md) | Finalized design for `cta_fit(recursive = TRUE)` (LORT): API, MDSA per-level MINDENOM, single seed / RNG streaming, right-then-left traversal, recursion guards, object structure, prediction/plot semantics, testing plan | Active design — implementation complete in `R/cta_ort.R` |
| [ORT_SELECTION_METHODS.md](ORT_SELECTION_METHODS.md) | Decision memo: three-role method stack (EO-CTA/MDSA → SDA staging → imbalanced recursive CTA); canon uncertainty notes; private rare-event test results; stale/invalid MINDENOM=117 case documented | Active design memo — **not implementation spec**; do not code from it without explicit approval |

**ORT_SELECTION_METHODS §5.1 note:** The MINDENOM=117 root member is explicitly
marked `INVALID — exposed as degenerate (2025-05)`. Sections §5.2–5.4 discuss
this result as a historical probe, not a current production outcome. The
degeneracy gate (2026-05-27) prevents this from recurring.

**Terminology note (2026-05-28):** `recursive = TRUE` is now named LORT internally.
"ORT" may appear in legacy code — the names `cta_ort`, `predict.cta_ort`,
`print.cta_ort`, `summary.cta_ort`, `plot.cta_ort`, `ort_plot_data()`, and
`cta_ort_node_table()` are **legacy compatibility names for LORT**.
They are retained for backward compatibility. New docs and agent instructions
should use LORT/SORT/GORT; do not introduce new bare-`ort` public names.
See `LORT_SORT_GORT_TAXONOMY.md §Legacy API Naming` for the full table and rules.

---

## SDA Docs

| File | Covers | Status |
|------|--------|--------|
| [SDA_AUTO_SDA_PLAN.md](SDA_AUTO_SDA_PLAN.md) | Comprehensive design for Sequential Discriminant Analysis and auto-SDA: canon anchor, SDA modes (legacy `unioda_max_ess` vs. MPE novometric `min_d`), object contract, prediction semantics, CTA/ORT interop, auto-SDA plan logic, generalized staged terminology, test plan, implementation slices | Active design — SDA-1 through SDA-4B complete; weighted/staged adjustment deferred |

---

## SORT / Staged CTA Workflow Docs

| File | Covers | Status |
|------|--------|--------|
| [STAGED_CTA_WORKFLOW_PLAN.md](STAGED_CTA_WORKFLOW_PLAN.md) | Lessons-ingested design for **SORT** (Sequentially Optimal Recursive Trees): FORCENODE (node-level) vs. staged EX CTA (workflow-level) differences, artifact inventory, model comparison table, proposed `staged_cta_workflow` object, implementation sequence, cleanup plan for data-raw artifacts | Active design — SCTA-0 (evidence capture) complete; SCTA-1+ deferred; **no implementation yet** |

---

## Vision and Architecture Doc

| File | Covers | Status |
|------|--------|--------|
| [ODACORE_VISION.md](ODACORE_VISION.md) | Package scope; public API; return value contract; fixture validation status; selection algorithm invariants; novometric axioms (Axiom 1–4); full phase map (Phase 0–4); invariants that must not regress; reference index | Living doc — authoritative for package scope and phase structure; **update when production state changes** |

---

## Stale / Archive Materials

These exist in the repository but are **not active roadmap**. Do not drive
implementation decisions from them.

| File | Status | Note |
|------|--------|------|
| `prompts/prod-roadmap-cta-mdsa-2026-05-19.md` | **Archive / stale** | Historical phase handoff prompt from 2026-05-19. Superseded by current CLAUDE.md, ODACORE_VISION.md, and this index. Do not treat as authoritative. |
| Any other files in `prompts/` | **Archive / historical** | Phase-handoff prompts from prior sessions. Not current roadmap. |

---

## Current Production Checkpoints

| Checkpoint | Status | Commit |
|-----------|--------|--------|
| UniODA / MultiODA parity (iris, synthetic) | complete | Phase 0 |
| CTA_DEMO MINDENOM=1 and 8 parity | complete | Phase 0 |
| Myeloma MINDENOM=1/30/56 parity | complete | Phase 0 |
| Weighted ordered scan + LOO STABLE gate | complete | 85459a4 |
| Root-only ENUMERATE stump phase | complete | a2e2a9d |
| CTA translation stack first production version | complete | f0328e2 |
| CTA graphics v1 (`cta_plot_data`, `plot.cta_tree`) | complete | f0328e2 |
| D statistic + endpoint denominators | complete | Phase 2D |
| MDSA descendant family (`cta_descendant_family`, `cta_family_table`) | complete | Phase 2E |
| SDA-4B novometric_min_d mode | complete | 2026-05-27 |
| Post-pruning degeneracy gate (CTA + ORT + family) | complete | fbaaa85 |
| CRAN check: 0 errors / 0 warnings / 1 note (clock) | current | fbaaa85 |

---

## Current Next-Slice Order

```
A. Docs index / coherence pass                      <- current (this file)
B. Synthetic degeneracy regression fixture          <- next hardening slice
C. CTA endpoint map + print clarity
D. Balance diagnostics v1 design
E. Graphics v3 plot-data contract
F. SDA -> CTA/ORT anchor
G. Staged CTA workflow implementation
H. Weighted/staged adjustment design
I. Vignettes / pkgdown
J. Release hardening
```

Slices D–J are ordered by dependency; do not skip ahead without completing
upstream prerequisites.

---

## External Vignette / Demo Candidate Note

**Credit-card fraud dataset:** A publicly available rare-event imbalanced dataset
(e.g., Kaggle credit card fraud) is a **future optional external demo candidate**
for rare-event degeneracy and balance workflows. It is:

- not package data
- not a test fixture
- must not be downloaded during CRAN checks or automated tests
- appropriate only as a user-facing vignette (Slice I) under an explicit
  network-skip guard, after the balance-diagnostics API exists (Slice D)

Do not reference it in tests or as a canon anchor.

---

## How to Use This Index

- **Finding canon behavior:** Read `ODA_CANON.md` and `CTA_CANON.md` first.
  `ODACORE_VISION.md` covers invariants and phase context.
- **Understanding the reporting stack:** `CTA_TRANSLATION_STACK.md` then
  `myeloma-cta-translation.md`.
- **LORT / recursive CTA:** `CTA_ORT_DESIGN.md` then `ORT_SELECTION_METHODS.md`
  (memo only, not spec).
- **SDA / auto-SDA:** `SDA_AUTO_SDA_PLAN.md`.
- **Staged CTA workflows:** `STAGED_CTA_WORKFLOW_PLAN.md`.
- **What not to read for current roadmap:** Anything in `prompts/`.
