# DOCS_INDEX.md

Navigation index for oda documentation.

**Start here:** [STATUS.md](STATUS.md) - current production state, what is implemented,
what is deferred.

---

## Canon Docs

These files define correct behavior. A divergence from a covered fixture or
from these specs is a **bug**, not an alternative interpretation.

| File | Covers |
|------|--------|
| [ODA_CANON.md](ODA_CANON.md) | `oda_fit()` / UniODA / MultiODA: degeneracy policy, objective, ordered/categorical/binary attributes, weights, missing values, MC, LOO |
| [CTA_CANON.md](CTA_CANON.md) | `cta_fit()` / `oda_cta_fit()`: ENUMERATE, pruning, LOO STABLE, MINDENOM, gold fixtures (myeloma + CTA_DEMO) |

---

## Method and Workflow Docs

| File | Covers |
|------|--------|
| [LORT_SORT_GORT_TAXONOMY.md](LORT_SORT_GORT_TAXONOMY.md) | LORT/SORT/GORT definitions, canon boundary, novometric canon hierarchy, output classification audit standard |
| [SDA_ANCHOR_CONTRACT.md](SDA_ANCHOR_CONTRACT.md) | `sda_anchor` object contract: schema, API (`as_sda_anchor`, `validate_sda_anchor`), task hooks |
| [FIT_OBJECT_EVIDENCE_CONTRACT.md](FIT_OBJECT_EVIDENCE_CONTRACT.md) | Fitted-object evidence wiring: ESS/WESS/D/confusion across public methods |

LORT distinctions:
- **LORT** (implemented): adjacent workflow-layer composition using canon CTA/MDSA components; greedy local min-D; not Chapter 12 canon.
- **SORT** (reserved): SDA-anchored sequential recursive CTA/MDSA; requires SDA source object; not implemented.
- **GORT** (future): global recursive search; not implemented.

---

## Reporting, Translation, and Balance Docs

| File | Covers |
|------|--------|
| [CTA_TRANSLATION_STACK.md](CTA_TRANSLATION_STACK.md) | Lean-fit principle, function map for all reporting/translation functions, pipeline overview |
| [COVARIATE_BALANCE_CONTRACT.md](COVARIATE_BALANCE_CONTRACT.md) | Balance diagnostics API: `oda_balance_table`, `smd_balance_table`, `cta_balance_table`, plot-data transforms |
| [GRAPHICS_V3.md](GRAPHICS_V3.md) | Graphics v3 function reference, ggplot2 dependency policy, no-fitting-inside-renderers rule |

---

## Implementation Audit Docs

| File | Covers |
|------|--------|
| [CTA_ORDERED_CUT_AUDIT.md](CTA_ORDERED_CUT_AUDIT.md) | Audit evidence: weighted ordered-cut selection rule and LOO STABLE gate; MPE.pdf anchors and myeloma V4/V15 empirical evidence |

---

## Worked Examples

| File | Covers |
|------|--------|
| [examples/myeloma-cta-translation.md](examples/myeloma-cta-translation.md) | Worked example: full reporting stack using myeloma MINDENOM=1/30/56; values computed from actual package runs |

---

## Current Production Checkpoints

| Checkpoint | Status |
|-----------|--------|
| UniODA / MultiODA parity (iris, synthetic) | complete |
| CTA_DEMO MINDENOM=1 and 8 parity | complete |
| Myeloma MINDENOM=1/30/56 parity | complete |
| Weighted ordered scan + LOO STABLE gate | complete |
| Root-only ENUMERATE stump phase | complete |
| CTA translation stack first production version | complete |
| CTA graphics v1/v3 (tree, balance, LORT, love plot) | complete |
| D statistic + endpoint denominators | complete |
| MDSA descendant family | complete |
| SDA-4B novometric_min_d mode | complete |
| Post-pruning degeneracy gate (CTA + ORT + family) | complete |
| LORT/SORT/GORT taxonomy | complete |
| Power/sample-size planning: `oda_power()` / `oda_sample_size()` (unit-weighted binary 2x2) | complete |
| ODA/SMD/CTA covariate balance tables + plot data | complete |
| SDA anchor object + task hooks (Slice O) | complete |
| Production tools: readiness checks + ODA/LORT propensity (Slice Q) | complete |
| Slice S: evidence wiring across public methods | complete |
| Slice T: live-fire rehearsal; public contract alignment | complete |
| Full fast suite: FAIL 0 / WARN 0 | current |
| CRAN check: 0 errors / 0 warnings / 1 note (clock) | current |

---

## Deferred Items

| Item | Note |
|------|------|
| SORT (`sort_fit`) | Reserved - requires SDA anchor; not implemented |
| GORT (`gort_fit`) | Reserved - future design only |
| Weighted SDA | Deferred to SDA-5 |
| SDA-derived propensity weights | Deferred; SDA produces stage order, not propensity strata |
| Multiclass / weighted / CTA power | Not yet scoped |
| Multiclass CTA | Future extension |

