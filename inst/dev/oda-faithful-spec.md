# odacore — Engine Spec and Readiness Assessment

## Current state: 96/96 tests passing. R CMD check: 0 errors.

---

## 1. Selection spec (MegaODA-faithful, proven from executable evidence)

```
For each cut position r (combn enumeration order):
  Level 2 — within this cut, select best segment assignment:
    PRIMARY = MAXSENS (mean PAC in priors-weighted space)
    → SECONDARY = SAMPLEREP (L1 distance in RAW count space)
    → FIRST ASSIGNMENT ENUMERATED

Level 1 — across cut positions, compare cut-level winners:
    PRIMARY
    → FIRST CUT POSITION ENUMERATED  (not SAMPLEREP)
```

Key proof points:
- iris V1: cuts 6.15 vs 6.25 have identical primary PAC and unique assignments.
  Raw SREP favors 6.25 (0.040 vs 0.093), but MegaODA chose 6.15 → first-cut wins.
- DAT5: at cut (1.5,2.5), assignments (2,1,3) vs (1,2,3) tie on primary →
  SAMPLEREP picks (2,1,3). Within-cut assignment tie, not across-cut.
- SAMPLEREP always uses raw counts regardless of priors_on setting.

---

## 2. API gaps to fix before CTA

### Gap 1: No unified entry point (blocking for CTA)
CTA nodes need a single `oda_fit()` that works for C=2 and C≥3 automatically.
Currently two separate functions must be called depending on class count.

### Gap 2: Parameter name `miss_codes` vs `missing_code`
`miss_codes` is canonical. `missing_code` now accepted as alias in unioda_core.
**Before CTA:** also add alias to multioda_core, and document canonical name.

### Gap 3: LOO parameter proliferation
Five LOO parameters on multioda_core will multiply when CTA adds tree-level LOO.
Consolidate into `loo_opts = list(grid_mode, priors_mode, ...)` before CTA.

---

## 3. CTA = chained ODA models

CTA differs from ODA in exactly three ways:
1. **Outer loop over attributes** (best-attribute selection per node)
2. **Stopping rules**: alpha_split, mindenom (min node n), max_depth, ess_min
3. **Tree structure** for the model object

The engine already provides everything needed at the node level.
CTA needs `oda_fit()` dispatcher + `oda_cta_node()` recursive builder.

### Must-fix list before CTA branch

1. `oda_fit()` — single public dispatcher (routes on C; ~30 lines)
2. `miss_codes` alias everywhere  
3. `loo_opts` list consolidation in multioda_core
4. Move internal functions off export list (oda_best_ordered_multiclass_partition,
   oda_loo_for_rule, oda_loo_multiclass, oda_loo_multiclass_ordered, harness utilities)
