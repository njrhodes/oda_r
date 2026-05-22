# CTA Reporting and Translation Stack

This file maps the on-demand reporting and translation functions available for
a fitted `cta_tree` object.  It is a navigation aid, not a behavior
specification.  For canonical engine behavior see `docs/CTA_CANON.md`; for a
fully worked numeric example see `docs/myeloma-cta-translation.md`.

---

## Lean-fit principle

`oda_cta_fit()` stores only what is needed to reproduce predictions and
compute reporting artifacts on demand:

- Tree nodes (split rules, child pointers, `n_obs`)
- `training_confusion` — C×C matrix from the training pass
- Per-leaf `class_counts_raw` and `class_counts_weighted` vectors

Nothing else is stored at fit time:

- No training `X` or `y`
- No row indices or observation membership
- No endpoint membership table
- No staging, propensity, or observation-weight tables
- No cached reporting artifacts of any kind

Every function in the pipeline below computes its output on demand from the
stored leaf counts or from caller-supplied `newdata`.  The `cta_tree` object
is never mutated by any reporting function.

---

## Pipeline overview

```
oda_cta_fit(X, y, w, ...)
  |
  +-- cta_endpoint_summary()       endpoint structure and denominators
  |
  +-- cta_endpoint_counts()        endpoint × class raw/weighted counts
  |     |
  |     +-- cta_staging_table()    target-class propensity, proportion, odds
  |     |
  |     +-- cta_propensity_weights()   stabilized weights per endpoint × class
  |
  +-- cta_assign_endpoints(newdata)    row → endpoint (traversal on demand)
  |     |
  |     +-- cta_observation_weights(newdata, y)   row → endpoint-level weight
  |
  +-- cta_confusion_table()        actual × predicted confusion matrix
  |
  +-- [cta_descendant_family()]
        |
        +-- cta_family_table()     MINDENOM family summary with D-statistic
```

All functions except `cta_assign_endpoints` and `cta_observation_weights`
operate solely on stored leaf counts — no `newdata` required.

---

## Function map

| Function | Row granularity | Purpose | Cached at fit time? |
|---|---|---|---|
| `cta_endpoint_summary()` | endpoint | Endpoint structure: `node_id`, path, majority-class prediction, `n_obs`, denominator, and reference sizes | No |
| `cta_endpoint_counts()` | endpoint × actual class | Raw and weighted class counts per leaf endpoint | No |
| `cta_staging_table()` | endpoint | Target-class proportion, propensity (odds-style), and adjusted values for perfectly predicted endpoints | No |
| `cta_propensity_weights()` | endpoint × actual class | Stabilized propensity weights (`n_s × Pr(Z=z) / n_{s,z}`) computed from raw counts; adjusted variant for perfectly predicted endpoints; `undefined_empirical` flag when class absent | No |
| `cta_assign_endpoints()` | observation | Traverses the fitted tree for each row of `newdata`; returns `endpoint_node_id` and sequential `endpoint_id`; supports `missing_action = "na"` (canonical) or `"majority"` | No |
| `cta_observation_weights()` | observation | Joins `cta_assign_endpoints()` output with `cta_propensity_weights()` on `(endpoint_id, actual_class)`; returns endpoint-level weight assigned to each observation | No |
| `cta_confusion_table()` | actual × predicted class | Confusion matrix with per-class sensitivity, specificity, PAC, and ESS from tree predictions | No |
| `cta_family_table()` | family member / MINDENOM | MINDENOM family summary: ESS, D-statistic, strata count, min-terminal-denom, and minimum-D selection flag | No |

---

## Perfect endpoint handling

`cta_staging_table()` and `cta_propensity_weights()` expose empirical values
directly from the stored leaf counts.  When an endpoint contains observations
from only one class, the empirical propensity for the absent class is undefined
(zero denominator).  Both functions flag these cases with a
`perfectly_predicted_endpoint` column (logical) and an `undefined_empirical`
column per endpoint × class row.

For perfectly predicted endpoints, adjusted values are available via the
Yarnold-Linden (2017) remedy: one hypothetical misclassified observation is
added to the absent-class profile.  The addition is applied uniformly within
each endpoint (all class rows for that endpoint share the same adjusted
denominator increment) to preserve the internal consistency of the weight
formula.

The adjustment is a reporting and translation artifact only.  It does not
alter the fitted tree, the stored leaf counts, or any prediction.  Callers
control whether adjusted or unadjusted values are used via the `adjusted`
argument.

---

## MDSA / family context

`cta_descendant_family()` fits a sequence of `oda_cta_fit()` calls with
increasing `mindenom` values, stepping each time by `next_mindenom =
min_terminal_denom + 1`.  The result is a `cta_family` object whose members
span from the most complex admissible tree to the no-tree terminal case.

`cta_family_table()` summarizes the family.  Each member's D-statistic
(ESS / maximum possible ESS, or WESS variant for weighted fits) is computed
alongside strata count and parsimony metrics.  The member with the minimum
D-statistic identifies the MINDENOM value that balances translational strength
against parsimony — the MDSA (Minimum Descriptive Sample Adequacy) selection
criterion.

The myeloma fixture illustrates the full family:

| MINDENOM | Status | Classified n | WESS |
|---|---|---|---|
| 1 | 2-level tree (V14 → V15) | 255 | 27.69% |
| 30 | Stump (V17 only) | 186 | 16.51% |
| 56 | No tree | — | — |

---

## Worked example

`docs/myeloma-cta-translation.md` walks through the full MINDENOM=1 and
MINDENOM=30 trees using actual computed values from the R package.  It covers
endpoint structure, staging table construction, propensity weight arithmetic,
confusion table recovery, and the lean-fit invariant check.

---

## Current limitations

- **No balance diagnostics.** Propensity weights are computed but there is no
  function yet to assess whether they improve covariate balance across
  endpoint-defined groups (standardized mean differences, variance ratios, etc.).
- **Raw-count propensity formula only.** `cta_propensity_weights()` uses raw
  class counts in the weight formula.  Weighted counts are stored in the tree
  (`class_counts_weighted`) but are not yet surfaced in the propensity
  calculation.
- **No downstream outcome model integration.** The observation-level weights
  from `cta_observation_weights()` are ready to pass to a weighted outcome
  model, but no wrapper or example for that step exists yet.
- **No package vignette.** `docs/myeloma-cta-translation.md` serves as the
  worked example; a formal `vignettes/` entry is not yet written.

---

## Next planned direction

The immediate next work is documentation synchronization and applied
translation hardening:

- Keep `docs/myeloma-cta-translation.md` aligned with the completed stack,
  including `cta_assign_endpoints()` and `cta_observation_weights()`.
- Add a small README pointer to this page in a separate README-only slice.
- Preserve the lean-fit invariant while extending examples.
- Defer balance diagnostics, downstream outcome modeling, and weighted
  propensity-weight variants until explicitly designed and canon-reviewed.

Balance diagnostics may be useful future utilities, but they are not part of
the current canon-backed CTA/MDSA propensity-weight translation path.
