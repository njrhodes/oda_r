# CTA Reporting and Translation Stack

This file maps the on-demand reporting and translation functions available for
a fitted `cta_tree` object.  It is a navigation aid, not a behavior
specification.  For canonical engine behavior see `docs/CTA_CANON.md`; for a
fully worked numeric example see `docs/examples/myeloma-cta-translation.md`.

---

## Lean-fit principle

`cta_fit()` (and the internal `oda_cta_fit()` engine it delegates to) stores
only what is needed to reproduce predictions and compute reporting artifacts
on demand:

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
cta_fit(X, y, w, ...)              # public entry point
  |  [delegates to oda_cta_fit()]   # internal engine
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
  |     |
  |     +-- cta_family_table()     MINDENOM family summary with D-statistic
  |
  +-- [output layer — derived, no refit, no training X/y]
        |
        +-- cta_plot_data()        tree diagram data (nodes/edges/endpoints)
        +-- plot.cta_tree()        native base-R tree diagram
```

Most reporting functions operate on stored tree/family summaries and leaf
counts.  `cta_assign_endpoints()` and `cta_observation_weights()` additionally
require caller-supplied `newdata`; `cta_confusion_table()` reads the stored
final-tree training confusion.

---

## Function map

| Function | Row granularity | Purpose | Cached at fit time? |
|---|---|---|---|
| `cta_endpoint_summary()` | endpoint | Endpoint structure: `node_id`, path, majority-class prediction, `n_obs`, denominator, and reference sizes | No |
| `cta_endpoint_counts()` | endpoint × actual class | Raw and weighted class counts per leaf endpoint | No |
| `cta_staging_table()` | endpoint | Target-class count, denominator, proportion, odds, perfect-endpoint flags, and adjusted reporting values | No |
| `cta_propensity_weights()` | endpoint × actual class | Stabilized propensity weights (`n_s × Pr(Z=z) / n_{s,z}`) computed from raw counts; adjusted variant for perfectly predicted endpoints; `undefined_empirical` flag when class absent | No |
| `cta_assign_endpoints()` | observation | Traverses the fitted tree for each row of `newdata`; returns `endpoint_node_id` and sequential `endpoint_id`; supports `missing_action = "na"` (canonical) or `"majority"` | No |
| `cta_observation_weights()` | observation | Joins `cta_assign_endpoints()` output with `cta_propensity_weights()` on `(endpoint_id, actual_class)`; returns endpoint-level weight assigned to each observation | No |
| `cta_confusion_table()` | actual × predicted class | Confusion matrix with per-class sensitivity, specificity, PAC, and ESS from tree predictions | No |
| `cta_family_table()` | family member / MINDENOM | MINDENOM family summary: ESS, D-statistic, strata count, min-terminal-denom, and minimum-D selection flag | No |
| `cta_plot_data()` | node / edge / endpoint | Derived tree diagram data: node geometry, edge labels, endpoint enrichment (target proportion, rank, color, label) when `target_class` is supplied | No |
| `plot.cta_tree()` | — (side effect: plot) | Native base-R tree diagram; structural or target-enriched mode; returns `invisible(pd)` | No |

---

## Perfect endpoint handling

`cta_staging_table()` exposes empirical values directly from the stored leaf
counts and flags perfectly predicted endpoints with a
`perfectly_predicted_endpoint` column (logical).

`cta_propensity_weights()` additionally marks each endpoint × class row with
an `undefined_empirical` flag (logical) when `class_n == 0` at that cell —
i.e., the absent class at a perfectly predicted endpoint.  The empirical
propensity weight is undefined (infinite) for those cells.

For perfectly predicted endpoints, adjusted values are available in
`cta_propensity_weights()` via the Yarnold-Linden (2017) remedy: one
hypothetical misclassified observation is added to the absent-class profile.
The adjusted denominator increment is applied uniformly within each endpoint
(all class rows for that endpoint receive the same increment) to preserve the
internal consistency of the weight formula.

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

`cta_family_table()` summarizes the family.  Each member's D-statistic is
reported as implemented by the family machinery alongside strata count and
parsimony metrics.  The member with the minimum D-statistic identifies the
MINDENOM value that balances translational strength against parsimony — the
MDSA (Minimum Descriptive Sample Adequacy) selection criterion.

The myeloma fixture illustrates the full family:

| MINDENOM | Status | Classified n | WESS |
|---|---|---|---|
| 1 | 2-level tree (V14 → V15) | 255 | 27.69% |
| 30 | Stump (V17 only) | 186 | 16.51% |
| 56 | No tree | — | — |

---

## Worked example

`docs/examples/myeloma-cta-translation.md` walks through the full MINDENOM=1 and
MINDENOM=30 trees using actual computed values from the R package.  It covers
endpoint structure, staging table construction, propensity weight arithmetic,
confusion table recovery, and the lean-fit invariant check.

---

## Current limitations

- **No balance diagnostics.** Propensity weights are computed but no function
  yet assesses covariate balance after weight application. Future
  balance-diagnostic work should be canon-reviewed against the Linden/Yarnold
  JEP covariate-balance paper in `docs/theory/`. That paper frames balance
  assessment as an ODA classification problem: after matching or weighting,
  observed pre-intervention covariates should have weak or non-significant
  ability to discriminate assignment/study groups. Conventional
  standardized-difference summaries may be useful complements, but they are not
  the canon center.
- **Raw-count propensity formula only.** `cta_propensity_weights()` uses raw
  class counts in the weight formula.  Weighted counts are stored in the tree
  (`class_counts_weighted`) but are not yet surfaced in the propensity
  calculation.
- **No downstream outcome model integration.** The observation-level weights
  from `cta_observation_weights()` are ready to pass to a weighted outcome
  model, but no wrapper or example for that step exists yet.
- **No formal package vignette yet.** `docs/examples/myeloma-cta-translation.md` serves
  as the worked example. Vignette conversion is planned for Slice 3; each
  vignette must be CRAN-safe by default (no slow MC/LOO/fixture fits at
  build time).

---

## Next planned direction

The first production version of the translation stack (through `cta_plot_data()`
and `plot.cta_tree()`) is complete. Ongoing and planned work:

**Vignette conversion (Slice 3):**
Convert in-repo analysis Rmd files to formal package vignettes. Every vignette
must use public APIs only. Any fit that invokes MC (`mc_iter`), LOO, or real
fixture data must use precomputed output or be guarded with `eval = FALSE` to
remain CRAN-safe by default. No slow canon fits should run during
`devtools::check()` or `R CMD build`.

**Graphics production polish (Phase 2I — remaining):**
- `plot.cta_family()` for MDSA descendant family comparison.
- Optional ggplot/grid renderer, ellipse/circle split nodes, endpoint
  rectangles, arrow-style edges, dynamic sizing, performance caption, legend.
- Optional Mermaid text export for Quarto/GitHub rendering.
  Mermaid is an export format, not the internal R graphics engine.

**Deferred design work (not yet canon-backed):**
- **Directional hypothesis support** — relevant to ODA/CTA rule reporting but
  specification is open. Tracked as issue #6; do not implement until canon
  target is written.
- **Balance diagnostics** — future work based on the ODA/CTA discrimination
  framing: after matching or weighting, pre-intervention covariates should
  have weak or non-significant ability to discriminate assignment/study groups.
  SMD-style summaries may be useful complements but are not the canon center.
  Canon review against the Linden/Yarnold JEP covariate-balance paper in
  `docs/theory/` is required before implementation.
- **Weighted propensity-weight variant** — `cta_propensity_weights()` uses raw
  observation counts exclusively. A `weighted = TRUE` path using `n_weighted`
  and weighted marginals remains deferred design.
- **Downstream outcome models** — propensity weights are ready to pass to a
  weighted outcome model; no wrapper or example for that step is provided yet.

The lean-fit invariant must be preserved throughout all future additions.
