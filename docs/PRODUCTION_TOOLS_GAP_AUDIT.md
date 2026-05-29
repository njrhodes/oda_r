# PRODUCTION_TOOLS_GAP_AUDIT.md

Production readiness and propensity tooling gap audit for odacore.

**Status:** Slice P — design / audit only. No code changed.
**Date:** 2026-05-29
**Slice:** P

---

## 1. Data Readiness / Preflight Functions

### Current coverage

None of the following exist as exported functions:

| Function | Status |
|----------|--------|
| `oda_readiness_check()` | MISSING |
| `oda_clean_missing_codes()` | MISSING |
| `oda_validate_group()` | MISSING |
| `oda_validate_weights()` | MISSING |
| `oda_infer_attr_types()` | MISSING |
| `oda_analysis_bundle()` | MISSING |

Internal equivalents exist inside fitting functions but are not exposed:

| Internal helper | Location | What it does |
|----------------|----------|--------------|
| `.validate_case_weights()` | `R/utils.R` | Validates weight length, finite, non-zero, non-negative |
| `.oda_clean_xy()` | `R/oda_fit.R` | Internal x/y cleaning with miss_codes |
| `.sda_validate_candidate_frame()` | `R/sda_core.R` | SDA candidate frame validation |
| `.sda_resolve_attr_types()` | `R/sda_core.R` | Attribute type inference |

### Gap assessment

There is no standalone preflight pipeline. A production operator cannot validate their
data, check missing codes, or confirm attribute types without running a fitting function.
This is a material readiness gap.

---

## 2. Propensity Weight Functions — Full Audit

### 2.1 cta_propensity_weights() — EXISTS and matches canon

**Exported:** yes (`NAMESPACE` line 43)
**Location:** `R/cta_s3.R`
**Signature:** `cta_propensity_weights(tree, target_class = NULL, adjusted = TRUE)`

**Yarnold/Linden 2017 formula — VERIFIED:**

| Component | Implementation |
|-----------|---------------|
| Empirical weight | `w = n_s * Pr(Z=z) / n_{z,s}` (endpoint n, marginal Pr, stratum count) |
| Undefined endpoint (`n_{z,s} = 0`) | `undefined_empirical = TRUE`; `Inf` before adjustment |
| One-hypothetical-misclassification adjustment | Adds 1 obs to absent class; recomputes `Pr^{adj}(Z=z)` and `n_{z,s}^{adj} = 1` |
| Adjusted weight reported separately | `adjusted_propensity_weight` column |
| Flags per row | `undefined_empirical`, `adjusted`, `perfectly_predicted_endpoint` |
| Empty tree handling | Returns empty 21-column data frame for no_tree |

**Return schema:** 21-column data frame — one row per observation, including:
`propensity_weight`, `undefined_empirical`, `adjusted`, `adjusted_propensity_weight`,
`adjusted_marginal_class_probability`, and full endpoint/stratum counts.

**Assessment:** Fully matches Yarnold/Linden 2017. No hardening needed.

### 2.2 lort_propensity_weights() — MISSING

**Status:** Not implemented.

**Terminal strata availability from cta_ort:**

| Accessor | What it surfaces |
|---------|-----------------|
| `cta_ort_node_table(ort)` | `ort_node_id`, `terminal`, `n`, `selected_endpoint_count`, `method`, `global_optimization`, `sda_anchored` |
| `ort_plot_data(ort)$strata` | `stratum_id`, `node_id`, `prop_class1` and endpoint-level class distributions |

Terminal strata metadata exists. A `lort_propensity_weights()` function is feasible:
it would traverse `ort$ort_settings` for method/metadata and apply the same endpoint-weight
logic as `cta_propensity_weights()` at the ORT terminal level.

**Required labels for any LORT propensity function:**
- `model_family = "lort"`
- `method = "lort"`
- `global_optimization = FALSE`
- `sda_anchored = FALSE`

### 2.3 oda_propensity_weights() — MISSING

**Status:** Not implemented.

ODA fits a single binary rule, producing two strata (left and right of cutpoint).
These rule strata can support propensity weights when the class variable is
group/treatment/exposure (not outcome). Weight formula:
`w = Pr(Z=z) / Pr(Z=z | rule_stratum)`.

This is simpler than CTA (two strata only) and does not require the
Yarnold/Linden one-hypothetical adjustment unless one stratum has zero counts for
a group, which should not occur with a valid ODA fit (both classes represented in
predictions by degeneracy gate).

### 2.4 sda_propensity_weights() — ABSENT (CORRECT)

SDA does not produce a full-sample propensity stratification. It produces
stage order and selected attributes for future SORT. SDA anchor task hook
explicitly prohibits `"propensity_weighting"` in `prohibited_downstream`.

This function must not be implemented.

### 2.5 Generic propensity functions — ABSENT (CORRECT)

The following do not exist and must not be created:

- `propensity_weights()`
- `propensity_scores()`
- `propensity_strata()`
- `propensity_balance_report()`

Per OBBP doctrine: use explicit odacore/model-family naming only.

### 2.6 SORT / GORT propensity functions — FUTURE ONLY

`sort_propensity_weights()` and `gort_propensity_weights()` are future capacity,
dependent on SORT/GORT implementation. Do not stub or implement now.

---

## 3. Agent-Readable Metadata Gaps

### Current task_hook coverage

| Object type | task_hook | Structured metadata | Node/strata accessible |
|-------------|-----------|--------------------|-----------------------|
| `oda_fit` | NO | Partial (rule, ess, pac, p_mc) | No |
| `cta_tree` | NO | Yes (full node tree) | Yes (cta_endpoint_table, cta_node_table) |
| `cta_ort` (LORT) | NO | Yes (ort_settings + node tree) | Yes (cta_ort_node_table, ort_plot_data) |
| `sda_fit` | NO | Yes (steps, settings, diagnostics) | N/A |
| **`sda_anchor`** | **YES** | Yes | N/A (anchor only) |
| `oda_balance_table` | NO | Yes (rows + meta) | No |
| `cta_balance_table` | NO | Yes (tree + tables) | Yes (via embedded tree) |

`sda_anchor` is the only object with structured task_hook metadata for agent/pipeline review.
Adding task hooks to ODA/CTA/LORT objects is possible but is not required for public go.
The OBBP Slice O scope required this audit; implementation of hooks on fitting objects
is optional and should be explicitly approved before Slice Q.

---

## 4. Proposed Slice Q Scope

Subject to review and explicit approval before starting:

### Required (blocking public go)

| Function | Why required |
|----------|-------------|
| `oda_readiness_check()` | No preflight exists; operators have no way to validate data before fitting |
| `oda_clean_missing_codes()` | Missing-code normalization embedded inside fitting only |
| `oda_validate_group()` | No standalone group/class validation |
| `oda_validate_weights()` | Weight validation is internal only via `.validate_case_weights()` |

### Likely required

| Function | Why |
|----------|-----|
| `oda_infer_attr_types()` | Attribute type inference is internal only; useful preflight |
| `oda_propensity_weights()` | ODA produces rule strata; propensity weighting is applicable and straightforward |
| `lort_propensity_weights()` | Terminal strata exposed; needed for LORT balance workflows |

### Not required (cta_propensity_weights already matches canon)

`cta_propensity_weights()` matches Yarnold/Linden 2017 including the one-hypothetical
adjustment. No hardening needed. Only regression-test confirmation required in Slice Q.

### Hard no (per OBBP)

- `sda_propensity_weights()` — SDA is not a propensity estimator; prohibited
- Generic `propensity_weights()` / `propensity_scores()` — use model-family names
- SORT/GORT propensity — future only; not implemented
- Agent framework — not in scope
- New datasets / private data
- Fraud / credit-card work

---

## 5. Confirmation

No code changed in this audit. Git status at audit completion: clean working tree.

---

## 6. Files Audited

- `NAMESPACE`
- `R/utils.R`, `R/oda_fit.R`, `R/oda_s3.R`
- `R/cta_core.R`, `R/cta_s3.R`
- `R/cta_ort.R`
- `R/sda_core.R`, `R/sda_interop.R`, `R/sda_anchor.R`
- `R/balance.R`
- `tests/testthat/test-cta-propensity-weights.R`
