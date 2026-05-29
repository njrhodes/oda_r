# COVARIATE_BALANCE_CONTRACT.md

Design memo for ODA/CTA covariate-balance analysis and Graphics v3A plot-data
contract.

Canon anchor: Linden & Yarnold (2016), "Using machine learning to assess
covariate balance in matching studies."

Status: **Design only.** No implementation in this slice.

---

## 1. Balance question and intended inputs

Covariate balance asks whether observed baseline covariates predict
treatment/exposure/study-group membership.  If no covariate, and no
combination of covariates, can predict group membership above chance with
statistical reliability, the study groups are considered balanced on the
observed covariates under the declared analytic constraints.

**Inputs:**

| Argument | Type | Description |
|---|---|---|
| `X` | `data.frame` | Baseline covariate frame. Columns are candidates for ODA/CTA. |
| `group` | integer vector | Binary group indicator (e.g., 0 = control, 1 = treated). This plays the role of `y` in all ODA/CTA calls. |
| `w` | numeric or `NULL` | Optional case weights (matching weights, IPW, etc.). Passed to all ODA/CTA calls. |
| `alpha` | numeric | Significance threshold. Default `0.05`. Used for multiplicity-adjusted significance flags. |
| `loo` | character/numeric | LOO mode. Default `"off"`. See §10. |
| `...` | — | Additional arguments forwarded to `oda_fit()` or `cta_fit()`. |

**Not inputs:** the outcome variable, treatment effect estimates, or any
post-baseline variables.  Balance analysis is strictly pre-outcome and
pre-adjustment.

---

## 2. Group variable vs. outcome variable

The `group` variable is the **class variable (y)** in every ODA and CTA call.
The covariates in `X` are the **attributes (x)**.

The outcome variable (clinical endpoint, response) is **not modeled** by any
balance function.  Balance analysis precedes outcome analysis; it characterizes
whether the study design produced comparable groups at baseline.

This distinction must be enforced in every doc, example, and error message:
do not pass the outcome as `group`.  Do not pass `group` as a covariate in `X`.

---

## 3. `oda_balance_table()` contract

### Purpose

Runs a univariate `oda_fit()` for each column of `X`, treating `group` as the
class variable.  Returns one row per covariate summarizing ODA-based balance
diagnostics.

### Signature (proposed)

```r
oda_balance_table(
  X,
  group,
  w            = NULL,
  alpha        = 0.05,
  loo          = "off",
  attr_types   = NULL,   # named character vector; NULL = auto-detect per column
  mc_iter      = 25000L,
  mc_seed      = NULL,
  ...
)
```

### Return value

A `data.frame` of class `c("oda_balance_table", "data.frame")`, one row per
covariate, with columns:

| Column | Type | Description |
|---|---|---|
| `attribute` | character | Column name from `X`. |
| `attr_type` | character | Attribute type used: `"ordered"`, `"categorical"`, `"binary"`. |
| `n_total` | integer | Total non-missing observations used in this fit. |
| `n_group_0` | integer | Observations in group 0. |
| `n_group_1` | integer | Observations in group 1. |
| `sensitivity` | numeric | Sensitivity of the ODA rule (proportion of group-1 correctly predicted). |
| `specificity` | numeric | Specificity of the ODA rule (proportion of group-0 correctly predicted). |
| `mean_pac` | numeric | Mean PAC = (sensitivity + specificity) / 2, in [0, 1]. |
| `ess` | numeric | ESS (%) when `w = NULL` or all weights equal. |
| `wess` | numeric | WESS (%) when case weights are active; `NA` otherwise. |
| `ess_display` | numeric | The operative measure: `wess` when weights active, else `ess`. |
| `p_mc` | numeric | Raw permutation p-value from `oda_fit()`. |
| `p_sidak` | numeric | Sidak-corrected p: `1 - (1 - p_mc)^k` where k = number of covariates. |
| `p_bonferroni` | numeric | Bonferroni-corrected p: `min(p_mc * k, 1)`. |
| `significant_raw` | logical | `p_mc < alpha`. |
| `significant_sidak` | logical | `p_sidak < alpha`. |
| `significant_bonferroni` | logical | `p_bonferroni < alpha`. |
| `rule_type` | character | Rule type from the winning ODA fit (e.g., `"ordered_cut"`, `"binary_map"`). |
| `rule_summary` | character | Human-readable rule string (e.g., `"age <= 55"`). |
| `loo_status` | character | LOO status: `"STABLE"`, `"PVALUE"`, `"REJECTED"`, or `NA` when `loo = "off"`. |
| `ess_loo` | numeric | LOO ESS/WESSL; `NA` when `loo = "off"`. |
| `has_weights` | logical | `TRUE` when case weights were active. |
| `fit_ok` | logical | `TRUE` when `oda_fit()` returned a valid (non-failed) result. |

**Multiplicity note:** `k` = number of covariates with a valid `p_mc` (not
`NA`).  Covariates where `fit_ok = FALSE` are excluded from k and flagged.

**Implementation constraint:** `oda_balance_table()` must call `oda_fit()` for
each covariate; it must not reimplement the ODA search or the MC p-value
procedure.

---

## 4. `cta_balance_table()` contract

### Purpose

Fits a single multivariate `cta_fit()` with `group` as the class variable and
all columns of `X` as candidate predictors.  Returns a structured summary of
the CTA balance result.

### Signature (proposed)

```r
cta_balance_table(
  X,
  group,
  w            = NULL,
  mindenom     = 1L,
  alpha        = 0.05,
  loo          = "off",
  mc_iter      = 5000L,
  mc_seed      = NULL,
  ...
)
```

### Return value

A list of class `c("cta_balance_table")` with fields:

| Field | Type | Description |
|---|---|---|
| `status` | character | `"valid_tree"`, `"stump"`, or `"no_tree"`. |
| `balance_interpretation` | character | `"discriminating"` (some combination predicts group) or `"no_discriminating_combinations"` (no_tree). See §11. |
| `root_attribute` | character | Root split variable; `NA` for no_tree. |
| `n_endpoints` | integer | Number of terminal endpoints; `NA` for no_tree. |
| `overall_ess` | numeric | Full-tree ESS (%); `NA` for no_tree. |
| `overall_wess` | numeric | Full-tree WESS (%); `NA` when weights not active or no_tree. |
| `ess_display` | numeric | Operative measure: `overall_wess` when weights active, else `overall_ess`. |
| `d_stat` | numeric | D statistic; `NA` for no_tree. |
| `mindenom` | integer | MINDENOM used. |
| `has_weights` | logical | `TRUE` when case weights were active. |
| `tree` | `cta_tree` | The raw fitted `cta_tree` object (for downstream `cta_endpoint_table()`, `cta_node_table()`, etc.). |
| `endpoint_table` | `data.frame` | `cta_endpoint_table(tree)` output; zero-row for no_tree. |
| `node_table` | `data.frame` | `cta_node_table(tree)` output. |

**Implementation constraint:** `cta_balance_table()` must call `cta_fit()` (or
`oda_cta_fit()`) once; it must not reimplement ENUMERATE or node-growth logic.

---

## 5. Conventional SMD companion table contract

### Purpose

Provides a conventional standardized mean difference (SMD) table as a
companion diagnostic.  SMD is not the odacore objective; it is included for
comparison with non-ODA balance reports and for plots that display both metrics
side-by-side.

### Signature (proposed)

```r
smd_balance_table(
  X,
  group,
  w       = NULL,   # when supplied: weighted means and weighted SMD
  digits  = 4L
)
```

### Return value

A `data.frame` of class `c("smd_balance_table", "data.frame")`, one row per
covariate, with columns:

| Column | Type | Description |
|---|---|---|
| `attribute` | character | Column name from `X`. |
| `mean_0` | numeric | Unweighted mean (or proportion) in group 0. |
| `sd_0` | numeric | Unweighted SD in group 0. |
| `mean_1` | numeric | Unweighted mean (or proportion) in group 1. |
| `sd_1` | numeric | Unweighted SD in group 1. |
| `smd` | numeric | Raw SMD: (mean_1 − mean_0) / pooled SD. |
| `abs_smd` | numeric | `|smd|`. |
| `balanced_020` | logical | `abs_smd < 0.20` (conventional rule of thumb). |
| `balanced_010` | logical | `abs_smd < 0.10` (stricter threshold). |
| `wmean_0` | numeric | Weighted mean in group 0; `NA` when `w = NULL`. |
| `wmean_1` | numeric | Weighted mean in group 1; `NA` when `w = NULL`. |
| `wsmd` | numeric | Weighted SMD (Rubin formula); `NA` when `w = NULL`. |
| `wabs_smd` | numeric | `|wsmd|`; `NA` when `w = NULL`. |
| `wbalanced_020` | logical | `wabs_smd < 0.20`; `NA` when `w = NULL`. |
| `wbalanced_010` | logical | `wabs_smd < 0.10`; `NA` when `w = NULL`. |

**SMD formula:** pooled SD = sqrt((var_0 * (n_0 - 1) + var_1 * (n_1 - 1)) /
(n_0 + n_1 - 2)).  For binary variables: pooled SD from proportions.
Weighted variant: follow Rubin (2001) / Austin (2009) formulas.

**Scope limit:** `smd_balance_table()` does not test significance, fit ODA, or
produce p-values.  It is a purely descriptive companion table.

---

## 6. `oda_balance_plot_data()` contract

### Purpose

Transforms an `oda_balance_table()` result (and optionally an
`smd_balance_table()` result) into a renderer-independent data structure
suitable for Graphics v3 plotting.  **Does not fit ODA.**

### Signature (proposed)

```r
oda_balance_plot_data(
  balance_table,          # required: oda_balance_table output
  smd_table    = NULL,    # optional: smd_balance_table output to join
  sort_by      = c("ess_display", "abs_smd", "p_mc", "name"),
  p_col        = c("p_mc", "p_sidak", "p_bonferroni"),
  alpha        = 0.05,
  digits       = 2L
)
```

### Return value

A list of class `c("oda_balance_plot_data")` with elements:

| Element | Type | Description |
|---|---|---|
| `rows` | `data.frame` | One row per covariate; see below. |
| `has_weights` | logical | Whether weights were active. |
| `ess_label` | character | `"WESS"` or `"ESS"` depending on weights. |
| `p_col_used` | character | Which p column was selected for plotting. |
| `alpha` | numeric | Significance threshold used. |
| `n_covariates` | integer | Total number of covariates. |
| `n_significant` | integer | Number significant on `p_col_used`. |
| `sort_by` | character | Sort method applied. |

`rows` columns:

| Column | Type | Description |
|---|---|---|
| `attribute` | character | Covariate name. |
| `attr_type` | character | Attribute type. |
| `ess_display` | numeric | ESS or WESS (operative measure). |
| `ess_display_bar` | numeric | `ess_display` clipped to [0, 100] for bar/dot display. |
| `p_plot` | numeric | The selected p column (`p_mc`, `p_sidak`, or `p_bonferroni`). |
| `significant` | logical | `p_plot < alpha`. |
| `significance_label` | character | `"*"` if significant, `""` otherwise. |
| `rule_summary` | character | Human-readable ODA rule. |
| `abs_smd` | numeric | From joined SMD table; `NA` if not supplied. |
| `wsmd_available` | logical | Whether weighted SMD is available. |
| `abs_smd_display` | numeric | `wabs_smd` if weights active and available, else `abs_smd`; `NA` if SMD table not supplied. |
| `rank` | integer | Display rank (sort order, 1 = top). |
| `fit_ok` | logical | Whether the ODA fit was valid. |

**Key invariant:** `oda_balance_plot_data()` reads from `balance_table` and
optionally joins `smd_table` by `attribute`.  It computes no new statistics.

---

## 7. `cta_balance_plot_data()` contract

### Purpose

Transforms a `cta_balance_table()` result into a renderer-independent data
structure for tree-diagram display with balance-specific annotations.  **Does
not fit CTA.**

### Signature (proposed)

```r
cta_balance_plot_data(
  cta_balance,            # required: cta_balance_table output
  target_class = 1L,      # group-1 (treated) as target for endpoint coloring
  digits       = 1L
)
```

### Return value

A list of class `c("cta_balance_plot_data")` with elements:

| Element | Type | Description |
|---|---|---|
| `status` | character | `"valid_tree"`, `"stump"`, or `"no_tree"`. |
| `balance_interpretation` | character | `"discriminating"` or `"no_discriminating_combinations"`. |
| `no_tree_message` | character | Human-readable no_tree annotation for renderers (see §11). `NA` when status ≠ `"no_tree"`. |
| `cta_pd` | list or `NULL` | `cta_plot_data()` output for the embedded tree; `NULL` for no_tree. |
| `ess_display` | numeric | Full-tree ESS or WESS (%); `NA` for no_tree. |
| `d_stat` | numeric | D statistic; `NA` for no_tree. |
| `has_weights` | logical | Whether case weights were active. |
| `ess_label` | character | `"WESS"` or `"ESS"`. |

**For valid tree / stump:** `cta_pd` is produced by calling `cta_plot_data()`
on the embedded `cta_balance$tree` with `target_class = target_class`.  All
existing `cta_plot_data()` semantics apply unchanged.

**For no_tree:** `cta_pd = NULL`; `no_tree_message` is populated (see §11).

**Key invariant:** `cta_balance_plot_data()` reads from `cta_balance_table()`
output; it does not refit CTA.

---

## 8. Multiplicity handling

Multiplicity arises from testing k covariates simultaneously in
`oda_balance_table()`.  Two corrections are always computed and returned:

### Sidak correction

```
p_sidak = 1 - (1 - p_mc)^k
```

Assumes independence across tests.  Conservative in the presence of correlated
covariates.

### Bonferroni correction

```
p_bonferroni = min(p_mc * k, 1)
```

Always valid; more conservative than Sidak when tests are positively
correlated.

### Implementation rules

- `k` counts only covariates with a valid (non-NA) `p_mc` from `oda_fit()`.
- Both corrections are always computed; the caller selects which column to
  use for significance decisions.
- `oda_balance_plot_data()` accepts a `p_col` argument to select which
  corrected p drives the `significant` flag.
- `alpha` is a display/flagging threshold; the raw and corrected p-values
  are always returned for user inspection regardless of alpha.
- No hierarchical testing, no FDR (Benjamini-Hochberg), and no step-down
  procedures are implemented in this version.  These are deferred.

---

## 9. Weights / WESS handling

Case weights affect both the ODA/CTA fitting and the SMD companion:

**ODA/CTA fitting:**
- When `w` is supplied, all `oda_fit()` and `cta_fit()` calls receive `w`.
- The resulting `ess_display` column uses `wess` when weights are active,
  `ess` otherwise.
- The `ess_label` metadata field (`"ESS"` vs `"WESS"`) must propagate into
  every downstream table and plot-data object.
- Never call weighted scores ESS; never call unit-weight scores WESS.

**SMD companion:**
- `smd_balance_table()` computes both unweighted and weighted SMD when `w`
  is supplied.
- The `abs_smd_display` column in `oda_balance_plot_data()` uses the weighted
  SMD when weights are active.

**LOO:**
- LOO uses the same weights as the fitting call.  LOO ESS/WESSL is returned
  in `ess_loo` following existing `oda_fit()` conventions.

---

## 10. LOO / generalizability as optional, not default

Balance is a property of the training sample under the declared constraints.
LOO cross-validation tests whether the ODA rule generalizes beyond that
sample — a secondary question for balance analysis.

**Default:** `loo = "off"` in all balance table functions.

**Available:** `loo = "stable"`, `"pvalue"`, or numeric threshold — forwarded
to every underlying `oda_fit()` call.  LOO results are stored in `loo_status`
and `ess_loo` columns.

**Caution:** LOO balance analysis is expensive (each fold refits the ODA model
on n−1 cases per covariate × fold) and is not needed for a standard balance
report.  Users who enable LOO must set `mc_iter` carefully to avoid excessive
runtimes.

**CTA balance:** LOO at the CTA level (LOO STABLE gate per node) follows the
existing `cta_fit()` `loo` parameter semantics unchanged.  Default remains
`"off"` for `cta_balance_table()`.

---

## 11. no_tree as favorable evidence of multivariable balance

When `cta_balance_table()` returns `status = "no_tree"`, the correct
interpretation is:

> No combination of baseline covariates in `X` was able to predict group
> membership at the declared significance level, LOO constraint, and minimum
> endpoint denominator.  This is **favorable evidence of multivariable
> covariate balance** under the declared analytic constraints.

This must be distinguished from model failure.  In outcome modeling, `no_tree`
means no discriminating structure was found, which may be a problem.  In
balance analysis, `no_tree` means the groups cannot be discriminated, which is
the goal.

**Required output field:** `balance_interpretation`:
- `"no_discriminating_combinations"` — when `no_tree`
- `"discriminating"` — when a valid tree or stump was found

**`no_tree_message` in `cta_balance_plot_data()`:**

```
"No combination of covariates predicted group membership
 (MINDENOM = {mindenom}, alpha = {alpha}).
 This is favorable evidence of multivariable balance
 under the declared constraints."
```

The renderer must display this message prominently when `status = "no_tree"`.
It must not display an empty or error state.

**Univariate note:** A `no_tree` result from `oda_balance_table()` at the
individual covariate level (i.e., `fit_ok = FALSE`) is different — it means
the individual covariate had no valid ODA solution (e.g., constant covariate,
all missing), not that balance was achieved.  These cases are flagged via
`fit_ok = FALSE`, not via `balance_interpretation`.

---

## 12. What is deferred

The following are explicitly out of scope for this contract and the initial
implementation:

| Deferred item | Reason / when |
|---|---|
| ggplot/grid renderer (`plot.oda_balance`, `plot.cta_balance`) | Graphics v3B — after plot-data contract is validated |
| Quarto/Mermaid export | Graphics v3C — after renderer exists |
| Causal-effect estimation (ATE, ATT, DR estimators) | Different problem; not odacore scope |
| Outcome modeling | Strictly post-balance; separate API |
| SORT/GORT recursive balance workflows | SORT/GORT are reserved; not implemented |
| Multiclass group variable (C > 2) | CTA balance currently binary-only; future extension |
| FDR / step-down multiplicity corrections | Deferred; Sidak and Bonferroni cover the initial need |
| `oda_balance_table()` LOO runtime optimization | Design only; no parallel execution or caching |
| Covariate adjustment layers | `w` covers weight-based adjustment; other methods deferred |
| Balance on matched subsets (caliper matching, NN matching) | Caller provides matched `X`, `group`, `w`; matching is out of scope |

---

## Implementation order (when slice opens)

```
1. smd_balance_table()                  ← pure arithmetic; no ODA/CTA calls
2. oda_balance_table()                  ← loops over oda_fit(); multiplicity columns
3. cta_balance_table()                  ← single cta_fit(); no_tree interpretation
4. oda_balance_plot_data()              ← reads oda_balance_table() + optional SMD
5. cta_balance_plot_data()              ← reads cta_balance_table() + cta_plot_data()
6. Tests (synthetic fixture + myeloma)  ← balance tables on known fixtures
7. Graphics v3B renderer                ← deferred; separate slice
```

Dependencies: steps 4–5 require steps 1–3 to exist.  Renderer (step 7) requires
steps 4–5 to be stable.

---

## Reference

Linden A, Yarnold PR (2016). Using machine learning to assess covariate
balance in matching studies.
*Journal of Evaluation in Clinical Practice*, 22(6), 861–867.
