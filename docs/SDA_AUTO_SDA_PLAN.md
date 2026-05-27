# SDA and Auto-SDA Planning Document

**Status:** Design / planning slice (SDA-0). Do not code from this document without
explicit approval of the targeted implementation slice.

**Boundary:** Extends `docs/ORT_SELECTION_METHODS.md`. Does not supersede it.

---

## Table of Contents

1. [Canon Anchor](#1-canon-anchor)
2. [SDA Modes — User-Facing Design](#2-sda-modes--user-facing-design)
3. [Object Contract](#3-object-contract)
4. [Prediction and Application Semantics](#4-prediction-and-application-semantics)
5. [CTA / ORT Interoperability](#5-cta--ort-interoperability)
6. [Auto-SDA](#6-auto-sda)
7. [Generalized Staged Adjustment Terminology](#7-generalized-staged-adjustment-terminology)
8. [Test Plan](#8-test-plan)
9. [Documentation Plan](#9-documentation-plan)
10. [Implementation Slices](#10-implementation-slices)
- [Weighted SDA Design Boundary](#weighted-sda-design-boundary)
11. [Open Design Questions](#11-open-design-questions)

---

## 1. Canon Anchor

### 1.1 What SDA is not

SDA is **not** CTA. SDA is **not** ORT. They share infrastructure
(UniODA, `oda_fit`, `cta_descendant_family`) but are conceptually distinct:

| Method | Role | Unit of search | Stopping criterion |
|--------|------|----------------|--------------------|
| UniODA / ODA | Single-attribute classification | One attribute | N/A — fits and returns |
| CTA / MDSA | Multi-attribute classification tree | All attributes jointly | MINDENOM, LOO gate, pruning |
| SDA | Staged attribute-set identification | One attribute per step | All resolved, Axiom 1, p > α, n < min |
| ORT (recursive CTA) | Imbalanced recursive extension | MDSA at each endpoint | Node guards (max_depth, min_n, etc.) |

SDA is an **upstream staged attribute-set identification and decomposition procedure**.
It identifies which attributes are needed and in what sequence, before or instead of
running MDSA on the full candidate pool.

### 1.2 MPE two-path synthesis

MPE Chapter 12 describes GO-CTA with multiple attributes as obtainable by two paths:

**Path A — MDSA without SDA:**
- Run `cta_descendant_family()` directly over all eligible candidate attributes.
- MDSA identifies the full descendant family.
- No pre-screening step.

**Path B — MDSA with SDA:**
- First run SDA to identify an optimal attribute subset.
- Apply MDSA to the SDA-selected reduced set.
- SDA result is a constrained starting point, not a ceiling on subsequent
  correction/rebalancing passes.

Both paths are valid. SDA-with-MDSA is preferred when the candidate pool is large,
when theoretical parsimony is required, or when temporally-staged confounding
adjustment is the goal.

### 1.3 SDA sequential mechanics

SDA proceeds as follows regardless of mode:

1. Evaluate each candidate attribute individually against the current working sample.
2. Select the best eligible attribute/rule/model by the declared mode criterion.
3. Apply the selected rule/model.
4. Remove correctly classified observations from the working sample.
5. Remove the selected attribute from the candidate set.
6. Repeat on the remaining unresolved observations.
7. Stop when a canon stopping condition is met.

**Canon stopping conditions (first met):**

- All observations are correctly classified.
- All observations in either class are correctly classified.
- Axiom 1 / statistical power is violated (no remaining attribute achieves MINDENOM).
- p > α: no remaining attribute achieves significance.
- `max_steps` reached (safety cap; not a canon condition, but a practical guard).
- `n_unresolved < min_n` (safety cap).

### 1.4 Two SDA modes

**Mode 1 — Legacy / pre-novometric SDA (`sda_mode = "unioda_max_ess"`):**

Per-step selection criterion: **maximum ESS** from iterative UniODA.

Mechanics:
- Run `oda_fit()` (or `oda_univariate_core()` directly) on each candidate attribute
  against the current working sample.
- Select the attribute/rule with maximum ESS among eligible candidates.
- Eligible = passes `min_n`, any stated `mindenom` guard, and `p ≤ α`.
- Apply selected rule; remove correctly classified obs; remove selected attribute.

This mode is a **legacy/approximate approximation**. It is acceptable for exploratory
analysis and pre-screening but must never be presented as MPE novometric SDA.

**Mode 2 — MPE novometric SDA (`sda_mode = "novometric_min_d"`):**

Per-step selection criterion: **minimum D** from per-attribute EO-CTA/MDSA.

MPE Chapter 12 mechanics:

*Step S = 1:*
- For each candidate attribute, run EO-CTA/MDSA (`cta_descendant_family()`) on
  that single attribute against the current working sample.
- Require: MINDENOM / Axiom 1 support (min endpoint n ≥ MINDENOM).
- Require: p ≤ α (MC significance gate).
- Among eligible candidates, select the GO-CTA model that achieves **minimum D**.

*Step S ≥ 2:*
- Delete observations correctly classified in step S from the dataset.
- Remove the selected attribute from the candidate set.
- Repeat EO-CTA/MDSA on each remaining attribute using the reduced sample.

*Termination:* First stopping condition above.

*Hand-off:* Attributes identified across all SDA steps form the constrained candidate
set passed to MDSA or to the recursive CTA extension. This is a starting point,
not a ban on subsequent correction passes (see §5).

### 1.5 Mode isolation requirement

Do not silently mix Mode 1 and Mode 2 within a single `sda_fit()` call.
They can share internal helpers (data removal, candidate tracking, step logging),
but the per-step selection rule must be explicitly governed by the declared `mode`.
When reporting or comparing SDA results, always state which mode was used.

---

## 2. SDA Modes — User-Facing Design

### 2.1 Public API

```r
sda_fit(
  X,                                          # data.frame of candidate attributes
  y,                                          # integer class vector
  mode        = c("novometric_min_d",         # default: MPE canon
                  "unioda_max_ess"),           # alternative: legacy/approximate
  attr_types  = NULL,                         # named chr vector, or NULL for auto-detect
  weights     = NULL,                         # case weights (NULL = unit)
  mindenom    = NULL,                         # integer; required for novometric_min_d
  mc_iter     = 5000L,
  mc_seed     = 42L,
  mc_stop     = 99.9,
  mc_stopup   = 20,
  alpha       = 0.05,
  loo         = "stable",
  max_steps   = NULL,                         # safety cap on SDA steps
  min_n       = NULL,                         # safety cap on working sample size
  min_class_n = NULL,                         # stop if either class n < this
  remove_correct = TRUE,                      # TRUE = canonical SDA; FALSE = diagnostic only
  collinearity   = c("skip", "warn", "allow"),
  verbose     = FALSE
)
```

**Notes:**

- `mindenom` is required when `mode = "novometric_min_d"` and no default is
  appropriate. It must be supplied or `sda_fit()` will error with a canon-grounded
  message (MINDENOM must be tied to statistical power; see §2.3).
- `mindenom` is ignored when `mode = "unioda_max_ess"` (UniODA does not run
  descendant families; no MINDENOM gate is applied). A warning is emitted if
  `mindenom` is supplied with `mode = "unioda_max_ess"`.
- `remove_correct = FALSE` is for diagnostic/dry-run use only. It runs the step
  logic without modifying the working sample. Produces a `stop_reason = "dry_run"`.
- `collinearity = "skip"` removes duplicate-column candidates silently. `"warn"` skips
  with a warning. `"allow"` retains duplicates (not recommended).

### 2.2 Separate wrappers: not recommended at this stage

`sda_unioda_fit()` and `sda_novometric_fit()` as separate public functions would
duplicate argument validation and object construction. Prefer `sda_fit(..., mode=)`.
If the two modes require substantially different argument sets in the future, revisit.

### 2.3 MINDENOM and power (novometric mode)

MINDENOM must be tied to statistical power (MPE):

- **Moderate effect:** N ≈ 32 cases per class per category → MINDENOM = 64 per endpoint.
- **Strong effect:** N ≈ 12 cases per class per category → MINDENOM = 24 per endpoint.

Do not default MINDENOM to 1 for novometric mode and claim the result satisfies power
requirements. The caller must declare MINDENOM explicitly for `mode = "novometric_min_d"`.

### 2.4 Internal helpers

Do not expose these publicly at this stage:

```r
.sda_step_unioda_max_ess(X, y, w, candidates, active_rows, settings)
  # Runs oda_fit() on each candidate; returns per-candidate ESS, p_mc, eligible flag.
  # Returns candidate_table + winner index.

.sda_step_novometric_min_d(X, y, w, candidates, active_rows, settings)
  # Runs cta_descendant_family() per candidate (single-attribute MDSA).
  # Returns per-candidate D, ESS, p_mc, mindenom, eligible flag.
  # Returns candidate_table + winner index (min-D among eligible).

.sda_remove_correctly_classified(y, y_pred, active_rows)
  # Returns: correct_indices, incorrect_indices (both global index vectors).

.sda_validate_candidate_frame(X, y, w, settings)
  # Checks dimensions, class count, weight validity, zero-variance columns.
  # Returns list(ok, warnings, candidate_names).
```

### 2.5 Argument guard matrix

| Argument | unioda_max_ess | novometric_min_d |
|----------|---------------|-----------------|
| `mindenom` | Warn + ignore | Required |
| `loo` | Passed to `oda_fit()` | Passed to `cta_descendant_family()` |
| `attr_types` | Passed to `oda_fit()` | Used in `cta_descendant_family()` call per attr |
| `weights` | **Error if non-NULL (SDA-1)** | **Error if non-NULL (SDA-1)** |
| `mc_stop`/`mc_stopup` | Passed to `oda_fit()` | Passed to `cta_descendant_family()` |
| `collinearity` | Applied before any step | Applied before any step |

`weights` is in the API signature to reserve the argument name, but supplying
a non-NULL value errors in SDA-1 with: `"Weighted SDA is not implemented yet;
use weights = NULL. Weighted SDA requires explicit weighted removal, WESS
labeling, and weighted candidate-table semantics."` See §Weighted SDA design boundary.

---

## 3. Object Contract

### 3.1 Returned class

```r
class(result) == c("sda_fit", "odacore_sda")
```

### 3.2 Top-level fields

```r
list(
  call                   = match.call(),
  mode                   = "novometric_min_d" | "unioda_max_ess",
  outcome_name           = character(1),  # or NA if y has no name
  class_levels           = integer(2),    # sorted unique class labels
  n_initial              = integer(1),    # total obs entering SDA
  n_final_unresolved     = integer(1),    # obs remaining at termination
  candidate_names_initial = character(),  # names of all candidates before step 1
  selected_attributes    = character(),   # names of attributes selected across all steps, in order
  steps                  = list(),        # list of step objects (see §3.3)
  unresolved_indices     = integer(),     # global row indices not classified at any step
  resolved_indices       = integer(),     # global row indices classified at some step
  settings               = list(),        # all sda_fit() arguments as supplied
  stop_reason            = character(1),  # see §3.4
  diagnostics            = list()         # warnings, collinearity flags, step counts
)
```

### 3.3 Step object

Each element of `steps` is a named list:

```r
list(
  step_id                = integer(1),
  attribute              = character(1),    # selected attribute name
  mode                   = character(1),    # mode at this step (always same as top-level)
  n_in                   = integer(1),      # working-sample size entering this step
  class_counts_in        = integer(2),      # class counts before removal
  model                  = list(),          # oda_fit or cta_descendant_family result for selected attr
  rule_summary           = list(),          # compact rule description (type, cut, direction, etc.)
  ess                    = numeric(1),
  d                      = numeric(1),      # NA for unioda_max_ess mode
  p_mc                   = numeric(1),
  ge_count               = integer(1),      # MC iterations achieving >= observed ESS
  iter_used              = integer(1),      # actual MC iterations run
  mindenom               = integer(1),      # NA for unioda_max_ess mode
  n_correct              = integer(1),
  n_incorrect            = integer(1),
  correct_indices_global  = integer(),      # row indices removed at this step
  incorrect_indices_global = integer(),     # row indices remaining after this step
  selected               = TRUE,
  reason                 = character(1),    # why this candidate won ("max_ess" | "min_d")
  candidate_table        = data.frame()     # see §3.5
)
```

### 3.4 Stop reason codes

| Code | Meaning |
|------|---------|
| `"all_resolved"` | All observations correctly classified |
| `"class_resolved"` | All observations in one class correctly classified |
| `"axiom1_violated"` | No candidate satisfies MINDENOM (novometric mode) |
| `"p_gate"` | No candidate achieves p ≤ α |
| `"max_steps"` | Safety cap on number of steps reached |
| `"min_n"` | Working sample too small to continue |
| `"min_class_n"` | A class has too few observations to continue |
| `"no_candidates"` | Candidate pool exhausted |
| `"dry_run"` | `remove_correct = FALSE`; no actual removal performed |

### 3.5 Candidate table (per step)

One row per candidate attribute evaluated at that step. This is the primary
auditability record — agents and reviewers must be able to see what was tried,
why something was rejected, and why the winner won.

```r
data.frame(
  attribute       = character(),   # candidate name
  status          = character(),   # "selected" | "eligible" | "ineligible" | "skipped"
  n               = integer(),     # working-sample n for this attribute (after missingness)
  class_counts    = list(),        # raw class counts as integer(2) vector
  ess             = numeric(),
  d               = numeric(),     # NA for unioda_max_ess
  p_mc            = numeric(),
  mindenom        = integer(),     # NA for unioda_max_ess
  min_ep_n        = integer(),     # minimum endpoint n achieved (novometric only)
  eligible        = logical(),
  ineligible_reason = character(), # NA if eligible; otherwise code from §3.6
  selected        = logical()
)
```

### 3.6 Ineligibility reason codes (candidate table)

| Code | Meaning |
|------|---------|
| `"p_gate"` | p > α |
| `"axiom1"` | Min endpoint n < MINDENOM |
| `"min_n"` | Working n below threshold |
| `"pure_node"` | Only one class present in working sample for this attribute |
| `"all_missing"` | All values are missing for this attribute in working sample |
| `"collinear"` | Excluded by collinearity check |
| `"no_model"` | `oda_fit()` / `cta_descendant_family()` returned ok = FALSE |
| `"zero_variance"` | Attribute is constant in working sample |

---

## 4. Prediction and Application Semantics

### 4.1 Function signature

```r
predict.sda_fit(
  object,
  newdata,
  type = c("class", "stage", "rule", "propensity", "weights", "trace"),
  ...
)
```

### 4.2 Sequential application logic

Prediction follows the **learned selected-step sequence** from the fitted object.
It does not re-scan `X`, does not re-evaluate candidates, and does not select a
"first attribute" from the newdata frame. The sequence is `object$steps[[1]]`,
`object$steps[[2]]`, …, in order — each step's selected attribute and rule are
fixed at fit time.

At step `s`:
1. Apply the selected rule/model from `object$steps[[s]]`.
2. If the step's rule classifies the observation (definite class prediction —
   not NA, not unresolved), record `class`, `stage_id = s`,
   `rule_id = object$steps[[s]]$attribute`, and exit the sequence for that
   observation.
3. If the step's rule does not classify the observation (missing attribute
   value, or observation falls outside the rule's covered region), advance
   to step `s+1`.
4. If no step classifies the observation, it is `resolved = FALSE`.

**Training-time vs. prediction-time analogy:**

During training, SDA removes correctly classified observations from the working
sample. For prediction, the analogous behavior is: once a step's selected rule
fires for an observation, that observation exits the sequence. Unresolved
observations continue.

This is a **sequential selected-step application**, not an ensemble vote and
not a first-column scan. The step sequence is the learned SDA structure; it is
fixed and ordered.

### 4.3 Output by type

| `type` | Returns |
|--------|---------|
| `"class"` | Integer vector of predicted class labels; `NA` for unresolved |
| `"stage"` | Integer vector of step_id at which classified; `NA` for unresolved |
| `"rule"` | Character vector of selected attribute name at classifying step; `NA` for unresolved |
| `"propensity"` | Numeric vector of propensity scores (stage-specific); requires step-level propensity support; deferred to SDA-5 |
| `"weights"` | Numeric vector of SDA-derived propensity weights; requires `sda_weights()`; deferred to SDA-5 |
| `"trace"` | Data frame: one row per observation × step; columns: obs_id, step_id, attribute, classified, class_pred |

`type = "trace"` is the agent/debug mode. It enables full replay of which step fired
for each observation and why.

**`type = "weights"` semantics:** This refers to **SDA-derived propensity weights**
computed from SDA strata after the fit. It does not mean the input `weights`
argument passed to `sda_fit()`. Do not overload. Input fitting weights are a
future SDA-5 concern; SDA-derived propensity weights are a distinct output concept.
In SDA-1, `type = "propensity"` and `type = "weights"` both error with a
"not yet implemented" message.

### 4.4 Missing attribute handling

If an observation's classifying-step attribute is missing:
- The rule is treated as non-classifying for that observation.
- Advance to the next step.
- This is path-local missingness — consistent with `predict.cta_tree()` canon.

Do not silently impute or majority-route missing values in prediction.

### 4.5 Wide newdata contract

`predict.sda_fit()` must accept `newdata` with additional columns beyond the
SDA candidate set. Column matching is name-based. Required columns are the
`selected_attributes` from the fitted object. Extra columns are ignored.
Missing required columns cause a clear error.

---

## 5. CTA / ORT Interoperability

### 5.1 Accessor helpers

```r
sda_selected_attributes(fit)
  # Returns character vector of selected attribute names in SDA step order.

sda_candidate_table(fit, step = NULL)
  # Returns candidate table for step (integer) or list of all steps.

sda_step_table(fit)
  # Returns a summary data.frame: one row per SDA step.
  # Columns: step_id, attribute, ess, d, p_mc, mindenom, n_in, n_correct, stop_reason.

as_cta_candidates(fit, X)
  # Returns X[, sda_selected_attributes(fit), drop = FALSE].
  # Subsets X to the SDA-selected attribute columns.

sda_to_cta_data(fit, X, y)
  # Returns list(X_cta = ..., y_cta = ...) where:
  #   X_cta: SDA-selected attributes only.
  #   y_cta: y for all observations (including those SDA resolved and unresolved).
  # Does NOT remove SDA-resolved observations — CTA sees the full sample with
  # a constrained candidate set. The resolved/unresolved distinction is recorded
  # in fit$resolved_indices / fit$unresolved_indices for downstream use.

sda_stage_table(fit)
  # Returns a data.frame of SDA strata:
  # stage_id, attribute, rule_summary, n_classified, class_counts, ess, d.

sda_propensity_table(fit, X, y)
  # Returns observation-level propensity estimates derived from SDA strata.
  # Requires predict.sda_fit(type = "stage") output.

sda_weights(fit, X, y, method = c("inverse", "stabilized", "overlap", "none"))
  # Returns numeric vector of observation-level weights.
  # method = "inverse": inverse propensity weights from SDA strata.
  # method = "stabilized": stabilized IPW.
  # method = "overlap": overlap/trimming weights.
  # method = "none": unit weights (diagnostic).
```

### 5.2 Use case 1: SDA → MDSA / CTA

Feed SDA-selected attributes as constrained candidate set for `cta_descendant_family()`
or `cta_fit()`.

```r
sda_result  <- sda_fit(X, y, mode = "novometric_min_d", mindenom = 24L, ...)
X_cta       <- as_cta_candidates(sda_result, X)
cta_result  <- cta_fit(X_cta, y, ...)
```

This is the Path B workflow from MPE Chapter 12. The MDSA / CTA search is constrained
to the SDA-identified attribute set. The full sample is passed to CTA — SDA resolution
does not restrict which observations CTA sees.

### 5.3 Use case 2: SDA → imbalanced recursive CTA extension (ORT)

SDA identifies the initial staged structure and constrained candidate set.
The ORT extension starts from that anchor, not from an unconstrained greedy
local min-D scan on the full candidate frame.

```r
sda_result  <- sda_fit(X, y, mode = "novometric_min_d", mindenom = 24L, ...)
X_ort       <- as_cta_candidates(sda_result, X)
ort_result  <- cta_fit(X_ort, y, recursive = TRUE, ...)
```

Later-stage or time-varying attributes may be added to `X_ort` when their
timing is valid and their purpose is explicit stratum correction/rebalancing.
They must be declared — do not silently absorb all columns.

### 5.4 Use case 3: SDA → staged propensity workflow

Convert SDA stages to propensity strata, compute observation-level weights,
and carry weights into later weighted CTA/EO-CTA stages.

```r
sda_result <- sda_fit(X_baseline, y_assign, mode = "novometric_min_d", ...)
wts        <- sda_weights(sda_result, X_baseline, y_assign, method = "stabilized")
cta_result <- cta_fit(X_outcome, y_outcome, weights = wts, ...)
```

This requires that `sda_weights()` is implemented before this workflow can be used.
See implementation slice SDA-5.

### 5.5 Use case 4: SDA → agents / auditable workflows

SDA objects must be interrogable by agents that need to know:
- What attributes were tried and in what order.
- Why a candidate was rejected.
- Which observations were resolved at each step.
- What the current step's working sample looks like.
- Whether SDA is complete or was interrupted.

The `candidate_table` per step (§3.5) and `stop_reason` (§3.4) are the
primary auditability fields. `sda_step_table()` and `sda_candidate_table()`
are the accessor interface.

For agent resumability, `settings` stores all `sda_fit()` arguments as supplied.
A workflow agent can reconstruct the call from `fit$call` or `fit$settings`.

**Requirement:** Do not emit hidden global side effects in `sda_fit()`.
No `options()` mutation, no global RNG reset, no implicit cache writes.
MC seed is set locally per step via `mc_seed` argument to `oda_fit()` /
`cta_descendant_family()`.

---

## 6. Auto-SDA

### 6.1 Purpose and boundary

SDA executes a declared candidate set and declared mode.

auto-SDA helps construct and validate the candidate set, timing set, and run plan
before `sda_fit()` is called. It does not silently choose scientifically meaningful
timing or causal roles without user declaration.

auto-SDA never runs `sda_fit()` by default. It proposes; the user approves.

### 6.2 Public function

```r
auto_sda_plan(
  data,
  outcome,
  candidates          = NULL,     # NULL: use all non-outcome columns
  exclude             = NULL,     # character vector of column names to force-exclude
  role_map            = NULL,     # named list: "assignment_mechanism" | "outcome" | "id" | "leakage"
  time_map            = NULL,     # named numeric/integer: column → time index
  stage_map           = NULL,     # named integer: column → stage assignment
  attr_types          = NULL,     # named chr: column → "ordered" | "categorical" | "binary"
  collinearity_threshold = 1.0,   # Kendall τ or exact-match threshold for exclusion
  min_n               = NULL,     # passed through to sda_fit()
  min_class_n         = NULL,
  mode                = c("novometric_min_d", "unioda_max_ess"),
  dry_run             = TRUE      # default: plan only, do not fit
)
```

### 6.3 Responsibilities

1. **Type inference:** Infer basic variable type for each candidate column
   (numeric ordered, integer binary, factor categorical, logical binary).
   Do not override a declared `attr_types` entry.

2. **Exclusion:**
   - Exclude the outcome column.
   - Exclude columns listed in `exclude`.
   - Exclude columns declared as `role_map = "id"` or `"leakage"`.
   - Exclude constants (zero variance).
   - Exclude all-missing columns.
   - Exclude character/non-numeric columns with no sensible coercion.
   - Record each exclusion with a reason code (see §6.5).

3. **Collinearity / duplicate flagging:**
   Flag (not auto-exclude) collinear candidate pairs. Emit as `warnings`.
   If `collinearity_threshold = 1.0` (default), flag only exact duplicates.

4. **Temporal ordering enforcement:**
   If `time_map` is supplied, flag any candidate with `time_map[col] > time_map[outcome]`
   as potential leakage. Emit as `warnings`; do not silently exclude.

5. **Role validation:**
   If `role_map` is supplied, validate declared roles against `time_map` if available.
   Flag `"assignment_mechanism"` columns that are post-exposure by timing as
   potential leakage.

6. **Candidate frame construction:**
   Return the validated candidate data frame in `proposed_call$X`.

7. **Run plan:**
   Construct a `proposed_call` that can be passed to `sda_fit()` via
   `do.call(sda_fit, plan$proposed_call)`.

8. **Reproducibility hash:**
   Store a hash of the input data dimensions, outcome name, candidate set,
   mode, and seed in `$metadata$hash`. This allows an agent to verify
   whether the plan is stale relative to the current data state.

### 6.4 Returned class

```r
class(result) == c("auto_sda_plan", "odacore_plan")
```

Fields:

```r
list(
  outcome           = character(1),
  candidate_names   = character(),     # final validated candidate set
  excluded_names    = character(),     # columns excluded from candidates
  exclusion_reasons = data.frame(),    # one row per excluded column; see §6.5
  role_map          = list(),          # as supplied or inferred
  time_map          = list(),          # as supplied
  attr_types        = character(),     # inferred + declared
  warnings          = character(),     # collinearity, leakage flags, etc.
  proposed_call     = list(),          # ready-to-pass args for sda_fit()
  dry_run           = logical(1),
  metadata          = list(
    n_rows     = integer(1),
    n_cols_in  = integer(1),
    n_cands    = integer(1),
    n_excluded = integer(1),
    hash       = character(1),
    timestamp  = character(1)
  )
)
```

### 6.5 Exclusion reason codes (auto-SDA)

| Code | Meaning |
|------|---------|
| `"is_outcome"` | Column is the declared outcome |
| `"force_excluded"` | Listed in `exclude` |
| `"role_id"` | Declared as `role_map = "id"` |
| `"role_leakage"` | Declared as `role_map = "leakage"` |
| `"zero_variance"` | Constant across all rows |
| `"all_missing"` | All values are NA or miss_code |
| `"post_exposure"` | `time_map` places column after exposure (flagged, not excluded — emit warning) |
| `"invalid_type"` | Non-coercible type (e.g., list column, raw) |
| `"collinear"` | Exact or near-exact duplicate of another candidate |

Note: `"post_exposure"` is emitted as a **warning**, not an exclusion by default.
Temporal validity is a scientific judgment; auto-SDA flags but does not decide.

### 6.6 Agent principle

auto-SDA proposes and validates. It does not silently decide causal validity,
temporal ordering, or role assignment. Scientific meaning must be declared by
the analyst, not inferred by the function.

When `dry_run = TRUE` (default), `auto_sda_plan()` returns the plan object
without calling `sda_fit()`. The caller inspects, modifies, and approves before
passing to `sda_fit()`.

To execute: `do.call(sda_fit, plan$proposed_call)`.

---

## 7. Generalized Staged Adjustment Terminology

SDA and all downstream tools use the following terminology. Do not substitute
medical-specific language (e.g., "treatment predictors") in public-facing code,
documentation, or error messages.

| Term | Definition |
|------|-----------|
| **Exposure / assignment / action** | The process state, intervention, allocation, policy, decision, behavior, or condition whose confounding structure is being adjusted |
| **Assignment propensity** | The probability or staged propensity that an observation receives an exposure, action, or assignment, conditional on eligible attributes at that stage |
| **Assignment-mechanism predictors** | Attributes used to model or rebalance the exposure/assignment process; not the same as the exposure variable itself and not the same as outcome predictors |
| **Outcome predictors** | Attributes used to model the final response; may overlap with assignment-mechanism predictors but role must be declared explicitly |
| **Baseline / pre-exposure stage** | Attributes measured before the exposure, assignment, or action; used to construct initial propensity strata; post-exposure variables must not enter this stage |
| **Correction / rebalancing stage** | Later-stage or time-varying attributes used within existing strata to refine balance or classification, subject to temporal validity |
| **Marginal / sequential confounding control** | Multi-stage adjustment where propensity or balance is updated across stages as new valid time-indexed information becomes available |
| **R→L recursive pass** | Recursive evaluation proceeds right branch then left branch |

**Core staging rule:**

> Start upstream. Rebalance downstream only when declared, temporally valid, and constrained.

SDA operationalizes this rule by identifying upstream attributes first (Step 1),
then proceeding to subsequent candidates on the reduced unresolved sample. Downstream
attributes enter only in explicitly declared correction passes — not through
uncontrolled all-column scans.

---

## 8. Test Plan

Tests must use only synthetic or publicly available data. Private PsA fixtures
remain untracked until a public-safe equivalent exists.

### 8.1 Mode = "unioda_max_ess" tests

| # | Test | Assertion |
|---|------|-----------|
| 1 | Known-winner selection | Selects the attribute with highest ESS at step 1 from a synthetic frame where the winner is unambiguous |
| 2 | Correct removal | After step 1, `fit$steps[[1]]$n_correct` observations are removed from the working sample at step 2 |
| 3 | Candidate removal | Selected attribute does not appear in step 2 candidates |
| 4 | Stop: no candidates | When one-attribute frame is used, stops after step 1 with `stop_reason = "no_candidates"` |
| 5 | Stop: min_n | Stops with `stop_reason = "min_n"` when working sample drops below threshold |
| 6 | Candidate table | Every step stores a `candidate_table` with one row per evaluated attribute |
| 7 | predict sequential | `predict.sda_fit(type="class")` classifies obs at correct step; unresolved obs return `NA` |
| 8 | Wide newdata | `predict()` works on newdata with extra columns; uses name-based matching |
| 9 | Missing column | `predict()` errors clearly when a required split attribute is absent from newdata |
| 10 | Weights accepted | `weights` argument passes through without error; ESS reflects weighting |
| 11 | Dry run | `remove_correct = FALSE` runs without modifying working sample; `stop_reason = "dry_run"` |
| 12 | Verbose mode | `verbose = TRUE` emits `[SDA]` progress lines; no warning/error |

### 8.2 Mode = "novometric_min_d" tests

| # | Test | Assertion |
|---|------|-----------|
| 1 | Per-attribute MDSA | `cta_descendant_family()` is called once per candidate per step (not UniODA) |
| 2 | min-D selection | Selects the candidate with lowest D among eligible |
| 3 | p gate | Candidate with p > α is marked ineligible; if all fail, `stop_reason = "p_gate"` |
| 4 | Axiom 1 gate | Candidate whose min endpoint n < MINDENOM is marked `ineligible_reason = "axiom1"` |
| 5 | Candidate table fields | Candidate table includes `d`, `p_mc`, `mindenom`, `min_ep_n` columns |
| 6 | Correct removal | Correctly classified obs removed after step 1; step 2 sees smaller working sample |
| 7 | Attribute removal | Selected attribute absent from step 2 candidate set |
| 8 | Stop: all resolved | When all obs classified, `stop_reason = "all_resolved"` |
| 9 | Stop: class resolved | When all obs in one class classified, `stop_reason = "class_resolved"` |
| 10 | CTA interop | `as_cta_candidates(fit, X)` returns exactly the SDA-selected columns; `cta_fit()` runs without error on that frame |
| 11 | mindenom required | `sda_fit(mode="novometric_min_d")` with no `mindenom` errors with canon-grounded message |
| 12 | mindenom warn on legacy | `sda_fit(mode="unioda_max_ess", mindenom=24L)` emits warning; runs without error |

### 8.3 Fixture tests

| Fixture | Type | Notes |
|---------|------|-------|
| Synthetic binary, 2 attrs | Public | Step 1 selects clear winner; step 2 selects second attr; stop on no_candidates |
| Synthetic rare-event | Public | Low prevalence; confirm Axiom 1 gate fires before p gate in novometric mode |
| Synthetic tie case | Public | Two attributes with equal ESS or equal D; confirm tie is broken by declared rule (enumeration order) |
| PsA private fixture | Private / untracked | Not used in package tests; referenced only in dev scripts |

---

## 9. Documentation Plan

### 9.1 Files to create (eventually, not in SDA-0)

| File | Slice | Content |
|------|-------|---------|
| `docs/SDA_DESIGN.md` | SDA-1 | Implementation notes, object contract, canon reference |
| `docs/AUTO_SDA_DESIGN.md` | SDA-3 | Auto-SDA design, exclusion codes, role/time map semantics |

### 9.2 Rd documentation (eventually)

| File | Slice |
|------|-------|
| `man/sda_fit.Rd` | SDA-1 |
| `man/predict.sda_fit.Rd` | SDA-1 |
| `man/sda_step_table.Rd` | SDA-2 |
| `man/sda_candidate_table.Rd` | SDA-2 |
| `man/auto_sda_plan.Rd` | SDA-3 |

### 9.3 Vignettes (eventually)

| File | Slice |
|------|-------|
| `vignettes/sda-staging.Rmd` | SDA-4 |
| `vignettes/staged-adjustment-workflow.Rmd` | SDA-5 |

### 9.4 This document

`docs/SDA_AUTO_SDA_PLAN.md` is the SDA-0 design output.
It is not implementation specification. Do not code from it without
explicit approval of the targeted slice.

---

## 10. Implementation Slices

### Slice SDA-0 — Design document (current)

Deliverable: `docs/SDA_AUTO_SDA_PLAN.md`.
No R code. No tests. No NAMESPACE changes.

### Slice SDA-1 — `sda_fit(mode="unioda_max_ess")` — minimal implementation

New files:
- `R/sda_core.R`: `sda_fit()`, `.sda_step_unioda_max_ess()`,
  `.sda_remove_correctly_classified()`, `.sda_validate_candidate_frame()`
- `R/sda_s3.R`: `predict.sda_fit()`, `print.sda_fit()`, `summary.sda_fit()`
- `tests/testthat/test-sda-unioda.R`: Tests from §8.1

Scope constraints:
- Only `mode = "unioda_max_ess"` implemented; `"novometric_min_d"` returns a
  clear "not yet implemented" error.
- `predict.sda_fit()`: only `type = "class"` and `type = "stage"` required in SDA-1.
- No propensity or weight outputs.
- No `auto_sda_plan()`.

NAMESPACE additions:
- `export(sda_fit)`
- `S3method(predict, sda_fit)`
- `S3method(print, sda_fit)`
- `S3method(summary, sda_fit)`

### Slice SDA-2 — CTA interop helpers

New functions (added to `R/sda_s3.R` or `R/sda_interop.R`):
- `sda_selected_attributes()`
- `sda_step_table()`
- `sda_candidate_table()`
- `as_cta_candidates()`
- `sda_to_cta_data()`

Tests:
- `tests/testthat/test-sda-interop.R`: Tests from §8.2 row 10 + end-to-end
  SDA → `cta_fit()` pipeline on synthetic data.

NAMESPACE additions:
- `export()` the five interop helpers.

### Slice SDA-3 — `auto_sda_plan()`

New file: `R/auto_sda.R`

Scope constraints:
- `dry_run = TRUE` by default; no fitting.
- Type inference, exclusion logic, collinearity flagging.
- `role_map`, `time_map`, `stage_map` accepted and stored but temporal-ordering
  enforcement is warning-only.
- Does not call `sda_fit()` unless `dry_run = FALSE` and user passes explicit
  `fit = TRUE` (TBD argument name).

Tests: `tests/testthat/test-auto-sda.R`

### Slice SDA-4 — `sda_fit(mode="novometric_min_d")`

Adds novometric mode to `sda_fit()`.
Requires per-attribute `cta_descendant_family()` calls; relies on existing
`cta_descendant_family()` infrastructure.

Tests: `tests/testthat/test-sda-novometric.R`: Tests from §8.2.

Requires synthetic canon anchors for novometric mode. Private PsA fixture
remains untracked.

### Slice SDA-5 — Staged adjustment workflow helpers

New functions:
- `sda_stage_table()`
- `sda_propensity_table()`
- `sda_weights()`
- `predict.sda_fit(type = "propensity")` and `type = "weights"`

These require careful design around propensity estimation from SDA strata
and stabilized weight computation. A separate design note is recommended
before implementing SDA-5.

---

---

## Weighted SDA Design Boundary

**SDA-1 is unweighted only.** Weighted SDA is deferred because the following
semantics must be specified and agreed before implementation. This section records
what that design must cover when the time comes.

### A. Weighted fitting

When `weights` are supplied, fitting calls (`oda_fit()` / `cta_descendant_family()`)
compute weighted accuracy. Accuracy must then be labeled **WESS**, not ESS.
The candidate table must be extended to carry both raw and weighted fields:

```
n_raw              integer   — observation count
class_counts_raw   integer(2) — raw class counts
weight_total       numeric   — sum of case weights for this attribute's valid obs
class_weights      numeric(2) — summed weights per class
ess                numeric   — ESS (unit-weight runs only; NA when weighted)
wess               numeric   — WESS (weighted runs only; NA when unweighted)
metric_label       character — "ESS" or "WESS" (never both non-NA)
```

Do not use a single `ess` column that silently holds WESS when weights are active.
That violates the package's WESS/ESS labeling convention.

### B. Removal mechanics

Training SDA removes observations by row identity, not by fractional weight units.

If row `i` is correctly classified, row `i` exits the unresolved set entirely.
Its case weight contributes to weighted summary counts but removal is row-based.
MPE describes SDA removal by observation count; no fractional-weight removal
semantics are described.

Record both raw and weighted correct/incorrect counts at each step:

```
n_correct_raw       integer
n_correct_weighted  numeric   — sum of weights for correctly classified obs
n_incorrect_raw     integer
n_incorrect_weighted numeric  — sum of weights for incorrectly classified obs
```

### C. Step object future fields

When weighted SDA is implemented, add to the step object:

```
weighted            = logical(1)
weight_total_in     = numeric(1)   — sum of weights entering this step
class_weights_in    = numeric(2)   — summed weights per class entering this step
n_correct_weighted  = numeric(1)
n_incorrect_weighted = numeric(1)
metric_label        = character(1) — "ESS" | "WESS"
```

### D. Candidate table future fields

When weighted SDA is implemented, extend the candidate table with:

```
n_raw              integer
weight_total       numeric
class_counts_raw   list        — integer(2) per row
class_weights      list        — numeric(2) per row
ess                numeric     — NA when weighted
wess               numeric     — NA when unweighted
metric_label       character   — "ESS" | "WESS"
selected_metric_value numeric  — the actual score used for selection
```

`selected_metric_value` disambiguates: for `unioda_max_ess` it is ESS or WESS;
for `novometric_min_d` it is D (D does not change label with weights, but
its ESS component reflects WESS when weights are active).

### E. Prediction semantics for weights

Input fitting weights do **not** automatically apply to newdata predictions.

`predict.sda_fit(type = "weights")` means **SDA-derived propensity weights**
computed from SDA strata — not the original fitting weights. Do not overload.

These are two distinct concepts:
- Input case weights: supplied by the caller at `sda_fit()` time; used in
  fitting only; not stored for prediction.
- SDA-derived propensity weights: computed post-fit from endpoint/stage
  propensity estimates; output via `predict(type="weights")` or `sda_weights()`.

### F. Three weight concept summary

| Concept | Source | SDA role | Status |
|---------|--------|----------|--------|
| Input case weights | Caller → `sda_fit(weights=)` | Used in per-step ODA/MDSA fitting | Deferred to SDA-5+ |
| SDA-derived propensity weights | Computed from SDA strata post-fit | Output via `sda_weights()` / `predict(type="weights")` | Deferred to SDA-5 |
| Downstream CTA case weights | Output of SDA-5 workflow | Consumed by weighted `cta_fit()` / `cta_descendant_family()` | Deferred to SDA-5 |

These must remain separate. Do not conflate fitting weights with propensity
weights or with downstream CTA weights.

### Implementation boundary

| Slice | Weight support |
|-------|---------------|
| SDA-1 | `weights != NULL` errors; unweighted only |
| SDA-2 | No change |
| SDA-3 | No change |
| SDA-4 | No change |
| SDA-5 | Full weighted SDA: input weights, WESS labeling, weighted removal counts, extended step/candidate objects, propensity weight outputs |

---

## 11. Open Design Questions

### Q1 — Prediction: first-classifying-rule vs. all-step ensemble

This document specifies first-classifying-rule semantics (§4.2): the first SDA
step that fires for an observation is authoritative. An alternative is to run all
steps and return a step-indexed classification table. The first-classifying-rule
approach matches how SDA is conceptualized in MPE (each step resolves a subset
and removes them). Confirm before implementing `predict.sda_fit()`.

### Q2 — Correctly-classified removal: strict vs. relaxed

At each SDA step, "correctly classified" means the predicted class matches the
true class. In training data, this is unambiguous. For predict-time application,
observations that would have been "correctly classified" are simply "classified"
— there is no ground truth. The `resolved = TRUE/FALSE` flag handles this without
needing a post-hoc removal semantics. Confirm this framing is sufficient.

### Q3 — MINDENOM specification for novometric mode

`mindenom` is a per-analysis power decision, not a function default. SDA-1
(unioda_max_ess) does not require it. SDA-4 (novometric_min_d) will require
either a mandatory argument or a principled default. Current proposal: error
if not supplied with a canon-grounded message. Confirm before SDA-4.

### Q4 — Tie-breaking in min-D selection (novometric mode)

When two candidates have equal D (rounded to the same tick), the tie-breaking
rule is: first enumerated (enumeration order = column order of X). Confirm this
is consistent with CTA.exe tie-breaking and document in SDA_DESIGN.md.

### Q5 — `sda_weights()` propensity method

Inverse propensity weights from SDA strata require a propensity model within
each SDA stratum. The simplest approach is empirical stratum proportions
(class frequency within stratum = propensity estimate). Stabilized and overlap
variants require a numerator model. Full design deferred to SDA-5; confirm
scope before implementing.

### Q6 — `auto_sda_plan()` and `sda_fit()` integration

The `proposed_call` field in `auto_sda_plan()` is designed to be passed directly
to `sda_fit()` via `do.call(sda_fit, plan$proposed_call)`. Confirm the argument
names match exactly between `auto_sda_plan()` output and `sda_fit()` signature
before finalizing either.

### Q7 — Agent resumability: partial SDA

If `sda_fit()` is interrupted mid-run (e.g., by a compute timeout or an agent
pause), can it be resumed from the last completed step? The current object design
stores completed steps in `$steps`. A `sda_resume()` function that takes an
interrupted `sda_fit` object and continues from step `length(steps) + 1` is
possible but out of scope for SDA-1. Flag for SDA-5.

### D1 — Case weights: design decision (not open question)

**Decision: SDA-1 is unweighted only.**

If `weights != NULL` is supplied to `sda_fit()`, error immediately:

```
"Weighted SDA is not implemented yet; use weights = NULL.
Weighted SDA requires explicit weighted removal, WESS labeling,
and weighted candidate-table semantics."
```

This is not deferred for convenience. Weighted SDA requires explicit design
across three distinct dimensions — see §Weighted SDA design boundary for the
full specification of what must be designed before implementation.
