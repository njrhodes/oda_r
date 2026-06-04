# SDA_ANCHOR_CONTRACT.md

Design memo: SDA anchor object for future staged / SORT-style workflows.

**Status:** Implemented  -  `R/sda_anchor.R`, `tests/testthat/test-sda-anchor.R` (Slice O).  
**Design date:** 2026-05-29  
**Implementation date:** 2026-05-29  
**Slice:** I (design) + O (implementation)

---

## 1. Purpose

`sda_fit()` produces a sequential attribute-selection history.  Future
staged workflows  -  specifically SORT (Sequentially Optimal Recursive Tree,
reserved)  -  need a stable, typed object that captures that history and
constrains downstream CTA/MDSA use.

This document defines the **`sda_anchor`** object: its required fields,
optional fields, constructor/converter API, and the rules for how a future
`sort_fit()` will consume it.

**What this document is:**
- A design contract for `sda_anchor`.
- An audit of the current `sda_fit` surface against downstream needs.
- An explicit-anchor design for non-`sda_fit`-derived workflows.
- A future SORT consumption sketch.

**What this document is not:**
- An implementation of SORT.
- An implementation of GORT.
- An implementation of any staged CTA workflow.
- A change to CTA, LORT, ODA, or SDA fitting behavior.

---

## 2. Canon Boundary

| Layer | Canon anchor | Current status |
|-------|-------------|----------------|
| CTA / MDSA | CTA.exe golden fixtures (myeloma, CTA_DEMO) | Implemented, passing |
| UniODA / MultiODA | MegaODA.exe golden fixtures (iris) | Implemented, passing |
| SDA novometric_min_d | MPE / Novometrics per-step min-D | Implemented, passing |
| LORT | Workflow layer over canon CTA/MDSA components | Implemented; not CTA.exe verified |
| SORT | Reserved  -  not yet implemented | Requires SDA anchor before starting |
| GORT | Reserved  -  not yet implemented | Requires formal design doc |

**LORT / SORT / GORT are oda workflow-layer names, not MPE Chapter 12
canonical method names.**  Do not claim they are canon-anchored against
CTA.exe unless explicit fixture parity tests exist.

**`sda_anchor` is not a fitting object.**  It does not call ODA, CTA, or
SDA.  It is a typed structural object that carries SDA selection history
and metadata for downstream consumption.

---

## 3. Current SDA Object Audit

### 3.1 Exported SDA functions

| Function | Class returned | Description |
|----------|---------------|-------------|
| `sda_fit(X, y, mode, ...)` | `c("sda_fit", "oda_sda")` | Core SDA fitting |
| `predict.sda_fit(object, newdata, type, ...)` | vector / data frame | Sequential prediction |
| `print.sda_fit(x, ...)` | invisible | Concise step summary |
| `summary.sda_fit(object, ...)` | `sda_fit_summary` | Structured summary |
| `sda_selected_attributes(fit)` | character vector | Selected names in step order |
| `sda_step_table(fit)` | data frame | One row per completed step |
| `sda_candidate_table(fit, step)` | data frame / list | Per-step candidate detail |
| `as_cta_candidates(fit, X)` | data frame | X restricted to SDA columns |
| `sda_to_cta_data(fit, X, y)` | named list | Downstream CTA/MDSA-ready data |
| `auto_sda_plan(data, outcome, ...)` | `c("auto_sda_plan", "oda_plan")` | Dry-run planning only |

### 3.2 `sda_fit` object structure

Top-level fields:

```r
$call                     # match.call()
$mode                     # "novometric_min_d" | "unioda_max_ess"
$outcome_name             # NA_character_  <- not stored (see gap section 9.2)
$class_levels             # integer vector  -  sorted unique class values
$n_initial                # integer  -  total observations
$n_final_unresolved       # integer  -  remaining after last step
$candidate_names_initial  # character vector  -  all eligible attributes at start
$selected_attributes      # character vector  -  selected in step order
$steps                    # list  -  one element per completed step
$unresolved_indices       # integer vector  -  row indices remaining
$resolved_indices         # integer vector  -  row indices classified
$settings                 # list (see section 3.3)
$stop_reason              # character  -  "no_candidates" | "min_n" | "min_class_n" |
                          #             "class_resolved" | "all_resolved" |
                          #             "p_gate" | "axiom1_violated" |
                          #             "max_steps" | "dry_run"
$diagnostics              # list(warnings, excluded data frame)
```

Settings sublist:

```r
$settings$mode            # character
$settings$mc_iter         # integer
$settings$mc_seed         # integer or NULL
$settings$mc_stop         # numeric (MC early-stop percentile, default 99.9)
$settings$mc_stopup       # integer (MC upper bound parameter)
$settings$alpha           # numeric
$settings$loo             # character ("off" | "stable" | "pvalue" | numeric)
$settings$max_steps       # integer or NULL
$settings$min_n           # integer or NULL
$settings$min_class_n     # integer or NULL
$settings$remove_correct  # logical
$settings$collinearity    # character ("skip" | "warn" | "allow")
$settings$mindenom        # integer or NULL (novometric mode only)
```

### 3.3 Per-step object structure (novometric_min_d mode)

Each `fit$steps[[i]]` contains:

```r
$step_id                  # integer
$attribute                # character  -  selected attribute name
$mode                     # "novometric_min_d"
$n_in                     # integer  -  unresolved count at step start
$class_counts_in          # named integer vector
$model                    # cta_family object (full EO-CTA/MDSA per attribute)
$min_d_idx                # integer  -  index of selected family member
$rule_summary             # list(type, mindenom, strata, ess, d)
$ess                      # numeric
$d                        # numeric
$p_mc                     # numeric
$ge_count                 # integer  -  MC iterations yielding >= observed ESS
$iter_used                # integer  -  MC iterations used
$mindenom                 # integer  -  MINDENOM used at this step
$n_correct                # integer  -  observations classified at this step
$n_incorrect              # integer
$correct_indices_global   # integer vector  -  global row indices now resolved
$incorrect_indices_global # integer vector
$selected                 # logical
$reason                   # character  -  why selected / not selected
$candidate_table          # data frame  -  all candidates evaluated at this step
```

Candidate table columns (from `sda_candidate_table()`):

```
attribute, status, eligible, ineligible_reason, n, mindenom,
min_terminal_denom, ess, d, p_mc, ge_count, iter_used, strata,
selected, tied_objective, selected_by_tie_break, class_counts (list col)
```

### 3.4 Current accessor outputs (sda_interop.R)

`sda_step_table()` columns:
```
step_id, attribute, n_in, n_correct, n_incorrect, ess, d, p_mc, mindenom
```

`sda_to_cta_data()` returns `list(X_cta, y_cta)`:
- `X_cta`  -  data frame restricted to `sda_selected_attributes()` columns
- `y_cta`  -  full-length integer class vector
- No observations removed (Path B: SDA identifies subset; CTA sees full sample)

### 3.5 Current SDA constraints relevant to anchor design

- **Weights:** `sda_fit()` currently raises an error when `weights` is
  non-NULL.  `weights_used` is therefore always `FALSE` in current results.
- **LOO:** `loo` is passed through to per-step CTA/ODA calls but SDA has no
  dedicated LOO gate of its own.  LOO behavior is whatever the underlying
  `cta_descendant_family()` / `oda_fit()` does per step.
- **Branch candidate map:** There is no branch-level attribute constraint in
  current SDA.  SDA selects globally from the remaining candidate pool.
- **outcome_name:** Not stored in `sda_fit`.  The outcome variable name is
  not part of the fitted object (gap  -  see section 9.2).

---

## 4. Downstream CTA / LORT / SORT Needs

### 4.1 What LORT needs from SDA (current  -  `as_cta_candidates` path)

Currently, `as_cta_candidates(sda, X)` returns X restricted to SDA-selected
columns.  `lort_fit()` then consumes that as its full candidate pool.

This is the **unanchored greedy LORT** pattern.  LORT knows nothing about
step order, which step selected each attribute, or which MINDENOM / D values
governed selection.  The anchor information is discarded.

### 4.2 What future SORT needs

SORT is SDA-anchored sequential CTA/MDSA.  It needs:

| Field | Source in sda_fit | Need |
|-------|-------------------|------|
| Selected attribute sequence | `$selected_attributes` | Ordered list for stage-by-stage recursion |
| Full candidate universe | `$candidate_names_initial` | Constrain which attributes may enter at each stage |
| Stage ESS / D / p per step | `$steps[[i]]$ess`, `$d`, `$p_mc` | Auditability and reporting |
| Rule summary per step | `$steps[[i]]$rule_summary` | Reproducibility, staging description |
| MINDENOM used per step | `$steps[[i]]$mindenom` | SORT must respect same or stricter MINDENOM |
| LOO mode | `$settings$loo` | SORT inherits LOO state |
| MC metadata | `$settings$mc_iter`, `$mc_seed` | Reproducibility |
| Alpha | `$settings$alpha` | Significance gating |
| Class levels | `$class_levels` | Group labels for confusion/staging tables |
| n_initial, n_final | `$n_initial`, `$n_final_unresolved` | Scope |
| Stop reason | `$stop_reason` | Understanding SDA termination |
| Removal history | `$steps[[i]]$correct_indices_global` | Which obs were resolved and when |
| Branch candidate map | _not available_ | SORT-specific; per-branch attribute constraint |
| Weights metadata | _not applicable yet_ | Future when sda_fit supports weights |

### 4.3 What explicit/manual anchor needs to provide

A user-declared anchor (without a real `sda_fit` run) must supply at minimum:
- `selected_attributes` in intended step order
- `candidate_universe` (or a reasonable default)
- `stage_table` with at least `stage_id` and `attribute`
- `anchor_type = "explicit"` to prevent false claims

Everything else is optional when the anchor is explicit.

---

## 5. Proposed `sda_anchor` Object Schema

```r
class(anchor) <- c("sda_anchor", "list")
```

### 5.1 Required fields (always present)

| Field | Type | Description |
|-------|------|-------------|
| `anchor_type` | character | `"sda_fit"` or `"explicit"` |
| `source_class` | character | `class()` of source object, or `NA_character_` |
| `source_call` | call or NULL | Captured `match.call()` from source, or NULL |
| `group_levels` | integer vector | Sorted unique class values |
| `selected_attributes` | character vector | SDA-selected attributes in step order |
| `candidate_universe` | character vector | Full eligible attribute pool at SDA start |
| `stage_table` | data frame | One row per SDA step (see section 5.2) |
| `n_initial` | integer | Total observations at SDA start |
| `n_final_unresolved` | integer | Observations unresolved after last step |
| `stop_reason` | character | SDA stop reason |
| `alpha` | numeric | Significance threshold used |
| `mc_iter` | integer | MC iterations per step |

### 5.2 `stage_table` columns

```
stage_id       integer     -  step number (1-based)
attribute      character   -  selected attribute name
rule_summary   character   -  human-readable rule (e.g. "V17 <= 0.5")
cutpoint       numeric     -  ordered cut value; NA for categorical / stump
direction      character   -  "0->1" / "1->0" / NA for categorical
ess            numeric     -  ESS (%) at this step
wess           numeric     -  WESS (%); NA when weights not active
d_stat         numeric     -  parsimony-adjusted D statistic; NA if unavailable
p_value        numeric     -  MC p-value
loo_status     character   -  LOO result string or NA
n_in           integer     -  unresolved count at start of step
n_correct      integer     -  observations classified at this step
n_removed      integer     -  = n_correct (alias for clarity)
stop_reason    character   -  NA unless this is the stopping step
```

All columns are present.  Unavailable values are NA, not missing.

### 5.3 Optional fields (NULL when not applicable)

| Field | Type | When present |
|-------|------|-------------|
| `weights_used` | logical | Always: FALSE until sda_fit supports weights |
| `weight_summary` | character | When weights were active (future) |
| `loo_mode` | character | Inherited from SDA settings |
| `mc_seed` | integer | When seed was set |
| `mindenom` | integer | novometric_min_d mode only |
| `removal_history` | data frame | When `remove_correct = TRUE` (see section 5.4) |
| `branch_candidate_map` | named list | SORT-specific; NULL until SORT is designed |
| `reproducibility_notes` | character | Free-text audit trail |
| `canon_notes` | character | Boundary statement for the anchor source |

### 5.4 `removal_history` data frame (optional)

When `remove_correct = TRUE` in the source `sda_fit`, observations are
removed from the unresolved pool at each step.  The removal history captures
this:

```
obs_index    integer    -  row index in the original X / y
removed_at   integer    -  step_id at which this observation was classified/removed
attribute    character  -  attribute that classified this observation
class_pred   integer    -  predicted class at removal
```

NULL when the source is an explicit anchor (removal tracking was not performed)
or when `remove_correct = FALSE`.

---

## 6. Proposed Constructor / Converter / Validator API

**These are designs.  No implementation in this slice.**

### 6.1 `as_sda_anchor` generic

```r
as_sda_anchor <- function(x, ...) UseMethod("as_sda_anchor")
```

#### `as_sda_anchor.sda_fit(x, include_removal_history = FALSE, ...)`

Extracts all required and available optional fields from a fitted `sda_fit`
object.  Sets `anchor_type = "sda_fit"`.

Key derivations:
- `stage_table$rule_summary`  -  from `x$steps[[i]]$rule_summary$type` and
  root-node attribute/cutpoint in the selected family member.
- `stage_table$cutpoint`  -  from `x$steps[[i]]$model$members[[min_d_idx]]$tree$nodes[[1]]$cut_value`
  when the rule is `ordered_cut`; NA otherwise.
- `removal_history`  -  built from `x$steps[[i]]$correct_indices_global` when
  `include_removal_history = TRUE`.
- `weights_used`  -  always FALSE (current sda_fit errors on weights).
- `canon_notes`  -  populated with mode label, whether novometric_min_d or
  unioda_max_ess, and a note that SORT has not been run.

#### `as_sda_anchor.data.frame(x, selected_attributes, group_levels, n_initial, ...)`

Explicit/manual path.  See section 7.

### 6.2 `sda_anchor()` direct constructor

```r
sda_anchor(
  anchor_type,
  group_levels,
  selected_attributes,
  candidate_universe,
  stage_table,
  n_initial,
  n_final_unresolved,
  stop_reason,
  alpha,
  mc_iter,
  # optional:
  source_class          = NULL,
  source_call           = NULL,
  weights_used          = FALSE,
  weight_summary        = NULL,
  loo_mode              = "off",
  mc_seed               = NULL,
  mindenom              = NULL,
  removal_history       = NULL,
  branch_candidate_map  = NULL,
  reproducibility_notes = NA_character_,
  canon_notes           = NA_character_
)
```

Returns a validated `sda_anchor` object.  Calls `validate_sda_anchor()`
internally before returning.

### 6.3 `validate_sda_anchor(anchor)`

Checks:
1. `inherits(anchor, "sda_anchor")`  -  class present.
2. Required fields all present and non-NULL.
3. `anchor_type` is `"sda_fit"` or `"explicit"`.
4. `selected_attributes` is a non-empty character vector.
5. All `selected_attributes` are in `candidate_universe`.
6. `stage_table` has required columns (`stage_id`, `attribute`, `n_in`,
   `n_correct`).
7. `nrow(stage_table) == length(selected_attributes)`.
8. `stage_table$attribute` matches `selected_attributes` in order.
9. `group_levels` is an integer vector with >= 2 distinct values.
10. `alpha` is a length-1 numeric in (0, 1).
11. If `removal_history` is not NULL, has required columns (`obs_index`,
    `removed_at`, `attribute`, `class_pred`).

Returns `anchor` invisibly on success.  Errors clearly on failure, naming
the failed check.

### 6.4 `print.sda_anchor(x, ...)`

Concise output:

```
SDA Anchor
  anchor_type : sda_fit  [or: explicit]
  n_initial   : 255
  n_unresolved: 0
  stop_reason : all_resolved
  steps       : 3  (V17, V15, V9)
  alpha       : 0.05
  mc_iter     : 5000
  mc_seed     : 42
  mindenom    : 30

Stage table:
  stage attribute  n_in n_correct  ESS   D     p
  1     V17         255      186   16.5  10.1  0.002
  2     V15          69       30   ...
  3     V9           39       ...
```

The word "SORT" must not appear as an implemented method.  At most a note:
`"[anchor ready for future SORT workflow  -  sort_fit() not yet implemented]"`.

### 6.5 `summary.sda_anchor(object, ...)`

Returns a list of class `"sda_anchor_summary"` with:
- All top-level fields minus `stage_table` and `removal_history`.
- A `stage_summary` data frame (same as `stage_table` but with added
  formatted columns for display).

---

## 7. Explicit / Manual Anchor

A researcher may want to declare a SORT-style attribute sequence without
having run `sda_fit()`.  The explicit anchor path supports this.

### 7.1 Design constraints

- `anchor_type` must be `"explicit"`.
- `print.sda_anchor` must label it as `"anchor_type: explicit (not SDA-derived)"`.
- `summary.sda_anchor` must include a warning line when `anchor_type = "explicit"`:
  `"WARNING: This anchor was not derived from sda_fit(). It is not canonical SDA output."`
- No claim of SDA canon, p-values, or ESS/D statistics unless the user provides
  them explicitly in `stage_table`.

### 7.2 Minimal explicit anchor example

```r
anchor <- sda_anchor(
  anchor_type         = "explicit",
  group_levels        = c(0L, 1L),
  selected_attributes = c("V17", "V15", "V9"),
  candidate_universe  = c("V4", "V9", "V11", "V14", "V15", "V16", "V17",
                          "V18", "V19"),
  stage_table = data.frame(
    stage_id      = 1:3,
    attribute     = c("V17", "V15", "V9"),
    rule_summary  = c("V17 <= 0.5", "V15 <= 0.5", "V9 <= 2.0"),
    cutpoint      = c(0.5, 0.5, 2.0),
    direction     = c("0->1", "0->1", "0->1"),
    ess           = c(16.5, NA_real_, NA_real_),
    wess          = c(NA_real_, NA_real_, NA_real_),
    d_stat        = c(10.1, NA_real_, NA_real_),
    p_value       = c(0.002, NA_real_, NA_real_),
    loo_status    = c(NA_character_, NA_character_, NA_character_),
    n_in          = c(255L, 69L, 39L),
    n_correct     = c(186L, 30L, NA_integer_),
    n_removed     = c(186L, 30L, NA_integer_),
    stop_reason   = c(NA_character_, NA_character_, "max_steps"),
    stringsAsFactors = FALSE
  ),
  n_initial        = 255L,
  n_final_unresolved = NA_integer_,   # unknown for explicit anchor
  stop_reason      = "explicit",
  alpha            = 0.05,
  mc_iter          = 5000L
)
```

`validate_sda_anchor(anchor)` must pass on a well-formed explicit anchor.

### 7.3 When to use explicit anchors

- Reproducing a prior published SDA result without the original `sda_fit`
  output.
- Researcher specifies an a-priori attribute order based on domain knowledge
  (not data-driven SDA).
- Comparison studies: running SORT on a fixed attribute sequence to compare
  against a data-driven SDA-anchored SORT result.

All three use cases require explicit labeling.  No results from an explicit
anchor should be reported as "SDA-derived" or "canonical SDA output."

---

## 8. Future SORT Consumption Sketch

**SORT is not implemented.  This section describes the intended API only.**

```r
# Step 1  -  Run SDA
sda <- sda_fit(X, y,
               mode     = "novometric_min_d",
               mindenom = 30L,
               mc_iter  = 5000L,
               mc_seed  = 42L)

# Step 2  -  Build anchor
anchor <- as_sda_anchor(sda)
validate_sda_anchor(anchor)      # errors if malformed

# Step 3  -  Run SORT  [NOT YET IMPLEMENTED]
# sort_fit() consumes:
#   anchor$selected_attributes    -  ordered attribute list
#   anchor$candidate_universe     -  eligible pool
#   anchor$stage_table            -  per-step rule history and MINDENOM
#   anchor$alpha, mc_iter, etc.   -  reproducibility metadata
sort_result <- sort_fit(X, y, anchor = anchor, ...)
```

### 8.1 How SORT would use the anchor

- `anchor$selected_attributes[1]` -> first split attribute (root).  SORT
  does not re-select the root; it is given by the anchor.
- `anchor$selected_attributes[2:k]` -> ordered candidate priorities for
  subsequent stages.
- `anchor$candidate_universe` -> upper bound on attributes SORT may consider
  at any branch.  SORT may not introduce attributes outside this universe.
- `anchor$stage_table` -> provides per-step MINDENOM context.  SORT must
  respect at least as strict a MINDENOM at each stage.
- `anchor$branch_candidate_map` (if present, future) -> maps recursive branch
  IDs to allowed attribute subsets.  When NULL, SORT defaults to
  `candidate_universe` minus already-used attributes.

### 8.2 What SORT must NOT do

- Re-run or modify SDA.
- Alter `lort_fit()` behavior.
- Alter `cta_fit()` behavior.
- Claim global optimization (it is sequential, not global).
- Use `degen = TRUE` (forbidden for CTA/LORT/SORT per degeneracy policy).
- Mutate the `sda_anchor` object.

### 8.3 SORT object metadata (future)

When SORT is eventually implemented:

```r
sort_result$sort_settings$method             == "sort"
sort_result$sort_settings$method_label       == "Sequentially Optimal Recursive Tree"
sort_result$sort_settings$sda_anchored       == TRUE
sort_result$sort_settings$global_optimization == FALSE
sort_result$sort_settings$anchor_type        == anchor$anchor_type  # "sda_fit" or "explicit"
```

---

## 9. Gap Analysis

### 9.1 Fields available in current `sda_fit`

These map directly to `sda_anchor` required fields with no additional work:

| sda_anchor field | Source |
|-----------------|--------|
| `anchor_type` | Derived: `"sda_fit"` |
| `source_class` | `class(sda)` |
| `source_call` | `sda$call` |
| `group_levels` | `sda$class_levels` |
| `selected_attributes` | `sda$selected_attributes` |
| `candidate_universe` | `sda$candidate_names_initial` |
| `n_initial` | `sda$n_initial` |
| `n_final_unresolved` | `sda$n_final_unresolved` |
| `stop_reason` | `sda$stop_reason` |
| `alpha` | `sda$settings$alpha` |
| `mc_iter` | `sda$settings$mc_iter` |
| `mc_seed` | `sda$settings$mc_seed` |
| `mindenom` | `sda$settings$mindenom` |
| `loo_mode` | `sda$settings$loo` |
| `weights_used` | Always FALSE (sda_fit errors on weights) |
| `stage_table$stage_id` | `step$step_id` |
| `stage_table$attribute` | `step$attribute` |
| `stage_table$ess` | `step$ess` |
| `stage_table$d_stat` | `step$d` |
| `stage_table$p_value` | `step$p_mc` |
| `stage_table$n_in` | `step$n_in` |
| `stage_table$n_correct` | `step$n_correct` |

### 9.2 Fields derivable from current object (require light extraction)

| sda_anchor field | How to derive |
|-----------------|---------------|
| `stage_table$rule_summary` | From `step$rule_summary$type` + root node of `step$model$members[[min_d_idx]]$tree` |
| `stage_table$cutpoint` | Root node `cut_value` when `rule_type == "ordered_cut"`; NA otherwise |
| `stage_table$direction` | Root node `direction` when `rule_type == "ordered_cut"`; NA otherwise |
| `stage_table$n_removed` | Same as `n_correct` (= observations removed from unresolved pool) |
| `removal_history` | From `step$correct_indices_global` at each step |
| `stage_table$wess` | From root node of selected family member when weights active |
| `stage_table$loo_status` | From selected family member LOO metadata |

### 9.3 Fields unavailable  -  desirable but require future work

| Field | Gap | Work required |
|-------|-----|---------------|
| `outcome_name` | `sda_fit$outcome_name` is always `NA_character_` | Store column name in sda_fit call |
| `weight_summary` | sda_fit errors on weights | Implement weighted SDA |
| `weights_used = TRUE` | Not possible today | Implement weighted SDA |
| `branch_candidate_map` | No concept of branch in SDA | SORT design work |
| Reproducibility hash | Not implemented | Future (SHA of call + settings + seed) |

### 9.4 Fields that should remain optional

- `removal_history`  -  large; useful for auditing but not required for SORT.
- `weight_summary`  -  not applicable until SDA supports weights.
- `branch_candidate_map`  -  SORT-specific; NULL is the correct default.
- `reproducibility_notes` / `canon_notes`  -  free-text audit fields.
- `stage_table$loo_status`  -  LOO in SDA context is per-step CTA LOO;
  should be surfaced when available but may be NA for most current fits.

### 9.5 Fields that require future implementation

| Field | Blocker |
|-------|---------|
| Weighted SDA -> anchor | `sda_fit()` currently errors on weights |
| `branch_candidate_map` | SORT design doc (not yet written) |
| LOO per step in stage_table | SDA does not currently surface per-step LOO status |

---

## 10. Future Tests

These tests are proposed for when `as_sda_anchor()` and `sda_anchor()` are
implemented.  **Do not write these tests in Slice I.**

| Test | Description |
|------|-------------|
| A1 | `as_sda_anchor(sda_fit_object)` returns class `sda_anchor` |
| A2 | All required fields present in result |
| A3 | `anchor$selected_attributes` order matches `sda$selected_attributes` |
| A4 | `stage_table` has correct number of rows and column names |
| A5 | `stage_table$attribute` matches `selected_attributes` in order |
| A6 | Explicit anchor with missing `selected_attributes` errors in `validate_sda_anchor` |
| A7 | Explicit anchor with `selected_attributes` outside `candidate_universe` errors |
| A8 | `validate_sda_anchor` catches wrong `anchor_type` value |
| A9 | `validate_sda_anchor` catches non-integer `group_levels` |
| A10 | `print.sda_anchor` output does not contain the word "SORT" as implemented |
| A11 | `summary.sda_anchor` includes `"explicit"` warning when anchor_type is explicit |
| A12 | CTA behavior unchanged: `cta_fit(X, y)` on the same data produces same result before and after sda_anchor landing |
| A13 | LORT behavior unchanged: `lort_fit(X, y)` produces same result |
| A14 | `sort_fit` is not exported (not implemented) |
| A15 | `gort_fit` is not exported (not implemented) |

---

## 11. Deferred Work

| Item | Why deferred |
|------|-------------|
| `as_sda_anchor.data.frame()` method | Use case not fully determined; explicit anchor path via `sda_anchor()` constructor is sufficient for now |
| Weighted SDA -> anchor | `sda_fit()` errors on weights; defer until weighted SDA is implemented |
| `branch_candidate_map` population | Requires SORT design doc; no SDA concept of branch-level candidates |
| `sort_fit()` implementation | Requires: (a) this anchor contract, (b) dedicated SORT design doc, (c) SDA source object or approved explicit anchor |
| `gort_fit()` implementation | Requires approved formal design doc; far future |
| LOO per step in `stage_table` | SDA currently does not surface per-step LOO status; add when SDA LOO is designed |
| Reproducibility hash | Future hardening; not needed for SORT anchor consumption |
| `sda_anchor` as formal R class | Will be implemented alongside `as_sda_anchor()` in a future implementation slice |

---

## 12. Confirmation: No SORT / GORT Implementation

This document is a design memo.  It contains:

- A field contract for the `sda_anchor` object.
- An API design for `as_sda_anchor()`, `sda_anchor()`, `validate_sda_anchor()`,
  `print.sda_anchor()`, `summary.sda_anchor()`.
- A SORT consumption sketch (future sketch only).
- A gap analysis and test plan.

It does **not** contain:

- Any implementation of SORT.
- Any implementation of GORT.
- Any change to `cta_fit()` behavior.
- Any change to `lort_fit()` behavior.
- Any change to `sda_fit()` behavior.
- Any change to the ODA / CTA optimization objective.

SORT and GORT remain reserved.  `sort_fit()` and `gort_fit()` are not
exported names and must not be introduced without the required preconditions
stated in `LORT_SORT_GORT_TAXONOMY.md section Agent Handoff Contract`.

---

## 13. Next Implementation Slice

When approved, the implementation slice for `sda_anchor` will:

1. Add `as_sda_anchor.sda_fit()` in `R/sda_interop.R` (or a new
   `R/sda_anchor.R`).
2. Add `sda_anchor()` direct constructor.
3. Add `validate_sda_anchor()`.
4. Add `print.sda_anchor()` and `summary.sda_anchor()`.
5. Export the new functions in NAMESPACE.
6. Add S3method registrations for `as_sda_anchor`.
7. Add tests A1-A15 (section 10) in `tests/testthat/test-sda-anchor.R`.
8. Update `SDA_AUTO_SDA_PLAN.md` and `DOCS_INDEX.md` to mark anchor
   contract as implemented.

Preconditions before starting that slice:
- This design doc (`SDA_ANCHOR_CONTRACT.md`) reviewed and approved.
- No pending changes to `sda_fit()` or `cta_descendant_family()` that would
  alter step object structure.
