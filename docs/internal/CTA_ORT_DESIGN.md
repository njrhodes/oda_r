# LORT (Locally Optimal Recursive Tree) — Design Document

> **Maintainer historical design memo. Not public API. Do not treat as implemented behavior.**

_Formerly: "CTA Optimal Recursive Tree (ORT)"_

**Status:** Finalized design. Implementation in R/cta_ort.R.
**Date:** 2026-05-26. **Updated:** 2026-05-28 (named LORT).
**Motivation:** Private rare-event staged workflow example (2026-05-25/26).

**Method taxonomy (2026-05-28):** `cta_fit(recursive = TRUE)` implements **LORT**
— greedy local min-D recursive CTA/MDSA. "ORT" remains as legacy S3 class name
(`cta_ort`). New code and docs use LORT/SORT/GORT. See
`docs/LORT_SORT_GORT_TAXONOMY.md` for the full agent handoff contract.

---

## 1. Problem

CTA.exe requires the analyst to manually chain models using `EX` filter commands.
Each sequential analysis is a separate program invocation:

```text
Run 1: full cohort  → MDSA family scan → min-D model → v21 splits
Run 2: EX v21>0.5   → MDSA family scan → min-D model → v19 splits
Run 3: EX v21<=0.5  → MDSA family scan → min-D model → v20 splits
Run 4: EX v21>0.5 v19<=0.5  → no-tree (terminal confirmed)
...
```

There is no automated stopping criterion, no composite prediction object, and no
single plot spanning the recursive structure. The analyst must track which strata
are terminal by inspecting each `MODEL*.TXT` output by hand.

`odacore` currently replicates this behavior faithfully but does not automate it.
`cta_fit(..., recursive = TRUE)` automates the recursion while preserving full
canon at each step.

---

## 2. Motivating Workflow: Private Rare-Event Staged Workflow Example

Dataset: 9,014 observations; class = rare target event (binary).
Exclusion criterion applied (n=8,898; 72 positive-class; 0.81%).
Attributes: v17, v19, v20, v21, v27 (all binary baseline covariates).

### Level 0 — Root

MDSA family scan on full cohort. Min-D member: MINDENOM=1354 stump on **v21**
(baseline covariate A), ESS=18.6%, p≈0, LOO STABLE.

```
v21 >  0.5  →  Attribute A present  (n=2189, positive-class=1.42%)  ← recurse RIGHT first
v21 <= 0.5  →  Attribute A absent   (n=6709, positive-class=0.61%)  ← recurse LEFT second
```

### Level 1 (right): Attribute A present stratum

MDSA family on n=2189. Min-D: MINDENOM=1 stump on **v20** (baseline covariate B),
ESS=19.5%, p=0.020, LOO STABLE.

```
v20 >  0.5  →  Attribute B present  (n=426,  positive-class=2.82%)  ← terminal (no-tree)
v20 <= 0.5  →  Attribute B absent   (n=1763, positive-class=1.08%)  ← terminal (no-tree)
```

### Level 1 (left): Attribute A absent stratum

MDSA family on n=6709. Min-D: MINDENOM=1 stump on **v19** (baseline covariate C),
ESS=5.9%, p=0.013, LOO STABLE.

```
v19 >  0.5  →  Attribute C present  (n=98,   positive-class=3.06%)  ← terminal (no-tree)
v19 <= 0.5  →  Attribute C absent   (n=6611, positive-class=0.58%)  ← terminal (no-tree)
```

### Terminal strata — ordered by positive-class proportion

| Stratum | Path | n | Positive-class rate |
|---------|------|---|---------------------|
| 1 (lowest)  | v21≤0.5 & v19≤0.5 | 6611 | 0.58% |
| 2           | v21>0.5 & v20≤0.5 | 1763 | 1.08% |
| 3           | v21>0.5 & v20>0.5 | 426  | 2.82% |
| 4 (highest) | v21≤0.5 & v19>0.5 | 98   | 3.06% |

Three binary predictors carve 8,898 observations into four strata (0.58%–3.06% positive-class rate)
via two MDSA-optimal stumps in sequence — currently requiring four separate `GO;`
commands with manual EX filters.

---

## 3. Proposed API

### 3.1 Entry point: `cta_fit(..., recursive = TRUE)`

No new wrapper function. `recursive` is an argument to the existing `cta_fit()`.
When `recursive = FALSE` (default), behavior is identical to current — no code
path changes. When `recursive = TRUE`, an internal engine handles the recursion
and returns a `cta_ort` / `cta_tree` composite object.

**Training X contract:** For recursive CTA, `X` is the declared candidate
predictor frame — pass only variables eligible for CTA search.  Prediction may
be performed on wider `newdata` (e.g., the full clinical table) as long as the
split variable names are present; internal fitting remains positional, and names
are used as a safe cross-reference when `ncol(newdata) != ncol(X)`.

```r
cta_fit(
  X, y, w = NULL,

  # --- existing args unchanged ---
  mindenom     = 1L,          # ignored when recursive=TRUE (MDSA determines per level)
  mc_iter      = 5000L,
  mc_seed      = 42L,         # single seed — set once, stream consumed per node
  alpha_split  = 0.05,
  prune_alpha  = 0.05,
  loo          = "stable",
  verbose      = FALSE,

  # --- recursive guards ---
  recursive    = FALSE,       # if TRUE, enables ORT recursion
  min_n        = 30L,         # minimum endpoint n to attempt recursion
  max_depth    = 8L,          # safety cap on recursion depth
  max_nodes    = 31L,         # safety cap on total ORT nodes

  # --- MDSA family scan budget (recursive-only) ---
  family_max_steps = 20L      # max CTA fits per-node MDSA scan; default = cta_descendant_family() default
)
```

### 3.2 MINDENOM: determined by MDSA at each level

When `recursive = TRUE`, `mindenom` is NOT a user parameter. At every node, the
full MDSA descendant family scan is run (MINDENOM = 1 to n), and the min-D member
is selected as the split for that level. This is the canonical novometric procedure.

`cta_descendant_family()` is the inner workhorse at each recursive node.

### 3.3 Single seed / single RNG stream

`mc_seed` initializes the RNG once at the start of `cta_fit(recursive=TRUE)`.
Child-node MDSA scans consume the RNG stream in deterministic right-then-left
traversal order; seeds are not reset per node.

Reproducibility: calling `cta_fit(recursive=TRUE, mc_seed=s)` twice with the
same `s` gives identical results. The stored `ort_settings$mc_seed` records
the top-level seed for auditing.

### 3.4 Traversal order: Right then Left

After the MDSA-optimal split is selected at a node:
1. Recurse into the **right** endpoint first (high-value / "positive" branch).
2. Then recurse into the **left** endpoint (low-value / "negative" branch).
3. Within each sub-branch, the same original attribute set is available.

### 3.5 Attribute passing: full set, natural self-elimination

The complete original attribute set is passed to every child node. No forced
exclusion of previously used attributes. Binary attributes used on an ancestor
path self-eliminate in practice — if the parent split on v21, the v21>0.5
sub-stratum has v21=1 for every observation; v21 adds no information and will
not be selected as a split. This is natural self-elimination, not a policy rule.

### 3.6 Recursion guards

Guard evaluation order (first match wins):

1. `max_nodes`: if total ORT nodes already exceed `max_nodes`, the current
   endpoint becomes terminal without a sub-fit. Stop reason: `"max_nodes"`.
2. `min_n`: if an endpoint has fewer than `min_n` observations, it becomes
   terminal immediately. Stop reason: `"min_n"`.
3. `max_depth`: if `depth >= max_depth`, the node becomes terminal.
   Stop reason: `"max_depth"`.
4. No-tree: if the MDSA family scan finds no admissible tree (all members
   are no-tree), the node is terminal. Stop reason: `"no_tree"`.

### 3.7 Compute budget: `family_max_steps`

Each ORT node runs a full `cta_descendant_family()` scan to select its min-D
winning model.  `cta_descendant_family()` steps from MINDENOM=1 upward, fitting
one `oda_cta_fit()` per step, until the chain produces a no-tree result or
`max_steps` is reached.

On large wide datasets (n >> 1000, p >> 20), each `oda_cta_fit()` call is
expensive.  Without a budget cap, a single ORT node scan may run up to 20 full
CTA fits — several minutes for diagnostic/iteration work.

`family_max_steps` controls this:

- **Default `20L`**: preserves `cta_descendant_family()` behavior.  Full MDSA
  family is evaluated at each node.
- **Smaller values** (e.g. `5L`): cap the per-node scan at 5 CTA fits.  Useful
  for diagnostic runs where speed matters more than full family coverage.
  The min-D winner is selected from whichever family members are evaluated.
  If the true min-D would appear at step 6+, it is missed.
- Stored in `ort_settings$family_max_steps` for auditing.
- Only used when `recursive = TRUE`; error if explicitly supplied with
  `recursive = FALSE`.

**Important:** `family_max_steps` bounds the MDSA search within each node.  It
does not change recursion semantics.  ORT always recurses only the min-D winner's
endpoints — right first, then left.  Family losers never become recursive branches,
regardless of `family_max_steps`.

- `max_depth` / `max_nodes`: bound recursive topology (how deep and wide the ORT grows).
- `family_max_steps`: bounds per-node MDSA search (how many CTA fits per node).

---

## 4. Object Structure

### 4.1 Class: dual-tag `cta_ort` / `cta_tree`

The returned object carries two class tags:

```r
class(obj) == c("cta_ort", "cta_tree")
```

- `inherits(obj, "cta_tree")` is TRUE → all existing `cta_tree` S3 methods
  (predict, print, summary, plot) dispatch and work on the root-level model
  without modification.
- `inherits(obj, "cta_ort")` is TRUE → recursive-aware methods dispatch for
  full composite-tree operations (predict through sub-trees, plot full composite).

Non-recursive `cta_fit()` returns `class = "cta_tree"` as today — unchanged.

### 4.2 ORT node structure

Each node in the composite tree is a list stored under `$ort_nodes`:

```r
list(
  node_id         = integer,    # unique across the entire composite tree
  depth           = integer,    # recursion depth (root = 0)

  path            = character,  # full conjunction string: "v21>0.5 AND v20<=0.5"
  path_conditions = character,  # character vector of per-level conditions

  n               = integer,    # n in this subset
  class_counts    = integer,    # named raw class counts (e.g. c("0"=18, "1"=2))

  # Traversal data -- stored for auditability
  level_mindenom  = integer,    # MINDENOM of the MDSA min-D model selected here
  level_ess       = numeric,    # ESS of the model selected here (NA if no-tree)
  level_d         = numeric,    # D-statistic of the model selected here (NA if no-tree)
  path_ess        = numeric,    # vector of ESS values from root to this node
  path_mindenom   = integer,    # vector of MINDENOM values from root to this node

  model           = cta_tree,   # the cta_tree (MDSA min-D member) fitted here; NULL if guard fired

  is_terminal     = logical,    # TRUE if no further recursion
  stop_reason     = character,  # "no_tree" | "min_n" | "max_depth" | "max_nodes"
  terminal_class  = integer,    # majority class in this subset (if terminal; NA otherwise)

  child_ids       = integer,    # vector of child ORT node_ids (length 0 if terminal)
  right_child_id  = integer,    # node_id of right child (NA if terminal)
  left_child_id   = integer     # node_id of left child (NA if terminal)
)
```

`mc_seed` is not stored per node.  The seed is set once in `.cta_ort_fit()` and
child MDSA scans consume the RNG stream in deterministic right-then-left traversal
order.  See Section 3.3.

### 4.3 Top-level `cta_ort` / `cta_tree` object fields

All existing `cta_tree` fields are preserved (root-level model is still accessible
as a standard `cta_tree`). Added fields for the recursive structure:

```r
$recursive       = TRUE
$ort_nodes       = list of ort_node   # flat indexed list, all nodes
$ort_root_id     = 1L                 # always 1
$strata          = data.frame         # flat terminal-strata table (see §4.4)
$n_strata        = integer
$ort_settings    = list(
  mc_seed          = integer,
  mc_iter          = integer,
  mc_stop          = numeric,
  mc_stopup        = numeric,
  alpha_split      = numeric,
  prune_alpha      = numeric,
  loo              = character,
  min_n            = integer,
  max_depth        = integer,
  max_nodes        = integer,
  family_max_steps = integer
)
```

### 4.4 Flat strata table

One row per terminal stratum, sorted by ascending target-class proportion
(Stage 1 = lowest risk, canonical staging order):

| Column | Type | Description |
|--------|------|-------------|
| `stratum_id` | integer | sequential, ascending by `prop_class1` |
| `node_id` | integer | ORT node_id of this terminal |
| `path` | character | full conjunction of conditions |
| `depth` | integer | recursion depth at termination |
| `n` | integer | n in this stratum |
| `class_counts` | list | raw counts per class |
| `prop_class1` | numeric | class 1 proportion |
| `odds_class1` | numeric | class 1 odds |
| `terminal_class` | integer | majority class prediction |
| `stop_reason` | character | why recursion stopped |
| `path_ess` | character | ESS values along path, e.g. `"18.6 → 5.9"` |
| `path_mindenom` | character | MINDENOM at each level, e.g. `"1354 → 1"` |

Sorted ascending by `prop_class1` so stratum 1 = lowest risk = Stage 1.

---

## 5. Prediction Semantics

### 5.1 `predict.cta_ort()`

S3 dispatch on `cta_ort` class. For each row of `newdata`:

1. Start at the root ORT node (node_id = 1).
2. Apply the node's `model` (a `cta_tree`) via `predict.cta_tree()` to route the
   observation to an endpoint (right or left).
3. Follow the corresponding child node (`right_child_id` or `left_child_id`).
4. Repeat until `is_terminal == TRUE`.
5. Return `terminal_class` plus optional `stratum_id` and `path`.

```r
predict.cta_ort(object, newdata,
                type = c("class", "stratum", "path", "all"),
                missing_action = c("na", "majority"))
```

`type = "all"` returns a data.frame:
`predicted_class`, `stratum_id`, `path`, `prop_class1`, `stop_reason`.

Path-local missing: if an observation is missing the split attribute at any node
on its actual traversal path and `missing_action = "na"`, it returns `NA_integer_`
for class. This extends the canonical `predict.cta_tree()` behavior to the
composite tree.

### 5.2 Lean prediction

The composite tree never stores training row indices. Routing is purely via
split rules in the fitted `cta_tree` models. The lean-fit invariant extends to
the ORT: no raw X/y stored at any node.

### 5.3 Strata/predict consistency check (internal)

At fit time, route all training observations through `predict.cta_ort()` and
verify that the n per stratum matches `strata$n`. This is a self-consistency check
computed once at the end of `cta_fit(recursive=TRUE)` and stored as
`$strata_check_passed = TRUE/FALSE`. It is not re-run on new data.

---

## 6. Plot Semantics

### 6.1 `plot.cta_ort()`

Renders the full composite tree using G1 base-R conventions:
- **Ellipses** for split nodes (each MDSA-optimal split at its depth level).
- **Rectangles** for terminal leaf nodes.
- **Directed arrows** for edges, labeled with the branch condition.
- Split nodes show: attribute name, ESS at this level, MINDENOM used, n.
- Terminal leaves show: n, majority class, and (if `target_class` supplied)
  target-class proportion, stratum_id.

```r
plot.cta_ort(object, target_class = NULL,
             class_labels = NULL, digits = 1,
             main = "ORT", show_caption = FALSE,
             cex = 0.75, ...)
```

### 6.2 Layout

The composite tree layout extends the current `cta_plot_data()` approach.
Terminal strata are placed at the lowest occupied depth, evenly spaced.
Internal nodes are centered over their children.

A new `ort_plot_data()` function (analogous to `cta_plot_data()`) provides the
renderer-independent data contract for the composite tree layout. `cta_plot_data()`
is unchanged.

**v1 scope:** correct for balanced or near-balanced trees up to depth 3 (a private
rare-event example). Deeper or highly unbalanced trees may need layout tuning deferred to
graphics v3.

### 6.3 Not in scope for plot v1

- Confusion panel.
- Per-stratum staging annotation beyond proportion label.
- ggplot/grid renderer.

---

## 7. Stop Reasons

| Reason | Condition | Sub-fit attempted? |
|--------|-----------|-------------------|
| `"no_tree"` | MDSA family scan finds no admissible tree in this stratum | Yes — confirmed no-tree |
| `"max_nodes"` | Total ORT nodes > `max_nodes` | No |
| `"min_n"` | Endpoint n < `min_n` | No |
| `"max_depth"` | `depth >= max_depth` | No |

`"no_tree"` is the canonical stopping criterion.
`"min_n"`, `"max_depth"`, and `"max_nodes"` are safety guards that may mask
further structure; both are stored and inspectable in the node and strata table.

---

## 8. Canon and Safety Boundaries

### 8.1 `cta_fit()` non-recursive path unchanged

`recursive = FALSE` (default) follows the identical code path as today.
Zero impact on existing behavior, tests, or fixtures.

### 8.2 CTA core unchanged

All ENUMERATE, LOO STABLE, PRUNE, and MINDENOM logic is delegated entirely to
`cta_descendant_family()` and `cta_fit()` at each node. The ORT engine is a
coordination wrapper — not a new search algorithm.

### 8.3 Lean-fit invariant extends to ORT

The `cta_ort` object stores fitted `cta_tree` models (each lean) plus the flat
strata table. No training rows, no per-stratum X/y, no cached membership vectors.

### 8.4 Each recursive child fit is an explicit modeling step

No global optimization across levels. No shared permutation budget. The composite
tree reflects the sequential analyst workflow: each node is an independent,
self-contained `cta_descendant_family()` run on its data subset.

### 8.5 Attribution

Results from `cta_fit(recursive=TRUE)` are "novometric optimal recursive trees"
or "sequential ODA-based recursive partitioning." Not CART, not random forest,
not any ensemble method. The node criterion at each level is canonical ODA ESS
via the MDSA family selection procedure.

---

## 9. Testing Plan

### 9.1 Synthetic recursive tree (CRAN-safe)

Construct a small dataset with known two-level structure:
- Level 0: attribute A splits cleanly (p < 0.05, LOO STABLE).
- Level 1 right (A=1): attribute B splits cleanly → two terminal leaves.
- Level 1 left (A=0): no admissible split → no-tree immediately.

Assert:
- 3 terminal strata (2 from level 1 right + 1 from level 1 left).
- `strata$n` sums to N.
- `strata$stop_reason`: two `"no_tree"` from right branch leaves, one
  `"no_tree"` from left branch.
- `predict(obj, X_train, type="stratum")` row counts match `strata$n` exactly
  (`strata_check_passed = TRUE`).
- `predict(obj, newdata)` routes correctly for held-out observations with known
  class.

### 9.2 No-tree root

If the root MDSA scan finds no admissible tree:
- Returns valid `cta_ort` / `cta_tree` with 1 terminal stratum.
- `strata` has 1 row, stop_reason = `"no_tree"`, n = N.
- `predict()` returns majority class for all observations.
- `plot()` renders "no tree (leaf only)" without error.

### 9.3 `min_n` guard

Construct data where one depth-1 endpoint has n < `min_n`. Assert:
- That endpoint is terminal, stop_reason = `"min_n"`, no sub-fit attempted.
- Other endpoint recurses normally.

### 9.4 `max_depth` guard

Set `max_depth = 1`. Assert all depth-1 endpoints are forced terminal
regardless of whether a sub-CTA would find a tree.

### 9.5 Strata/predict consistency (all synthetic tests)

For every fitted `cta_ort` in the test suite:
`sum(strata$n) == N_train` and per-stratum predict counts match `strata$n`.

### 9.6 Private stress test (manual, data private)

Private scratch data is untracked and must NOT be committed.
The private rare-event workflow is a manual validation, not an automated fixture.
If a synthetic stand-in can reproduce the same four-stratum structure it may
be added as a smoke-tier fixture.

---

## 10. Out of Scope

- **Multiclass CTA:** ORT v1 targets binary class only. Multiclass ORT is a
  future extension contingent on multiclass CTA.
- **Global optimal tree search:** The recursion is greedy (locally optimal at
  each node by MDSA canon). Global search is not canon.
- **CART / random forest equivalence:** Not claimed, not intended.
- **Graphics v3 / ggplot renderer:** Deferred.
- **Automated MDSA at every level with user-visible family tables per node:**
  The per-node MDSA scan is internal. Family tables for individual nodes are
  accessible via `$ort_nodes[[i]]$model` but not printed by default.
- **Weighted ORT canon parity:** Case weights flow through to each sub-fit but
  weighted ORT is not tested against CTA.exe canon in v1.
- **`min_target_n` / `min_class_n` guard (backlog):** In rare-event clinical
  classification, total `min_n` is insufficient to replicate conservative
  manual stopping rules. A future guard like `min_class_n = NULL` would
  terminate recursion when any terminal child stratum contains fewer than
  `min_class_n` cases of the target class. Not in v1; add when a use case
  requires it.

---

## Open Questions (resolved)

| # | Decision |
|---|----------|
| MINDENOM policy | MDSA family scan at every node; min-D member selected; no user MINDENOM param for recursive mode |
| Seed | Single `mc_seed`, same value at all levels |
| Traversal order | Right then Left; full attribute set passed at each level |
| Entry point | `recursive = TRUE` argument in `cta_fit()` — no separate wrapper |
| S3 class | Dual-tag: `c("cta_ort", "cta_tree")` — backward compatible with all existing `cta_tree` methods |
| Traversal data stored | `level_mindenom`, `level_ess`, `level_d`, `path_ess` per node |
| Strata sort order | Ascending `prop_class1` (Stage 1 = lowest risk) |
| File placement | New file `R/cta_ort.R` |
| Attribute reuse | Full set passed; binary attributes self-eliminate naturally |
| `max_nodes` guard | `max_nodes = 31L`; counts total nodes allocated; terminal immediately if exceeded |
