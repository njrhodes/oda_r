# CTA Canon

This file defines canonical `cta_fit()` / `oda_cta_fit()` behavior for odacore.

CTA chains UniODA models across recursively defined sample strata. Each split node is a UniODA model.

## Core Definitions

- CTA = Classification Tree Analysis.
- HO-CTA = Hierarchically Optimal CTA.
- EO-CTA = Enumerated Optimal CTA.
- ESS = Effect Strength for Sensitivity.
- WESS = weighted ESS.
- If weights are declared, tree selection uses WESS.
- Reported confusion matrices in gold fixtures are raw unit-count confusion matrices unless explicitly stated otherwise.

## Degenerate Solutions

CTA never produces degenerate trees. A CTA split is **degenerate** if its rule maps all observations at the node to a single predicted class. Such a split is invalid and must be rejected at the node-level candidate gate — not post-hoc.

Rules:

- A candidate attribute that yields only one predicted class label across the node is **ineligible**, regardless of its ESS/WESS value.
- A CTA tree in which all terminal endpoints predict the same class is not a valid tree. It signals that no admissible split existed; the correct result is `no_tree`.
- `degen = TRUE` does not exist for `cta_fit()` or `oda_cta_fit()`. It is a UniODA/MultiODA-only option.
- ORT (`cta_fit(recursive = TRUE)`) follows the same rule at every recursive MDSA level.

Do not treat an all-same-class CTA or ORT result as canonical output. It is a model failure that should surface as `no_tree` or a rejected candidate, not as a valid tree with ESS > 0.

---

## Node-Level Construction

At each candidate split node:

1. Evaluate candidate attributes using UniODA.
2. A candidate attribute is eligible only if it passes:
   - MC alpha rule
   - LOO rule, if requested
   - minimum denominator / endpoint constraints
   - non-degeneracy: predicted labels must cover both classes (binary) or all C classes (multiclass)
3. Compute local ESS/WESS for eligible candidates.
4. Select according to the active CTA mode.

LOO and alpha are branch/node-level gates. They are not applied only to the final tree.

## HO-CTA

HO-CTA is greedy hierarchical construction.

At each node:

1. Screen all attributes.
2. Keep only candidates passing MC alpha + LOO.
3. Select the candidate with highest local ESS/WESS.
4. Split observations.
5. Recurse on child branches.
6. Stop when no eligible split remains or depth / denominator constraints stop growth.

HO-CTA may later be pruned.

## EO-CTA / ENUMERATE

Plain `ENUMERATE` is not root-only enumeration.

Canonical `ENUMERATE` evaluates valid combinations in the top three nodes:

```text
root
left child of root
right child of root
```

Operationally:

```text
best_tree = NULL
best_wess = -Inf

for each valid root candidate A:
  split by A

  for each valid left-child candidate B, including leaf option if no valid split:
    for each valid right-child candidate C, including leaf option if no valid split:
      construct candidate tree using A/B/C
      grow remaining allowed portions according to CTA rules
      apply pruning
      compute full-tree ESS/WESS
      if full-tree objective improves best_so_far:
        retain this candidate tree
```

Final output is the retained best tree.

Do not implement plain `ENUMERATE` as:

* greedy root selection
* root-only ranking
* first passing candidate
* missing-adjusted proxy ranking
* “effective WESS” root proxy
* final-tree-only LOO filtering

## ENUMERATE ROOT

`ENUMERATE ROOT` is distinct from plain `ENUMERATE`.

`ENUMERATE ROOT` enumerates root candidates only. Deeper nodes follow HO-CTA logic.

Do not conflate `ENUMERATE` and `ENUMERATE ROOT`.

## Pruning

Canonical CTA pruning includes:

1. sequentially rejective Sidak-Bonferroni pruning
2. maximum-accuracy pruning

`PRUNE 0.05` means pruning must be applied using alpha 0.05.

Tree growth followed by no pruning is not canonical when `PRUNE` is specified.

If a tree is structurally too deep and has higher training ESS than the gold fixture, suspect missing or incorrect pruning before suspecting RNG.

## Weights

If weights exist:

* node and tree objective uses WESS
* final objective uses weighted confusion / weighted class sensitivities
* raw confusion is still reported as raw unit counts for gold fixture parity

## Missing Values

Missing values include:

* `NA`
* declared missing codes, for example `-9`

Candidate split search uses observations non-missing for that candidate attribute.

CTA retains observations when possible. An observation missing an attribute needed for its path does not proceed down that branch and is classified according to the applicable node/branch rule.

The `OBS` shown for a node may refer to the number of observations non-missing for that node’s attribute, not necessarily the full sample size reaching that logical tree position.

## LOO

LOO applies to every attribute in the tree.

Canonical syntax:

```text
LOO {pvalue | STABLE};
```

Rules:

* `LOO STABLE`: allow only attributes whose LOO ESS/WESS equals training ESS/WESS
* `LOO pvalue`: allow only attributes whose LOO p-value is less than or equal to `pvalue`

LOO is a filter/gate, never the objective.

## Monte Carlo

Canonical command shape:

```text
MC ITER <max_iter> CUTOFF <alpha> STOP <confidence>;
```

Rules:

* `ITER` is maximum iterations
* `CUTOFF` is the alpha boundary
* `STOP` is the confidence threshold
* stop early once the alpha decision is sufficiently confident
* do not force all iterations if STOP has already resolved the decision

## Public Wrapper Policy

Public users should generally call:

```r
cta_fit(X, y, ...)
```

Current package policy:

* `cta_fit()` supports binary class variables only
* if `y` has more than two classes, fail clearly:

```text
cta_fit currently supports binary class variables only
```

Internal engine functions may exist but should not be the primary user interface.

## Current Gold Fixtures

### Myeloma

Command includes:

```text
ENUMERATE;
LOO STABLE;
PRUNE 0.05;
WEIGHT V2;
```

Expected selected tree:

```text
root = V14
node 2 = V15
raw confusion = [[146, 40], [36, 33]]
WESS = 27.69
```

### CTA_DEMO

Command includes:

```text
MC ITER 25000 CUTOFF 0.05 STOP 99.9;
MINDENOM 1;
PRUNE 0.05;
ENUMERATE;
LOO STABLE;
```

This fixture is used to validate EO-CTA structure and pruning behavior.

**Current status:** Covered CTA_DEMO fixture tests currently pass. Preserve
the documented gold expectations; do not change expected values without
rerunning the matching CTA.exe fixture.

## Non-Canonical Future Extensions

Extensions are allowed only behind explicit new options.

Examples:

* deeper enumeration
* beam search
* global search
* parallel enumeration
* robustness-weighted selection

Do not alter canonical CTA behavior to implement research extensions.

## Performance note:
CTA.exe is expected to be much faster than the R reference implementation. However, large slowdowns usually indicate repeated MC/LOO fits. Canonical behavior should be implemented first; later optimization should use memoization of node × attribute fits and early MC stopping, without changing results.
