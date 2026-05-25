# ODA Canon

This file defines canonical `oda_fit()` behavior for odacore.

## Core Definitions

- ODA = Optimal Data Analysis.
- UniODA is the atomic model used by CTA.
- ESS = Effect Strength for Sensitivity.
- WESS = weighted ESS.
- If weights are declared, model selection uses WESS.
- Reported confusion matrices may be raw counts or weighted counts, but the objective must be explicit.

## Objective

For binary class problems, UniODA searches for the rule that maximizes ESS/WESS.

For multiclass problems, UniODA searches for the rule set that maximizes mean class sensitivity, then transforms it to ESS/WESS.

ESS is computed from Mean PAC:

```text
ESS = (MeanPAC - chancePAC) / (100 - chancePAC) * 100
```

For two classes:

```text
ESS = (MeanPAC - 50) / 50 * 100
```

## Ordered Attributes

For ordered attributes:

1. Sort observations by attribute value.
2. Identify valid cutpoints.
3. Evaluate valid direction(s).
4. Compute class sensitivities.
5. Select the cutpoint + direction maximizing ESS/WESS.

For binary numeric attributes `{0,1}`, canonical CTA display may normalize the rule to:

```text
<= 0.5 -> class 0
>  0.5 -> class 1
```

when prediction behavior is equivalent.

## Categorical Attributes

For categorical attributes:

1. Evaluate valid category assignments.
2. Compute class sensitivities.
3. Select assignment maximizing ESS/WESS.

## Weights

If weights are supplied:

* training objective uses weighted class sensitivities
* WESS is reported as the objective metric
* raw confusion may still be used for fixture comparison and human-readable output

Do not mix raw confusion and weighted objective.

## Missing Values

Missing values include:

* `NA`
* declared missing codes, for example `-9`

For a candidate attribute, missing observations are excluded from that attribute’s UniODA split search.

CTA may retain observations with missing data if their path through the tree does not require the missing attribute.

## Monte Carlo

Canonical command shape:

```text
MC ITER <max_iter> CUTOFF <alpha> STOP <confidence>;
```

Rules:

* `ITER` is a maximum iteration count, not a required count.
* `CUTOFF` is the p-value decision boundary.
* `STOP` is the confidence threshold for early stopping.
* MC stops when any stopping criterion is reached:

  * max iterations reached
  * confidence is high enough that `p <= CUTOFF`
  * confidence is high enough that `p > CUTOFF`

Do not force all MC iterations when the STOP criterion already resolves the alpha decision.

## LOO

Canonical syntax:

```text
LOO {pvalue | STABLE};
```

LOO is a filter, not the optimization objective.

Modes:

* `loo = "off"`: no LOO filter
* `loo = "stable"`: allow only attributes whose LOO ESS/WESS equals training ESS/WESS
* `loo = numeric`: allow only attributes whose LOO p-value is less than or equal to the supplied threshold

For CTA, LOO applies to every attribute used in the tree, not only to the final tree confusion.

## DIRECTION (MPE Chapter 2 binary ordered — Phase 6A/6B)

Canonical MegaODA Appendix A command shape:

```text
DIRECTION < 0 1   -- high values predict class 1 ("greater" / "0->1")
DIRECTION > 0 1   -- low values predict class 1  ("less"  / "1->0")
```

R `direction` parameter values and mapping for binary ordered/binary:

| `direction`  | Meaning                                              | Internal      |
|--------------|------------------------------------------------------|---------------|
| `"both"`     | Non-directional (canonical default; both sides)      | `"off"`       |
| `"off"`      | Backward-compatible synonym for `"both"`             | `"off"`       |
| `"greater"`  | High attribute values predict class 1 (Chapter 2 >) | `"0->1"`      |
| `"less"`     | Low attribute values predict class 1 (Chapter 2 <)  | `"1->0"`      |

Implementation scope (Phase 6A/6B — complete):
- Candidate filter: only directions matching the constraint are evaluated.
- MC permutation: directional constraint applied identically in each permutation refit.
- LOO fold refits: directional constraint applied identically in each fold refit.
- LOO Fisher `alternative`: `"two.sided"` for `"both"`/`"off"`; `"greater"` for directional.

Binary class errors: `"ascending"` and `"descending"` are rejected for binary class
(Chapter 4 semantics require multiclass). Use `"greater"` or `"less"` instead.

Categorical attributes + Chapter 2 direction: `direction` in `{"greater","less"}` with
categorical attribute returns `ok = FALSE` with `reason = "direction_not_supported_for_categorical"`.
Use `direction_map` for categorical directional hypotheses.

## DIRECTION (MPE Chapter 4 categorical/table DIRECTIONAL — Phase 6C)

Two additional direction values for multiclass and categorical identity-map DIRECTIONAL:

| `direction`     | Meaning                                                    | Scope              |
|-----------------|------------------------------------------------------------|--------------------|
| `"ascending"`   | Segment/level s → class s (identity / ascending order)    | multiclass ordered + categorical k=C |
| `"descending"`  | Segment/level s → class C+1-s (reverse order)             | multiclass ordered + categorical k=C |

**Multiclass ordered (`direction = "ascending"` or `"descending"`):**
- Constrains the segment-to-class assignment in `oda_best_ordered_multiclass_partition()`.
- For each cut position, only the ascending (s→s) or descending (s→C+1-s) assignment
  is evaluated; no other assignments are considered.
- MC permutation: same direction constraint applied per permutation refit.
- Covers MPE Chapter 4 stability analysis (Bowker Table 4.1 DIRECTIONAL < 1 2 3 4).

**Multiclass categorical, k = C (`direction = "ascending"`):**
- When L == C, auto-creates identity `direction_map` = setNames(1..C, attribute_levels).
- When L != C and no `direction_map` supplied: error (must supply `direction_map`).
- Covers MPE Chapter 4 political affiliation and protein convergent validity.

**`direction_map` parameter (categorical fixed-partition DIRECTIONAL):**
- Named integer vector: names = attribute levels (character), values = predicted class labels.
- For binary class: values 0/1 in coded space (recoded by `oda_fit()` from original labels).
- For multiclass: values 1..C in internal encoded space.
- All attribute levels must be covered exactly once; at least two distinct classes.
- Bypasses the partition search entirely; evaluates only the specified mapping.
- MC permutation: each permutation evaluates the SAME fixed mapping on permuted y labels.
- LOO: with fixed categorical rule, predictions are determined by x alone → trivially stable.
- Covers MPE Chapter 4 gully erosion example (fixed partition, k ≠ C).

**Compatible combinations:**
- `direction = "both"` (default) + `direction_map` = nondirectional fit with a fixed partition.
- `direction = "ascending"` + categorical + L==C = auto identity map (Chapter 4 ascending).
- `direction = "ascending"` + categorical + L≠C + `direction_map` = custom fixed partition.
- `direction = "ascending"` + ordered + multiclass = constrained segment-to-class (Chapter 4).
- `direction = "greater"/"less"` + binary ordered = Chapter 2 directional (existing).

## Public Wrapper Policy

Public users should generally call:

```r
oda_fit(x, y, ...)
```

`oda_fit()` should auto-detect the number of classes:

* two classes -> binary UniODA path
* more than two classes -> multiclass ODA path, if implemented

Internal helper functions may exist but should not be required for normal use.
