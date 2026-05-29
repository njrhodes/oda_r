# FIT_OBJECT_EVIDENCE_CONTRACT.md

Evidence fields carried by each fitted-object class and the public methods that
consume them.  This document is the Slice S audit deliverable.

**Invariants enforced here:**
- No refitting in reporting, graphics, or bootstrap methods.
- No fabricated evidence (values not derivable from stored fit fields).
- Evidence gaps are documented; absent evidence is exposed as `NA`, not silently
  dropped or invented.

---

## 1. Fitted-object evidence map

### 1.1 `oda_fit` (binary: `oda_fit_binary`)

| Field | Type | Notes |
|-------|------|-------|
| `ok` | logical | Fit succeeded |
| `reason` | character | Reason when `ok = FALSE` |
| `rule` | list | Rule object (`type`, `cut_value`, `direction`, `label_0`, `label_1`) |
| `ess` | numeric | ESS (%) on training data |
| `pac` | numeric | Priors-weighted PAC (%) |
| `p_mc` | numeric or NULL | Monte Carlo p-value |
| `loo` | list or NULL | LOO results (`allowed`, `confusion`, `ess_loo`, `p_value`) |
| `n_eff` | integer | Effective N (non-missing, non-excluded) |
| `attr_type` | character | `"ordered"`, `"binary"`, or `"categorical"` |
| `engine` | character | `"binary"` |
| `has_weights` | logical | Any weight ≠ 1 |
| `priors_on` | logical | Priors weighting active |
| `miss_codes` | numeric or NULL | Miss-code values |
| `confusion` | list | `TP`, `TN`, `FP`, `FN`, `sensitivity`, `specificity`, `mean_pac` |
| `confusion_wt` | list | Same fields, weighted |

**2×2 confusion matrix layout** (actual × predicted, 0-indexed classes):

```
         pred=0  pred=1
actual=0  [TN]    [FP]
actual=1  [FN]    [TP]
```

### 1.2 `oda_fit` (multiclass: `oda_fit_multiclass`)

| Field | Type | Notes |
|-------|------|-------|
| `ok` | logical | |
| `rule` | list | `type`, `cut_values`, `seg_classes` for ordered; `cat_map` for nominal |
| `ess` | numeric | Mean PAC-based ESS |
| `mean_pac` | numeric | |
| `pac_by_class` | numeric vector | Per-class PAC |
| `p_mc` | numeric or NULL | |
| `loo` | list or NULL | |
| `n_eff` | integer | |
| `attr_type` | character | |
| `engine` | character | `"multiclass"` |
| `has_weights` | logical | |
| `priors_on` | logical | |
| `miss_codes` | numeric or NULL | |
| `confusion` | matrix | C×C integer, raw counts |
| `confusion_wt` | matrix | C×C weighted counts |

**Note:** NOVOboot (`novo_boot_ci`) requires a 2×2 confusion.  It is not
applicable to multiclass ODA.  `novo_boot_ci.oda_fit` errors with a clear
message if `inherits(fit, "oda_fit_binary")` is FALSE.

### 1.3 `cta_tree`

| Field | Type | Notes |
|-------|------|-------|
| `nodes` | named list | Per-node: `node_id`, `depth`, `leaf`, `attribute`, `n_obs`, `majority_class`, `child_ids`, `ess`, `loo_ess`, `loo_p`, `loo_status`, `p_mc`, `class_counts_raw`, `class_counts_weighted` |
| `root_id` | integer | |
| `no_tree` | logical | TRUE when no admissible split |
| `training_confusion` | matrix or NULL | C×C integer, actual × predicted; NULL for no_tree |
| `overall_ess` | numeric | WESS when weights active, ESS otherwise; NA for no_tree |
| `has_weights` | logical | |
| `n` | integer | Training N |
| `C` | integer | Number of classes |
| `attr_names` | character | Attribute names |
| `miss_codes` | numeric or NULL | |
| `priors_on` | logical | |
| `alpha_split` | numeric | |
| `mindenom` | integer | |
| `prune_alpha` | numeric | |
| `max_depth` | integer | |
| `loo` | character | LOO mode string |

**`training_confusion` convention:** rows = actual class (0-indexed), cols =
predicted class.  For binary CTA: `[1,1]`=TN, `[1,2]`=FP, `[2,1]`=FN,
`[2,2]`=TP.  Classified observations only (path-local missingness may exclude
some obs from the root-path).

**`overall_ess` / D-statistic:** `cta_d_stat(tree)` = `100 / (ess/s) - s` where
`s = cta_strata(tree)`.  Returns `NA_real_` for no_tree or when `ess <= 0` or
`s < 2`.

### 1.4 `cta_ort` (LORT)

Inherits `cta_tree`, adds:

| Field | Type | Notes |
|-------|------|-------|
| `ort_nodes` | named list | LORT node tree |
| `ort_root_id` | integer | |
| `strata` | data.frame | Per-terminal: `stratum_id`, `node_id`, `terminal_class`, `class_counts` (list col of named int vector), `prop_class1`, `n`, `n_weighted` |
| `n_strata` | integer | |
| `ort_settings` | list | `method="lort"`, `global_optimization=FALSE`, `sda_anchored=FALSE` |
| `recursive` | logical | Always TRUE |

**`training_confusion` for LORT:** The `training_confusion` on a `cta_ort` is
the ROOT model confusion only, not the full-LORT confusion.  For NOVOboot
purposes the correct confusion for LORT is derived from `$strata$class_counts`
by summing over all terminal strata.  See §3.3 below.

---

## 2. Public methods and their evidence consumption

### 2.1 Accessors — confirmed correct, no changes needed

| Method / Function | Object | Evidence consumed | Status |
|---|---|---|---|
| `oda_confusion(fit)` | `oda_fit` | `$confusion`, `$confusion_wt`, `$loo$confusion` | ✓ correct |
| `oda_metrics(fit)` | `oda_fit` | `$ess`, `$pac`, `$confusion`, `$p_mc`, `$loo` | ✓ correct |
| `oda_d_stat(fit)` | `oda_fit` | `$ess`, `$rule$type`, strata count | ✓ correct |
| `summary.oda_fit()` | `oda_fit` | `$ess`, `$confusion`, `$p_mc`, `$loo`, `$rule` | ✓ correct |
| `print.oda_fit()` | `oda_fit` | delegates to summary | ✓ correct |
| `summary.cta_tree()` | `cta_tree` | `$overall_ess`, `$training_confusion`, `$no_tree` | ✓ correct |
| `cta_confusion_table()` | `cta_tree` | `$training_confusion` (stored) | ✓ correct |
| `cta_d_stat()` | `cta_tree` | `$overall_ess`, `cta_strata()` | ✓ correct |
| `cta_family_table()` | `cta_family` | `$overall_ess`, `$d` per member | ✓ correct |

### 2.2 Plot-data extractors — gaps identified and fixed

| Function | Object | Was missing | Fix |
|---|---|---|---|
| `cta_plot_data()` | `cta_tree` | `overall_ess`, `ess_label`, `d`, `model_label`, `training_n` | Added in Slice S |
| `ort_plot_data()` | `cta_ort` | `overall_ess`, `ess_label`, `d`, `model_label`, `training_n` | Added in Slice S |

### 2.3 Graphics renderers — `show_metrics` added

| Function | Parameter added | Behavior |
|---|---|---|
| `plot_cta_tree()` | `show_metrics = FALSE` | When TRUE: appends `"ESS: xx.xx% | D: x.xx"` subtitle |
| `plot_lort_tree()` | `show_metrics = FALSE` | Same; uses LORT evidence from `ort_plot_data()` |

### 2.4 NOVOboot — S3 dispatch added

| Method | Evidence path | Notes |
|---|---|---|
| `novo_boot_ci.default(x, ...)` | `x` = 2×2 matrix (existing behavior) | |
| `novo_boot_ci.oda_fit(fit, ...)` | `fit$confusion` → build 2×2 | Binary only; errors for multiclass |
| `novo_boot_ci.cta_tree(tree, ...)` | `tree$training_confusion` | Must be 2×2; errors for no_tree or multiclass |
| `novo_boot_ci.cta_ort(ort, ...)` | Sum `ort$strata$class_counts` by `terminal_class` | Full-LORT confusion; errors when strata absent |

---

## 3. Implementation notes

### 3.1 `novo_boot_ci.oda_fit`

```r
conf <- fit$confusion       # list: TP, TN, FP, FN
m <- matrix(c(conf$TN, conf$FP, conf$FN, conf$TP),
            nrow = 2L, byrow = TRUE)
novo_boot_ci.default(m, ...)
```

Requires `inherits(fit, "oda_fit_binary")`.

### 3.2 `novo_boot_ci.cta_tree`

```r
m <- tree$training_confusion    # stored 2×2 integer matrix
novo_boot_ci.default(m, ...)
```

Errors when `no_tree = TRUE` (NULL confusion) or when confusion is not 2×2.

### 3.3 `novo_boot_ci.cta_ort`

Accumulate from `strata` df:

```r
st <- ort$strata
conf <- matrix(0L, 2L, 2L)
for each row in st:
  p <- terminal_class + 1L          # predicted col (1-indexed)
  cc <- class_counts                # named int vector, e.g. c("0"=18, "1"=2)
  for class_name in names(cc):
    a <- as.integer(class_name) + 1L  # actual row (1-indexed)
    conf[a, p] <- conf[a, p] + cc[class_name]
novo_boot_ci.default(conf, ...)
```

Errors when `strata` is absent, empty, or coverage < 2 obs per cell.

### 3.4 Evidence fields in `cta_plot_data()` / `ort_plot_data()`

Added to the returned list (NA when no_tree or value unavailable):

| Field | Source | No-tree value |
|---|---|---|
| `overall_ess` | `tree$overall_ess` | `NA_real_` |
| `ess_label` | `"WESS"` if `has_weights`, else `"ESS"` | `NA_character_` |
| `d` | `cta_d_stat(tree)` | `NA_real_` |
| `model_label` | `"CTA"` or `"LORT"` | `"CTA"` / `"LORT"` |
| `training_n` | `tree$n` | `NA_integer_` |

For `ort_plot_data()`, `model_label = "LORT"` always.

### 3.5 `show_metrics` subtitle construction

When `show_metrics = TRUE` and `pd$overall_ess` is not NA:

```
metrics_line <- sprintf("%s: %.2f%%  |  D: %s",
                        pd$ess_label,
                        pd$overall_ess,
                        if (is.finite(pd$d)) sprintf("%.4f", pd$d) else "NA")
```

If a `subtitle` argument was also supplied by the caller:
```
effective_subtitle <- paste0(subtitle, "\n", metrics_line)
```
Otherwise `effective_subtitle <- metrics_line`.

---

## 4. Deferred gaps

| Gap | Reason deferred |
|-----|----------------|
| NOVOboot for multiclass ODA | NOVOboot requires 2×2 confusion; no correct extension exists for C>2 without redesign |
| LOO-confusion NOVOboot | LOO confusion stored per-fold only; aggregation semantics not defined |
| `show_confusion` renderer option | Out of scope for Slice S; add in Phase 2A if requested |
| `model_evidence()` public generic | Adds abstraction for one-time use; not justified |
| `.evidence_*()` internal helpers | Same; adds indirection without benefit |

---

## 5. Tests

Test file: `tests/testthat/test-fit-evidence.R`

Covers:
- `novo_boot_ci.default` — matrix path unchanged
- `novo_boot_ci.oda_fit` — binary only; errors on multiclass
- `novo_boot_ci.cta_tree` — correct confusion; errors on no_tree
- `novo_boot_ci.cta_ort` — full-LORT confusion from strata
- `cta_plot_data()` evidence fields present, correct type, NA for no_tree
- `ort_plot_data()` evidence fields present, correct type
- `plot_cta_tree(show_metrics = TRUE)` — subtitle contains ESS string
- `plot_lort_tree(show_metrics = TRUE)` — subtitle contains ESS string
- No-fabrication: evidence fields absent from `oda_fit` multiclass `novo_boot_ci` call
- No-refitting: `novo_boot_ci` methods do not call `oda_fit`, `cta_fit`, `lort_fit`
- Reserved names: `sda_propensity_weights`, `sort_propensity_weights` do not exist
