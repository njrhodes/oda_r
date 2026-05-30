# Myeloma CTA: Translation and Reporting Walkthrough

> **Worked fixture example.** All values computed from actual `odacore` package runs
> on the myeloma fixture (`tests/testthat/fixtures/myeloma/`). Not active roadmap.

This document is a workflow walkthrough showing how the `odacore`
CTA reporting stack maps from a fitted tree to endpoint summaries, staging
tables, propensity weights, confusion table, and family-level interpretation.
All output shown below is computed from the myeloma fixture included in this
repository (`tests/testthat/fixtures/myeloma/data.txt`).

---

## 1. Purpose

The reporting stack follows a strict separation:

- **Fitting** (`oda_cta_fit`) — produces a lean `cta_tree` object that stores
  only the tree topology, per-leaf majority class, raw classification
  confusion (`training_confusion`), and per-leaf class-count vectors
  (`class_counts_raw`, `class_counts_weighted`). No tables, no row indices,
  no staging or propensity objects are stored at fit time.
- **Reporting** (`cta_endpoint_summary`, `cta_endpoint_counts`,
  `cta_staging_table`, `cta_propensity_weights`, `cta_confusion_table`) —
  pure on-demand views over the stored leaf counts. No refitting, no
  prediction, and no training-data recomputation.
- **Translation** (`cta_assign_endpoints`, `cta_observation_weights`) —
  on-demand traversal and weight assignment requiring caller-supplied
  `newdata` and class labels.

This walkthrough covers the myeloma MINDENOM family: MINDENOM=1 (full
two-split tree), MINDENOM=30 (stump), and MINDENOM=56 (no-tree terminal).

---

## 2. Fixture Setup

```r
# Data path is relative to repo root
raw <- read.table("tests/testthat/fixtures/myeloma/data.txt",
                  header = FALSE, sep = "")
names(raw) <- paste0("V", seq_len(ncol(raw)))
df  <- raw[raw$V2 != 0, ]          # EX V2=0 — exclude zero-weight rows

attrs <- c("V4","V9","V11","V12","V14","V15","V16","V17","V18","V19")
y     <- as.integer(df$V1)
w     <- df$V2
```

```r
fit_myeloma <- function(mindenom) {
  cta_fit(
    X           = df[, attrs],
    y           = y,
    w           = w,
    mindenom    = mindenom,
    priors_on   = TRUE,
    miss_codes  = -9,
    alpha_split = 0.05,
    prune_alpha = 0.05,
    max_depth   = 20L,
    mc_iter     = 5000L,
    mc_target   = 0.05,
    mc_stop     = 99.9,
    mc_stopup   = 99.9,
    loo         = "stable",
    attr_names  = attrs
  )
}

t1  <- fit_myeloma(1L)
t30 <- fit_myeloma(30L)
t56 <- fit_myeloma(56L)
```

Key fixture characteristics:
- Binary class: 0 = alive, 1 = deceased.
- Case weights active (`w = V2`); scores are WESS, not ESS.
- `miss_codes = -9` — value −9 is treated as missing on any attribute.
  Path-local missingness is canonical: an observation is excluded (returned
  as `NA`) when its split attribute is missing on its actual traversal path.
- 255 classified observations at MINDENOM=1 (69 obs with V17 = −9 are
  excluded at MINDENOM=30/56 root, which splits on V17).

---

## 3. Endpoint Structure

### MINDENOM=1 — three-endpoint tree (V14 → V15)

```r
cta_endpoint_summary(t1)
#   endpoint_id endpoint_node_id                  path depth terminal_prediction
# 1           1                3 V14<=0.5 AND V15<=0.5     3                   0
# 2           2                4  V14<=0.5 AND V15>0.5     3                   1
# 3           3                5               V14>0.5     2                   1
#   n_obs n_weighted denominator
# 1   182    6477.61         182
# 2    29     875.97          29
# 3    44    1461.48          44
```

`endpoint_id` is sequential by node order. `depth` reflects the number of
splits traversed to reach the leaf. `n_obs` is the raw observation count;
`n_weighted` is the sum of case weights at the leaf.

```r
cta_endpoint_counts(t1)
#   endpoint_id endpoint_node_id                  path terminal_prediction class
# 1           1                3 V14<=0.5 AND V15<=0.5                   0     0
# 2           1                3 V14<=0.5 AND V15<=0.5                   0     1
# 3           2                4  V14<=0.5 AND V15>0.5                   1     0
# 4           2                4  V14<=0.5 AND V15>0.5                   1     1
# 5           3                5               V14>0.5                   1     0
# 6           3                5               V14>0.5                   1     1
#   n_raw n_weighted
# 1   146    5704.71
# 2    36     772.90
# 3    14     552.09
# 4    15     323.88
# 5    26    1030.64
# 6    18     430.84
```

One row per endpoint × actual class. The per-leaf `class_counts_raw` vectors
are the source; no training-data indices are stored or accessed.

### MINDENOM=30 — stump (V17 only)

```r
cta_endpoint_summary(t30)
#   endpoint_id endpoint_node_id     path depth terminal_prediction n_obs
# 1           1                2 V17<=0.5     2                   0   131
# 2           2                3  V17>0.5     2                   1    55
#   n_weighted denominator

cta_endpoint_counts(t30)
#   endpoint_id endpoint_node_id     path terminal_prediction class n_raw n_weighted
# 1           1                2 V17<=0.5                   0     0   101    3824.40
# 2           1                2 V17<=0.5                   0     1    30     624.36
# 3           2                3  V17>0.5                   1     0    34    1226.83
# 4           2                3  V17>0.5                   1     1    21     430.31
```

Classified n = 186 (69 obs with V17 = −9 are path-locally excluded;
255 − 186 = 69).

### MINDENOM=56 — no-tree terminal

```r
isTRUE(t56$no_tree)  # TRUE

nrow(cta_endpoint_summary(t56))    # 0
nrow(cta_endpoint_counts(t56))     # 0
```

No-tree fits return zero-row data frames with correct column structure from
all reporting functions. This is the canonical termination state.

---

## 4. Staging Table

The staging table orders endpoints from lowest to highest target-class
propensity (ascending). Stage 1 is the lowest-risk stratum; the highest
stage is the highest-risk stratum.

### MINDENOM=1

```r
cta_staging_table(t1, target_class = 1L)[,
  c("stage","path","terminal_prediction","target_n",
    "denominator","target_proportion","perfectly_predicted")]
#   stage                  path terminal_prediction target_n denominator
# 1     1 V14<=0.5 AND V15<=0.5                   0       36         182
# 2     2               V14>0.5                   1       18          44
# 3     3  V14<=0.5 AND V15>0.5                   1       15          29
#   target_proportion perfectly_predicted
# 1         0.1978022               FALSE
# 2         0.4090909               FALSE
# 3         0.5172414               FALSE
```

Interpretation:
- Stage 1 (`V14<=0.5 AND V15<=0.5`, predicted class 0): 36 of 182 obs (19.8%)
  are class 1 (deceased). Lowest target-class propensity.
- Stage 2 (`V14>0.5`, predicted class 1): 18 of 44 obs (40.9%).
- Stage 3 (`V14<=0.5 AND V15>0.5`, predicted class 1): 15 of 29 obs (51.7%).
  Highest target-class propensity.

No endpoint is perfectly predicted at MINDENOM=1 (`perfectly_predicted =
FALSE` for all), so the `adjusted` flag is `FALSE` and adjusted values equal
empirical values.

### MINDENOM=30

```r
cta_staging_table(t30, target_class = 1L)[,
  c("stage","path","terminal_prediction","target_n",
    "denominator","target_proportion")]
#   stage     path terminal_prediction target_n denominator target_proportion
# 1     1 V17<=0.5                   0       30         131         0.2290076
# 2     2  V17>0.5                   1       21          55         0.3818182
```

Stump: two endpoints. The V17<=0.5 endpoint predicts class 0 but 30/131
(22.9%) are actually class 1 (deceased). The V17>0.5 endpoint predicts class
1 with 21/55 (38.2%) actually class 1.

**Perfect endpoints and the one-hypothetical-observation adjustment:**

When an endpoint contains only one class (`target_n == 0` or
`non_target_n == 0`), the empirical odds and target proportion cannot be
ordered. `cta_staging_table()` applies the Yarnold-Linden (2017) remedy by
default (`adjust_perfect = TRUE`): one hypothetical misclassified observation
is added to the absent-class profile at that endpoint, and to the global
marginal totals, so all endpoints can be ordered and compared. The `adjusted`
and `perfectly_predicted` columns flag which endpoints received the
correction. The myeloma fixtures at MINDENOM=1 and MINDENOM=30 do not contain
perfectly predicted endpoints.

---

## 5. Propensity Weights

Stabilized propensity-style weights at the endpoint × actual-class level,
following Yarnold and Linden (2017):

```
weight_{s,z} = n_s × Pr(Z=z) / n_{s,z}
```

where `n_s` is the endpoint denominator, `Pr(Z=z)` is the marginal class
probability across all classified observations, and `n_{s,z}` is the raw
count of actual class `z` at endpoint `s`.

**Scope:** `cta_propensity_weights()` returns endpoint × actual-class weights —
one row per endpoint × class combination. `cta_assign_endpoints()` maps
caller-supplied rows to endpoint IDs on demand. `cta_observation_weights()`
joins endpoint assignment to the endpoint × class weights and returns one row
per observation. No endpoint membership is stored at fit time; assignment and
weight lookup are performed on explicit function call and remain lean-fit
compatible.

The endpoint-level weights (this section) are the building block. The
observation-level join (`cta_observation_weights()`) is the downstream step.

### MINDENOM=1

```r
cta_propensity_weights(t1, target_class = 1L)[,
  c("endpoint_id","class","class_n","endpoint_n",
    "marginal_class_probability","propensity_weight",
    "undefined_empirical","adjusted_propensity_weight")]
#   endpoint_id class class_n endpoint_n marginal_class_probability
# 1           1     0     146        182                  0.7294118
# 2           1     1      36        182                  0.2705882
# 3           2     0      14         29                  0.7294118
# 4           2     1      15         29                  0.2705882
# 5           3     0      26         44                  0.7294118
# 6           3     1      18         44                  0.2705882
#   propensity_weight undefined_empirical adjusted_propensity_weight
# 1         0.9092667               FALSE                  0.9092667
# 2         1.3679739               FALSE                  1.3679739
# 3         1.5109244               FALSE                  1.5109244
# 4         0.5231373               FALSE                  0.5231373
# 5         1.2343891               FALSE                  1.2343891
# 6         0.6614379               FALSE                  0.6614379
```

Marginal class probabilities: class 0 = 186/255 ≈ 0.729; class 1 = 69/255 ≈
0.271. These are computed from the pooled leaf counts and are the same for
every row within a class.

**Algebraic identity:** Within each endpoint, the weighted mean of
`propensity_weight` (weights = `class_n`, divided by `endpoint_n`) equals
exactly 1.0. For endpoint 1: `(146 × 0.9093 + 36 × 1.3680) / 182 = 1.0`.

When `undefined_empirical = TRUE` (absent class at a perfect endpoint), the
empirical `propensity_weight = Inf`. With `adjusted = TRUE` (default), the
one-hypothetical-observation remedy makes all adjusted weights finite. The
myeloma MINDENOM=1 and MINDENOM=30 fixtures have no absent classes.

### MINDENOM=30

```r
cta_propensity_weights(t30, target_class = 1L)[,
  c("endpoint_id","class","class_n","endpoint_n",
    "marginal_class_probability","propensity_weight")]
#   endpoint_id class class_n endpoint_n marginal_class_probability
# 1           1     0     101        131                  0.7258065
# 2           1     1      30        131                  0.2741935
# 3           2     0      34         55                  0.7258065
# 4           2     1      21         55                  0.2741935
#   propensity_weight
# 1         0.9413925
# 2         1.1973118
# 3         1.1740987
# 4         0.7181260
```

Marginal class probabilities shift because the classified sample at MINDENOM=30
is n=186 (not 255): class 0 = 135/186 ≈ 0.726; class 1 = 51/186 ≈ 0.274.

### MINDENOM=56

```r
nrow(cta_propensity_weights(t56, target_class = 1L))  # 0
```

Zero-row data frame with correct schema — no tree, no endpoints, no weights.

---

## 5.1. Endpoint Assignment

`cta_assign_endpoints()` traverses the fitted tree for each row of `newdata`
and returns the endpoint each observation reaches.  It requires only the fitted
`cta_tree` and a data frame with the model attributes — no training data is
stored or accessed.  Path-local missingness is canonical: any observation whose
split attribute is missing on its actual traversal path receives `NA` for
`endpoint_node_id` and `endpoint_id`.

### MINDENOM=1 — first 5 training rows

```r
# df[1:5, attrs] has V14 in {0, 0, 1, 1, 0} and V15 in {0, 0, 0, 0, 0}.
# t1 root is V14 (<=0.5 → node 2, >0.5 → node 5); node 2 splits on V15.
cta_assign_endpoints(t1, newdata = df[1:5, attrs])
#   row_id endpoint_node_id endpoint_id
# 1      1                3           1
# 2      2                3           1
# 3      3                5           3
# 4      4                5           3
# 5      5                3           1
```

Rows 1, 2, 5 (V14 ≤ 0.5 and V15 ≤ 0.5) land at endpoint 1 (node 3).
Rows 3 and 4 (V14 > 0.5) land at endpoint 3 (node 5).
Row 5 has V17 = −9 but that attribute is not used in the t1 tree; the
V17 missingness does not affect assignment.

### MINDENOM=30 — path-local NA when root attribute is missing

The MINDENOM=30 stump roots at V17.  Rows with V17 = −9 cannot be traversed
and receive `NA` endpoints.

```r
# First 3 rows have V17 in {1, 0, 0}; next 3 have V17 = -9.
demo_rows <- c(which(df$V17 != -9)[1:3], which(df$V17 == -9)[1:3])

cta_assign_endpoints(t30, newdata = df[demo_rows, attrs])
#   row_id endpoint_node_id endpoint_id
# 1      1                3           2
# 2      2                2           1
# 3      3                2           1
# 4      4               NA          NA
# 5      5               NA          NA
# 6      6               NA          NA
```

Rows 4–6 have V17 = −9 and are excluded with `NA` — canonical path-local
missingness.  These 69 obs with V17 = −9 are why MINDENOM=30 classifies
n = 186, not n = 255.

---

## 5.2. Observation-Level Weights

`cta_observation_weights()` calls `cta_assign_endpoints()` internally, joins
the resulting `endpoint_id` and `actual_class` to the endpoint × class weights
from `cta_propensity_weights()`, and returns one row per observation.
Observations that receive `NA` endpoint (path-local missing on root) carry
`NA` weights and `assigned = FALSE`.

### MINDENOM=1 — first 5 training rows

```r
cta_observation_weights(t1, newdata = df[1:5, attrs], y = y[1:5])[,
  c("row_id","actual_class","endpoint_id",
    "propensity_weight","adjusted_propensity_weight","assigned")]
#   row_id actual_class endpoint_id propensity_weight adjusted_propensity_weight
# 1      1            0           1         0.9092667                  0.9092667
# 2      2            0           1         0.9092667                  0.9092667
# 3      3            0           3         1.2343891                  1.2343891
# 4      4            1           3         0.6614379                  0.6614379
# 5      5            0           1         0.9092667                  0.9092667
#   assigned
# 1     TRUE
# 2     TRUE
# 3     TRUE
# 4     TRUE
# 5     TRUE
```

Each observation receives the endpoint × actual-class propensity weight for
its endpoint.  Rows 1, 2, 5 (actual class 0, endpoint 1) all receive
weight ≈ 0.909.  Row 4 (actual class 1, endpoint 3) receives weight ≈ 0.661.
`adjusted` is `FALSE` for all rows here because no endpoint is perfectly
predicted at MINDENOM=1.

### MINDENOM=30 — path-local NA propagates to weights

```r
cta_observation_weights(t30, newdata = df[demo_rows, attrs],
                         y = y[demo_rows])[,
  c("row_id","actual_class","endpoint_id","propensity_weight","assigned")]
#   row_id actual_class endpoint_id propensity_weight assigned
# 1      1            0           2         1.1740987     TRUE
# 2      2            0           1         0.9413925     TRUE
# 3      3            0           1         0.9413925     TRUE
# 4      4            0          NA                NA    FALSE
# 5      5            1          NA                NA    FALSE
# 6      6            1          NA                NA    FALSE
```

Rows 4–6 (V17 = −9) cannot be assigned to an endpoint; their weights are
`NA` and `assigned = FALSE`.

---

## 6. Confusion Table

The confusion table covers the full selected (pruned) tree applied to all
classified training observations.

### MINDENOM=1

```r
cta_confusion_table(t1)
#   actual predicted   n
# 1      0         0 146
# 2      0         1  40
# 3      1         0  36
# 4      1         1  33
```

Total classified n = 255. [[146, 40], [36, 33]] — rows are actual class,
columns are predicted class. OVERALL ESS = 26.32%; WEIGHTED ESS = 27.69%.

### MINDENOM=30

```r
cta_confusion_table(t30)
#   actual predicted   n
# 1      0         0 101
# 2      0         1  34
# 3      1         0  30
# 4      1         1  21
```

Total classified n = 186 (69 obs with V17 = −9 excluded path-locally;
255 − 186 = 69).
[[101, 34], [30, 21]]. OVERALL ESS = 15.99%; WEIGHTED ESS = 16.51%.

The confusion table is stored as `training_confusion` on the `cta_tree`
object at fit time and read back by `cta_confusion_table()`. It is the only
aggregate table cached on the fit object.

---

## 7. Family-Level Interpretation

The MDSA descendant family traces the tree from MINDENOM=1 through the
novometric stepping rule (next MINDENOM = min terminal denominator + 1),
terminating when a no-tree fit is reached.

```r
fam <- cta_descendant_family(
  X              = df[, attrs],
  y              = y,
  w              = w,
  start_mindenom = 1L,
  max_steps      = 20L,
  priors_on      = TRUE,
  miss_codes     = -9,
  alpha_split    = 0.05,
  prune_alpha    = 0.05,
  max_depth      = 20L,
  mc_iter        = 5000L,
  mc_target      = 0.05,
  mc_stop        = 99.9,
  mc_stopup      = 99.9,
  loo            = "stable",
  attr_names     = attrs
)

cta_family_table(fam)
#   index mindenom     status no_tree strata min_terminal_denom next_mindenom
# 1     1        1 valid_tree   FALSE      3                 29            30
# 2     2       30      stump   FALSE      2                 55            56
# 3     3       56    no_tree    TRUE     NA                 NA            NA
#   overall_ess has_weights        d selected_min_d
# 1    27.68636        TRUE  7.83566           TRUE
# 2    16.51269        TRUE 10.11190          FALSE
# 3          NA        TRUE       NA          FALSE
```

**Reading the family table:**

- `overall_ess`: For weighted fits this is WESS (weighted ESS); column name
  in `cta_family_table` is `overall_ess` regardless. Never call weighted
  scores generic ESS when quoting CTA output.
- `d`: D statistic = `100 / (ESS / strata) - strata`. Lower D is more
  parsimonious relative to translational strength. MINDENOM=1 (D = 7.84)
  is selected by minimum D.
- `strata`: Number of terminal endpoints. MINDENOM=1 has 3 strata (ordered
  risk groups); MINDENOM=30 has 2.
- `selected_min_d`: The feasible (non-no-tree) member with minimum D.
  MINDENOM=1 is selected here.
- The MINDENOM=56 member (`no_tree = TRUE`) terminates the family.

**D and ESS are distinct diagnostics.** ESS measures translational strength;
D measures parsimony. A model with high ESS and high strata count may have
worse D than a simpler model with fewer strata.

The staging table and propensity weights for the selected member
(`fam$members[[fam$min_d_idx]]$tree`) are accessed using the same
`cta_staging_table()` and `cta_propensity_weights()` functions shown above.

---

## 8. Lean-Fit Invariant

The `cta_tree` object stores only:

1. **`training_confusion`** — a C × C integer matrix of raw classification
   counts for the final selected tree (read by `cta_confusion_table()`).
2. **Per-leaf class-count vectors** — `class_counts_raw` (integer) and
   `class_counts_weighted` (double), each of length C (number of classes),
   stored directly on every terminal leaf node.

Nothing else related to reporting is cached:

- No endpoint summary tables.
- No staging tables.
- No propensity weight tables.
- No training observation indices or row identifiers.
- No predicted-label vectors for training observations.

All reporting functions (`cta_endpoint_summary`, `cta_endpoint_counts`,
`cta_staging_table`, `cta_propensity_weights`) read the leaf-count vectors
and compute everything on demand.  The translation functions
(`cta_assign_endpoints`, `cta_observation_weights`) additionally traverse the
tree on demand using caller-supplied `newdata`.  The fit object does not grow
when new reporting or translation functions are added.

---

## 9. Current Limitations

**Implemented since this doc was first written:**
- Observation-level propensity weights: `cta_assign_endpoints()` and
  `cta_observation_weights()` assign endpoint-level weights to individual
  observations on demand.

**Not yet implemented (deferred):**
- **Weighted propensity weights.** `cta_propensity_weights()` uses raw
  observation counts (`n_raw`) exclusively. A `weighted = TRUE` path using
  `n_weighted` and weighted marginals is not yet implemented.
- **Balance diagnostics.** Standardized mean differences or other
  covariate-balance checks after propensity-weight application are not
  provided. Future balance-diagnostic work should be canon-reviewed against
  the Linden/Yarnold JEP covariate-balance paper stored in `docs/theory/`.
- **Downstream outcome models.** Propensity-weighted outcome regression or
  marginal structural model fitting is not part of this package.

---

## References

Yarnold PR, Linden A (2017). Computing propensity score weights for CTA
models involving perfectly predicted endpoints. *Optimal Data Analysis*,
**6**, 43–46.

Yarnold PR, Soltysik RC (2016). *Maximizing Predictive Accuracy*. ODA Books.
