# SORT / Staged CTA Workflow Plan (SCTA)

_Created: 2026-05-27. Updated 2026-05-28: method renamed SORT._

**Method name (updated 2026-05-28):** The workflow described in this document
is now named **SORT (Sequentially Optimal Recursive Trees)**. "Staged CTA" may
remain as a generic descriptive phrase, but **SORT** is the specific method name.

SORT begins from an SDA result or an explicitly declared sequential anchor (the
SDA-discovered root/sequence), then recurses R-to-L using EO-CTA/MDSA within
the allowed branch-specific candidate universe. See
`docs/LORT_SORT_GORT_TAXONOMY.md` for the full LORT/SORT/GORT contract.

---

## 1. Problem Statement

**FORCENODE is a node-level constraint; staged EX CTA is a workflow-level decomposition.**

These are not the same operation. Conflating them leads to systematically incorrect FP/TN
allocation even when overall ESS appears similar. This document records the empirical evidence
and derives the design requirements for a future `staged_cta_workflow` object.

### 1.1 What FORCENODE does

`FORCENODE 1 V37` forces V37 as the root split. CTA.exe then grows the best subtrees
below that root. MINDENOM is applied **globally** across all branches.

Consequence: FORCENODE cannot anchor a split when global MINDENOM exceeds the forced root's
smaller endpoint denominator. Example: FORCENODE V37 with MINDENOM=117 fails because V37's
smaller endpoint (V37>0.5) has n=116 < 117. CTA.exe emits Warning 4 and returns no tree.

### 1.2 What staged EX CTA does

Staged EX CTA separates the data into branches using `EX` filters, then fits independent
CTA models on each branch with branch-specific:
- candidate attribute sets
- MINDENOM values
- terminal no-tree decisions

This preserves the branch filter at every stage, allowing branch 1 to have MINDENOM=47
while branch 2 (no tree) is declared terminal, all while keeping the integrated confusion
table interpretable.

### 1.3 The single-global-MINDENOM limitation

A single global MINDENOM cannot simultaneously satisfy:
- branch-specific power requirements (each branch has a different n)
- branch-specific complexity targets (some branches warrant deep trees; others warrant stumps)
- explicit no-tree terminal endpoints (CTA.exe has no per-branch stop mechanism)

---

## 2. Data-Raw Artifact Summary

The following scratch artifacts exist at `data-raw/` (private, gitignored, not committed).
They document an exploratory CTA analysis on a private dataset with 9,014 observations,
binary class V1 (class 0: ~8,915 obs; class 1: ~99 obs), attributes V17 V19 V20 V21 V27 V37,
MISSING ALL (-9), LOO STABLE, ENUMERATE (where noted).

**Do not commit data-raw files. Do not copy raw data into any tracked file.**
Only model summaries and confusion counts are recorded here.

### Artifact inventory

| File | Type | Purpose | Delete after ingestion? |
|---|---|---|---|
| `CTA.EXE` | executable | CTA.exe binary used to run analyses | yes |
| `data.csv` | private data | Full dataset in CSV format | yes |
| `data.txt` | private data | Full dataset in CTA.exe TEXT format | yes |
| `cta0_1.pgm` | CTA program | Unrestricted EO, MD=1, FORCENODE commented out | yes |
| `cta0_1_FN1.pgm` | CTA program | FORCENODE V37, MD=1, no ENUMERATE | yes |
| `cta0_99_FN1.pgm` | CTA program | FORCENODE V37, MD=99, no ENUMERATE | yes |
| `cta0_117_FN1.pgm` | CTA program | FORCENODE V37, MD=117, no ENUMERATE | yes |
| `cta1_1.pgm` | CTA program | Staged EX: V37=0 branch (EX V37=0), MD=1 | yes |
| `cta2_1.pgm` | CTA program | Staged EX: V37=1 branch (EX V37=1), MD=1 | yes |
| `cta2_47.pgm` | CTA program | Staged EX: V37=1 branch, MD=47 | yes |
| `cta2_183.pgm` | CTA program | Staged EX: V37=1 branch, MD=183 | yes |
| `cta2_427.pgm` | CTA program | Staged EX: V37=1 branch, MD=427 | yes |
| `cta2_1354.pgm` | CTA program | Staged EX: V37=1 branch, MD=1354 | yes |
| `cta2_2190.pgm` | CTA program | Staged EX: V37=1 branch, MD=2190 | yes |
| `cta3_1.pgm` | CTA program | Staged EX: V37=1 AND V21>0.5 branch, MD=1 | yes |
| `cta4_1.pgm` | CTA program | Staged EX: V37=1 AND V21<=0.5 branch, MD=1 | yes |
| `cta5_1.pgm` | CTA program | Staged EX: V37=1 AND V21=1 AND V19=1 deep branch, MD=1 | yes |
| `MODEL0.TXT` | CTA output | Unrestricted EO results | yes |
| `MODEL0_FN1.TXT` | CTA output | FORCENODE V37, MD=1 results | yes |
| `MODEL0_FN1_99.TXT` | CTA output | FORCENODE V37, MD=99 results | yes |
| `MODEL0_FN1_117.TXT` | CTA output | FORCENODE V37, MD=117 — no tree | yes |
| `MODEL1.TXT` | CTA output | V37=0 branch results — no tree | yes |
| `MODEL2_1.TXT` | CTA output | V37=1 branch, MD=1 results | yes |
| `MODEL2_47.TXT` | CTA output | V37=1 branch, MD=47 — selected model | yes |
| `MODEL2_183.TXT` | CTA output | V37=1 branch, MD=183 results | yes |
| `MODEL2_427.TXT` | CTA output | V37=1 branch, MD=427 results | yes |
| `MODEL2_1354.TXT` | CTA output | V37=1 branch, MD=1354 results | yes |
| `MODEL2_2190.TXT` | CTA output | V37=1 branch, MD=2190 — no tree | yes |
| `MODEL3.TXT` | CTA output | Sub-branch: V37=1, V21=0, MD=1 | yes |
| `MODEL4.TXT` | CTA output | Sub-branch: V37=1, V21=1, MD=1 | yes |
| `MODEL5.TXT` | CTA output | Deep endpoint: V37=1, V21=0, V19=0, MD=1, attrs V17/V27 — no tree | yes |
| `README.md` | existing doc | Fixture provenance notes (keep) | no |
| `tmp_ort_blackbox_psa_test.R` | private probe | ORT blackbox acceptance gate on private data | yes |
| `tmp_ort_enumerate_abc_test.R` | private probe | ORT enumerate test | yes |
| `tmp_ort_manual_equivalence.R` | private probe | ORT manual equivalence gate | yes |
| `tmp_root_family_audit.R` | private probe | Root family audit probe | yes |

---

## 3. Model Comparison Table

All runs: LOO STABLE, MC ITER=5000 CUTOFF=0.05 STOP=99.9, PRUNE=0.05, MISSING ALL(-9).
Class: V1. Attributes unless noted: V17 V19 V21 V27 V37.

**Integrated confusion is reported; individual-branch confusion counts are not personally
identifiable. Raw data values are not reproduced.**

| Model | MINDENOM | EX filter | Tree | TN | FP | FN | TP | ESS | D | S | Notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Unrestricted EO (MODEL0, Enumerated) | 1 | none | V27→V37→V21→V19 | 6583 | 2332 | 38 | 61 | 35.46% | 3.64 | 2 | automated EO does not match staged EX error allocation |
| FORCENODE V37, MD=1 (MODEL0_FN1, Pruned) | 1 | none (FORCENODE) | V37→V21→V19 | 6573 | 2342 | 38 | 61 | 35.35% | 3.66 | 2 | similar ESS, FP=2342 vs staged FP=598 |
| FORCENODE V37, MD=99 (MODEL0_FN1_99) | 99 | none (FORCENODE) | V37→V21 | 6668 | 2247 | 41 | 58 | 33.38% | 3.99 | 2 | higher MD caps complexity but FP still high |
| FORCENODE V37, MD=117 (MODEL0_FN1_117) | 117 | none (FORCENODE) | no tree | — | — | — | — | — | — | — | FORCENODE fails: V37>0.5 branch n=116 < MD=117 |
| **Manual staged EX (integrated)** | 47 (V37=0 branch) | V37 split root | V21→V19/V20 + leaf | **8317** | **598** | **57** | **42** | **35.72%** | **8.999** | **5** | branch-specific EX workflow preserves FP control |

### Branch detail for manual staged EX result

| Stage | Branch filter | n (obs) | Attributes | MINDENOM | Result | Branch confusion |
|---|---|---|---|---|---|---|
| Root split | V37 ≤ 0.5 | 8898 | V21 V20 V19 V17 V27 | 47 | MODEL2_47 unpruned tree | TN=8317, FP=509, FN=57, TP=15 |
| Root split | V37 > 0.5 | 116 | V21 V20 V19 V17 V27 | 1 | MODEL1: no tree → leaf class 1 | FP=89, TP=27 |
| **Integrated** | — | **9014** | — | — | — | **TN=8317, FP=598, FN=57, TP=42** |

### MINDENOM sensitivity for V37=0 branch (EX V37=1)

| MINDENOM | V37=0 branch tree | Unpruned depth | Key structure |
|---|---|---|---|
| 1 | MODEL2_1 | 4 splits | V21→V19/V20→V19 (deepest) |
| 47 | **MODEL2_47** | 3 splits | V21→V19/V20 (**selected**) |
| 183 | MODEL2_183 | 2 splits | V21→V20 |
| 427 | MODEL2_427 | 1 split | V21 stump only |
| 1354 | MODEL2_1354 | 1 split | V21 stump only |
| 2190 | MODEL2_2190 | no tree | V21 right endpoint n=2189 < MD=2190 |

---

## 4. Methodological Lessons

### Lesson 1: Node-level vs. workflow-level

**FORCENODE is a node-level constraint; staged EX CTA is a workflow-level decomposition.**

FORCENODE tells CTA.exe which variable must be the root split. It does not:
- change how MINDENOM is enforced on the forced split's endpoints
- allow different MINDENOM values on the two branches below the forced root
- allow branch-specific candidate attribute sets
- permit "no tree" as an explicit terminal decision for one branch while continuing growth on another

Staged EX CTA does all of these by running independent CTA.exe invocations with explicit EX
filter stacks. Each branch is a separate model with its own MINDENOM, candidate set, and
terminal decision.

### Lesson 2: FORCENODE is gated by global MINDENOM

A forced root split requires **both** endpoints to have n >= MINDENOM. If the forced
variable's smaller endpoint falls below global MINDENOM, CTA.exe emits Warning 4 and
returns no tree — regardless of how much predictive signal the forced variable has.

Empirical proof:
- V37's smaller endpoint (V37>0.5): n = 116
- MINDENOM = 99: FORCENODE succeeds (116 >= 99)
- MINDENOM = 117: FORCENODE fails (116 < 117) → Warning 4, no tree

This means: **global MINDENOM constrains which roots can be forced**. A practitioner who
sets MINDENOM based on the full sample n (e.g., MPE power formula using total N) may
inadvertently make FORCENODE inadmissible for splits with small endpoints.

### Lesson 3: Similar ESS hides radically different error allocation

FORCENODE V37 MD=1 and manual staged EX produce nearly identical ESS (35.35% vs 35.72%).
However their FP counts differ by a factor of nearly 4:

| Model | FP | TP | ESS |
|---|---|---|---|
| FORCENODE V37 MD=1 | 2342 | 61 | 35.35% |
| Manual staged EX (MD=47) | 598 | 42 | 35.72% |

The staged EX approach trades 19 TP (61→42) for 1,744 fewer FP (2342→598). Whether
this trade-off is clinically or scientifically acceptable is a design decision, but it
**cannot be seen from ESS alone**. Error allocation must always be reported alongside ESS.

### Lesson 4: Single global MINDENOM is too blunt for staged workflows

The V37=0 branch (n=8898) supports a MINDENOM of up to ~1354 before losing the root split.
The V37=1 branch (n=116) supports a maximum MINDENOM of 116.

A global MINDENOM that satisfies the V37=1 branch power requirement (e.g., MINDENOM <= 116)
will under-constrain the V37=0 branch (where MINDENOM=47 was selected as the appropriate
level for the branch-specific analysis).

Branch-specific MINDENOM must be set per-branch, not globally.

### Lesson 5: No-tree is a legitimate terminal decision

MODEL1 (V37=0 branch) produces no tree with any of the candidate attributes. This is a
valid finding: the V37=0 subgroup (high-risk by root split) has no further structured
predictive signal from the available attributes. The leaf predicting class 1 for all
V37=0 observations (TP=27, FP=89) is the correct terminal action.

A staged CTA workflow object must explicitly represent "no tree" as a terminal stage result,
not as an error or exception. The integrated confusion table must account for all leaf
contributions, including no-tree terminals.

### Lesson 6: EX filter stacking produces branch-specific subsamples

CTA.exe EX filters are OR-combined exclusions. Multiple EX conditions in one program
file exclude a row if it satisfies any of the conditions:

```
EX V37=1 V21>0.5   -- keep: V37=0 AND V21<=0.5
EX V37=1 V21<=0.5  -- keep: V37=0 AND V21>0.5
EX V37=1 V21=1 V19=1  -- keep: V37=0 AND V21=0 AND V19=0
```

A staged CTA workflow object must represent each branch filter as an explicit predicate
(or list of predicates) that can be evaluated programmatically against new data.

---

## 5. Future Staged CTA Object Design

### 5.1 Proposed object: `staged_cta_workflow`

```r
staged_cta_workflow <- list(
  workflow_id          = character(1),   # unique run identifier
  mode                 = character(1),   # "manual" | "auto_ex" | future
  n_stages             = integer(1),
  stages               = list(           # one element per branch/stage
    list(
      stage_id          = integer(1),
      parent_stage_id   = integer(1),    # NA for root stages
      branch_label      = character(1),  # human-readable, e.g., "V37>0.5"
      branch_filter     = list(),        # programmatic predicates
      class_variable    = character(1),
      candidate_attributes = character(),
      fit_settings      = list(          # mindenom, mc_iter, alpha, loo, etc.
        mindenom        = integer(1),
        mc_iter         = integer(1),
        alpha_split     = numeric(1),
        prune_alpha     = numeric(1),
        loo             = character(1),
        enumerate       = logical(1)
      ),
      terminal_action   = character(1),  # "tree" | "leaf_class_0" | "leaf_class_1" | "no_tree"
      model             = NULL,          # cta_tree or NULL
      leaf_table        = NULL,          # data.frame: endpoint x class
      confusion_contribution = NULL,     # C x C matrix for this branch
      n_branch          = integer(1),    # n obs in this branch
      n_classified      = integer(1)     # n obs classified by this stage
    )
    # ... one entry per stage
  ),
  integrated_confusion  = NULL,          # C x C matrix summed across all stages
  integrated_ess        = numeric(1),
  integrated_d          = numeric(1),
  integrated_strata     = integer(1),
  settings              = list()         # global settings
)
```

### 5.2 Proposed functions

| Function | Purpose |
|---|---|
| `staged_cta_plan()` | Dry-run validator: check branch filters, MINDENOM feasibility, candidate sets |
| `validate_staged_cta_plan()` | Validate plan object without fitting |
| `fit_staged_cta()` | Fit each stage, build integrated confusion and leaf table |
| `staged_cta_confusion()` | Extract integrated confusion from fitted workflow |
| `staged_cta_leaf_table()` | Extract integrated leaf table across all branches |
| `staged_cta_propensity_table()` | Propensity weights per integrated endpoint |
| `predict.staged_cta_workflow()` | Predict new observations by routing through branch filters |

### 5.3 Branch filter representation

Each branch filter must be a list of atomic predicates:

```r
branch_filter <- list(
  list(variable = "V37", op = ">",  value = 0.5),   # V37=1 branch
  list(variable = "V21", op = "<=", value = 0.5)    # V21=0 sub-branch
)
```

Evaluation: keep row if ALL predicates are satisfied.

The EX filter in CTA.exe is the complement: a row is excluded if any EX condition
is satisfied. The branch filter is the logical AND of negated EX conditions, which is
equivalent to the set of conditions the row must satisfy to be included.

---

## 6. Implementation Sequence

### SCTA-0 (complete): Design doc only

This document. No code. Lessons ingested from data-raw. Cleanup plan prepared.

### SCTA-1: Integrated confusion helper from explicit terminal strata

Implement `staged_cta_confusion()` that takes a named list of (branch_filter, confusion)
pairs and returns the integrated C×C confusion matrix, ESS, and D statistic.

No CTA fitting required. Pure aggregation utility.

### SCTA-2: `staged_cta_plan()` dry-run validator

Validate a proposed staged workflow:
- Check all branch filters are non-overlapping and exhaustive
- Check MINDENOM feasibility per branch (warn if forced split endpoint < MINDENOM)
- Check candidate attribute availability per branch
- Return a plan object with proposed settings and diagnostics

### SCTA-3: `fit_staged_cta()` explicit branch workflow

Fit the CTA model for each stage:
- Apply branch filter to training data
- Call `cta_fit()` (or `cta_descendant_family()`) with branch-specific settings
- Record terminal_action: tree | leaf_class_0 | leaf_class_1 | no_tree
- Accumulate integrated confusion

### SCTA-4: SDA / auto_sda / CTA interop

Wire `staged_cta_plan()` / `fit_staged_cta()` into the SDA→CTA handoff:
- SDA identifies attribute subsets and candidate EX filters
- Staged CTA receives branch-specific candidate sets and MINDENOM recommendations
- Integrated confusion feeds back into model comparison (ESS, D, error allocation)

---

## 7. Cleanup Plan

The following files may be deleted after user approval. Do not delete until explicitly
authorized.

**Exact deletion commands (do not run — awaiting approval):**

```bash
# In repo root:
rm -f data-raw/CTA.EXE
rm -f data-raw/data.csv data-raw/data.txt
rm -f data-raw/*.pgm
rm -f data-raw/MODEL*.TXT
rm -f tmp_ort_blackbox_psa_test.R
rm -f tmp_ort_enumerate_abc_test.R
rm -f tmp_ort_manual_equivalence.R
rm -f tmp_root_family_audit.R
```

`data-raw/README.md` should be retained — it documents the fixture provenance for
the tracked test fixtures (CTA_DEMO, myeloma) and contains no private data.
