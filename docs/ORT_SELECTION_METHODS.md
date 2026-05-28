# LORT / SORT / GORT Selection Methods — Design Memo

**Status:** Decision memo. Not implementation specification.
Do not code from this document without explicit approval.
Do not proliferate methods. One coherent staged build strategy.

---

## 1. Method Stack Overview

Three sequential roles. Not alternatives to each other — each feeds the next.

| Role | Method | Status |
|------|--------|--------|
| 1. Classification | Enumerated optimal CTA / MDSA | Implemented (canon-faithful) |
| 2. Attribute-set staging | SDA — identifies initial staged structure and candidate set | Not implemented |
| 3. Imbalanced recursive extension | Recursive CTA seeded from SDA structure; R→L with staged correction | Partial; greedy min-D only; canon not verified |

Do not add a fourth role without a canonical anchor in MPE.

### Concise staging rule

> Start upstream. Rebalance downstream only when declared, temporally valid,
> and constrained.

Upstream predictors are the **first** place to start — they anchor baseline
propensity and protect against leakage. They are not the only variables that
can ever enter the analysis. Later and time-varying attributes may rebalance
or refine strata within a recursive correction pass when their timing is valid
and their purpose is explicit.

### Terminology

**Exposure / assignment / action:** The process state, intervention, allocation,
policy, decision, behavior, or condition whose confounding structure is being
adjusted. Not limited to medicine — applies to any longitudinal, staged, or
policy/process setting requiring marginal confounding adjustment.

**Assignment propensity:** The probability or staged propensity that an
observation receives an exposure, action, or assignment, conditional on eligible
attributes at that stage.

**Assignment-mechanism predictors:** Attributes used to model or rebalance the
exposure/assignment process. These are **not** the same as the exposure variable
itself, and not the same as outcome predictors.

**Outcome predictors:** Attributes used to model the final response or outcome.
May overlap with assignment-mechanism predictors, but the role must be declared
explicitly.

**Baseline / pre-exposure stage:** Attributes measured before the exposure,
assignment, or action. Used to construct initial propensity strata. Temporal
ordering must be preserved; post-exposure variables must not enter this stage.

**Correction / rebalancing stage:** Later-stage or time-varying attributes used
*within* existing strata to refine balance or classification, subject to temporal
validity. Structured correction passes — not uncontrolled all-column scans.

**Marginal / sequential confounding control:** Multi-stage adjustment where
propensity or balance is updated across stages as new valid time-indexed
information becomes available.

**R→L recursive pass:** Recursive evaluation proceeds right branch then left
branch, using the declared stage-eligible variable set at each depth.

---

## 2. Enumerated Optimal CTA / MDSA

Default first step for problems without severe class imbalance.

`oda_cta_fit()` / `cta_descendant_family()` / ENUMERATE pipeline.
Canon-faithful against CTA.exe fixtures.

### Objective: max ESS vs. min D

Declare an objective explicitly. These two objectives are separable and can
diverge — both are meaningful.

- **Max ESS** — translational prediction accuracy. Use when the goal is
  accurate classification of new cases.
- **Min D** — theoretical parsimony. The most parsimonious representation
  of the classification structure, achieved when the D-minimizing model in the
  MDSA family is selected.

```
D = S × (100 / ESS − 1)
  = 100 / (ESS / S) − S
```

MINDENOM does not directly enter D. It governs admissibility and family
generation. Different MINDENOM values produce trees with different ESS and S;
D varies as a consequence.

Report both when they diverge. Do not conflate them.

### MINDENOM and Axiom 1 / power

MINDENOM must be tied to statistical power, not convenience. MPE guidance:

- For a **moderate effect**: N ≈ 32 cases per class per category required.
  → MINDENOM = 64 per endpoint.
- For a **strong effect**: N ≈ 12 cases per class per category required.
  → MINDENOM = 24 per endpoint.

Do not set MINDENOM to 1 and claim the result satisfies power requirements
unless the sample clearly exceeds these thresholds.

### Optimal pruning: bounded complete-model evaluation

After a tree is already grown, optimal pruning enumerates sub-branches to
identify the best complete-model candidate:

1. Start from an already-built tree.
2. Identify all sub-branches of each emanating branch.
3. Enumerate unique complete-model combinations.
4. Compute integrated confusion tables for each complete candidate.
5. Select max-ESS candidate for translational accuracy.
6. Compare min-D candidate for theoretical parsimony.

In the influenza CART example from the optimal-pruning paper: L1-R2 maximizes
ESS = 33.54%; L1-R1 minimizes D = 5.77. These diverge explicitly. The paper
reports both by purpose.

This is bounded complete-model pruning/evaluation after the tree exists. It is
not SDA and not the current greedy ORT. It belongs here as an eventual
`cta_optimal_prune()` helper; not a current build priority.

---

## 3. SDA Staging

SDA (Sequential Discriminant Analysis) is the **upstream attribute-set
identification procedure** used before MDSA/GO-CTA. It is not a variant of
CTA or ORT — it identifies the initial staged structure and candidate set that
MDSA receives.

MPE Chapter 12 canon: GO-CTA with multiple attributes can be obtained two ways:

1. **MDSA without SDA**: identify the full descendant family directly across all
   candidate attributes.
2. **MDSA with SDA**: first use SDA to identify the optimal attribute subset,
   then apply MDSA to that reduced set.

SDA is a disciplined attribute-set reduction and decomposition step. It is not
a permanent ban on subsequent structured modeling. After SDA establishes initial
strata, correction/rebalancing passes using later or time-varying attributes
remain valid when their timing is explicitly declared and their purpose is
stratum refinement, not uncontrolled global scanning.

### Two SDA modes: legacy vs. MPE novometric

These modes differ in the per-step selection rule. The package must distinguish
them explicitly rather than silently mixing them.

---

**Mode 1 — Legacy / pre-novometric SDA (iterative UniODA)**

Used before novometric theory. Per-step selection rule: **maximum ESS** from
iterative UniODA analyses.

Mechanics:
1. Run UniODA on each candidate attribute.
2. Select the attribute/rule with maximum ESS.
3. Apply the selected rule.
4. Remove correctly classified observations.
5. Repeat on the remaining unresolved sample until stopping conditions are met.

---

**Mode 2 — MPE novometric SDA (EO-CTA/MDSA per attribute)**

Novometric SDA uses D as the per-step selection criterion because per-attribute
descendant families and multicategorical parses can have more than two strata.
Minimum D is the theoretical parsimony measure appropriate for these richer
model spaces.

MPE mechanics (Chapter 12):

**Step S = 1:**
- Apply EO-CTA and MDSA to each individual attribute separately.
- Select the GO-CTA model that:
  (a) satisfies the minimum denominator criterion (MINDENOM/Axiom 1),
  (b) achieves p < 0.05 (MC significance),
  (c) achieves the **minimum D statistic** among eligible candidates.

**Step S = 2 and beyond:**
- Delete observations correctly classified in step S from the total dataset.
- Remove the selected attribute from the candidate attribute set (ATTR command).
- Repeat EO-CTA/MDSA on each remaining attribute using the reduced sample.

**Termination** (first condition met):
- All observations are correctly classified.
- All observations in either class are correctly classified.
- Axiom 1 / statistical power is violated.
- p > 0.05 (no remaining attribute achieves significance).

**Hand-off:** Upon termination, retain the attributes identified across prior
SDA steps. Apply MDSA to this reduced SDA-selected attribute set. This is the
initial constrained candidate set — it is the starting point for recursive CTA,
not a ceiling on the full analysis. Subsequent correction/rebalancing passes
(see §4) may admit additional declared variables when temporally valid.

---

### First implementation note

A first implementation using `sda_mode = "unioda_max_ess"` (Mode 1) is
acceptable **only if** it is explicitly documented as a legacy/iterative UniODA
staging approximation and not presented as MPE novometric SDA.

### Private rare-event example status

A private staged analysis begins with upstream/pre-exposure baseline propensity staging.
A baseline indicator attribute may have been selected by univariate max ESS at step 1 (Mode 1).

MPE novometric SDA (Mode 2) would require verifying per-attribute EO-CTA/MDSA
with MINDENOM, p, and min-D gates before declaring step-1 selection. That
verification has not been performed.

**Do not call the private staged tree "MPE SDA canon."**
State instead: *the private example supports the need for SDA/staging; exact canon mode
remains to be implemented and validated.*

After upstream staging establishes initial propensity strata, later-stage or
time-varying attributes may be used within strata to rebalance or refine
classification — provided their timing is valid and their role is explicitly
declared as assignment-mechanism correction, not outcome prediction.

The divergence between greedy min-D ORT and the private SDA staging tree is expected and
correct. They are outputs of different methods with different objectives.

**Do not make `cta_fit(recursive = TRUE)` chase the SDA tree.** SDA needs
its own API.

---

## 4. Imbalanced Recursive CTA / ORT-Type Extension

For problems with **severe class imbalance**, a recursive R→L ORT-type method
can be employed as a third-role tool, after SDA has identified the upstream
attribute set.

### Current state

`cta_fit(recursive = TRUE)` implements greedy recursive min-D CTA/MDSA:

```
fit_node(rows, depth):
  1. check guards (max_nodes → min_n → max_depth → pure no_tree)
  2. run cta_descendant_family()
  3. select fam$min_d_idx  (local min-D winner)
  4. recurse only the selected winner's endpoints, right then left
  5. discard family losers
```

This is useful infrastructure. It is **not** claimed as canonical ORT until
verified against CTA.exe program files.

### Intended approach

1. **Seed / root from SDA-selected upstream structure.** Do not begin with
   unconstrained greedy local min-D on the full candidate frame. SDA identifies
   the initial staged structure and constrained candidate set. The recursive
   CTA extension starts from that anchor.

2. **Proceed R→L recursively** using the stage-eligible variable set at each
   depth:
   - Baseline/pre-exposure assignment-mechanism predictors not yet exhausted by
     prior SDA steps.
   - Later-stage or time-varying attributes — attributes that become valid at
     a subsequent process state — **if** their timing is valid and their
     purpose is explicit stratum correction/rebalancing.
   - Do not run uncontrolled all-column scans. Do not violate temporal/causal
     ordering. Do not confuse assignment-mechanism predictors with outcome
     predictors.

3. **Declare the recursive objective explicitly**: max ESS, min D, or another
   canon-supported target. Do not let unconstrained greedy local decisions
   masquerade as globally optimal ORT.

The correction/rebalancing passes are inside this third-role method — not a
fourth method.

### Generalized staged workflow

This framework applies to any longitudinal, staged, or policy/process setting
requiring marginal confounding adjustment — not only medicine.

1. **Baseline / pre-exposure stage:** Attributes measured before the exposure,
   assignment, or action. Build initial propensity or baseline strata.
2. **Recursive correction / rebalancing stage:** Refine strata using declared
   stage-valid assignment-mechanism predictors. Preserve temporal ordering.
   No uncontrolled all-column scans.
3. **Later-stage / time-varying assignment stage:** Define the assignment or
   process state at that stage. Model assignment propensity using eligible
   predictors. Update weights, strata, or balance diagnostics.
4. **Outcome / response stage:** Model the final outcome using the adjusted
   staged structure. Keep assignment-mechanism modeling distinct from outcome
   modeling.
5. **Final object:** Auditable staged ODA workflow for marginal or sequential
   confounding control.

### Canon uncertainty

- Greedy recursive min-D ORT has not been verified against any CTA.exe program
  file.
- Do not claim ORT is globally optimal until verified.
- Do not claim current ORT reproduces SDA staging without an explicit
  SDA-anchored root and constrained attribute set.

---

## 5. What the Recent Tests Showed

### 5.1 Root family audit — 6-variable private rare-event frame

n = 9,014 | p = 6 | target_n = 99 | outcome = `target_event`

Candidates: `baseline_signal_A`, `baseline_covariate_E`, `baseline_covariate_D`,
`baseline_covariate_C`, `baseline_covariate_B`, `baseline_covariate_F`.

Root family (seed = 42, mc_iter = 25,000, max_steps = 20):

| MINDENOM | root_split | S | ESS | D | min_term_denom |
|---|---|---|---|---|---|
| 1 | baseline_covariate_E | 12 | 38.28% | 19.3464 | 3 |
| 4 | baseline_covariate_E | 11 | 38.14% | 17.8441 | 16 |
| 17 | baseline_covariate_E | 6 | 37.38% | 10.0497 | 39 |
| 40 | baseline_signal_A | 4 | 35.35% | 7.3168 | 116 |
| **117** | **baseline_covariate_D** | **3** | **30.51%** | **6.8333** | 206 |
| 207 | baseline_covariate_B | 3 | 25.28% | 8.8654 | 940 |
| 941 | baseline_covariate_C | 3 | 25.28% | 8.8654 | 1390 |
| 1391 | baseline_covariate_B | 2 | 22.55% | 6.8690 | 2269 |
| 2270 | baseline_covariate_E | 2 | 11.24% | 15.7941 | 2640 |
| 2641 | no_tree | — | — | — | — |

**Previously reported min-D ORT winner: MINDENOM = 117, root = `baseline_covariate_D`,
D = 6.8333.**

**INVALID — exposed as degenerate (2025-05):** The MINDENOM = 117 S = 3 member
was subsequently found to produce an all-class-0 terminal prediction set after
Sidak-Bonferroni pruning collapsed the class-1-predicting branch to a
majority_class = 0 leaf (the dataset has severe class imbalance).  A tree whose
terminal predictions cover only one class is degenerate and must not be accepted
as a valid CTA or ORT result.  With the degeneracy gate added to the ENUMERATE
candidate selection loop, this member is now rejected before entering the
best-tree competition.  `min_d_idx` must not point to a degenerate member;
if all members at or above MINDENOM = 117 are also degenerate, the correct
output is no_tree at that recursive node.

Key constraint: at MINDENOM = 117, `baseline_signal_A` is ineligible — its
positive-class branch has n = 116 < 117. At MINDENOM = 40 it is eligible but
D = 7.317 > 6.833, so it loses under min-D selection.

### 5.2 Black-box ORT test result

Current greedy ORT (root = `baseline_covariate_D`, D = 6.8333) failed all
structural checks for the private SDA staging target tree. Wide-newdata prediction passed.

**Conclusion:** Implementation predicts safely from wide newdata. Greedy
min-D ORT is not and should not be expected to reproduce the SDA staging tree.

### 5.3 Manual final-tree distinction

The private manual staged tree is a **five-stratum SDA/staging final tree**:

| TN | FP | FN | TP | Specificity | Sensitivity | ESS | D | S |
|---|---|---|---|---|---|---|---|---|
| 8,317 | 598 | 57 | 42 | 93.292% | 42.424% | 35.716% | 8.999 | 5 |

D check: `5 × (100 / 35.716 − 1) ≈ 8.999` ✓

This is **not** the MINDENOM = 40 CTA family member:

| | MINDENOM = 40 family member | Private manual tree |
|---|---|---|
| Object | Single CTA/MDSA family member | Final five-stratum SDA staging tree |
| S | 4 | 5 |
| ESS | ~35.35% | 35.716% |
| D | ~7.317 | 8.999 |

Do not compare these as if they are the same tree.

### 5.4 A×B×C design probe — tested, not explanatory

Private probe: enumerate root candidates A, evaluate B/C endpoints before
committing, score by D_product, D_sum, ESS_min, ESS_hmean.

| Criterion | Best A |
|---|---|
| min D_product | `baseline_covariate_D` |
| min D_sum | `baseline_covariate_D` |
| max ESS_min | `baseline_covariate_B` |
| max ESS_hmean | `baseline_covariate_B` |
| Target | `baseline_signal_A` — not recovered |

For `baseline_signal_A`: all four endpoints returned no-tree. Not competitive
under any global score.

**Conclusion:** A×B×C lookahead does not explain or recover the private SDA staging
tree. Do not implement ORT-level A×B×C.

---

## 6. Commit Plan

### Commit 1 — Safe recursive CTA infrastructure

Method-neutral, canon-safe changes only:

- `R/cta_ort.R`: `mc_stop`/`mc_stopup` threading, `family_max_steps` budget
- `R/oda_fit.R`: `mc_stop`/`mc_stopup`/`family_max_steps` plumbing
- `R/cta_s3.R`: wide-newdata safe column routing
- `docs/CTA_ORT_DESIGN.md`: training X contract, §3.7 (`family_max_steps`),
  `mc_stop`/`mc_stopup` signature entries
- `man/cta_fit.Rd`: `X` item expanded, `family_max_steps` Rd item
- `tests/testthat/test-cta-ort.R`: Tests 31–37 (mc threading, family_max_steps),
  T38–T44 (wide-newdata routing)

### Commit 2 — Design memo

- `docs/ORT_SELECTION_METHODS.md` (this file)

### Deferred — do not commit until resolved

- `min_target_n` / `target_class` implementation and Tests 22–30: not yet
  canon-verified; entangled with unresolved ORT objective question
- SDA implementation (either mode)
- ORT-level A×B×C (tested; not explanatory; do not implement)
- Global/lookahead ORT (unresolved; no canonical anchor)
- Private tmp scripts (deleted, gitignored — see Checkpoint 9 cleanup)
- Any claim that current greedy min-D ORT reproduces the private SDA staging tree
