# CTA Ordered-Cut Eligibility Audit  -  Myeloma MINDENOM=1, Tree 1, Nodes 2 and 4

**Date:** 2026-05-11  
**Status:** Confirmed audit findings. MPE.pdf anchors are canon; fixture-derived CTA.exe parity findings are empirical canon for the named fixtures.

---

## Canon Anchors (MPE.pdf)

1. **ESS formula.** ESS is computed from Mean PAC across classes:
   ```
   ESS = (Mean PAC - C*) / (100 - C*) x 100
   ```
   For binary class (C* = 50): `ESS = 2 x MeanPAC - 100 = (spec + sens - 1) x 100`.

2. **Weighted ESS (WESS).** When WEIGHT is active, class-specific PAC is
   computed from priors-adjusted (weighted) confusion. WESS uses weighted
   class-specific PAC; ESS uses raw-count PAC.

3. **LOO STABLE is law (MPE.pdf verbatim):**
   > "LOO STABLE allows only attributes with LOO ESS equal to the ESS for
   > that attribute."

   **In weighted output:** LOO STABLE requires WESSL = WESS.  
   Signif T alone is insufficient. Candidate eligibility requires **both**
   Signif T **and** ESSL/WESSL = ESS/WESS.

4. **LOO is true re-fit/re-scan LOO (MPE.pdf verbatim):**
   > "for each observation, hold it out, obtain an ODA model on the remaining
   > observations, classify the held-out observation, then tabulate LOO
   > performance."

   The LOO model is obtained via generic ODA (`oda_univariate_core()`,
   global ESS maximum), not the CTA-specific ordered scan.

---

## Confirmed: CTA Ordered-Cut Eligibility Rule

**Evidence:** Myeloma MINDENOM=1, Tree 1, Node 4 (V4) and Node 2 (V4).

> **CTA.exe selects the rightmost ordered-cut candidate for which the
> class-1 priors-adjusted PAC on the right branch exceeds 0.5.**

Formally, direction "0->1":
```
sens_wa(j) = 1 - cumsum(adj_w1)[j]     (class-1 adj weight on right)

Eligible cuts:  { j : sens_wa(j) > 0.5 }
Selected cut:   argmax_j { eligible j }  (rightmost = largest j = highest cut value)
```

Equivalently: **maximise spec_wa subject to sens_wa > 0.5**, because spec_wa
is monotonically non-decreasing with j.

### Node 4 Evidence

**Node path:** V17 <= 0.5 AND V15 <= 0.5.  
**n = 113**, n0 = 92, n1 = 21.  
**Unique V4 values:** 109; candidate cuts: 108.

Key rows (direction "0->1", priors-adjusted):

| j  | cut    | nL  | nR | spec_wa | sens_wa | WESS%    |
|----|--------|-----|----|---------|---------|----------|
| 61 | 85.70  |  65 | 48 | 0.6445  | 0.8254  | 46.9929% | <- generic ODA peak |
| 83 | 309.50 |  96 | 17 | 0.8118  | 0.5206  | 33.2384% |
| 84 | 359.80 |  97 | 16 | 0.8284  | 0.5206  | 34.9000% | <- CTA canonical |
| 85 | 402.45 |  98 | 15 | 0.8284  | 0.4172  | 24.5557% | <- sens drops below 0.5 |

j = 84 (cut = 359.80) is the rightmost cut with sens_wa > 0.5.  
WESS = **34.90%**  -  matches CTA.exe fixture exactly.  
CTA.exe reports: **V4 Node 4 Tree 1 Signif F WESS 34.90%** (not significant; not selected).

### Node 2 Evidence

**Node path:** V17 <= 0.5.  
**n = 131**, n0 = 101, n1 = 30.

Key rows around threshold:

| j   | cut    | spec_wa | sens_wa | WESS%    |
|-----|--------|---------|---------|----------|
|  ~  | ~84.85 | ~0.55   | ~0.80   | 42.54%   | <- generic ODA peak |
| 95  | 309.50 | 0.7987  | 0.5148  | 31.3458% |
| 96  | 327.70 | 0.8138  | 0.5148  | 32.8579% |
| 97  | 339.10 | 0.8227  | 0.5148  | 33.7514% |
| 98  | 371.20 | 0.8341  | 0.5148  | 34.8886% | <- CTA canonical |
| 99  | 402.45 | 0.8341  | 0.4388  | 27.2920% | <- sens drops below 0.5 |

j = 98 (cut = 371.20) is the rightmost cut with sens_wa > 0.5.  
WESS = **34.89%**  -  matches CTA.exe fixture exactly.  
CTA.exe reports: **V4 Node 2 Tree 1 Signif T WESS 34.89%** (significant; but not selected  -  see LOO STABLE below).

### Contrast with generic ODA

| Criterion | R `oda_univariate_core` | CTA.exe |
|-----------|-------------------------|---------|
| Admissibility | Any class-0 left + any class-1 right | sens_wa (class-1 right PAC) > 0.5 |
| Selection | Global ESS maximum | Rightmost admissible (max spec_wa within gate) |
| Node 4 result | cut = 85.70, WESS = 46.99% | cut = 359.80, WESS = 34.90% |
| Node 2 result | cut ~ 84.85, WESS ~ 42.54% | cut = 371.20, WESS = 34.89% |

### Scope of confirmed evidence

The "rightmost sens_wa > 0.5" rule is confirmed for:
- Myeloma dataset, MINDENOM=1, WEIGHT V2 active
- Attribute V4 (continuous ordered), direction "0->1"
- Node 4 (V17 <= 0.5 AND V15 <= 0.5) and Node 2 (V17 <= 0.5)

Not yet verified against:
- CTA-demo (no WEIGHT command); see regression note below
- Other myeloma attributes (V9, V11, V12, V14, V16, V18, V19)
- MINDENOM=30 or MINDENOM=56 myeloma runs
- Direction "1->0" patterns

---

## Confirmed: CTA LOO STABLE Law at Node 2

**Evidence:** Myeloma MINDENOM=1, Tree 1, Node 2, candidates V4 and V15.

### Weighted confusion matrices

| Attr | Cut    | Dir | spec_wa | sens_wa | WESS     |
|------|--------|-----|---------|---------|----------|
| V4   | 371.20 | 0->1 | 0.8341  | 0.5148  | 34.8886% |
| V15  | 0.50   | 0->1 | 0.9100  | 0.2656  | 17.5654% |

Both match CTA.exe fixture values exactly (34.89% and 17.57%).

### LOO STABLE analysis

Canon: LOO STABLE allows only attributes with LOO ESS equal to ESS. In
weighted output this is WESSL equal to WESS.

The LOO model is obtained via generic ODA re-fit (`oda_univariate_core()`)
on n-1 observations per fold. For non-uniform case weights,
`oda_loo_ordered_cut_counts()` returns NULL and the per-fold full-refit path
executes.

| Attr | Full-model rule      | Cut    | WESS     | LOO (generic ODA) | WESSL    | \|delta\|   | STABLE |
|------|----------------------|--------|----------|-------------------|----------|---------|--------|
| V4   | CTA scan (sens>0.5)  | 371.20 | 34.8886% | Generic ODA refit | 38.5206% | 3.6320% | **NO** |
| V15  | Generic ODA          | 0.50   | 17.5654% | Generic ODA refit | 17.5654% | 0.0000% | **YES** |

**V4 is UNSTABLE** because the CTA-canonical cut (371.20) is not the generic
ODA optimum (84.85). Generic ODA per-fold LOO finds cut ~ 84.85 each fold,
producing WESSL = 38.52% != WESS = 34.89%. The 3.63 pp discrepancy marks V4
UNSTABLE and excludes it from selection.

**V15 is STABLE** because V15 is binary (values 0/1; only possible cut = 0.5).
Both the CTA scan and generic ODA resolve to the same cut (0.5) in the same
direction (0->1). Generic ODA LOO gives WESSL = WESS = 17.57% exactly.

Note: `.cta_ordered_scan()` applied to V15 selects direction "1->0" with
WESS = -17.57% (because the CTA eligibility gate requires sens_wa > 0.5 in
the chosen direction, and "0->1" at cut=0.5 gives sens_wa = 0.27 < 0.5). The
fixture value of +17.57% is the generic ODA result. CTA.exe uses generic ODA
for binary attributes at this node, not the CTA-specific ordered scan.

**Selection outcome:** Both V4 and V15 are Signif T. V4 is UNSTABLE ->
excluded. V15 is Signif T and STABLE -> selected as the highest eligible WESS
candidate. Confirmed against MODEL1.TXT: `V15 WESS 17.57% STABLE WESSL 17.57%`.

### LOO STABLE in CTA node selection

A candidate passes the node-selection gate only if:
1. **Signif T**  -  MC p-value < threshold (default 0.05)
2. **STABLE**  -  WESSL = WESS (generic ODA LOO performance = full-model performance)

Among all candidates passing both gates, the one with highest WESS is selected.

---

## Regression Note and Open Questions

### Regression discovered

Wiring `.cta_ordered_scan()` universally into `.fast_screen()` and
`.full_fit_one()` caused cta_demo fixture tests to fail: V4 selected as root
(ESS ~ 35.4%) instead of canonical V2 (ESS = 52.63%). CTA.exe does not apply
the "rightmost sens_wa > 0.5" rule to V2 in cta_demo.

`.cta_ordered_scan()` and `.cta_mc_ordered()` remain implemented as isolated
helpers in `cta_core.R` but are **not** wired into `.fast_screen()` or
`.full_fit_one()`.

### Open questions

1. What rule does CTA.exe apply to V2 in cta_demo? Does it use generic ODA?
2. Is the "rightmost sens_wa > 0.5" rule conditioned on WEIGHT being active,
   the attribute type, the node position, or something else?
3. Does cta_demo CTA.exe output show reported cuts for V2-V6 at root?

### Do not patch until resolved

- Do not integrate `.cta_ordered_scan()` into `.fast_screen()` or
  `.full_fit_one()` until the rule is confirmed against cta_demo.
- Do not patch ENUMERATE.
- Do not patch scoring.
- Do not wire `.cta_ordered_scan()` globally without the LOO STABLE gate
  (WESSL = WESS check using generic ODA LOO) implemented in the same patch.

---

## Known EXE LOO Limitations

The leave-one-out (LOO) results printed by CTA.exe and MegaODA.exe do not
always equal a true leave-one-out refit, i.e. a model refit on the other n-1
observations and then used to classify the held-out one. oda_r performs a
true refit. When oda_r's LOO disagrees with an EXE, check it against a manual
LOO run in the EXE itself before treating it as an oda_r bug (issue #25).

### CTA.exe: weighted LOO

**Observed.** With a `WEIGHT` variable, CTA.exe can report a node as
`LOO UNSTAB`, with a weighted LOO ESS near 0%. On the same data and weights:

- oda_r reports the LOO as stable (LOO table equals the training table);
- MegaODA.exe (`WEIGHT`, `LOO`) gives the same LOO table, LOO ESS and Fisher
  p as oda_r.

**Pattern.** CTA.exe's LOO ESS falls further below the true value as the
weights become more uneven, while the training cut, leaves, WESS and MC p
stay the same. With no weights, or nearly equal weights, CTA.exe's LOO agrees
with oda_r.

**Ruled out:**

- weight scale (rescaling the weights to sum to n fails the same way);
- the weight's column position in the data file;
- invalid weight values (none missing, non-finite or <= 0);
- the objective used to refit each fold (count priors, weighted ESS,
  unweighted);
- re-selecting the attribute within each fold (node-level LOO).

The cause inside CTA.exe is unknown.

**Consequence.** With uneven weights, CTA.exe `LOO STABLE` can reject splits
whose LOO is actually stable. CTA.exe remains a valid reference for cut
selection, WESS and MC p.

### MegaODA.exe: built-in LOO vs. its own refits

**Observed.** For an unweighted, non-directional UniODA model, MegaODA's
built-in `LOO` table differed from a manual LOO built by running MegaODA once
per fold on n-1 observations (LOO off). The manual LOO matched oda_r fold by
fold: same cut, same direction, same ESS, same held-out prediction. The
built-in table differed by one held-out observation.

**Probable cause (unconfirmed).** The held-out observation was one of the two
values the full-sample cut lies between, which suggests the built-in LOO does
not fully refit in that case.

### Recommended practice

- Weighted LOO: compare oda_r with MegaODA (`WEIGHT`, `LOO`), not with
  CTA.exe `LOO STABLE`.
- Any LOO disagreement: build a manual LOO in the EXE (one run per fold on
  n-1 rows, LOO off, then classify the held-out case with that run's rule)
  and compare oda_r fold by fold.
- CTA.exe parity for cut selection, WESS and MC p is unaffected.

### Weighted UniODA in oda_r (0.1.5)

Weighted fits with priors on select the cut by the count-priors objective,
`(W1/n1) * sens + (W0/n0) * spec`, as MegaODA.exe and CTA.exe do; the
reported ESS is still the weighted ESS. The MC permutes (class, weight) pairs
together and uses the same objective. The myeloma Node 4 V4 fixture
(`UNIV4.OUT`, `UNIV4MC.OUT`) pins cut, WESS, MC p and LOO to MegaODA output
(`tests/testthat/test-unioda-weighted-megaoda.R`).

---

## Current Code Status

- `.cta_ordered_scan()` and `.cta_mc_ordered()`  -  isolated helpers in
  `R/cta_core.R`. Not wired into `.fast_screen()` or `.full_fit_one()`.
- `tests/testthat/test-cta-ordered-scan.R`  -  unit tests scoped to myeloma V4
  audit evidence at Node 4 and Node 2. Does not assert universal integration.
- All existing CTA fixture tests pass.
- Generic `oda_univariate_core()`: weighted fits with priors on use count
  priors (0.1.5); unweighted fits are unchanged.
