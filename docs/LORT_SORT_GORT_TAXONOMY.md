# LORT / SORT / GORT Method Taxonomy

**Status:** Active agent handoff contract. Read before touching recursive CTA code.
**Supersedes:** Generic "ORT" and "staged CTA" as primary method names.

**Canon boundary:** LORT, SORT, and GORT are odacore workflow-layer names, not
canonical MegaODA/CTA.exe commands or MPE Chapter 12 method names. They are valid
only insofar as their components and claims are explicitly anchored to novometric
canon: Axiom 1 statistical sample/power, Axiom 2 SDA when applicable, Axiom 3
MDSA descendant-family construction, D/ESS objective clarity, and Axiom 4
reproducibility.

- **LORT** is non-canonical workflow composition using canon CTA/MDSA components.
  It is greedy-local, not globally verified, and has no direct MegaODA/CTA.exe
  equivalent. It is valid insofar as each node-level CTA/MDSA fit is canon-anchored.
- **SORT** is a reserved SDA-anchored workflow, not yet implemented and not a canon
  method. It will be valid only after an explicit SDA-anchor and canon-constrained
  design is approved.
- **GORT** is future design only. No canon anchor is claimed or implied.
- **Chapter 12 / Novometrics** (MPE) is the governing anchor for all component
  claims. Do not describe LORT, SORT, or GORT as "canonical" without citing the
  specific canon-anchored component that supports the claim.

---

## Glossary

| Name | Full Name | Status |
|------|-----------|--------|
| **LORT** | Locally Optimal Recursive Tree | **Implemented** — current `recursive = TRUE` behavior |
| **SORT** | Sequentially Optimal Recursive Tree | **Reserved** — not yet implemented |
| **GORT** | Globally Optimal Recursive Tree | **Reserved** — future design/research only |

---

## LORT — Locally Optimal Recursive Tree

**What it is:**
LORT is the automatic greedy recursive CTA/MDSA workflow currently implemented
in `cta_fit(recursive = TRUE)`.

**Algorithm at each LORT node:**
1. Check guards (`max_nodes` → `min_n` → `max_depth` → pure-node).
2. Run `cta_descendant_family()` on the node-local sample.
3. Select the locally min-D complete CTA family member (`fam$min_d_idx`).
4. Embed that complete selected CTA member in the recursive tree.
5. Recurse under each terminal endpoint of the selected member (right first).
6. Stop when guards exhaust or no admissible tree found.

**"Optimal" scope:** Local to the node-level descendant-family scan only. LORT
does not evaluate downstream recursive consequences before selecting a root
member. There is no lookahead, no beam search, and no global comparison of
alternative recursive configurations.

**Object metadata (confirmed):**
```r
ort$ort_settings$method              == "lort"
ort$ort_settings$method_label        == "Locally Optimal Recursive Tree"
ort$ort_settings$selection_scope     == "local_node"
ort$ort_settings$objective           == "min_d_within_node_family"
ort$ort_settings$recursive_selection == "greedy_local_min_d"
ort$ort_settings$global_lookahead    == FALSE
ort$ort_settings$global_optimization == FALSE
ort$ort_settings$sda_anchored        == FALSE
ort$ort_settings$sort_compatible     == TRUE
```

**Print header (confirmed):**
```
Locally Optimal Recursive Tree (LORT)
  selection: greedy local min-D per recursive node
  global optimization: no
  SDA anchored: no
```

**Depth audit result (2026-05-28):** LORT produces genuine multi-level recursion.
On a 1,200-observation nested synthetic dataset, LORT reached ORT depth 3 with
4 non-terminal nodes embedding independent CTA members on disjoint data subsets.
strata_check_passed = TRUE in all tested cases.

---

## SORT — Sequentially Optimal Recursive Tree

**What it is:**
SORT is the reserved name for the SDA-anchored sequential recursive CTA/MDSA
workflow. It is not yet implemented.

**How it differs from LORT:**
- SORT begins from an SDA result or explicitly declared sequential anchor.
- The SDA pass supplies the first/root classifying rule, the selected attribute
  set, and the valid candidate universe for downstream recursion.
- SORT recurses R-to-L / endpoint-order using EO-CTA/MDSA within the
  permitted branch-specific attribute universe.
- SORT is not unconstrained: the SDA anchor constrains which attributes may
  enter at which stage.
- "Staged CTA" as a generic phrase may still be used descriptively; **SORT** is
  the specific method name when the procedure is SDA-anchored.

**Relationship to manual staged CTA.exe workflow:**
SORT is the package-level equivalent of the manual MegaODA/CTA.exe workflow:
`SDA → EX filter → CTA on each branch → integrate confusion`. The manual process
is now described as "the SORT pattern" to distinguish it from LORT.

**Relationship to HO-CTA / EO-CTA:**
SORT is spiritually related to higher-order CTA because it builds branchwise CTA
structure. It is not identical to HO-CTA unless proven by canon. Describe SORT
as SDA-anchored sequential recursive CTA/MDSA using EO-CTA/MDSA within branches.
Do not claim full HO-CTA canon status unless source evidence supports it.

**Entry point (reserved, not yet exported):** `sort_fit()` / `cta_sort_fit()`

**Required to begin a SORT task:**
1. An SDA result object or explicit SDA-equivalent anchor.
2. The SDA-selected attribute set and allowed candidate universe.
3. Branch/endpoint traversal rule (e.g., R-to-L, endpoint order).
4. EO-CTA/MDSA settings per branch.
5. Integrated terminal-strata and confusion reporting contract.
6. Proof that LORT behavior remains unchanged.

---

## GORT — Globally Optimal Recursive Tree

**What it is:**
GORT is a future design concept for global recursive search over all possible
recursive configurations. It is not implemented and is a reserved namespace only.

**What GORT would do (design concept only):**
- Evaluate candidate family members at a node together with downstream recursive
  consequences before selecting a root member.
- Compare alternative root choices after simulating recursion.
- Select the configuration that achieves a declared global objective:
  min global D, max global ESS, balanced objective, or constrained objective.
- Requires budget controls: beam width, candidate cap, `family_max_steps`,
  depth limit, node cap, memoization/cache policy, timeout.

**Entry point (reserved, not yet exported):** `gort_fit()` / `cta_gort_fit()`

---

## Method Relationships

```
LORT: automatic greedy local min-D at each node
       ↓
SORT: SDA-anchored; LORT-like recursion but constrained by SDA-selected set
       ↓ (extends further)
GORT: global search over recursive configurations before committing to root
```

- SORT is not LORT (it has a sequential anchor and constrained universe).
- SORT is not GORT (it does not globally search all recursive configurations).
- LORT can run without SDA. SORT requires SDA or an equivalent anchor.
- GORT subsumes LORT/SORT conceptually but is exponentially more expensive.

---

## Reserved Names (do not export or stub without explicit approval)

**SORT:**
- `sort_fit()`, `cta_sort_fit()`, `sort_plan()`, `sort_from_sda()`
- `sort_control()`, `sort_tree_table()`, `sort_endpoint_table()`

**GORT:**
- `gort_fit()`, `cta_gort_fit()`, `gort_plan()`, `gort_search_control()`
- `gort_compare_objective()`

---

## Agent Handoff Contract

### What agents may modify in LORT

- `print.cta_ort` / `summary.cta_ort` / `print.cta_ort_summary` output labels.
- `ort_settings` metadata fields.
- `cta_ort_node_table()` columns (additive only; do not remove existing columns).
- LORT test coverage (additive only).
- LORT docs (clarifications only; do not change the algorithm).

### What agents may NOT do without approval

- Change LORT selection logic (`fam$min_d_idx` local selection, right-first recursion).
- Add lookahead, SDA anchoring, or global optimization to `recursive = TRUE`.
- Rename S3 classes (`cta_ort`, `cta_tree`) or break existing dispatch.
- Export SORT/GORT names.
- Implement any part of SORT without an SDA source object and explicit approval.
- Implement any part of GORT without a formal design doc, objective definition,
  and compute budget spec.
- Change `ort_settings$method` from `"lort"` for current `recursive = TRUE` fits.

### GORT reserved namespace

Do not implement GORT mechanics inside any current function. Do not add
`global_search = TRUE` as a hidden flag inside LORT code. GORT is a clean future
API, not a parameter to the current recursion.

---

## Suggested Future SORT Task Template

```text
Task: implement SORT — SDA-anchored sequential recursive CTA/MDSA

Preconditions:
- SDA implementation complete and validated (SDA-1 through SDA-4B).
- An SDA result object is available with selected attributes, sequential anchor,
  and branch order.
- LORT tests passing (to prove LORT behavior unchanged by SORT addition).

Deliverables:
1. sort_fit(sda_result, X, y, ...) entry point.
2. SDA anchor validation (require non-NULL sda_result or equivalent).
3. Allowed candidate universe enforcement per branch.
4. R-to-L traversal using EO-CTA/MDSA within allowed attribute set.
5. Integrated terminal-strata table and confusion reporting.
6. predict.cta_sort() with path routing through SDA-anchored branches.
7. cta_sort_node_table() and cta_sort_endpoint_table().
8. Tests: SDA anchor required; LORT unchanged; SORT metadata stored.

SORT metadata:
  method = "sort"
  method_label = "Sequentially Optimal Recursive Tree"
  selection_scope = "sda_anchored"
  sda_anchored = TRUE
  global_optimization = FALSE
```

---

## Suggested Future GORT Task Template

```text
Task: implement GORT — globally optimal recursive search

Preconditions (all required before any code):
- Formal design doc: docs/GORT_DESIGN.md with objective definition,
  search strategy (beam/DP/exhaustive), compute budget spec,
  comparison metric, and proof LORT/SORT unchanged.
- Explicit approval of design doc.
- LORT tests passing.
- SORT tests passing (if SORT implemented).

Deliverables:
1. gort_search_control() budget spec object.
2. gort_fit(X, y, search_control, ...) entry point.
3. Comparison infrastructure: gort_compare_objective().
4. Budget controls: beam_width, candidate_cap, depth_limit, node_cap, timeout.
5. Tests: budget respected; LORT/SORT unchanged; GORT metadata stored.

GORT metadata:
  method = "gort"
  method_label = "Globally Optimal Recursive Tree"
  global_optimization = TRUE
  global_scope = declared_objective
```

---

## Legacy API Naming

The following public and S3 names contain `ort` rather than `lort`. They are
retained for backward compatibility and refer to the implemented **LORT** method
unless explicitly documented otherwise.

| Name | Type | Status |
|------|------|--------|
| S3 class `cta_ort` | class tag on `cta_fit(recursive=TRUE)` result | Legacy compat — means LORT |
| `predict.cta_ort()` | S3 method | Legacy compat — dispatches for LORT objects |
| `print.cta_ort()` | S3 method | Legacy compat — dispatches for LORT objects |
| `summary.cta_ort()` | S3 method | Legacy compat — dispatches for LORT objects |
| `plot.cta_ort()` | S3 method | Legacy compat — dispatches for LORT objects |
| `ort_plot_data()` | exported function | Legacy compat — computes layout for LORT plots |
| `cta_ort_node_table()` | exported function | Legacy compat — node table for LORT objects |
| file `R/cta_ort.R` | source file | Legacy compat — implements LORT engine |

**Rules:**
- Do not introduce new exported or public names that contain bare `ort` (without
  `lort`/`sort`/`gort` prefix or explicit legacy-compat annotation).
- Future LORT-named aliases (`lort_fit`, `predict.lort`, etc.) should be added in
  a deliberate compatibility slice, not this one.
- If a new `cta_ort` dispatch method is needed, annotate it explicitly as
  "legacy-compat name for LORT" in the roxygen title.
- Agents reading `cta_ort` in code should understand it means LORT unless the
  object's `ort_settings$method` says otherwise.

---

## Canon Boundary

| Component | Canon status |
|-----------|--------------|
| CTA/MDSA, ESS/WESS, MINDENOM, PRUNE, ENUMERATE, LOO, MC, D | Canon-anchored against CTA.exe fixtures |
| SDA (legacy `unioda_max_ess` / novometric `min_d`) | Canon-anchored where implemented per MPE |
| LORT | Package workflow built from canon-anchored pieces; not verified against CTA.exe |
| SORT | Reserved; will be SDA-anchored package workflow |
| GORT | Future design notation; no canon anchor claimed |

Do not describe LORT, SORT, or GORT as "canon" unless the specific mechanics
are verified against a CTA.exe program file or MPE source document.
