# WIRING_AUDIT.md

Reconciliation document for the odacore public release hygiene pass.
Records what was audited, what was found, and what disposition was applied.

---

## 1. Export and Public Surface Classification

All exported symbols inspected against canon contracts.

### ODA surface -- KEEP

| Symbol | Status |
|--------|--------|
| `oda_fit` | KEEP - primary dispatcher, CTA node entry point |
| `oda_univariate_core` | KEEP - canon UniODA engine |
| `oda_multiclass_unioda_core` | KEEP - canon MultiODA engine |
| `oda_rule_predict` | KEEP - canonical rule application |
| `oda_rule_predict_multiclass` | KEEP - multiclass variant |
| `oda_confusion_binary`, `oda_confusion_multiclass` | KEEP - metrics |
| `oda_mean_pac`, `oda_ess_from_meanpac`, `oda_ess_from_mean` | KEEP - metrics |
| `oda_predictions`, `oda_confusion`, `oda_metrics`, `oda_d_stat` | KEEP - S3 accessors |
| S3: predict/print/summary for oda_fit | KEEP |

### CTA surface -- KEEP

| Symbol | Status |
|--------|--------|
| `cta_fit`, `oda_cta_fit` | KEEP - public entry + compat alias |
| `lort_fit` | KEEP - LORT public API (adjacent composition, not Chapter 12) |
| `cta_node_table`, `cta_strata`, `cta_endpoint_*` | KEEP - translation stack |
| `cta_staging_table`, `cta_propensity_weights` | KEEP - translation stack |
| `cta_assign_endpoints`, `cta_observation_weights` | KEEP - translation stack |
| `cta_confusion_table` | KEEP - confusion from training_confusion |
| `cta_descendant_family`, `cta_family_table` | KEEP - MDSA/novometric |
| `cta_d_stat`, `cta_min_terminal_denom` | KEEP - metrics |
| `cta_plot_data`, `ort_plot_data` | KEEP - graphics contract |
| `cta_ort_node_table` | KEEP - LORT node accessor |
| S3: predict/print/summary/plot for cta_tree, cta_ort, cta_family | KEEP |

### SDA surface -- KEEP

| Symbol | Status |
|--------|--------|
| `sda_fit` | KEEP - SDA procedure entry |
| `sda_selected_attributes`, `sda_step_table`, `sda_candidate_table` | KEEP - accessors |
| `as_cta_candidates`, `sda_to_cta_data` | KEEP - interop bridge |
| `auto_sda_plan` | KEEP - planning utility |
| `sda_anchor`, `as_sda_anchor`, `validate_sda_anchor` | KEEP - anchor contract |
| S3: predict/print/summary for sda_fit, sda_anchor | KEEP |

### Graphics and balance surface -- KEEP

| Symbol | Status |
|--------|--------|
| `plot_cta_tree`, `plot_lort_tree` | KEEP |
| `plot_oda_balance`, `plot_smd_balance`, `plot_balance_love`, `plot_cta_balance` | KEEP |
| `oda_balance_table`, `smd_balance_table` | KEEP |
| `oda_balance_plot_data`, `cta_balance_plot_data` | KEEP |

### Bootstrap and matrix surface -- KEEP

| Symbol | Status |
|--------|--------|
| `novo_boot_ci` (default, oda_fit, cta_tree, cta_ort) | KEEP — all six dispatch paths wired |
| `as_confusion_matrix` | KEEP — tidy df → 2x2 matrix bridge for `novo_boot_ci.default` |

### Production tools surface -- KEEP

| Symbol | Status |
|--------|--------|
| `oda_readiness_check`, `oda_clean_missing_codes` | KEEP |
| `oda_validate_group`, `oda_validate_weights`, `oda_infer_attr_types` | KEEP |
| `oda_propensity_weights`, `lort_propensity_weights` | KEEP |

### Harness utilities -- KEEP (test-only)

`same_int_mat`, `pac_overall`, `loo_match_status` exported for test scripts; not part of
production inference API. Documented in NAMESPACE.

---

## 2. Evidence Wiring Status

### ODA evidence
- UniODA: covered against MegaODA.exe gold outputs (iris, unioda, tie-breaking, synthetic multiclass tests).
- Weighted priors: `oda_apply_priors`, `oda_apply_priors_multiclass` -- production.
- LOO: true refit-per-fold; weighted categorical LOO forbidden.
- MC p-value: Fisher randomization with Clopper-Pearson early stopping -- production.

### CTA evidence
- CTA_DEMO: `cta_demo/` fixture (MINDENOM=1, MINDENOM=8); parity tests passing.
- Myeloma: `myeloma/` fixture (MINDENOM=1/30/56, WEIGHT V2); all anchors verified.
- Weighted ordered scan: `.cta_ordered_scan` / LOO STABLE gate -- production (commit 85459a4).
- Root-only stump ENUMERATE phase -- production (commit a2e2a9d).
- Path-local missingness (`missing_action = "na"`) -- canonical; `"majority"` -- legacy compat.

### LORT evidence
- `lort_fit(recursive = TRUE)`: greedy local min-D at each node; no lookahead; no SDA anchor.
- Labeled `method = "lort"`, `global_optimization = FALSE`, `sda_anchored = FALSE`.
- Adjacent workflow-layer composition; NOT Chapter 12 SORT canon.

### SDA evidence
- `sda_fit`: Structural Decomposition Analysis procedure; traverses attribute space by class.
- Two modes: `novometric_min_d` (MPE-canon MDSA via `cta_descendant_family`), `unioda_max_ess`.
- SDA-1 (unweighted) complete; SDA-5 (weighted) deferred.
- SDA produces stage order, not propensity scores.

---

## 3. Graphics Bridge Status

- `cta_plot_data()` is the data-layer contract; `plot_cta_tree()` renders via ggplot2/grid.
- `ort_plot_data()` + `plot_lort_tree()` -- LORT variant.
- Balance plots: `plot_oda_balance`, `plot_smd_balance`, `plot_balance_love`, `plot_cta_balance`.
- Graphics v3 design documented in `docs/GRAPHICS_V3.md`.

---

## 4. Bootstrap Bridge Status

Six dispatch paths wired; each return object carries `source_type`, `source_id`, `weighted` provenance.

| Path | Evidence source | `source_type` |
|------|----------------|--------------|
| A. `novo_boot_ci(matrix)` | 2x2 matrix | `"matrix"` |
| B. `novo_boot_ci(oda_fit)` | `fit$confusion` (raw counts) | `"oda_fit"` |
| C. `novo_boot_ci(cta_tree)` | `tree$training_confusion` | `"cta_tree"` |
| D. `novo_boot_ci(cta_tree, node_id)` | leaf `class_counts_raw/weighted` | `"cta_tree_node"` |
| E. `novo_boot_ci(cta_ort)` | sum of `strata$class_counts` | `"cta_ort"` |
| F. `novo_boot_ci(cta_ort, stratum_id)` | single-stratum `class_counts` | `"cta_ort_stratum"` |

`as_confusion_matrix(cta_confusion_table(tree))` bridges the tidy df output to path A.
CI extraction uses lean stored objects; no fit-time artifact storage required.

---

## 5. Propensity Bridge Status

- CTA propensity: `cta_propensity_weights()` -- CTA-level stabilized endpoint weights; production.
- ODA propensity: `oda_propensity_weights()` -- production tools layer.
- LORT propensity: `lort_propensity_weights()` -- production tools layer.
- SDA propensity: deferred to SDA-5; forbidden in SDA-1. The `sda_anchor` metadata object
  enforces this boundary: `prohibited_downstream` includes `"propensity_weighting"`.

---

## 6. SDA / SORT / GORT Boundary

| Name | Status | Notes |
|------|--------|-------|
| SDA | Implemented (SDA-1) | Structural Decomposition Analysis; procedure, not a model |
| SORT | Reserved | Requires SDA anchor; reserved entry `sort_fit()`; not implemented |
| GORT | Reserved | Future design only; reserved entry `gort_fit()`; not implemented |
| LORT | Implemented | Adjacent composition; `method = "lort"`; not Chapter 12 canon |

Verification: `grep -rn "SORT.*implemented|GORT.*implemented" R/ docs/` returns no false positives.
All SORT/GORT references in docs are explicitly labeled "reserved" or "not implemented."

---

## 7. Agent / Task-Hook Boundary

- `sda_anchor` stores a `task_hook` metadata field (language or NULL); metadata-only.
- `prohibited_downstream` list in `sda_anchor` enforces contract: rejects `"propensity_weighting"`,
  `"fraud_demo"`, and similar out-of-canon uses.
- `source_call` (`match.call()` expression) is stored in sda_fit and novo_boot_ci objects.
  Print/summary methods do not expose `$call` or `$source_call` by default.
  Verified: no `cat`/`print`/`paste`/`sprintf` of `source_call` in print.sda_fit,
  summary.sda_fit, or print.sda_fit_summary.

---

## 8. Metadata and Privacy Hygiene

### Maintainer files

| File | Status |
|------|--------|
| `CLAUDE.md` | Added to `.gitignore`; not tracked |
| `docs/internal/PROD_CHECKPOINT.md` | Added to `.gitignore`; not tracked |
| `prompts/*` | Covered by `.gitignore` glob |

### Private data scan

- Scan for local filesystem paths in R/, man/, vignettes/: no hits outside
  `system.file`, temp functions, and test fixture paths.
- `data-raw/build_article_artifacts.R` referenced in `vignettes/articles/myeloma-cta.Rmd`
  inside an HTML comment (`<!-- TODO ... -->`); not rendered, not a privacy leak.
- `fraud_demo` in `sda_anchor.R` `prohibited_downstream`: example prohibited task type;
  no actual fraudulent data; acceptable as metadata.

### Author / institutional identifiers

- No email addresses or institutional identifiers found outside DESCRIPTION.
- DESCRIPTION contains only the declared package maintainer entry; appropriate for a public package.

### MegaODA.exe / CTA.exe binaries

- Not committed. `.gitignore` covers EXE files. Fixture outputs (.TXT, .pgm) are retained
  as canon text records; no binary executables in tree.

---

## 9. Fixture Provenance

| Fixture set | Source | Privacy |
|------------|--------|---------|
| `myeloma/` | survminer R package, GEO ID GSE4581 | Public; no PHI |
| `cta_demo/` | Synthetic CTA.exe demonstration data | No real-world subjects |
| `vignettes/` | Published MPE example datasets | Published; PMID 6643432 for protein set |

README files added:
- `tests/testthat/fixtures/README.md` -- top-level provenance index.
- `tests/testthat/fixtures/myeloma/README.md` -- myeloma-specific provenance.

Fixture size: myeloma/ 197K, cta_demo/ 182K, vignettes/ 65K -- within acceptable range.

---

## 10. Public Documentation Hygiene

### ASCII hygiene

All non-ASCII characters replaced with ASCII equivalents across:
- `R/*.R` source files (14 files cleaned)
- `man/*.Rd` generated documentation
- `vignettes/*.Rmd` (5 CRAN vignettes: gully, migraine, odacore, protein, refugee)
- `vignettes/articles/*.Rmd` (myeloma-cta, validation-tiers)
- `docs/*.md` (STATUS, DOCS_INDEX, LORT_SORT_GORT_TAXONOMY, ODA_CANON, CTA_CANON,
  CTA_TRANSLATION_STACK, COVARIATE_BALANCE_CONTRACT, GRAPHICS_V3,
  FIT_OBJECT_EVIDENCE_CONTRACT, SDA_ANCHOR_CONTRACT, CTA_ORDERED_CUT_AUDIT,
  examples/myeloma-cta-translation)
- `README.md`

Substitutions applied: section symbol -> "section ", arrows -> "->", em/en dash -> "-",
checkmarks -> "OK", smart quotes -> plain quotes, Greek letters (alpha/delta) -> spelled out,
box-drawing -> "+--"/"\--", not-equal -> "!=", triple-equals -> "===", ellipsis -> "...",
subscript arrows -> "v".

### SDA terminology

All occurrences of "Sequential Discriminant Analysis" replaced with "Structural Decomposition
Analysis" in README, docs, R, man, vignettes. Verified: no stale hits outside `docs/internal/`.

All "SDA model", "fitted SDA model", "Fit a... SDA model" replaced with "SDA procedure result"
or "Run a... procedure". Verified: no stale hits.

### Canon boundary language

- SDA: described as a procedure that traverses the attribute space by class; not a model or fit.
- LORT: described as adjacent workflow-layer composition; not Chapter 12 SORT canon.
- SORT/GORT: all references labeled reserved or not implemented.

---

## 11. Live-Fire Validation Results

### Targeted smoke (live-fire + fit-evidence + graphics + novo + sda-anchor)

```
FAIL 0 / WARN 0 / SKIP 8 / PASS 550
```

### Full fast suite (ODACORE_TEST_TIER=fast)

```
FAIL 0 / WARN 0 / SKIP 165 / PASS 1838
```

### devtools::check(vignettes=FALSE)

```
0 errors / 0 warnings / 1 note (clock skew -- environment artifact, not a package defect)
```

All canon fixtures passing: myeloma MINDENOM=1/30/56, CTA_DEMO MINDENOM=1/8.

---

## Summary

| Area | Finding | Disposition |
|------|---------|-------------|
| Export surface | All symbols classified | KEEP or DEFER-HIDE documented |
| ODA/CTA/LORT/SDA evidence | All canon anchors passing | KEEP |
| SORT/GORT | Not implemented | Reserved entries only |
| SDA propensity | Deferred to SDA-5 | Prohibited in SDA-1 by sda_anchor |
| Graphics bridge | cta_plot_data + renderers wired | KEEP |
| Bootstrap bridge | All 6 paths wired; source_type/source_id/weighted on all returns | KEEP |
| as_confusion_matrix | tidy df bridge for novo_boot_ci.default | KEEP |
| Private paths | None found in R/man/vignettes | Clean |
| Maintainer files | gitignored | Not tracked |
| Fixture provenance | README files added; sources documented | Clean |
| ASCII hygiene | Full pass across R/man/vignettes/docs | Clean |
| SDA terminology | "Structural Decomposition Analysis" everywhere | Fixed |
| source_call privacy | Not printed by default | Clean |
| Validation | 0 errors, 0 warnings | Pass |
