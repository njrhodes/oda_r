# oda 0.1.1

## User-facing fixes

* Fixed `oda_fit(mcarlo = N)` so numeric `mcarlo` values enable Monte Carlo and set the iteration cap.
* Preserved rule, confusion, training metrics, and LOO details when a fitted model does not pass a LOO gate.
* Improved `print()` and `summary()` for ODA fits so rule-bearing models show the learned rule, training evidence, and LOO evidence.
* Enabled unweighted multiclass / multicategorical LOO reporting when `loo = "on"` is requested.
* Clarified that multicategorical LOO reports LOO confusion and ESS but does not report a LOO p-value because there is no canon-aligned C x C Fisher p-value.
* Clarified p-value directionality: MC permutation p-values are two-tailed;
  binary LOO Fisher p-values are one-tailed (MPE hold-out canon). Output
  labels now reflect the distinction.
* Fixed nominal multiclass rule display.
* Fixed multiclass summary output so class-level LOO PAC is visible.
* Fixed PAC percentage formatting to avoid double-scaling percent values.
* Added `myeloma` and `cta_demo` as package datasets for public examples.
* Re-enabled top-level vignette builds.
* Updated protein and migraine vignettes to reflect corrected LOO, rule, and p-value behavior.
* Kept `oda_sample_size()` runnable example meaningful while reducing CRAN example runtime.

## Tests and release hygiene

* Added fixture tests for directional categorical LOO: binary fixed `direction_map`
  (LOO ESS, Fisher p one-tailed, p < 0.001) and multiclass `direction = "ascending"`
  (LOO ESS, confusion, no LOO Fisher p, print states "not reported").
* Corrected directional-oda article: multiclass categorical LOO is supported with
  `loo = "on"`; clarified that MC p and LOO p are separate calculations.

* Fixed binary and categorical rule display - was showing `<categorical/binary rule>` placeholder. Now shows actual level-to-class mappings, e.g. `{low} --> 0 | {high} --> 1`.
* Added `cta_confusion_matrix(tree)` convenience wrapper that returns the 2x2 integer training confusion matrix directly from a `cta_tree` (previously required two-step `cta_confusion_table()` + `as_confusion_matrix()` call).
* Added `attr_names` length guard in `oda_cta_fit()` - error when supplied names do not match number of attributes.
* Improved `as_cta_candidates()` error message when `X` argument is omitted.
* Added regression tests for multiclass / multicategorical LOO reporting.
* Added tests for summary/print LOO evidence.
* Added scope guardrail tests for LORT lean-fit invariants and SORT/GORT export absence.
* Added CONTRIBUTING.md checklist for recursive CTA scope and release checks.
