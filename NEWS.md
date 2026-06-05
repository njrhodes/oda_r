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

* Added regression tests for multiclass / multicategorical LOO reporting.
* Added tests for summary/print LOO evidence.
* Added scope guardrail tests for LORT lean-fit invariants and SORT/GORT export absence.
* Added CONTRIBUTING.md checklist for recursive CTA scope and release checks.
