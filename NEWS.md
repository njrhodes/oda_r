# oda 0.1.4

## Bug fixes

* Weighted CTA Monte Carlo p-values now match CTA.exe: each permutation shuffles
  the attribute relative to the observations (class and case weight stay
  together) and the test statistic is the count-priors objective that selects
  the cut, not WESS. Reproduces CTA.exe significance for all 10 weighted myeloma
  root attributes. Unweighted CTA is unchanged.
* LOO refits each fold with the same optimizer that chose the training rule.
  UniODA: tied folds are refit (MAXSENS / SAMPLEREP) instead of taking the
  first-identified cut, and directional rules always refit (the shortcut
  ignored the direction); matches MegaODA LOO tables. CTA: node LOO refits each
  fold with the CTA scan (`.cta_loo_scan()`) instead of generic ODA; matches
  CTA.exe (cta_demo MODEL8 node 6 STABLE; myeloma node 2 V4 UNSTABLE).
* CTA fits every numeric attribute with the CTA scan, with or without case
  weights (uniform weights reduce its objective to ESS), so each node uses one
  optimizer for the cut, the MC p and LOO. Categorical attributes still use
  generic ODA.

# oda 0.1.3

## Bug fixes

* Weighted ordered CTA splits now match CTA.exe. `.cta_ordered_scan()` and
  `.cta_fast_scan_perm()` choose the cut that maximises count-priors weighted
  accuracy (sum of correct case weights divided by the raw class count) among
  admissible cuts, instead of the outermost cut with sensitivity > 0.5. This
  reproduces all 118 attribute WESS reports in the CTA.exe myeloma MINDENOM=1
  output (#19).
* `mindenom` is now a gate on the selected cut, as in CTA.exe: a split whose
  smaller side has fewer than `mindenom` observations is rejected rather than
  replaced by a different cut. This stops ESS rising as `mindenom` increases.
* Reported ESS is always the weighted ESS (class weight totals rescaled to 1),
  including with `priors_on = FALSE`, so it stays within [-100, 100] and does
  not depend on the scale of the weights (#20).
* The Monte Carlo permutation scan shares one implementation with the observed
  scan; the `mindenom` gate is not applied to permutations.
* Weighted numeric binary (0/1) attributes now use the same CTA scan, so their
  direction (and WESS, which can be negative) matches CTA.exe.
* The weighted CTA path now drops observations missing on the attribute before
  the scan, Monte Carlo permutations and LOO. Previously LOO scored missing rows
  through the rule, which understated LOO WESS for attributes with missing values.

# oda 0.1.2

## CRAN resubmission

* Added missing return-value documentation for SDA print and summary methods, as requested by CRAN.
* Replaced inappropriate `\dontrun{}` examples with `\donttest{}` or runnable examples, as requested by CRAN.
* Audited exported functions and registered S3 methods for missing return-value documentation.

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
