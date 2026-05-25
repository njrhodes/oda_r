# Production Checkpoint

**HEAD:** `661920930349c99900fac45b78ee3bd2e43e7bb0`
**Date:** 2026-05-25

---

## Local Validation Gates

| Suite | Pass | Skip | Fail |
|---|---|---|---|
| Targeted fast | 99 | 10 | 0 |
| Fast tier (full) | 981 | 163 | 0 |
| `devtools::check(vignettes=FALSE)` | 0 errors / 0 warnings / 1 timestamp note | — | — |
| Smoke targeted | 821 | 15 | 0 |

---

## CI Status

GitHub Actions currently unavailable / quota exhausted until approximately June 18.

The most recent scheduled run (2026-05-25, run ID 26402945139) showed failure across all 4 matrix jobs (ubuntu/windows/macos release + ubuntu devel) before any steps executed — `steps: []`, jobs completed in ~2 seconds. This is an infrastructure/quota failure, not a package failure. The same HEAD (`f503339`) passed CI cleanly on 2026-05-16 as a push-triggered run.

**Local gates are the source of truth for this checkpoint.**

---

## Major Shipped Items

- **Public API:** `oda_fit()`, `cta_fit()` dispatchers; `predict`/`summary`/`print` S3 methods
- **DIRECTIONAL resolved:** Chapter 2 binary ordered (`"greater"`/`"less"`) and Chapter 4 categorical directional ODA fully implemented and canon-tested
- **CTA/MDSA/reporting stack:** `cta_endpoint_summary()`, `cta_endpoint_counts()`, `cta_staging_table()`, `cta_propensity_weights()`, `cta_assign_endpoints()`, `cta_observation_weights()`, `cta_confusion_table()`, `cta_family_table()` — all on-demand; lean-fit invariant maintained
- **NOVOboot fixed-confusion CI:** `novo_boot_ci()` with myeloma regression anchor
- **pkgdown / article scaffold:** `_pkgdown.yml`, gateway vignette `vignettes/odacore.Rmd`, 10 article skeletons in `vignettes/articles/`
- **G1 native CTA graphics polish:** `plot.cta_tree()` renders split nodes as ellipses, terminal nodes as rectangles, edges as directed arrows; adds `border_col`, `text_col`, `arrow_col`, `show_caption` parameters; `cta_plot_data()` is the renderer-independent data contract

---

## Known Future Work

- **Graphics v3:** ggplot2/grid renderer variant behind an explicit option
- **Full article content:** Slices 3B–3D (ODA articles, CTA articles, graphics + NOVOboot + validation tiers)
- **Precomputed article artifacts:** `inst/extdata/articles/*.rds` + `data-raw/build_article_artifacts.R` (Slice 3E)
- **pkgdown CI:** Wire `pkgdown` GitHub Actions deploy once Actions quota is restored
- **Balance diagnostics:** Deferred post-translation-stack
- **Weighted propensity variant:** `cta_propensity_weights()` adjusted variant stabilization
- **Full release hardening:** Full canon tier (`ODACORE_TEST_TIER=full`), multiclass CTA, SDA/auto_SDA scope
