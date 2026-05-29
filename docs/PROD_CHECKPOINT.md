# PROD_CHECKPOINT.md

Production checkpoint after Graphics v3 (balance diagnostics + ggplot2 renderers)
and docs/export polish.

HEAD: pending Graphics v3D docs commit
Date: 2026-05-29

---

## Commits included since last checkpoint (89605b7)

| Commit | Message |
|--------|---------|
| 6b5bd06 | feat(balance): add univariate ODA balance tables and plot data |
| 88f84c6 | feat(balance): add CTA balance table and plot data |
| 9ccfdbd | feat(graphics): add ggplot CTA and LORT tree renderers |
| 1838ae7 | feat(graphics): add ggplot balance renderers |
| current | docs(graphics): document Graphics v3 ggplot renderers |

---

## Features landed in this checkpoint

### Balance diagnostics (v3B1/B2)

- `oda_balance_table(group, X, w, ...)` — univariate ODA balance: one ODA
  fit per covariate; returns ESS/WESS, p_mc, Sidak/Bonferroni corrections,
  rule summary, fit_ok flag.
- `smd_balance_table(group, X, w)` — conventional SMD companion (no p-values;
  descriptive arithmetic only).
- `oda_balance_plot_data(balance_table, smd_table, ...)` — pure transform;
  no fitting; renderer-ready rows with rank ordering.
- `cta_balance_table(group, X, w, mindenom, ...)` — single CTA fit for
  multivariate balance; status `"valid_tree"` / `"stump"` / `"no_tree"` /
  `"fit_error"`; `no_tree` = favorable balance evidence.
- `cta_balance_plot_data(cta_balance, ...)` — pure transform; populates
  `no_tree_message`; `cta_pd` is NULL for no_tree.

### Graphics v3C1 — tree renderers

- `plot_cta_tree(x, color_by, ...)` — ggplot2 CTA tree diagram.
- `plot_lort_tree(x, color_by, ...)` — ggplot2 LORT tree diagram.
- Both accept `cta_tree` / `cta_ort` directly or pre-computed plot-data.
- `color_by`: `"target_rate"` (gradient), `"prediction"` (discrete), `"none"`.
- `no_tree` CTA → message panel (not error).
- ggplot2 in Suggests; `.require_ggplot2()` guard; base methods unchanged.

### Graphics v3C2 — balance renderers

- `plot_oda_balance(x, ...)` — horizontal ESS/WESS lollipop dot plot; pure
  renderer; only what is in the pre-computed plot-data is shown.
- `plot_smd_balance(x, ref_010, ref_020, ...)` — horizontal |SMD| dot plot
  with optional threshold reference lines.
- `plot_balance_love(x, ...)` — direct alias for `plot_smd_balance`.
- `plot_cta_balance(x, color_by, ...)` — dispatches to message panel
  (no_tree/fit_error) or `plot_cta_tree` (valid tree/stump).
- Both `plot_oda_balance` and `plot_cta_balance` accept the parent table type
  via convenience coercion; never call the fitting function inside the renderer.

### Graphics v3D — docs/export polish

- `docs/GRAPHICS_V3.md` — full function reference, dependency policy,
  no-fitting rule, examples for all 6 functions, `no_tree` interpretation
  guide, `ggsave` usage.
- README — Graphics v3 section with tree diagram examples and balance plot
  examples; architecture table updated.
- DOCS_INDEX.md — Graphics v3 section added; checkpoints and next-slice order
  updated.
- ODACORE_VISION.md — Public API expanded with v3 graphics and balance
  functions; validation numbers updated.

---

## Validation results

| Command | Result |
|---------|--------|
| Full fast suite (`ODACORE_TEST_TIER=fast`) | FAIL 0 / WARN 0 / SKIP 165 / PASS 1566 |
| Targeted graphics+balance (`filter='graphics\|balance'`) | FAIL 0 / WARN 0 / SKIP 0 / PASS 110 |
| `devtools::check(vignettes=FALSE)` | 0 errors / 0 warnings / 1 NOTE (clock — known Windows environment issue) |

---

## Production invariants confirmed at this checkpoint

**Pure renderer contract:**
- All six Graphics v3 functions are pure renderers. None call fitting
  functions (`oda_fit`, `cta_fit`, `oda_balance_table`, `cta_balance_table`).
- The no-fitting contract is enforced by T9/B8 formal-arg inspection tests.

**ggplot2 dependency policy:**
- ggplot2 is in `Suggests` only. Package installs and loads without it.
- All v3 functions guard via `.require_ggplot2()` and error clearly.
- Base `plot.cta_tree` / `plot.cta_ort` methods are unchanged and do not
  require ggplot2.

**Balance semantics:**
- ODA balance table: ESS/WESS per covariate; p_mc via Monte Carlo Fisher
  randomization; Sidak/Bonferroni multiplicity correction uses only k = valid
  (non-NA) p_mc values.
- SMD companion: descriptive arithmetic; no p-values; no inference.
- CTA balance `no_tree` = favorable evidence of multivariable balance under
  declared constraints. It is not proof of universal balance and is not a
  model failure.
- Outcome variable is strictly out of scope for all balance functions.

**LORT / SORT / GORT:**
- LORT implemented; `method = "lort"`, `global_optimization = FALSE`,
  `sda_anchored = FALSE`.
- SORT and GORT remain reserved and not implemented.

**CTA / LORT output contract:**
- Post-pruning degeneracy gate active: no tree where all terminals predict
  the same class is returned.
- `degen = TRUE` does not exist for CTA or LORT.

---

## Previous checkpoint

See git history for 89605b7 (2026-05-28) — LOO semantics hardening.
Fast suite at that checkpoint: FAIL 0 / WARN 0 / SKIP 165 / PASS 1439.

---

## Known open items at this checkpoint

- Slice I: SDA → CTA/ORT anchor.
- Slice J: Staged CTA workflow implementation.
- Slice K: Weighted/staged adjustment design.
- Slice L: Vignettes / pkgdown.
- Slice M: Release hardening.
- Phase 2A: ODA S3 class + predict + summary/accessors — not yet started.
- Phase 2C: LOO p-value design — design notes not yet written.
- Phase 2G: Model comparison — not yet started.
