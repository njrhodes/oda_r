# PROD_CHECKPOINT.md

Production checkpoint after LOO semantics hardening, LORT taxonomy landing, and
CRAN-ability fix for test-cta-ort.R.

HEAD: 89605b7  
Date: 2026-05-28

---

## Commits included in this checkpoint

| Commit | Message |
|--------|---------|
| 9b542b4 | fix(loo): honor numeric LOO gates and report pvalue mode |
| 34d0fd2 | docs(loo): align LOO gate semantics with implemented behavior |
| 89605b7 | test(lort): cache recursive fixture for CRAN-safe timing |

Slices landed before this checkpoint (not listed above): LORT/SORT/GORT taxonomy
(LORT_SORT_GORT_TAXONOMY.md), synthetic no_tree / non-degenerate-output regression
tests (Slice B), post-pruning degeneracy gate (fbaaa85), CTA translation stack
first production version (f0328e2), all Phase 2 D/E/I components.

---

## Validation results

| Command | Result |
|---------|--------|
| `git status --short --untracked-files=all` | Clean — no modified or untracked files |
| Full fast suite (`ODACORE_TEST_TIER=fast`) | FAIL 0 / WARN 0 / SKIP 165 / PASS 1439 |
| Targeted CTA/LOO/LORT (`filter='cta\|cta-ort\|cta-loo-gate\|oda-fit'`) | FAIL 0 / WARN 0 / SKIP 142 / PASS 841 |
| Smoke canon fixtures (`ODACORE_TEST_TIER=smoke`, myeloma/cta-demo/loo/confusion) | FAIL 0 / WARN 0 / SKIP 15 / PASS 98 |
| `devtools::check(vignettes=FALSE)` | 0 errors / 0 warnings / 1 NOTE (clock — known Windows environment issue, not a package defect) |

---

## Production invariants confirmed at this checkpoint

**LORT / SORT / GORT:**
- LORT (Locally Optimal Recursive Tree) is implemented as a non-canonical
  workflow layer using canon CTA/MDSA components. It is not a new ODA engine.
  It does not add lookahead, SDA anchoring, or global search. Object metadata:
  `ort_settings$method == "lort"`, `global_optimization == FALSE`,
  `sda_anchored == FALSE`.
- SORT (Sequentially Optimal Recursive Tree) remains reserved and not
  implemented. Do not implement SORT inside `recursive = TRUE`.
- GORT (Globally Optimal Recursive Tree) remains reserved and not implemented.
  Requires approved design doc before any work.

**LOO gate semantics:**
- LOO STABLE and LOO PVALUE/numeric are distinct gates. They must not be
  conflated.
- `loo = "stable"`: binary only; accept when |WESSL − WESS| ≤ 0.01 pp; reports
  `loo_status = "STABLE"`.
- `loo = "pvalue"` or `loo = numeric`: binary only; accept when LOO Fisher p
  is **strictly less than** the threshold (default 0.05); reports
  `loo_status = "PVALUE"`. Numeric LOO is user-declared p-value gating, not
  STABLE.
- The threshold comparison is strict (`p < threshold`). The old doc phrase
  "less than or equal to" was wrong and has been corrected.
- `loo = "on"` is a synonym for `"pvalue"` when used with the multiclass engine.

**CTA / LORT output contract:**
- CTA and LORT public outputs are either a valid non-degenerate tree or
  `no_tree` / a terminal-guarded node (LORT recursion guards: `min_n`,
  `max_depth`, `max_nodes`).
- Post-pruning degeneracy gate (fbaaa85): any candidate where all terminal
  predictions cover fewer than C classes after pruning is rejected in the
  ENUMERATE loop before entering best-tree competition.
- `degen = TRUE` does not exist for CTA or LORT. It is a UniODA/MultiODA-only
  option.
- A CTA or LORT result in which all endpoints predict the same class is a model
  failure (`no_tree`), not a valid tree.

**CRAN timing:**
- The CRAN timing issue in `test-cta-ort.R` (~96 sec on default tier) was
  resolved by module-level fixture caching (37 repeated identical deterministic
  fits → 1 fit). Coverage is fully retained. The reproducibility test retains
  one independent fresh refit. Variant-argument tests continue to call
  `cta_fit()` directly.
- Current timing: ~7.4 sec (CRAN-safe).

---

## Known open items at this checkpoint

- Phase 2F (full status taxonomy on `cta_tree` S3: `degenerate` status field,
  ODA failed/degenerate policy) — partially complete; CTA/ORT/family degeneracy
  gate is production, ODA S3 status field is not.
- Phase 2A (ODA S3 class + predict + summary/accessors) — not yet started.
- Phase 2C (LOO p-value design) — design notes not yet written.
- Phase 2G (model comparison) — not yet started.
- Phase 2H (auto_SDA) — blocked on 2A + 2E.
- Phase 2J (vignettes) — blocked on 2A + 2B.
- Graphics v3 (mature ggplot/grid tree diagram) — not yet started.
- Balance diagnostics v1 — not yet designed.

---

## What is NOT in this checkpoint

- No production code changes. All three commits in this checkpoint are
  tests-only or docs-only.
- No new fixtures.
- No coverage deletions.
- No behavior changes to any fitting or inference function.
