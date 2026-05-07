# Copilot Instructions for odacore

This repository is a pure-R reimplementation of the MegaODA / CTA classification engine. All algorithmic behavior must match canonical MegaODA.exe / CTA.exe output exactly. The following rules apply when reviewing or suggesting changes.

## Behavioral compatibility

- Preserve exact parity with legacy MegaODA.exe and CTA.exe outputs. Deviations from canonical output are bugs, not style choices.
- Do not suggest changes that alter candidate enumeration order, tie-breaking priority (MAXSENS → SAMPLEREP → FIRST IDENTIFIED), or LOO refit semantics without explicit justification against a canonical fixture.

## Gold fixtures and tests

- Do not weaken gold fixture assertions unless the sole purpose is consolidating duplicate stochastic fits into a single cached fit. Removing assertions without consolidation is a regression.
- All gold expectations come from verified MegaODA.exe / CTA.exe output. Do not adjust expected values without re-running the canonical executable.

## LOO STABLE

- LOO STABLE is a correctness gate, not a heuristic filter. The per-fold refit must replicate what the full engine would produce on n-1 observations.
- The algebraic ordered-cut LOO helper must match the full-refit result exactly. Any deviation is a bug.
- Do not short-circuit LOO by reusing the training rule. The previous admissibility shortcut was a known bug.

## eval_order and UniODA/CTA distinction

- Standalone `oda_univariate_core()` defaults to `eval_order = "mc_then_loo"`. This is the public API default and must not change.
- CTA tree-building uses `eval_order = "loo_then_mc"` internally to reject LOO-unstable ordered-cut candidates before spending MC iterations. This is a performance optimization, not a behavioral change.
- Do not conflate these two modes when reviewing parameter defaults or forwarding logic.

## Non-determinism and seeds

- Do not suggest `mc_seed` as the first fix for CTA gold parity failures. Investigate harness/settings/order issues first (duplicate fits, wrong STOPUP, wrong eval_order).
- CTA.exe is deterministic without a seed. If a test is flaky, the cause is almost always two independent stochastic fits of the same fixture, not RNG variance in a single fit.

## Candidate accumulation

- Candidate ordering in the ordered-cut inner loop is semantically meaningful (TERTIARY tie-break = first identified). Do not reorder, deduplicate, or batch candidates in ways that change enumeration order.
- Column names and types in the candidate data frame are used downstream by name. Do not rename or drop columns without auditing all callsites.
