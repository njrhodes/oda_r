# Copilot Instructions for oda

This repository is a pure-R reimplementation of the MegaODA / CTA classification engine. All algorithmic behavior must match canonical MegaODA.exe / CTA.exe output exactly for covered fixtures. The following rules apply when reviewing or suggesting changes.

## 1. MegaODA / CTA parity is the primary correctness standard

- Do not "improve" algorithms if doing so breaks golden fixture parity.
- CTA.exe / MegaODA.exe outputs are oracle references for all covered fixtures.
- Deviations from canonical output are bugs, not style choices.
- Do not suggest changes that alter candidate enumeration order, tie-breaking priority (MAXSENS → SAMPLEREP → FIRST IDENTIFIED), or LOO refit semantics without explicit justification against a canonical fixture.

## 2. LOO STABLE is a correctness gate

- LOO STABLE is a correctness gate, not a heuristic filter.
- Under LOO STABLE, Signif T alone is insufficient to accept a candidate split.
- Candidate eligibility requires ESSL/WESSL = ESS/WESS (|delta| ≤ 0.01 pp).
- The per-fold refit must replicate what the full engine would produce on n−1 observations. Do not short-circuit LOO by reusing the training rule.
- Weighted CTA output uses WESS/WESSL terminology. Never call weighted scores ESS/ESSL, and never call unit-weight scores WESS/WESSL.

## 3. Candidate ordering preservation

- Do not reorder candidate enumeration, attribute order, cut order, segment assignment order, or tie-breaking order unless a fixture trace proves the current order is wrong.
- "First identified" is a real tertiary tie-breaker. Reordering or deduplicating candidates changes outcomes.
- Column names and types in candidate data structures are used downstream by name. Do not rename or drop fields without auditing all call sites.

## 4. Gold fixture assertion policy

- Add or update tests against exact fixture anchors before changing any parity logic.
- Never change fixture expectations to match R output unless CTA.exe / MegaODA.exe evidence justifies the change.
- Use raw confusion matrices and exact/rounded ESS/WESS anchors from fixture files.
- Do not weaken gold fixture assertions unless the sole purpose is consolidating duplicate stochastic fits into a single cached fit. Removing assertions without consolidation is a regression.

## 5. eval_order: UniODA / CTA distinction

- UniODA / MultiODA candidate search order and CTA tree/root enumeration order are distinct concepts.
- Do not apply UniODA eval_order assumptions to CTA ENUMERATE unless CTA fixture traces prove the correspondence.
- Standalone `oda_univariate_core()` defaults to `eval_order = "mc_then_loo"`. This is the public API default and must not change.
- CTA tree-building uses `eval_order = "loo_then_mc"` internally (reject LOO-unstable candidates before spending MC iterations). This is a performance optimization, not a behavioral change. Do not conflate these two modes.

## 6. Seed / nondeterminism policy

- MC randomness is for significance testing only; it is not a license to ignore fixture parity.
- When debugging nondeterminism, record: mc_iter, mc_seed, node path, candidate attribute, ESS/WESS, p-value, and Signif flag.
- Do not "fix" stochastic mismatches by changing deterministic scoring rules.
- Do not suggest mc_seed as the first fix for CTA gold parity failures. Investigate harness settings, candidate order, and STOPUP first.
- Prefer fixture-anchored deterministic tests wherever possible. Stochastic tests are a last resort.

## 7. Missingness and ENUMERATE

- `predict.cta_tree()` supports `missing_action = "majority"` (default, majority-fallback) and `"na"` (canonical path-local missingness, returns NA_integer_).
- Do not globally apply path-local scoring to expanded ENUMERATE candidates. Expanded candidates are scored with majority-fallback `.predict_all()`.
- Root-only ENUMERATE stump candidates (CTA.exe Trees 5–7 pattern) are scored path-locally: observations missing the root attribute are excluded from ESS/WESS, not routed to majority class.
- Expanded ENUMERATE candidate behavior is fixture-locked. Do not patch the expanded loop to use path-local scoring — this was attempted twice and reverted because it incorrectly displaces correct expanded trees for MINDENOM=1.

## 8. MINDENOM

- MINDENOM is raw child-node row-count admissibility (unweighted, no priors adjustment).
- A candidate split is admissible only if both resulting child nodes have at least MINDENOM raw observations.
- Do not use weighted counts, priors-adjusted counts, or classified-only counts for MINDENOM enforcement.
- Do not touch MINDENOM enforcement logic without a fixture that demonstrates a failure.

## 9. Scope discipline

- Do not touch generic ODA (`oda_univariate_core`, `oda_multiclass_unioda_core`, `oda_fit`) when fixing CTA bugs unless the bug is proven to originate in those functions.
- Do not touch the weighted ordered scan / LOO STABLE gate (`.cta_ordered_scan`, `.cta_mc_ordered`, `.cta_full_fit_ordered`, `.full_fit_one`) unless a node-selection regression test explicitly names those functions as the cause.
- Do not modify fixtures, README correctness claims, or CLAUDE.md canon summaries casually. Each has a specific evidence standard.
- Do not add features, refactor code, or make "improvements" beyond what a specific failing test requires.

## 10. Extension vs. parity

- Binary CTA and UniODA/MultiODA behavior are parity-driven against MegaODA.exe / CTA.exe.
- Multiclass CTA has no gold executable benchmark. Any multiclass CTA implementation must be explicitly labeled as an `oda` extension, not MegaODA/CTA parity. Fixture strategy for multiclass CTA must be synthetic/property-based, not gold parity.
