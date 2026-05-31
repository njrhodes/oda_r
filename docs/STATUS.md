# odacore - Current Production Status

`odacore` is a pure-R reimplementation of the MegaODA.exe / CTA.exe statistical
classification engine. Where fixture coverage exists, correctness is defined by
the Windows executables; any divergence on a covered fixture is a bug in odacore.

---

## What is implemented

### Canonical core (fixture-validated)

| Component | Entry point | Notes |
|-----------|-------------|-------|
| Binary ODA | `oda_fit()` | C = 2; ordered, categorical, binary attributes |
| Multiclass ODA | `oda_fit()` | C >= 3; priors weighting, multiclass ordered/categorical |
| CTA | `cta_fit()` | ENUMERATE, LOO STABLE, PRUNE, MINDENOM; myeloma + CTA_DEMO fixtures |
| SDA | `sda_fit()` | Novometric `min_d` and legacy `unioda_max_ess` modes; SDA-1 through SDA-4B complete |
| MDSA descendant family | `cta_descendant_family()` / `cta_family_table()` | D parsimony, family comparison |
| Novometric bootstrap | `novo_boot_ci()` | S3 generic; ODA / CTA / LORT dispatch |
| Power / sample-size planning | `oda_power()` / `oda_sample_size()` | Unit-weighted binary 2×2 only; Fisher exact isomorphism (Rhodes 2020) |

### Adjacent workflow-layer (uses canon components)

| Component | Entry point | Canon components used |
|-----------|-------------|----------------------|
| LORT | `lort_fit()` | CTA/MDSA at each node; greedy local min-D; not Chapter 12 canon |

### Reporting and translation (on-demand; lean-fit invariant)

- `cta_confusion_table()`, `cta_staging_table()`, `cta_propensity_weights()`
- `cta_endpoint_summary()`, `cta_assign_endpoints()`, `cta_observation_weights()`
- `cta_family_table()`, `cta_node_table()`, `cta_strata()`
- `oda_balance_table()`, `smd_balance_table()`, `cta_balance_table()`
- `oda_propensity_weights()`, `lort_propensity_weights()`
- `sda_anchor`, `as_sda_anchor()`, `validate_sda_anchor()`

### Graphics (ggplot2; Suggests)

- `plot_cta_tree()`, `plot_lort_tree()`
- `plot_oda_balance()`, `plot_smd_balance()`, `plot_balance_love()`, `plot_cta_balance()`

---

## What is deferred / not implemented

| Item | Status |
|------|--------|
| SORT (`sort_fit`) | Reserved - requires SDA anchor; not implemented |
| GORT (`gort_fit`) | Reserved - future design only; not implemented |
| Weighted SDA | Deferred to SDA-5 |
| SDA-derived propensity weights | Deferred; SDA produces stage order, not propensity strata |
| Multiclass ODA power / sample-size | Deferred; not yet scoped |
| Multiclass CTA | Future extension |
| Generic propensity API | Not planned |
| Fraud / credit-card demos | Out of scope |

---

## Key invariants

- Covered fixture divergence = bug.
- `degen = TRUE` does not exist for CTA or LORT.
- No fit-time storage of training X/y in `cta_tree` (lean-fit invariant).
- LORT is adjacent workflow; not Chapter 12 canon.
- SDA does not estimate propensity scores.
- ESS = unweighted; WESS = case-weighted. Do not mix labels.

---

## Fixture coverage

| Dataset | MINDENOM variants | Status |
|---------|------------------|--------|
| CTA_DEMO | 1, 8 | Passing |
| Myeloma | 1, 30, 56 | Passing |
| Iris (multiclass ODA) | - | Passing |

---

## Test tiers

| `ODACORE_TEST_TIER` | Scope |
|---------------------|-------|
| unset / `cran` | CRAN-safe unit tests only |
| `fast` | Developer loop; skips slow canon fixtures |
| `smoke` | Production gate; required before CTA/MDSA/graphics commits |
| `full` | Release gate; required before release tags |
