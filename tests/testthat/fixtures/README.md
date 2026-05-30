# Fixture provenance

This directory contains canonical test fixtures used for odacore parity and regression testing.
Fixtures are validated against MegaODA.exe / CTA.exe golden outputs.

## myeloma/

Public survminer/GEO fixture, transformed for CTA parity tests.

- Source: `myeloma` dataset from the survminer R package (GEO ID GSE4581).
- Public data; no PHI; no private institutional data.
- Used for weighted CTA, MINDENOM (1/30/56), LOO STABLE, missing-code handling,
  endpoint counts, and no-tree behavior.
- See `myeloma/README.md` for full provenance.

## cta_demo/

Synthetic demonstration fixture distributed with CTA.exe.

- Used for unweighted CTA parity (MINDENOM = 1 and 8).
- No confidential data.

## vignettes/

Vignette article fixtures (Example-1 through Example-4 from MPE Chapter 2 and 4).

- These are published-example datasets from Yarnold & Soltysik (2005/2016).
- Used for ODA/MultiODA parity tests (UniODA, categorical, directional).

### Protein fixture (Chapter 3 multiclass example)

Protein fixture provenance: derived from a previously published dataset;
publication PMID 6643432. Used only for fixture parity/regression testing.
