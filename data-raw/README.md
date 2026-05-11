# data-raw/

This directory is reserved for raw data preparation scripts. None exist yet —
all fixture data lives directly under `tests/testthat/fixtures/`.

## Fixture provenance

### UniODA / MultiODA gold values (`tests/testthat/test-iris.R`, `test-unioda.R`, etc.)

Confusion matrices and ESS values were extracted verbatim from MegaODA.exe
output (2025-04-06):

    MCARLO ITER 25000; LOO; GO;
    Priors ON | Degen OFF | Primary MAXSENS | Secondary SAMPLEREP

### CTA gold fixtures (`tests/testthat/fixtures/`)

CTA fixture outputs were produced by CTA.exe and are stored verbatim under
`tests/testthat/fixtures/<dataset>/`:

| Dataset  | MINDENOM | Key output file          | Notes                              |
|----------|----------|--------------------------|------------------------------------|
| cta_demo | 1        | MODEL1.TXT               | No WEIGHT command → ESS            |
| cta_demo | 8        | MODEL8.TXT               | mc_iter=25000 required for parity  |
| myeloma  | 1        | MODEL1.TXT               | WEIGHT V2 active → WESS            |
| myeloma  | 30       | MODEL30.TXT              | V17 stump, path-local n=186        |
| myeloma  | 56       | MODEL56.TXT              | No tree (all child sizes < 56)     |

CTA.exe command files (`.pgm`) and raw console output (`.txt`) are stored
alongside each MODEL*.TXT for full reproducibility. Do not adjust fixture
expectations without re-running the corresponding CTA.exe command file.
