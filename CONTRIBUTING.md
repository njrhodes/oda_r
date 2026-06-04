# Contributing to oda

Thank you for your interest in contributing. This document covers environment
setup on Linux, macOS, and Windows, running tests, and the minimal smoke check
required before opening a pull request.

---

## Requirements

- R >= 4.1.0
- The following R packages (all in `Suggests`): `devtools`, `testthat`, `pak`
- For graphics: `ggplot2 >= 3.4.0`, `patchwork >= 1.1.0`
- Git

---

## Clone

```bash
git clone git@github.com:njrhodes/oda_r.git
cd oda_r
```

---

## Install dependencies

### macOS / Linux

```r
install.packages("pak")
pak::local_install_deps(dependencies = TRUE)
```

### Windows

`pak` is preferred; fall back to `devtools` if needed:

```r
install.packages("pak")
pak::local_install_deps(dependencies = TRUE)

# Alternative if pak fails:
# install.packages("devtools")
# devtools::install_deps(dependencies = TRUE)
```

---

## Load the package

```r
devtools::load_all()
```

---

## Running tests

Tests are gated by the `ODA_TEST_TIER` environment variable.

| Tier | Scope | When to use |
|------|-------|-------------|
| unset / `cran` | CRAN-safe unit tests only | Default |
| `fast` | Developer loop; skips slow canon fixtures | Most development |
| `smoke` | Production gate; includes myeloma/CTA fixtures | Before any commit touching CTA/MDSA/graphics |
| `full` | Release gate; all canon tests | Before release tags only |

### macOS / Linux

```bash
# Fast dev loop
ODA_TEST_TIER=fast Rscript --vanilla -e "devtools::test(reporter='progress')"

# Full smoke gate
ODA_TEST_TIER=smoke Rscript --vanilla -e "devtools::test(reporter='progress')"
```

### Windows (PowerShell)

```powershell
# Fast dev loop
$env:ODA_TEST_TIER = "fast"
Rscript --vanilla -e "devtools::test(reporter='progress')"

# Full smoke gate
$env:ODA_TEST_TIER = "smoke"
Rscript --vanilla -e "devtools::test(reporter='progress')"
```

### Windows (Command Prompt)

```cmd
set ODA_TEST_TIER=fast
Rscript --vanilla -e "devtools::test(reporter='progress')"
```

### Windows (R console)

```r
Sys.setenv(ODA_TEST_TIER = "fast")
devtools::test(reporter = "progress")
```

---

## Minimal pre-PR smoke check

Run this before opening any pull request:

```bash
# macOS / Linux
ODA_TEST_TIER=fast Rscript --vanilla -e "devtools::test(reporter='progress')"
Rscript --vanilla -e "devtools::check(vignettes=FALSE)"
```

```powershell
# Windows PowerShell
$env:ODA_TEST_TIER = "fast"
Rscript --vanilla -e "devtools::test(reporter='progress')"
Rscript --vanilla -e "devtools::check(vignettes=FALSE)"
```

Expected result:

```
[ FAIL 0 | WARN 0 | SKIP 34 | PASS 3000+ ]
0 errors / 0 warnings / at most 1 note
```

The 34 skips are intentional: slow canon fixtures (myeloma, CTA_DEMO, graphics)
that require `ODA_TEST_TIER=smoke` or `full` to run. The timestamp note is a
filesystem artifact on cloud-synced or network-mounted paths (e.g., Dropbox)
and is not a package failure.

---

## Pull request checklist

- [ ] `ODA_TEST_TIER=fast` suite: FAIL 0 / WARN 0
- [ ] `devtools::check(vignettes=FALSE)`: 0 errors / 0 warnings
- [ ] If touching CTA, MDSA, reporting, or graphics: `ODA_TEST_TIER=smoke` suite passes
- [ ] No non-ASCII characters in R source, tests, or docs
- [ ] No fit-time storage added to `cta_tree` (lean-fit invariant)
- [ ] New public functions exported in `NAMESPACE` (hand-maintained)

---

## Architecture notes

See `CLAUDE.md` for engine internals, canon fixture locations, and terminology
rules. See `docs/STATUS.md` for what is implemented and what is deferred.
