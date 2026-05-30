# GRAPHICS_V3.md

Graphics v3  -  direct ggplot2 renderers for CTA trees, LORT trees, and
covariate balance diagnostics.

**Status:** v3C1 (tree renderers) and v3C2 (balance renderers) landed.

---

## Overview

Graphics v3 adds six direct ggplot2 rendering functions as an alternative
to the base-R `plot.cta_tree()` / `plot.cta_ort()` methods.  All functions:

- Return `ggplot` objects (printable, modifiable, `ggsave()`-compatible).
- Require `ggplot2 >= 3.4.0` (in `Suggests`; a clear error is raised if absent).
- Accept only pre-fitted objects or pre-computed plot-data  -  no model fitting
  inside renderers.
- Coexist with base `plot.*` methods; neither replaces the other.

---

## Function reference

### Tree renderers (v3C1)

| Function | Input | Description |
|----------|-------|-------------|
| `plot_cta_tree(x, ...)` | `cta_tree` or `cta_plot_data()` output | CTA tree diagram |
| `plot_lort_tree(x, ...)` | `cta_ort` or `ort_plot_data()` output | LORT tree diagram |

Both accept a `color_by` argument:

| Value | Leaf fill |
|-------|-----------|
| `"target_rate"` (default) | Continuous gradient by target-class proportion (0-100%) |
| `"prediction"` | Discrete fill by predicted class |
| `"none"` | White fill (no legend) |

### Balance renderers (v3C2)

| Function | Input | Description |
|----------|-------|-------------|
| `plot_oda_balance(x, ...)` | `oda_balance_plot_data` or `oda_balance_table` | ODA ESS/WESS dot plot |
| `plot_smd_balance(x, ...)` | `smd_balance_table` | Absolute SMD dot plot with threshold lines |
| `plot_balance_love(x, ...)` | `smd_balance_table` | Love plot (direct alias for `plot_smd_balance`) |
| `plot_cta_balance(x, ...)` | `cta_balance_plot_data` or `cta_balance_table` | CTA balance tree or no-tree message |

---

## ggplot2 dependency policy

ggplot2 is listed in `Suggests`, not `Imports`.  The package installs and
loads without ggplot2.  All Graphics v3 functions call an internal
`.require_ggplot2()` guard that raises a clear error if ggplot2 is absent:

```r
plot_cta_tree(tree)
# Error: Install ggplot2 to use Graphics v3 plotting functions:
#   install.packages("ggplot2")
```

Base-R plot methods (`plot.cta_tree`, `plot.cta_ort`) do not require ggplot2
and remain fully functional regardless.

---

## No fitting inside renderers

All six functions are pure renderers.  They read pre-computed objects and
build ggplot layers.  The only computation is coordinate adjustment and
aesthetic mapping.

**Coercion rules (convenience only  -  no fitting):**

- `plot_oda_balance(x)`  -  accepts `oda_balance_plot_data` (preferred) or
  `oda_balance_table` (coerced via `oda_balance_plot_data()` before rendering;
  never calls `oda_balance_table()` inside the renderer).  SMD is not
  recomputed; only what is already in the plot-data is shown.
- `plot_cta_balance(x)`  -  accepts `cta_balance_plot_data` (preferred) or
  `cta_balance_table` (coerced via `cta_balance_plot_data()`; never calls
  `cta_balance_table()` inside the renderer).

Attempting to pass raw data (`group`, `X`, `y`) directly to any renderer
will error.

---

## Examples

### CTA tree

```r
library(odacore)

X <- data.frame(
  A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
  B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
)
y <- c(rep(0L, 40), rep(1L, 20))

tree <- cta_fit(X, y, mindenom = 5L, mc_iter = 200L, mc_seed = 42L, loo = "off")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  # Default: leaf fill = target-rate gradient
  p <- plot_cta_tree(tree, color_by = "target_rate")
  print(p)

  # Discrete prediction fill
  p2 <- plot_cta_tree(tree, color_by = "prediction")

  # Save to file
  ggplot2::ggsave("cta_tree.png", p, width = 8, height = 5, dpi = 300)
}
```

### LORT tree

```r
lort <- lort_fit(X, y, mc_iter = 200L, mc_seed = 42L, loo = "off", min_n = 5L)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- plot_lort_tree(lort, color_by = "prediction")
  print(p)

  ggplot2::ggsave("lort_tree.png", p, width = 8, height = 5, dpi = 300)
}
```

### ODA balance plot

`plot_oda_balance()` plots ESS or WESS per covariate.  Point color reflects
significance status.  Only what is in the pre-computed `oda_balance_plot_data`
object is rendered; no ODA models are fit inside the renderer.

```r
group <- c(rep(0L, 30), rep(1L, 30))
X_bl  <- data.frame(
  age   = c(rep(40L, 30), rep(60L, 30)),   # imbalanced
  score = rnorm(60, 50, 10)                 # balanced
)

bt <- oda_balance_table(group, X_bl, mcarlo = TRUE, mc_iter = 500L)
pd <- oda_balance_plot_data(bt)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- plot_oda_balance(pd)
  print(p)

  ggplot2::ggsave("oda_balance.png", p, width = 7, height = 4, dpi = 300)
}
```

### SMD balance / Love plot

`plot_smd_balance()` plots |SMD| per covariate with reference lines at 0.10
(and optionally 0.20).  `plot_balance_love()` is a direct alias.

```r
smd <- smd_balance_table(group, X_bl)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  # SMD dot plot with both reference lines
  p_smd  <- plot_smd_balance(smd, ref_010 = TRUE, ref_020 = TRUE)
  print(p_smd)

  # Love plot (identical output, conventional name)
  p_love <- plot_balance_love(smd)
  print(p_love)

  ggplot2::ggsave("love_plot.png", p_love, width = 7, height = 4, dpi = 300)
}
```

### CTA balance

`plot_cta_balance()` renders a tree diagram when a discriminating tree was
found, or a message panel when `status = "no_tree"`.

```r
ct  <- cta_balance_table(group, X_bl, mindenom = 5L,
                          mc_iter = 500L, mc_seed = 42L)
cpd <- cta_balance_plot_data(ct)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- plot_cta_balance(cpd)
  print(p)

  ggplot2::ggsave("cta_balance.png", p, width = 8, height = 5, dpi = 300)
}
```

---

## Interpreting CTA balance no_tree results

When `status = "no_tree"`, `plot_cta_balance()` renders a message panel
that reads:

> No combination of covariates predicted group membership
> (MINDENOM = \<n\>, alpha = \<alpha\>).
> This is favorable evidence of multivariable balance
> under the declared constraints.

**What this means:**

- No combination of baseline covariates in `X` predicted group membership
  above the declared significance threshold (`alpha_split`, default 0.05) with
  the required minimum endpoint sample size (`mindenom`).
- In balance analysis, the inability to discriminate groups **is the desired
  result**.  It is favorable evidence of multivariable covariate balance under
  the declared analytic constraints.

**What this does NOT mean:**

- It does not prove universal balance.  Unobserved confounders, covariates
  not included in `X`, or constraints that are too lenient may allow imbalance
  to exist undetected.
- It is not a model failure.  `no_tree` in the balance context is the
  equivalent of "pass" on the balance diagnostic, not a fitting error.

---

## Base plot methods (unchanged)

The base-R `plot.cta_tree()` and `plot.cta_ort()` dispatch methods are
unchanged.  Graphics v3 does not modify or replace them.

```r
plot(tree)                     # base-R
plot(tree, target_class = 1L)  # base-R with target-class coloring

p <- plot_cta_tree(tree)       # ggplot2 v3  -  same tree, ggplot output
```

Both approaches read the same `cta_plot_data()` derived layout internally.

---

## Saving plots

All Graphics v3 functions return standard `ggplot` objects.  Use
`ggplot2::ggsave()`:

```r
# Raster output
ggplot2::ggsave("tree.png",    p, width = 8, height = 5, dpi = 300)

# Vector output
ggplot2::ggsave("tree.pdf",    p, width = 10, height = 7)
ggplot2::ggsave("tree.svg",    p, width = 8,  height = 5)
```

Since the return value is a ggplot, any standard ggplot2 customization layer
can be added before saving:

```r
p <- plot_cta_tree(tree) +
  ggplot2::labs(caption = "Source: odacore v0.1.0")
```

---

## Implementation notes

- `R/graphics_v3.R` contains all six renderers plus internal helpers
  (`.require_ggplot2()`, `.render_tree_gg()`, `.balance_pal()`, etc.).
- `R/balance.R` contains the balance table / plot-data pipeline.
- All six functions use `ggplot2::` namespace prefix; ggplot2 is never attached.
- `utils::globalVariables()` suppresses R CMD CHECK notes for `aes()` column
  name lookups.
- The NAMESPACE is hand-maintained; exports are added explicitly.
