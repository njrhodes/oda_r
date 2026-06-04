# GRAPHICS_V3.md

Graphics v3  -  direct ggplot2 renderers for CTA trees, LORT trees, and
covariate balance diagnostics.

**Status:** v3C1 (tree renderers) and v3C2 (balance renderers) landed.

---

## Graphics v4 visual remediation targets

Visual reference: ODA Journal articles and MPE book (Yarnold & Soltysik 2016)
chapter figures. Target: figures that could appear in a journal submission
without post-processing.

### CTA tree - publication standard

**Split nodes** (internal):
- Content: attribute name, cut value, n at node, ESS or WESS
- Label example: `x3\n<= 0.5  >0.5\nn=96, ESS=64.2%`
- Border: solid, medium weight; fill: light blue

**Terminal nodes** (endpoints/leaves):
- Content: Stage ID, predicted class, n classified, target-class rate%
- Label example: `Stage 1\nClass 0\nn=24, 17%`
- Fill: gradient by target rate (low = light, high = dark red); or discrete by prediction

**Edges**:
- Label: **short only** - just the branch condition at THIS split, not the full path
- Correct: `<= 0.5` or `> 0.5`
- Incorrect: `x3 <= 0.5` (attribute already shown in parent) or compound `x3>0.5 AND x1>25` (path belongs in companion table)
- The full root-to-leaf path logic is the job of a companion table (e.g., `cta_endpoint_table()`), not the diagram edges

**Title**:
- Format: `CTA - n=N, K endpoints, ESS=XX.X%, D=X.XXXX`
- The D-statistic and global ESS/WESS belong in the title, not buried in nodes

**v3 failures corrected in v4:**
1. Fixed `half_h = 0.28` - cramped boxes for 4-line terminal labels -> replaced with adaptive height from max line count
2. No auto-title with ESS/D/n -> added default auto-title
3. Edge labels included attribute name prefix (`V14<=0.5`) -> now stripped to `<=0.5` by default
4. No stage IDs on terminal nodes -> `cta_plot_data()` already adds Stage N; renderer exposes it

### LORT tree - publication standard

LORT is a **sequence of CTA sub-trees**, one per stratum. It is NOT a single
composite CTA tree. The visual must reflect this:

**Layout requirements:**
- Horizontal depth bands (shaded alternating stripes) separate LORT levels
- Terminal LORT nodes (final strata with no further sub-tree) have a dashed border to distinguish them from CTA internal-split nodes
- Stage number is prominent on each terminal stratum block
- Depth-level header annotations ("Level 1", "Level 2", ...) at the left margin

**Edge labels:**
- Each edge connects a stratum block to its child stratum (through a CTA endpoint)
- Label: the branch condition at the CTA endpoint split only - NOT the full path
- `path_conditions` is a vector; the last element may be a compound string from multi-split CTA sub-trees; strip to the final `AND` segment

**v3 failures corrected in v4:**
1. Used same `.render_tree_gg()` as CTA - LORT looked like a malformed CTA plot -> now uses LORT-specific visual differentiation (depth bands + dashed terminal borders)
2. Edge labels were compound paths (`V17>0.5 AND V15>0.5`) from multi-level sub-trees -> fixed in `ort_plot_data()` to extract final segment

### Balance plots - publication standard

**SMD lollipop / Love plot** (`plot_smd_balance`, `plot_balance_love`):
- Had the right diagnostic intent in v3 but failed public visual review
- v4 targets: typography hierarchy (attribute labels larger, axis label smaller), threshold line must carry an inline annotation (`"|SMD| = 0.10"`) not rely on the legend alone, legend must be compact and right-positioned (not dominating), point size increased for journal print, balanced/imbalanced color contrast improved

**ODA ESS balance** (`plot_oda_balance`):
- Dot plot of ESS per attribute, colored by significance
- v4: same as v3, no major change needed

**CTA balance no-tree** (`plot_cta_balance`, `status = "no_tree"`):
- v3 failure: `.gg_message_panel()` produced a blank void canvas with floating text - looked like an error screen or software crash
- v4 target: styled result card with:
  - Bordered box panel (not void canvas)
  - Clear positive header: "No Discriminating Tree Found"
  - "This is favorable evidence of multivariable covariate balance"
  - Footer with constraint parameters (MINDENOM, alpha)
  - Color: success green tones, NOT error red

### Figure + table pairing rule (v4 mandate)

Every public graphic in the practitioner guide must be paired with a companion
table. The figure answers "what did the model find?" visually; the table provides
the quantitative evidence. Never overload a single plot.

| Figure | Companion table |
|--------|-----------------|
| `plot_cta_tree()` | `cta_endpoint_table()` or `cta_staging_table()` |
| `plot_lort_tree()` | `cta_ort_node_table()` |
| `plot_oda_balance()` | `oda_balance_table()` output |
| `plot_smd_balance()` | `smd_balance_table()` output |
| `plot_cta_balance()` | `cta_balance_table()` output |

### Deterministic guide dataset (v4 mandate)

Replace rnorm-based random `dat` in guide with a fully deterministic 4-block dataset:

```r
# 96 obs x 4 blocks x 24, no randomness, exact 4-stage CTA tree structure
# x1 (binary, root split): x1=15 for obs 1-24 and 49-72, x1=35 for obs 25-48 and 73-96
# x3 (binary, child split): x3=0 for obs 1-48, x3=1 for obs 49-96
# y (binary outcome): controlled target rates per (x1,x3) block
# x2 (ordered, non-informative): 6-cycle, breaks row-order correlation
dat_guide <- local({
  x3    <- c(rep(0L, 48), rep(1L, 48))
  x1    <- c(rep(15L, 24), rep(35L, 24), rep(15L, 24), rep(35L, 24))
  y     <- c(
    rep(c(0L, 1L), times = c(22L, 2L)),   # Stage 1: x1=15,x3=0 ->  2/24 =  8.3%
    rep(c(0L, 1L), times = c(14L,10L)),   # Stage 2: x1=15,x3=1 -> 10/24 = 41.7%
    rep(c(0L, 1L), times = c(10L,14L)),   # Stage 3: x1=35,x3=0 -> 14/24 = 58.3%
    rep(c(0L, 1L), times = c(2L, 22L))    # Stage 4: x1=35,x3=1 -> 22/24 = 91.7%
  )
  x2    <- rep(c(25L, 30L, 35L, 40L, 45L, 50L), length.out = 96L)
  treat <- rep(c(0L, 1L), times = 48L)
  data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, treat = treat)
})
```

Expected CTA tree: root split on x1 (binary: {15} vs {35}), then x3 split (cut 0.5)
in each x1 branch, yielding 4 terminal endpoints (Stages 1-4) with rates
8%/42%/58%/92%, ESS ~= 50%, D = 4.0.

Note: x1 and x3 are both binary attributes. CTA uses binary_map rule for each.
x2 uses a 6-cycle pattern to break row-order correlation with y.

---

## Graphics v4 design decision

Deep and recursive tree structures (LORT, MDSA family) are displayed by indexed
sub-tree inspection and optional multipanel overview, not by forcing all logic into
one composite plot.  `plot_lort_tree(x, index=k)` renders the CTA sub-tree at LORT
node k.  `plot_cta_family(family, index=k)` renders the k-th family member.
`show_all=TRUE` returns a named list of ggplot objects.

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

### Tree renderers (v3C1 / v4)

| Function | Input | Description |
|----------|-------|-------------|
| `plot_cta_tree(x, ...)` | `cta_tree` or `cta_plot_data()` output | CTA tree diagram |
| `plot_lort_tree(x, index, show_all, ...)` | `cta_ort` | LORT indexed sub-tree |
| `plot_cta_family(family, index, show_all, ...)` | `cta_family` | MDSA family member |

All three accept a `color_by` argument:

| Value | Leaf fill |
|-------|-----------|
| `"none"` (default) | White fill, no legend (B/W publication default) |
| `"target_rate"` | Continuous gradient by target-class proportion (0-100%) |
| `"prediction"` | Discrete fill by predicted class |

`plot_lort_tree()` and `plot_cta_family()` render one member at a time (`index`)
or return a named list of ggplot objects (`show_all = TRUE`).

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
library(oda)

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
  ggplot2::labs(caption = "Source: oda v0.1.0")
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
