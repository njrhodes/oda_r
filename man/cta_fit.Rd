\name{cta_fit}
\alias{cta_fit}
\title{Fit a Classification Tree Analysis (CTA) model (public wrapper)}
\description{
Public entry point for CTA. Currently supports binary (two-class) outcome
variables only.

When \code{recursive = FALSE} (default), validates the class variable and
delegates to \code{\link{oda_cta_fit}}.  When \code{recursive = TRUE},
runs the Locally Optimal Recursive Tree (LORT) engine: at each endpoint a full MDSA
family scan (\code{\link{cta_descendant_family}}) is performed, the min-D
member is selected, and recursion continues until no further structure is
found or a compute guard fires.  Returns a dual-tagged
\code{cta_ort} / \code{cta_tree} object.
}
\usage{
cta_fit(X, y, verbose = FALSE,
        recursive        = FALSE,
        min_n            = 30L,
        max_depth        = 8L,
        max_nodes        = 31L,
        family_max_steps = 20L,
        ...)
}
\arguments{
  \item{X}{Data frame or matrix of attribute columns.  For recursive CTA,
    \code{X} is the declared candidate predictor frame; pass only variables
    eligible for CTA search.  Prediction may be performed on wider
    \code{newdata} as long as the split variable names are present.}
  \item{y}{Integer class variable vector. Must have exactly two distinct values.}
  \item{verbose}{Logical; if \code{TRUE}, emit \code{[CTA]} and \code{[ORT]}
    progress messages. Default \code{FALSE}.}
  \item{recursive}{Logical; if \code{TRUE}, run the Locally Optimal Recursive Tree (LORT)
    engine.  Default \code{FALSE}.  Cannot be combined with an explicit
    \code{mindenom} argument (error).}
  \item{min_n}{Integer; minimum endpoint n to attempt recursion.  Endpoints
    smaller than \code{min_n} become terminal with stop reason
    \code{"min_n"}.  Default \code{30L}.  Only used when
    \code{recursive = TRUE}.}
  \item{max_depth}{Integer; safety cap on recursion depth.  Nodes at
    \code{depth >= max_depth} become terminal with stop reason
    \code{"max_depth"}.  Default \code{8L}.  Only used when
    \code{recursive = TRUE}.}
  \item{max_nodes}{Integer; safety cap on total ORT nodes allocated.  When
    the node count exceeds \code{max_nodes} the current endpoint becomes
    terminal with stop reason \code{"max_nodes"}.  Default \code{31L}.
    Only used when \code{recursive = TRUE}.}
  \item{family_max_steps}{Integer or \code{NULL} (default \code{20L}); maximum
    number of MDSA family members evaluated at each recursive node.  This bounds
    the per-node \code{\link{cta_descendant_family}} scan (its \code{max_steps}
    argument).  The default of \code{20L} preserves
    \code{cta_descendant_family()} behavior.  Smaller values reduce the per-node
    compute budget at the cost of possibly missing the true min-D member if the
    family is long.  Stored in \code{ort_settings$family_max_steps}.  Only used
    when \code{recursive = TRUE}; error if explicitly supplied with
    \code{recursive = FALSE}.}
  \item{...}{Additional arguments passed to \code{\link{oda_cta_fit}}
    (non-recursive) or used as ORT settings (recursive).  Supported ORT
    settings via \code{...}: \code{w}, \code{mc_seed}, \code{mc_iter},
    \code{alpha_split}, \code{prune_alpha}, \code{loo}.  When
    \code{recursive = TRUE}, \code{mc_seed} initializes the RNG once at
    ORT start; child-node MDSA scans consume the stream in deterministic
    right-then-left traversal order without resetting the seed.}
}
\value{
Non-recursive: a \code{cta_tree} object.

Recursive: a dual-tagged \code{cta_ort} / \code{cta_tree} object.
All existing \code{cta_tree} S3 methods (\code{predict}, \code{print},
\code{summary}, \code{plot}) operate on the root-level model.
\code{cta_ort}-aware methods (\code{predict.cta_ort},
\code{print.cta_ort}, \code{summary.cta_ort}, \code{plot.cta_ort}) operate
on the full composite tree.  Use \code{predict(obj, newdata, type="all")}
to retrieve stratum assignments.
}
\seealso{
\code{\link{oda_fit}}, \code{\link{cta_descendant_family}},
\code{\link{cta_node_table}}, \code{\link{cta_staging_table}},
\code{\link{plot.cta_tree}}, \code{\link{plot.cta_ort}},
\code{\link{ort_plot_data}}
}
\note{
\code{oda_cta_fit()} is the internal engine name; \code{cta_fit()} is the
preferred public entry point for non-recursive CTA.  Both are exported and
functionally equivalent for non-recursive use.

\code{cta_fit(..., recursive = TRUE)} is a legacy-compatible interface for
the LORT workflow layer.  Prefer \code{\link{lort_fit}()} for new code.
SORT and GORT are reserved and not implemented.
}
\examples{
# Small synthetic two-class example (non-recursive)
X <- data.frame(
  x1 = c(1, 2, 3, 4, 5, 6, 7, 8),
  x2 = c(0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L)
)
y <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L)

tree <- cta_fit(X, y,
  priors_on   = TRUE,
  mindenom    = 1L,
  mc_iter     = 500L,
  mc_seed     = 42L,
  loo         = "off",
  attr_names  = c("x1", "x2")
)
print(tree)

# Recursive ORT — two-level synthetic dataset
X2 <- data.frame(
  A = c(rep(0, 20), rep(1, 20), rep(1, 20)),
  B = c(rep(0, 20), rep(0, 20), rep(1, 20))
)
y2 <- c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
ort <- cta_fit(X2, y2, recursive = TRUE,
               mc_iter = 100L, mc_seed = 42L, loo = "off",
               min_n = 5L)
print(ort)
}
