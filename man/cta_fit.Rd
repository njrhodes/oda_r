\name{cta_fit}
\alias{cta_fit}
\title{Fit a Classification Tree Analysis (CTA) model (public wrapper)}
\description{
Public entry point for CTA. Currently supports binary (two-class) outcome
variables only.

When \code{recursive = FALSE} (default), validates the class variable and
delegates to \code{\link{oda_cta_fit}}.  When \code{recursive = TRUE},
runs the Optimal Recursive Tree (ORT) engine: at each endpoint a full MDSA
family scan (\code{\link{cta_descendant_family}}) is performed, the min-D
member is selected, and recursion continues until no further structure is
found or a compute guard fires.  Returns a dual-tagged
\code{cta_ort} / \code{cta_tree} object.
}
\usage{
cta_fit(X, y, verbose = FALSE,
        recursive = FALSE,
        min_n     = 30L,
        max_depth = 8L,
        max_nodes = 31L,
        ...)
}
\arguments{
  \item{X}{Data frame or matrix of attribute columns.}
  \item{y}{Integer class variable vector. Must have exactly two distinct values.}
  \item{verbose}{Logical; if \code{TRUE}, emit \code{[CTA]} and \code{[ORT]}
    progress messages. Default \code{FALSE}.}
  \item{recursive}{Logical; if \code{TRUE}, run the Optimal Recursive Tree
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
  \item{...}{Additional arguments passed to \code{\link{oda_cta_fit}}
    (non-recursive) or used as ORT settings (recursive).  Supported ORT
    settings via \code{...}: \code{w}, \code{mc_seed}, \code{mc_iter},
    \code{alpha_split}, \code{prune_alpha}, \code{loo}.}
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
preferred public entry point.  Both are exported and functionally equivalent
for non-recursive use.
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
