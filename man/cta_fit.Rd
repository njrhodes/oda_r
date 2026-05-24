\name{cta_fit}
\alias{cta_fit}
\title{Fit a Classification Tree Analysis (CTA) model (public wrapper)}
\description{
Public entry point for CTA. Currently supports binary (two-class) outcome
variables only. Validates the class variable and delegates to
\code{\link{oda_cta_fit}}.
}
\usage{
cta_fit(X, y, verbose = FALSE, ...)
}
\arguments{
  \item{X}{Data frame or matrix of attribute columns.}
  \item{y}{Integer class variable vector. Must have exactly two distinct values.}
  \item{verbose}{Logical; if \code{TRUE}, emit \code{[CTA]} progress messages
    via \code{message()} at each major stage (fast screen, MC+LOO fit,
    ENUMERATE scoring, final selection). Default \code{FALSE}.}
  \item{...}{Additional arguments passed to \code{\link{oda_cta_fit}}.}
}
\value{
A \code{cta_tree} object. Use \code{predict()} to classify new observations,
\code{print()} for a node summary, and \code{plot()} for a tree diagram.
}
\seealso{
\code{\link{oda_fit}}, \code{\link{cta_descendant_family}},
\code{\link{cta_node_table}}, \code{\link{cta_staging_table}},
\code{\link{plot.cta_tree}}
}
\note{
\code{oda_cta_fit()} is the internal engine name; \code{cta_fit()} is the
preferred public entry point.  Both are exported and functionally equivalent.
}
\examples{
# Small synthetic two-class example (no MC, no LOO for speed)
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
}
