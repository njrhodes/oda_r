\name{cta_plot_data}
\alias{cta_plot_data}
\title{Extract layout data for plotting a CTA tree}
\description{
Returns a pure-data list describing tree topology and layout coordinates.
No graphics are produced.  Use this as input to \code{\link{plot.cta_tree}}
or to custom rendering code.

Layout algorithm: leaves receive sequential integer x-positions in
depth-first (left-to-right) order; internal nodes are centred over their
children.  \code{y = -depth} so the root sits at the top.
}
\usage{
cta_plot_data(tree)
}
\arguments{
  \item{tree}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
}
\value{
A list with elements:
\describe{
  \item{\code{nodes}}{A \code{data.frame} with one row per node and columns
    \code{node_id} (integer), \code{parent_id} (integer), \code{depth}
    (integer), \code{x} (numeric), \code{y} (numeric), \code{leaf}
    (logical), \code{attribute} (character; \code{NA} for leaves),
    \code{n_obs} (integer), \code{majority_class} (integer), \code{ess}
    (numeric; \code{NA} for leaves), \code{label} (character multi-line
    display text).}
  \item{\code{edges}}{A \code{data.frame} with one row per parent-to-child
    edge and columns \code{from_node_id} (integer), \code{to_node_id}
    (integer), \code{x0}, \code{y0}, \code{x1}, \code{y1} (numeric
    centre-to-centre coordinates), \code{label} (character branch condition,
    e.g. \code{"V14<=0.5"}).}
  \item{\code{no_tree}}{Logical; \code{TRUE} for leaf-only fits.}
  \item{\code{has_weights}}{Logical; \code{TRUE} when case weights are
    active.}
}
}
\seealso{
\code{\link{plot.cta_tree}}, \code{\link{oda_cta_fit}}
}
\examples{
data(mtcars)
X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
y    <- as.integer(mtcars$am)
tree <- suppressMessages(
  oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L,
              loo = "off")
)
pd <- cta_plot_data(tree)
head(pd$nodes)
pd$edges
}
