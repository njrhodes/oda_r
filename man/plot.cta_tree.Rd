\name{plot.cta_tree}
\alias{plot.cta_tree}
\title{Plot a fitted CTA tree}
\description{
Produces a base-R tree diagram.  Calls \code{\link{cta_plot_data}} for
layout; uses only base graphics — no external dependencies.

Split (internal) nodes show the split attribute, node-level ESS or WESS,
and observation count.  Leaf (terminal) nodes show the majority-class
prediction and observation count.  Edge labels show the branch condition
(e.g. \code{"V14<=0.5"}).
}
\usage{
\method{plot}{cta_tree}(x, main = "CTA Tree",
     node_col_split = "#D9EAF7", node_col_leaf = "#D9F7E6",
     edge_col = "grey40", cex = 0.75, \dots)
}
\arguments{
  \item{x}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
  \item{main}{Character plot title.  Default \code{"CTA Tree"}.}
  \item{node_col_split}{Fill colour for split (internal) nodes.
    Default \code{"#D9EAF7"} (light blue).}
  \item{node_col_leaf}{Fill colour for leaf (terminal) nodes.
    Default \code{"#D9F7E6"} (light green).}
  \item{edge_col}{Colour for edge lines.  Default \code{"grey40"}.}
  \item{cex}{Text expansion factor for node labels.  Default \code{0.75}.}
  \item{\dots}{Unused; included for S3 compatibility.}
}
\value{
\code{invisible(x)}.
}
\seealso{
\code{\link{cta_plot_data}}, \code{\link{oda_cta_fit}}
}
\examples{
data(mtcars)
X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
y    <- as.integer(mtcars$am)
tree <- suppressMessages(
  oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L,
              loo = "off")
)
plot(tree)
}
