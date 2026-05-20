\name{cta_endpoint_table}
\alias{cta_endpoint_table}
\title{Terminal-leaf report table for a fitted CTA tree}
\description{
Returns one row per terminal leaf (endpoint) of a \code{cta_tree}.  All
values are read directly from stored node fields; no refitting or prediction
is performed.

Per-endpoint class-count and target-class summaries are deferred to the
future CTA confusion/endpoint-count table because leaf class counts are not
currently stored at fit time.
}
\usage{cta_endpoint_table(tree)}
\arguments{
  \item{tree}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
}
\value{
A \code{data.frame} with one row per terminal leaf and columns:
\describe{
  \item{\code{node_id}}{Integer node identifier.}
  \item{\code{parent_id}}{Integer parent node identifier
    (\code{0} for the root leaf in a no-tree fit).}
  \item{\code{depth}}{Integer depth from root (root = 1).}
  \item{\code{path}}{Character; AND-joined branch labels from root to this
    leaf (e.g. \code{"V14<=0.5 AND V15>0.5"}).  \code{"/"} for the root
    leaf.}
  \item{\code{majority_class}}{Integer class label assigned to this leaf.}
  \item{\code{n_obs}}{Integer raw observation count at this leaf.}
  \item{\code{n_weighted}}{Numeric weighted observation count
    (\code{NA} when weights not active).}
  \item{\code{ess}}{Numeric ESS stored on this leaf (\code{NA_real_} —
    ESS is a split-node metric and is always \code{NA} for leaf nodes).}
  \item{\code{ess_weighted}}{Numeric WESS stored on this leaf
    (\code{NA_real_}).}
  \item{\code{loo_status}}{Character LOO status (\code{NA} for leaf nodes).}
  \item{\code{loo_ess}}{Numeric LOO ESS (\code{NA} for leaf nodes).}
}
For a no-tree fit the returned data frame has zero rows but the correct
column structure.
}
\seealso{
  \code{\link{oda_cta_fit}}, \code{\link{summary.cta_tree}},
  \code{\link{cta_strata}}, \code{\link{cta_endpoint_denominators}}
}
\examples{
data(mtcars)
X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
y    <- as.integer(mtcars$am)
tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
cta_endpoint_table(tree)
}
