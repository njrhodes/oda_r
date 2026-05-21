\name{cta_endpoint_summary}
\alias{cta_endpoint_summary}
\title{Endpoint reporting summary for a fitted CTA tree}
\description{
Returns one row per terminal leaf (endpoint) with stable endpoint identifiers
and stored node fields suitable for downstream reporting.  All values are read
directly from stored node fields; no refitting, no prediction, and no
recomputation of tree metrics is performed.

\strong{Scope:} This function reports structural endpoint fields only.
It does \emph{not} include endpoint class counts, target-class proportions,
event rates, odds, or staging order.  Staging-table and event-rate reporting
require per-endpoint class counts, which are not currently stored at fit time.
That capability is deferred until fit-time endpoint class-count storage is
designed.
}
\usage{cta_endpoint_summary(tree)}
\arguments{
  \item{tree}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
}
\value{
A \code{data.frame} with one row per terminal leaf and columns:
\describe{
  \item{\code{endpoint_id}}{Integer sequential index 1..n in node order.}
  \item{\code{endpoint_node_id}}{Integer tree node identifier for this leaf,
    corresponding to \code{node_id} in \code{\link{cta_endpoint_table}}.}
  \item{\code{path}}{Character; AND-joined branch labels from root to this
    leaf (e.g. \code{"V14<=0.5 AND V15>0.5"}).}
  \item{\code{depth}}{Integer depth from root (root = 1).}
  \item{\code{terminal_prediction}}{Integer class label assigned to this
    endpoint (stored leaf \code{majority_class}).}
  \item{\code{n_obs}}{Integer raw observation count at this endpoint.}
  \item{\code{n_weighted}}{Numeric weighted observation count.  Equals
    \code{n_obs} when case weights are not active (not \code{NA}).}
  \item{\code{denominator}}{Integer endpoint denominator (equal to
    \code{n_obs}); included to align with MPE/MDSA terminology.}
}
For a no-tree fit the returned data frame has zero rows but the correct
column structure and types.
}
\seealso{
  \code{\link{oda_cta_fit}}, \code{\link{cta_endpoint_table}},
  \code{\link{cta_strata}}, \code{\link{cta_endpoint_denominators}}
}
\examples{
data(mtcars)
X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
y    <- as.integer(mtcars$am)
tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
cta_endpoint_summary(tree)
}
