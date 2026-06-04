\name{cta_endpoint_counts}
\alias{cta_endpoint_counts}
\title{Per-endpoint class count table for a fitted CTA tree}
\description{
Returns one row per terminal endpoint (leaf) per actual class, read
directly from stored leaf node fields.  No refitting, no prediction,
and no recomputation from training data is performed.

Class counts are stored at fit time by \code{\link{oda_cta_fit}} on
every terminal leaf.  Row order within each endpoint follows the order
of \code{names(leaf$class_counts_raw)}, which is ascending by class
label.  Endpoints are ordered by \code{node_id}, matching
\code{\link{cta_endpoint_summary}}.

\strong{Scope:} This function exposes stored raw and weighted class
counts only.  It does \emph{not} include target-class proportions,
event rates, odds, or staging order.  Staging-table and event-rate
summaries are available via \code{\link{cta_staging_table}}.

If any terminal leaf is missing the stored class counts (i.e., the
\code{cta_tree} was fitted by an earlier version of oda that did
not store endpoint counts), the function stops with a clear error.
}
\usage{cta_endpoint_counts(tree)}
\arguments{
  \item{tree}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
}
\value{
A \code{data.frame} with one row per terminal endpoint per actual class
and columns:
\describe{
  \item{\code{endpoint_id}}{Integer sequential endpoint index 1..n in
    node order, matching \code{\link{cta_endpoint_summary}}.}
  \item{\code{endpoint_node_id}}{Integer tree node identifier for this
    endpoint leaf.}
  \item{\code{path}}{Character; AND-joined branch labels from root to
    this leaf (e.g. \code{"V14<=0.5 AND V15>0.5"}).}
  \item{\code{terminal_prediction}}{Integer class label assigned to
    this endpoint (stored leaf \code{majority_class}).}
  \item{\code{class}}{Character; actual class label for this row
    (e.g. \code{"0"}, \code{"1"}).}
  \item{\code{n_raw}}{Integer raw count of observations of this actual
    class reaching this endpoint.}
  \item{\code{n_weighted}}{Numeric weighted total for this actual class
    reaching this endpoint.  Equals \code{n_raw} when case weights are
    not active.}
}
For a no-tree fit the returned data frame has zero rows but the correct
column structure and types.
}
\seealso{
  \code{\link{oda_cta_fit}}, \code{\link{cta_endpoint_summary}},
  \code{\link{cta_confusion_table}}, \code{\link{cta_endpoint_table}},
  \code{\link{cta_staging_table}}, \code{\link{cta_propensity_weights}}
}
\examples{
data(mtcars)
X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
y    <- as.integer(mtcars$am)
tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
cta_endpoint_counts(tree)
}
