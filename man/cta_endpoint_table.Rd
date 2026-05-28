\name{cta_endpoint_table}
\alias{cta_endpoint_table}
\title{Canonical terminal endpoint map for a fitted CTA tree}
\description{
Returns one row per terminal leaf (endpoint) of a \code{cta_tree}.  All
values are read directly from stored node fields; no refitting or prediction
is performed.  This is the canonical endpoint map for reporting, translation,
ORT, and staged workflows.

Leaf class counts are stored on every terminal node at fit time
(\code{class_counts_raw}, \code{class_counts_weighted}).  \code{target_n}
and \code{target_prop} are derived from the stored counts.

ESS, WESS, p, LOO status, LOO ESS/WESSL, and LOOp are canonical split-node
report metrics (see \code{\link{cta_node_table}}).  Terminal endpoints are
connected to those metrics through their parent split-node lineage.  The
\code{parent_split_*} columns expose the immediate parent split's canonical
metrics for auditability.  They are not recomputed ESS at the leaf.
}
\usage{cta_endpoint_table(tree, target_class = NULL)}
\arguments{
  \item{tree}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
  \item{target_class}{Integer class label to use as the target (positive)
    class for \code{target_n} and \code{target_prop}.  When \code{NULL}
    (default), the function auto-detects: for binary trees with classes
    \code{0} and \code{1}, class \code{1} is used; otherwise \code{target_n}
    and \code{target_prop} are \code{NA}.}
}
\value{
A \code{data.frame} with one row per terminal leaf and columns:
\describe{
  \item{\code{endpoint_id}}{Integer sequential endpoint index 1..n.}
  \item{\code{leaf_node_id}}{Integer tree node identifier for this leaf.}
  \item{\code{terminal_marker}}{Character \code{"*"} on every row.}
  \item{\code{terminal}}{Logical \code{TRUE} on every row.}
  \item{\code{depth}}{Integer depth from root (root = 1).}
  \item{\code{parent_split_node_id}}{Integer parent split node identifier.}
  \item{\code{path}}{Character; AND-joined branch labels from root to this
    leaf (e.g. \code{"V14<=0.5 AND V15>0.5"}).}
  \item{\code{n}}{Integer raw observation count at this endpoint.}
  \item{\code{class_counts_raw}}{List column; each element is a named
    integer vector of raw per-class counts, or \code{NULL}.}
  \item{\code{class_counts_weighted}}{List column; each element is a named
    numeric vector of weighted per-class counts, or \code{NULL}.}
  \item{\code{predicted_class}}{Integer class label assigned to this
    endpoint (stored leaf majority class).}
  \item{\code{target_n}}{Integer count of \code{target_class} observations
    at this endpoint (\code{NA} when not resolvable).}
  \item{\code{target_prop}}{Numeric proportion \code{target_n / n}
    (\code{NA} when not resolvable).}
  \item{\code{parent_split_attribute}}{Attribute name of the parent split.}
  \item{\code{parent_split_ess}}{ESS of the parent split node.}
  \item{\code{parent_split_wess}}{WESS of the parent split node.}
  \item{\code{parent_split_loo_status}}{LOO status of the parent split node.}
  \item{\code{parent_split_loo_ess}}{LOO ESS/WESSL of the parent split node.}
  \item{\code{parent_split_p_mc}}{MC p-value of the parent split node.}
}
For a no-tree fit the returned data frame has zero rows but the correct
column structure and types.
}
\seealso{
  \code{\link{oda_cta_fit}}, \code{\link{cta_node_table}},
  \code{\link{summary.cta_tree}}, \code{\link{cta_strata}},
  \code{\link{cta_endpoint_denominators}}, \code{\link{cta_endpoint_summary}},
  \code{\link{cta_endpoint_counts}}
}
\examples{
data(mtcars)
X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
y    <- as.integer(mtcars$am)
tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
cta_endpoint_table(tree)
}
