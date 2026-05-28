\name{cta_ort_node_table}
\alias{cta_ort_node_table}
\title{Node-level summary table for a fitted LORT (legacy name: cta_ort)}
\description{
Returns one row per LORT node from a \code{cta_ort} (LORT) object.  Each row
exposes the embedded CTA member selected at that node (MINDENOM, ESS, D,
root attribute, split/leaf counts, endpoint count) plus the LORT method
taxonomy metadata (\code{method}, \code{selection_scope},
\code{global_optimization}, \code{sda_anchored}).

Terminal nodes have \code{NA} for all selected-model columns.  Non-terminal
nodes have \code{NA} for \code{stop_reason} and non-empty \code{child_ids}.

\emph{Naming note:} The function name \code{cta_ort_node_table} and the class
\code{cta_ort} are legacy compatibility names for the implemented LORT method.
They are retained for backward compatibility; new code and documentation should
refer to the method as LORT.
}
\usage{cta_ort_node_table(object)}
\arguments{
  \item{object}{A \code{cta_ort} from \code{\link{cta_fit}}\code{(..., recursive = TRUE)}.}
}
\value{
A \code{data.frame} with one row per ORT node and columns:
\describe{
  \item{\code{ort_node_id}}{Integer ORT node identifier.}
  \item{\code{parent_ort_node_id}}{Integer parent ORT node id; \code{NA} for root.}
  \item{\code{depth}}{Integer recursion depth (root = 0).}
  \item{\code{n}}{Integer observations at this ORT node.}
  \item{\code{class_counts}}{Character; named class counts, e.g. \code{"0=60, 1=40"}.}
  \item{\code{terminal}}{Logical; \code{TRUE} for terminal leaf ORT nodes.}
  \item{\code{stop_reason}}{Character stop reason for terminal nodes;
    \code{NA} for non-terminal.}
  \item{\code{selected_mindenom}}{Integer MINDENOM of the embedded CTA member.}
  \item{\code{selected_ess}}{Numeric ESS of the embedded CTA member (\%).}
  \item{\code{selected_d}}{Numeric D-statistic of the embedded CTA member.}
  \item{\code{selected_root_attribute}}{Character root attribute of the embedded
    CTA member.}
  \item{\code{selected_tree_nodes}}{Integer split-node count in the embedded CTA
    member.}
  \item{\code{selected_tree_leaves}}{Integer leaf count in the embedded CTA
    member.}
  \item{\code{selected_endpoint_count}}{Integer endpoint (terminal leaf) count of
    the embedded CTA member; equals number of ORT child nodes.}
  \item{\code{child_ids}}{Character comma-separated child ORT node ids;
    empty string for terminal nodes.}
  \item{\code{method}}{Character; always \code{"lort"} for current fits.}
  \item{\code{selection_scope}}{Character; always \code{"local_node"} for LORT.}
  \item{\code{global_optimization}}{Logical; always \code{FALSE} for LORT.}
  \item{\code{sda_anchored}}{Logical; always \code{FALSE} for LORT.}
}
}
\seealso{
  \code{\link{cta_fit}}, \code{\link{predict.cta_ort}},
  \code{\link{summary.cta_ort}}
}
\examples{
X <- data.frame(A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
                B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20)))
y <- c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
ort <- cta_fit(X, y, recursive = TRUE, mc_iter = 100L, mc_seed = 42L,
               loo = "off", min_n = 5L)
cta_ort_node_table(ort)
}
