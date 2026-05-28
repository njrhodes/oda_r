\name{cta_node_table}
\alias{cta_node_table}
\title{Canonical CTA node report table}
\description{
Returns a data frame with one row per node, mirroring the CTA.exe-style
node report (ATTRIBUTE, NODE, LEV, OBS, p, ESS/WESS, LOO, WESSL/LOO ESS,
LOOp, TYP, MODEL columns).

Split nodes carry canonical split metrics (ESS, WESS, p, LOO status, LOO
ESS/WESSL, LOOp) and a MODEL field with branch strings and terminal-leaf
\code{*} markers.  Leaf rows have \code{NA} for all split metrics and MODEL.
}
\usage{cta_node_table(tree)}
\arguments{
  \item{tree}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
}
\value{
Data frame with columns:
\describe{
  \item{\code{node_id}}{Integer node identifier.}
  \item{\code{parent_id}}{Integer parent node identifier (0 for root).}
  \item{\code{level}}{Integer level from root (root = 1); alias for
    \code{depth}.}
  \item{\code{depth}}{Integer depth from root (root = 1).}
  \item{\code{leaf}}{Logical; \code{TRUE} for terminal leaf nodes.}
  \item{\code{attribute}}{Character attribute name (\code{NA} for leaves).}
  \item{\code{attr_type}}{Character attribute type (\code{NA} for leaves).}
  \item{\code{n_obs}}{Integer observation count at this node.}
  \item{\code{n_weighted}}{Numeric weighted observation count.}
  \item{\code{p_mc}}{Numeric Monte Carlo p-value (\code{NA} for leaves).}
  \item{\code{ess}}{Numeric ESS at this split (\code{NA} for leaves).}
  \item{\code{ess_weighted}}{Numeric WESS at this split (\code{NA} for
    leaves); equals \code{ess} when case weights are not active.}
  \item{\code{loo_status}}{Character LOO status, e.g. \code{"STABLE"}
    (\code{NA} for leaves).}
  \item{\code{loo_ess}}{Numeric LOO ESS/WESSL (\code{NA} for leaves).}
  \item{\code{loo_p}}{Numeric LOO p-value (\code{NA} for leaves).}
  \item{\code{model}}{Character CTA.exe-style branch string with
    terminal-leaf \code{*} markers, e.g.
    \code{"<=0.5-->0,101/131,77.10\%*; >0.5-->1,21/55,38.18\%*"}.
    \code{NA} for leaf nodes.}
}
}
\seealso{
  \code{\link{oda_cta_fit}}, \code{\link{cta_endpoint_table}},
  \code{\link{summary.cta_tree}}
}
