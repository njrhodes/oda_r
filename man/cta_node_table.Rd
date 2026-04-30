\name{cta_node_table}
\alias{cta_node_table}
\title{Extract CTA node summary table}
\description{
Returns a data frame with one row per node, mirroring the MegaODA node
table output (ATTRIBUTE, NODE, LEV, OBS, p, ESS, LOO columns).
}
\usage{cta_node_table(tree)}
\arguments{
  \item{tree}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
}
\value{
Data frame with columns: node_id, parent_id, depth, attribute, n_obs,
n_weighted, p_mc, ess, loo_status, loo_ess, leaf, attr_type.
}
