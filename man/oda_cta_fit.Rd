\name{oda_cta_fit}
\alias{oda_cta_fit}
\title{Fit a Classification Tree Analysis (CTA) model}
\description{
Builds a classification tree by recursively applying ODA at each node.
At each split, all attributes are evaluated and the attribute with the
highest ESS passing the significance threshold is selected. Matches
MegaODA CTA behaviour including MINDENOM, PRUNE, ENUMERATE, LOO STABLE,
and WEIGHT parameters.
}
\usage{
oda_cta_fit(X, y, w = NULL, priors_on = TRUE, miss_codes = NULL,
  alpha_split = 0.05, mindenom = 5L, prune_alpha = 1.0,
  max_depth = 10L, ess_min = 0,
  mc_iter = 25000L, mc_target = 0.05, mc_stop = 99.9, mc_stopup = 20,
  mc_seed = NULL, loo = "off", attr_names = NULL, K_segments = NULL,
  verbose = FALSE)
}
\arguments{
  \item{X}{Data frame or matrix of attribute columns.}
  \item{y}{Class variable vector.}
  \item{w}{Optional numeric case weights (MegaODA WEIGHT). Same length as y.}
  \item{priors_on}{Use prior-odds weighting at each node. Default TRUE.}
  \item{miss_codes}{Numeric vector of missing-value codes (MegaODA MISSING).}
  \item{alpha_split}{Significance threshold to split a node (MegaODA MC CUTOFF).
    Default 0.05.}
  \item{mindenom}{Minimum weighted node size to attempt a split (MegaODA
    MINDENOM). Default 5.}
  \item{prune_alpha}{Branches with p >= prune_alpha are not grown (MegaODA
    PRUNE). Default 1.0 = no pruning (unpruned tree).}
  \item{max_depth}{Maximum tree depth. Default 10.}
  \item{ess_min}{Minimum ESS required to split. Default 0.}
  \item{mc_iter}{Maximum MC iterations per node. Default 25000.}
  \item{mc_target, mc_stop, mc_stopup}{MC stopping parameters.}
  \item{mc_seed}{Base RNG seed; each node uses mc_seed + node_id * 1000 + attr_j.}
  \item{loo}{LOO mode per node: \code{"off"} (default), \code{"stable"}
    (MegaODA LOO STABLE), or \code{"pvalue"}.}
  \item{attr_names}{Attribute names. Defaults to column names of X.}
  \item{K_segments}{Segments for multiclass ordered splits. Default = C.}
  \item{verbose}{Logical; if \code{TRUE}, emit \code{[CTA]} progress messages
    via \code{message()} at each major stage. Default \code{FALSE}.}
}
\value{
An object of class \code{cta_tree} containing:
\describe{
  \item{\code{nodes}}{Named list of node objects, each with fields:
    node_id, parent_id, depth, n_obs, n_weighted, attribute, rule,
    ess, p_mc, loo_status, loo_ess, confusion, child_ids, split_labels,
    majority_class, leaf.}
  \item{\code{root_id}}{Integer ID of the root node.}
  \item{\code{n_nodes}}{Total number of nodes grown.}
}
Use \code{\link{predict.cta_tree}} to classify new data and
\code{\link{cta_node_table}} to extract the node summary table.
}
\examples{
## Binary CTA on mtcars
data(mtcars)
mt <- mtcars
X  <- mt[, c("cyl","disp","hp","wt")]
y  <- as.integer(mt$am)
tree <- oda_cta_fit(X, y, alpha_split = 0.05, mindenom = 5L,
                    mc_iter = 500L, mc_seed = 42L)
print(tree)
preds <- predict(tree, X)
mean(preds == y)   # training accuracy
}
\seealso{\code{\link{predict.cta_tree}}, \code{\link{cta_node_table}}}
