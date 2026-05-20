\name{summary.cta_tree}
\alias{summary.cta_tree}
\title{Summarize a fitted CTA tree}
\description{
Returns a structured list with class \code{"cta_tree_summary"} capturing
tree-level metadata.  All fields are read directly from stored objects;
no refitting or prediction is performed.
}
\usage{\method{summary}{cta_tree}(object, ...)}
\arguments{
  \item{object}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
  \item{...}{Unused.}
}
\value{
A list of class \code{"cta_tree_summary"} with fields:
\describe{
  \item{\code{status}}{Character: \code{"valid_tree"}, \code{"stump"},
    or \code{"no_tree"}.}
  \item{\code{no_tree}}{Logical; \code{TRUE} for leaf-only fits.}
  \item{\code{root_attribute}}{Character attribute name at the root split;
    \code{NA_character_} for no-tree fits.}
  \item{\code{n_nodes}}{Total number of nodes including leaves.}
  \item{\code{n_splits}}{Number of non-leaf (split) nodes.}
  \item{\code{n_leaves}}{Number of terminal leaf endpoints (= \code{strata}).}
  \item{\code{strata}}{Alias for \code{n_leaves}; \code{NA_integer_} for
    no-tree fits.}
  \item{\code{overall_ess}}{WESS when weights are active, ESS otherwise;
    \code{NA_real_} when absent.}
  \item{\code{d}}{D statistic (\code{NA_real_} for no-tree or ESS \eqn{\le} 0).}
  \item{\code{min_terminal_denom}}{Smallest leaf \code{n_obs};
    \code{NA_integer_} for no-tree fits.}
  \item{\code{endpoint_denominators}}{Named integer vector of leaf
    \code{n_obs}; \code{integer(0)} for no-tree fits.}
  \item{\code{has_weights}}{Logical; \code{TRUE} when case weights are active.}
  \item{\code{mindenom}}{MINDENOM used when fitting.}
  \item{\code{alpha_split}}{Significance threshold used when fitting.}
  \item{\code{prune_alpha}}{Pruning threshold used when fitting.}
  \item{\code{loo}}{LOO mode string used when fitting.}
}
}
\seealso{\code{\link{oda_cta_fit}}, \code{\link{cta_node_table}},
  \code{\link{cta_strata}}, \code{\link{cta_d_stat}},
  \code{\link{print.cta_tree_summary}}}
\examples{
data(mtcars)
X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
y    <- as.integer(mtcars$am)
tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
s    <- summary(tree)
print(s)
}
