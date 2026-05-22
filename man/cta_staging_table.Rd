\name{cta_staging_table}
\alias{cta_staging_table}
\title{Staging table for a fitted CTA tree}
\description{
Returns one row per terminal endpoint ordered by ascending target-class
propensity (lowest to highest risk stratum).  Empirical counts,
proportions, and odds are computed from the stored leaf class counts.
When an endpoint is perfectly predicted (100 percent one class), the
empirical odds and proportion are undefined; the \code{adjust_perfect}
option adds one hypothetical misclassified observation to the undefined
profile so all endpoints can be ranked and compared — a canon remedy
anchored in Yarnold and Linden (2017).

\strong{Scope:} The two-class case is handled automatically when
\code{target_class = NULL} (defaults to the numerically larger class
label, typically 1).  For trees with three or more classes
\code{target_class} must be supplied explicitly.
}
\usage{cta_staging_table(tree, target_class = NULL, weighted = FALSE,
                   adjust_perfect = TRUE)}
\arguments{
  \item{tree}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
  \item{target_class}{Integer (or coercible); the class label treated as
    the target (positive / high-risk) class.  \code{NULL} (default) uses
    the numerically largest class label for binary trees, and stops for
    trees with three or more classes.}
  \item{weighted}{Logical.  \code{FALSE} (default) uses raw observation
    counts; \code{TRUE} uses case-weighted counts.}
  \item{adjust_perfect}{Logical.  \code{TRUE} (default) applies the
    one-hypothetical-misclassification adjustment to perfectly predicted
    endpoints so that all endpoints can be ordered by propensity.}
}
\value{
A \code{data.frame} with one row per terminal endpoint, ordered by
ascending target-class propensity (lowest to highest risk stratum),
with columns:
\describe{
  \item{\code{stage}}{Integer rank 1..n, ascending by target proportion.}
  \item{\code{endpoint_id}}{Integer sequential endpoint index, matching
    \code{\link{cta_endpoint_summary}}.}
  \item{\code{endpoint_node_id}}{Integer tree node identifier.}
  \item{\code{path}}{Character; AND-joined branch labels from root.}
  \item{\code{terminal_prediction}}{Integer majority-class prediction.}
  \item{\code{target_class}}{Integer; the target class used for this
    table.}
  \item{\code{target_n}}{Numeric; raw (or weighted) count of target-class
    observations at this endpoint.}
  \item{\code{denominator}}{Numeric; total raw (or weighted) observations
    at this endpoint.}
  \item{\code{target_proportion}}{Numeric; empirical target-class
    proportion (\code{target_n / denominator}).}
  \item{\code{non_target_n}}{Numeric; denominator minus target_n.}
  \item{\code{odds}}{Numeric; empirical odds
    (\code{target_n / non_target_n}); \code{NA} when
    \code{perfectly_predicted} is \code{TRUE}.}
  \item{\code{perfectly_predicted}}{Logical; \code{TRUE} when the
    endpoint is 100 percent one class (\code{target_n == 0} or
    \code{non_target_n == 0}).}
  \item{\code{adjusted}}{Logical; \code{TRUE} when the
    one-hypothetical-misclassification adjustment has been applied.
    Always \code{FALSE} when \code{adjust_perfect = FALSE}.}
  \item{\code{adjusted_target_n}}{Numeric; target_n after adjustment.
    Equal to \code{target_n} when \code{adjusted} is \code{FALSE}.}
  \item{\code{adjusted_denominator}}{Numeric; denominator after
    adjustment.}
  \item{\code{adjusted_target_proportion}}{Numeric; adjusted proportion.}
  \item{\code{adjusted_non_target_n}}{Numeric; adjusted non-target
    count.}
  \item{\code{adjusted_odds}}{Numeric; adjusted odds.}
  \item{\code{weighted}}{Logical; the value of the \code{weighted}
    argument.}
  \item{\code{n_obs}}{Integer; raw observation count at this endpoint
    (from \code{\link{cta_endpoint_summary}}).}
  \item{\code{n_weighted}}{Numeric; weighted observation count.}
}
For a no-tree fit the returned data frame has zero rows but the correct
column structure and types.
}
\references{
Yarnold PR, Linden A (2017). Computing propensity score weights for CTA
models involving perfectly predicted endpoints. \emph{Optimal Data
Analysis}, \strong{6}, 43-46.
}
\seealso{
  \code{\link{oda_cta_fit}}, \code{\link{cta_endpoint_summary}},
  \code{\link{cta_endpoint_counts}}, \code{\link{cta_propensity_weights}}
}
\examples{
data(mtcars)
X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
y    <- as.integer(mtcars$am)
tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
cta_staging_table(tree)
}
