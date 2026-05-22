\name{cta_propensity_weights}
\alias{cta_propensity_weights}
\title{Endpoint-level propensity-score weights for a fitted CTA tree}
\description{
Returns one row per terminal endpoint per actual class, containing the
CTA-derived stabilized propensity-style weights described in Yarnold and
Linden (2017).  All values are computed on demand from the stored leaf
class counts; no refitting, no prediction, and no training-data
recomputation is performed.

\strong{Formula:} For endpoint \eqn{s} and actual class \eqn{z},
\deqn{w_{s,z} = \frac{n_s \cdot \Pr(Z=z)}{n_{s,z}}}
where \eqn{n_s} is the endpoint denominator, \eqn{n_{s,z}} is the raw
count of class \eqn{z} observations at endpoint \eqn{s}, and
\eqn{\Pr(Z=z)} is the marginal class probability across the full
classified analytic sample.

\strong{Perfect endpoints:} When \eqn{n_{s,z} = 0} for some class, the
empirical weight is undefined (\code{Inf}).  When \code{adjusted = TRUE}
(default), one hypothetical misclassified observation is added to the
absent class profile — and to the global marginal totals — so that all
endpoint × class cells yield finite adjusted weights.  This is the canon
remedy from Yarnold and Linden (2017).

\strong{Scope:} Raw observation counts (\code{n_raw}) are used exclusively.
The function does not return observation-level weights; those require
endpoint membership per training observation, which is not stored on the
fitted tree.
}
\usage{cta_propensity_weights(tree, target_class = NULL, adjusted = TRUE)}
\arguments{
  \item{tree}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
  \item{target_class}{Integer (or coercible); annotation column only —
    does not filter output rows.  \code{NULL} (default) uses the
    numerically largest class label for binary trees, and stops for
    trees with three or more classes.}
  \item{adjusted}{Logical.  \code{TRUE} (default) applies the
    one-hypothetical-misclassification adjustment so that all cells
    yield finite adjusted weights.  \code{FALSE} leaves undefined
    weights as \code{Inf} and adjusted columns equal to empirical.}
}
\value{
A \code{data.frame} with one row per terminal endpoint per actual class,
with columns:
\describe{
  \item{\code{endpoint_id}}{Integer sequential endpoint index.}
  \item{\code{endpoint_node_id}}{Integer tree node identifier.}
  \item{\code{path}}{Character; AND-joined branch labels from root.}
  \item{\code{terminal_prediction}}{Integer majority-class prediction.}
  \item{\code{class}}{Character; actual class label for this row.}
  \item{\code{target_class}}{Integer; design-annotation class label.}
  \item{\code{class_n}}{Integer; raw count of this class at this
    endpoint (empirical \eqn{n_{s,z}}).}
  \item{\code{endpoint_n}}{Integer; total raw observations at this
    endpoint (empirical \eqn{n_s}).}
  \item{\code{marginal_class_n}}{Integer; total raw observations of this
    class across all endpoints (empirical \eqn{N_z}).}
  \item{\code{marginal_total_n}}{Integer; total classified observations
    across all endpoints (empirical \eqn{N}).}
  \item{\code{marginal_class_probability}}{Numeric; empirical marginal
    class probability \eqn{\Pr(Z=z) = N_z / N}.}
  \item{\code{propensity_weight}}{Numeric; empirical stabilized weight
    \eqn{n_s \cdot \Pr(Z=z) / n_{s,z}}.  \code{Inf} when
    \code{class_n == 0}.}
  \item{\code{undefined_empirical}}{Logical; \code{TRUE} when
    \code{class_n == 0}.}
  \item{\code{perfectly_predicted_endpoint}}{Logical; \code{TRUE} when
    any class has \code{class_n == 0} at this endpoint.}
  \item{\code{adjusted}}{Logical; \code{TRUE} when the
    one-hypothetical-observation adjustment was applied to this row.}
  \item{\code{adjusted_class_n}}{Numeric; \code{class_n + 1} where
    \code{adjusted}, otherwise \code{class_n}.}
  \item{\code{adjusted_endpoint_n}}{Numeric; endpoint denominator
    after adjustment.}
  \item{\code{adjusted_marginal_class_n}}{Numeric; global class count
    after all hypothetical additions.}
  \item{\code{adjusted_marginal_total_n}}{Numeric; global total after
    all hypothetical additions.}
  \item{\code{adjusted_marginal_class_probability}}{Numeric; adjusted
    marginal class probability.}
  \item{\code{adjusted_propensity_weight}}{Numeric; adjusted weight.
    Finite whenever \code{adjusted_class_n > 0}.}
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
  \code{\link{oda_cta_fit}}, \code{\link{cta_endpoint_counts}},
  \code{\link{cta_staging_table}}
}
\examples{
data(mtcars)
X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
y    <- as.integer(mtcars$am)
tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
cta_propensity_weights(tree)
}
