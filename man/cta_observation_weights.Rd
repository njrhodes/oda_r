\name{cta_observation_weights}
\alias{cta_observation_weights}
\title{Assign per-observation CTA propensity weights}
\description{
Convenience wrapper that calls \code{\link{cta_assign_endpoints}} and
\code{\link{cta_propensity_weights}} and returns a joined observation-level
data frame.  The \code{cta_tree} object is not mutated; all computation is
on demand.

\strong{Column order requirement:} \code{newdata} must have the same attribute
column order as the \code{X} matrix passed to \code{\link{oda_cta_fit}}.
Traversal uses the stored integer column positions (\code{attr_col}) from the
fit, not column names.

\strong{Unroutable observations:} Observations with \code{NA} endpoint
(missing root split attribute under \code{missing_action = "na"}) or
\code{NA} class label receive \code{assigned = FALSE} and \code{NA} for all
weight columns.  The output always contains \code{nrow(newdata)} rows.

\strong{Unmatched classified observations:} When a non-\code{NA} endpoint
observation's class is not present in the propensity weight table (e.g.,
a class unseen at fit time), a warning is issued and \code{assigned = FALSE}.
}
\usage{
cta_observation_weights(tree, newdata, y, target_class = NULL,
                        adjusted = TRUE,
                        missing_action = c("na", "majority"))
}
\arguments{
  \item{tree}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
  \item{newdata}{A \code{data.frame} (or coercible object) with the same
    column order as the training \code{X}.}
  \item{y}{Class labels for each row of \code{newdata}.  Any type coercible
    to character; length must equal \code{nrow(newdata)}.}
  \item{target_class}{Passed to \code{\link{cta_propensity_weights}} as an
    annotation parameter.  Identifies which class is treated as the design
    target (high-risk class) for the \code{target_class} output column; it does
    \emph{not} filter the endpoint \eqn{\times} class rows used for the join.
    Each observation is matched to its own \code{actual_class} regardless of
    this value.  \code{NULL} (default) lets \code{cta_propensity_weights}
    resolve the target class automatically (numerically largest class for binary
    trees; must be supplied explicitly for trees with three or more classes).}
  \item{adjusted}{Logical; passed to \code{\link{cta_propensity_weights}}.
    Default \code{TRUE}.}
  \item{missing_action}{Character; one of \code{"na"} (default) or
    \code{"majority"}.  Passed to \code{\link{cta_assign_endpoints}}.}
}
\value{
A \code{data.frame} with \code{nrow(newdata)} rows and columns:
\describe{
  \item{\code{row_id}}{Integer; positional row index (1 to
    \code{nrow(newdata)}).}
  \item{\code{actual_class}}{Character; class label from \code{y},
    coerced to character.}
  \item{\code{endpoint_node_id}}{Integer; node ID of the terminal leaf reached
    by traversal, or \code{NA_integer_} when unroutable.}
  \item{\code{endpoint_id}}{Integer; sequential endpoint index matching
    \code{\link{cta_endpoint_summary}}, or \code{NA_integer_}.}
  \item{\code{target_class}}{Integer; resolved design target class annotation
    from \code{\link{cta_propensity_weights}}, or \code{NA_integer_} when
    unassigned.}
  \item{\code{propensity_weight}}{Numeric; unadjusted propensity weight for the
    observation's endpoint--class cell, or \code{NA} when unassigned.}
  \item{\code{adjusted_propensity_weight}}{Numeric; adjusted propensity weight
    (Yarnold-Linden correction for perfectly predicted endpoints), or \code{NA}
    when unassigned.}
  \item{\code{undefined_empirical}}{Logical; \code{TRUE} when the endpoint--class
    cell has zero observed frequency, or \code{NA} when unassigned.}
  \item{\code{perfectly_predicted_endpoint}}{Logical; \code{TRUE} when all
    observations at the endpoint belong to one class, or \code{NA} when
    unassigned.}
  \item{\code{adjusted}}{Logical; \code{TRUE} when the adjusted weight was
    applied at this endpoint, or \code{NA} when unassigned.}
  \item{\code{assigned}}{Logical; \code{TRUE} when a propensity weight was
    successfully matched for this observation.}
}
}
\details{
No observation-level data are stored in the \code{cta_tree} object at fit
time.  This function performs traversal and weight lookup on demand.

\strong{No-tree fits:} When the tree has no splits (leaf-only), all rows have
\code{endpoint_id = NA_integer_} and \code{assigned = FALSE}.

\strong{Join semantics:} The join key is
\code{paste(endpoint_id, actual_class)}.  Each observation is matched to the
propensity weight row whose \code{class} equals its \code{actual_class}.
The \code{target_class} parameter annotates all rows with the resolved design
target class but does not affect which rows participate in the join.
}
\seealso{
  \code{\link{cta_assign_endpoints}}, \code{\link{cta_propensity_weights}},
  \code{\link{oda_cta_fit}}, \code{\link{cta_endpoint_summary}}
}
\examples{
data(mtcars)
X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
y    <- as.integer(mtcars$am)
tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
ow   <- cta_observation_weights(tree, X, y)
head(ow)
}
