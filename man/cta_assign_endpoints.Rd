\name{cta_assign_endpoints}
\alias{cta_assign_endpoints}
\title{Assign observations to CTA terminal endpoints}
\description{
Traverses the fitted \code{cta_tree} for each row of \code{newdata} and
returns the terminal leaf reached, expressed as both its stored node
identifier (\code{endpoint_node_id}) and its sequential endpoint index
(\code{endpoint_id}) matching \code{\link{cta_endpoint_summary}}.

No endpoint membership is stored at fit time.  This function performs the
traversal on demand so the \code{cta_tree} object remains lean.  The
returned \code{endpoint_id} can be joined with the output of
\code{\link{cta_propensity_weights}} to assign endpoint-level stabilized
weights to individual observations.

\strong{Column order requirement:} \code{newdata} must have the same
attribute column order as the \code{X} matrix passed to
\code{\link{oda_cta_fit}}.  Traversal uses the stored integer column
positions (\code{attr_col}) from the fit, not column names.  If both
\code{names(newdata)} and \code{tree$attr_names} are non-NULL, a warning is
issued when they disagree at the split attribute positions.

\strong{Missingness:}
\describe{
  \item{\code{"na"} (default)}{Canonical path-local behaviour: when a split
    attribute value is \code{NA} or a stored miss-code on the observation's
    actual traversal path, the row returns \code{NA} for both output
    columns.  This matches the canonical \code{missing_action = "na"}
    semantics of \code{\link[=predict.cta_tree]{predict}}.}
  \item{\code{"majority"}}{Routes the observation to the child subtree with
    the larger \code{n_obs}, then continues traversal to a terminal leaf.
    Ties are resolved by selecting the first child.}
}
}
\usage{cta_assign_endpoints(tree, newdata, missing_action = c("na", "majority"))}
\arguments{
  \item{tree}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
  \item{newdata}{A \code{data.frame} (or coercible object) with the same
    column order as the training \code{X} supplied to
    \code{\link{oda_cta_fit}}.}
  \item{missing_action}{Character; one of \code{"na"} (default) or
    \code{"majority"}.  See Description.}
}
\value{
A \code{data.frame} with one row per row of \code{newdata} and columns:
\describe{
  \item{\code{row_id}}{Integer; positional row index in \code{newdata}
    (1 to \code{nrow(newdata)}).}
  \item{\code{endpoint_node_id}}{Integer; \code{node_id} of the terminal
    leaf reached by traversal.  \code{NA_integer_} when the observation
    cannot be routed to a terminal leaf (missing split attribute with
    \code{missing_action = "na"}, or no-tree fit).}
  \item{\code{endpoint_id}}{Integer; sequential endpoint index matching
    \code{\link{cta_endpoint_summary}}.  \code{NA_integer_} under the
    same conditions as \code{endpoint_node_id}.}
}
For no-tree fits all rows have \code{endpoint_node_id = NA_integer_} and
\code{endpoint_id = NA_integer_}.
}
\details{
\strong{Observation-level propensity weights (workflow sketch):}

\preformatted{
ep  <- cta_assign_endpoints(tree, X_train, missing_action = "na")
pw  <- cta_propensity_weights(tree, target_class = 1L, adjusted = TRUE)

# One row per classified training observation with its weight:
obs <- merge(
  data.frame(row_id = seq_len(nrow(X_train)),
             class  = as.character(y_train)),
  merge(ep, pw[, c("endpoint_id", "class", "adjusted_propensity_weight")],
        by = "endpoint_id"),
  by = c("row_id", "class")
)
# Rows with NA endpoint_id (missing root attribute) drop naturally.
}

Observation-level propensity weight expansion is intentionally left to the
caller so that the \code{cta_tree} object stores no observation indices.
}
\seealso{
  \code{\link{oda_cta_fit}}, \code{\link{cta_endpoint_summary}},
  \code{\link{cta_propensity_weights}}, \code{\link{predict.cta_tree}}
}
\examples{
data(mtcars)
X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
y    <- as.integer(mtcars$am)
tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
ep   <- cta_assign_endpoints(tree, X)
head(ep)
}
