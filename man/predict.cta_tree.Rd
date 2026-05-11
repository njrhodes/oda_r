\name{predict.cta_tree}
\alias{predict.cta_tree}
\title{Classify new observations using a CTA tree}
\description{
Applies a fitted \code{cta_tree} to new data by routing each observation
through the tree until it reaches a leaf node.
}
\usage{\method{predict}{cta_tree}(object, newdata,
  missing_action = c("majority", "na"), ...)}
\arguments{
  \item{object}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
  \item{newdata}{Data frame or matrix with the same columns as training X.}
  \item{missing_action}{How to handle observations whose split attribute is
    missing on their traversal path. \code{"majority"} (default) routes the
    observation to the current node's majority class. \code{"na"} returns
    \code{NA_integer_} for that observation (canonical path-local
    missingness).}
  \item{...}{Unused.}
}
\value{Integer vector of predicted class labels, length \code{nrow(newdata)}.
  When \code{missing_action = "na"}, observations missing a split attribute
  on their path receive \code{NA_integer_}.}
