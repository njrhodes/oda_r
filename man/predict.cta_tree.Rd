\name{predict.cta_tree}
\alias{predict.cta_tree}
\title{Classify new observations using a CTA tree}
\description{
Applies a fitted \code{cta_tree} to new data by routing each observation
through the tree until it reaches a leaf node.
}
\usage{\method{predict}{cta_tree}(object, newdata, ...)}
\arguments{
  \item{object}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
  \item{newdata}{Data frame or matrix with the same columns as training X.}
  \item{...}{Unused.}
}
\value{Integer vector of predicted class labels.}
