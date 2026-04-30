\name{oda_rule_predict}
\alias{oda_rule_predict}
\title{Apply a binary ODA rule to new data}
\description{
Predict class labels (0 or 1) for new attribute values using a fitted
binary ODA rule.
}
\usage{oda_rule_predict(x, rule)}
\arguments{
  \item{x}{Numeric or character attribute values.}
  \item{rule}{A rule list returned in \code{fit$rule} from
    \code{\link{oda_univariate_core}} or \code{\link{oda_fit}}.}
}
\value{Integer vector of predicted class labels (0 or 1).}
