\name{oda_rule_predict_multiclass}
\alias{oda_rule_predict_multiclass}
\title{Apply a multiclass ODA rule to new data}
\description{
Predict class labels for new attribute values using a fitted multiclass
ODA rule.
}
\usage{
oda_rule_predict_multiclass(x, rule,
  boundary = c("megaoda_halfopen","right_closed"))
}
\arguments{
  \item{x}{Numeric attribute values.}
  \item{rule}{A rule list from \code{\link{oda_multiclass_unioda_core}}
    or \code{\link{oda_fit}}.}
  \item{boundary}{Boundary convention. Default \code{"megaoda_halfopen"}
    matches MegaODA.exe.}
}
\value{Integer vector of predicted class labels.}
