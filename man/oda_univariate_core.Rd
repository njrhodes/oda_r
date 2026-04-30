\name{oda_univariate_core}
\alias{oda_univariate_core}
\title{Fit a univariate binary-class ODA model}
\description{
Low-level engine for binary-class Optimal Data Analysis. Handles ordered,
categorical, and binary attributes with optional prior-odds weighting,
Monte Carlo p-value, and leave-one-out validity analysis. Most users should
call \code{\link{oda_fit}} instead.
}
\usage{
oda_univariate_core(x, y, w = NULL,
  attr_type = c("auto","ordered","categorical","binary"),
  priors_on = TRUE, primary = NULL, secondary = NULL,
  miss_codes = NULL, missing_code = NULL,
  loo = c("off","stable","pvalue"), loo_alpha = 0.05,
  mcarlo = TRUE, mc_iter = 25000L, mc_target = 0.05,
  mc_stop = 99.9, mc_stopup = NA_real_, mc_adjust = FALSE,
  mc_seed = NULL, chance_model = c("class","attribute"))
}
\arguments{
  \item{x}{Attribute values.}
  \item{y}{Binary class labels, coercible to 0/1 integers.}
  \item{w}{Optional numeric case weights.}
  \item{attr_type}{Attribute type. \code{"auto"} detects from data.}
  \item{priors_on}{If \code{TRUE}, use inverse-frequency weighting.}
  \item{primary}{Primary tie-break heuristic. \code{NULL} = default by
    \code{priors_on}.}
  \item{secondary}{Secondary tie-break. \code{NULL} = \code{"samplerep"}.}
  \item{miss_codes}{Additional missing-value codes.}
  \item{missing_code}{Scalar alias for \code{miss_codes}.}
  \item{loo}{\code{"off"}, \code{"stable"}, or \code{"pvalue"}.}
  \item{loo_alpha}{Alpha threshold for \code{loo = "pvalue"}.}
  \item{mcarlo}{Run Monte Carlo p-value?}
  \item{mc_iter}{Maximum MC iterations.}
  \item{mc_target}{Significance threshold.}
  \item{mc_stop}{Confidence level (percent) for STOP early stopping.}
  \item{mc_stopup}{Confidence level (percent) for STOPUP.}
  \item{mc_adjust}{Legacy parameter; unused.}
  \item{mc_seed}{RNG seed.}
  \item{chance_model}{\code{"class"} (1/2) or \code{"attribute"}
    (1/k_attr) baseline.}
}
\value{
Named list. Key fields: \code{ok}, \code{rule}, \code{confusion} (list
with integer counts \code{TP}, \code{TN}, \code{FP}, \code{FN} and rate
fields \code{sensitivity}, \code{specificity} as proportions in [0,1]),
\code{ess}, \code{pac}, \code{p_mc}, \code{loo}, \code{n_eff}.
}
\seealso{\code{\link{oda_fit}}, \code{\link{oda_multiclass_unioda_core}}}
