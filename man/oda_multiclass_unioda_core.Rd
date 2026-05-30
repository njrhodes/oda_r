\name{oda_multiclass_unioda_core}
\alias{oda_multiclass_unioda_core}
\title{Fit a univariate multiclass ODA model}
\description{
Low-level engine for multiclass (C >= 3) Optimal Data Analysis. Handles
ordered and categorical attributes. Most users should call
\code{\link{oda_fit}} instead.
}
\usage{
oda_multiclass_unioda_core(x, y, w = NULL,
  attr_type = c("auto","ordered","categorical","binary"),
  priors_on = TRUE, miss_codes = NULL, missing_code = NULL,
  K_segments = NULL, degen = FALSE,
  mcarlo = TRUE, mc_iter = 25000L, mc_target = 0.05,
  mc_stop = 99.9, mc_stopup = NA, mc_adjust = FALSE, mc_seed = NULL,
  loo = c("off","on"),
  boundary_mode = c("megaoda_halfopen","right_closed"),
  loo_opts = list(),
  direction = "off", direction_map = NULL)
}
\arguments{
  \item{x}{Attribute values (numeric or factor).}
  \item{y}{Integer class labels (will be re-coded to 1..C internally).}
  \item{w}{Optional numeric case weights.}
  \item{attr_type}{Attribute type.}
  \item{priors_on}{Inverse-frequency weighting.}
  \item{miss_codes}{Additional missing codes (scalar or vector).}
  \item{missing_code}{Alias for \code{miss_codes}.}
  \item{K_segments}{Number of segments; default = C.}
  \item{degen}{Allow degenerate solutions?}
  \item{mcarlo}{Run Monte Carlo p-value?}
  \item{mc_iter, mc_target, mc_stop, mc_stopup, mc_adjust, mc_seed}{MC parameters.}
  \item{loo}{\code{"off"} or \code{"on"}.}
  \item{boundary_mode}{Boundary convention for ordered cut values.}
  \item{loo_opts}{Named list of LOO options passed to the LOO engine.}
  \item{direction}{Directional constraint (MPE Chapter 4). \code{"ascending"}
    forces segment s to class s for ordered multiclass; auto-creates identity
    \code{direction_map} for categorical when L == C. \code{"descending"} uses
    reverse assignment. Default \code{"off"} (nondirectional).}
  \item{direction_map}{Named integer vector for categorical fixed-partition
    DIRECTIONAL. Names are attribute levels; values are class labels 1..C.
    When supplied, bypasses the partition search and evaluates only the
    specified mapping. Default \code{NULL}.}
}
\value{
Named list. Key fields: \code{ok}, \code{rule} (with \code{cut_values}
and \code{seg_classes}), \code{confusion} (weighted count matrix),
\code{pac}, \code{mean_pac}, \code{ess_pac}, \code{p_mc}, \code{loo},
\code{n_eff}.

\strong{Note on confusion matrix:} \code{confusion} contains weighted
counts (priors-adjusted when \code{priors_on = TRUE}). For raw integer
counts use \code{loo$confusion_raw}.
}
\seealso{\code{\link{oda_fit}}}
