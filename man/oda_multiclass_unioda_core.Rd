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
  priors_on = TRUE, miss_codes = NULL,
  K_segments = NULL, degen = FALSE,
  mcarlo = TRUE, mc_iter = 25000L, mc_target = 0.05,
  mc_stop = 99.9, mc_stopup = 20, mc_adjust = FALSE, mc_seed = NULL,
  loo = c("off","on"),
  cut_value_mode = c("midpoint","lower","upper"),
  loo_grid_mode = c("fixed","refit"),
  boundary_mode = c("megaoda_halfopen","right_closed"),
  loo_use_samplerep = FALSE, loo_return_folds = FALSE,
  loo_priors_mode = c("fold","global"))
}
\arguments{
  \item{x}{Attribute values (numeric or factor).}
  \item{y}{Integer class labels (will be re-coded to 1..C internally).}
  \item{w}{Optional numeric case weights.}
  \item{attr_type}{Attribute type.}
  \item{priors_on}{Inverse-frequency weighting.}
  \item{miss_codes}{Additional missing codes.}
  \item{K_segments}{Number of segments; default = C.}
  \item{degen}{Allow degenerate solutions?}
  \item{mcarlo}{Run Monte Carlo p-value?}
  \item{mc_iter, mc_target, mc_stop, mc_stopup, mc_adjust, mc_seed}{MC parameters.}
  \item{loo}{\code{"off"} or \code{"on"}.}
  \item{cut_value_mode}{\code{"midpoint"} (default), \code{"lower"}, or
    \code{"upper"}.}
  \item{loo_grid_mode}{\code{"fixed"} (global grid) or \code{"refit"}
    (true per-fold rebuild, MegaODA-faithful).}
  \item{boundary_mode}{Boundary convention.}
  \item{loo_use_samplerep, loo_return_folds, loo_priors_mode}{LOO options.}
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
