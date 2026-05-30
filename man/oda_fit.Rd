\name{oda_fit}
\alias{oda_fit}
\title{Fit an ODA model}
\description{
Unified entry point for Optimal Data Analysis. Dispatches to the
binary-class engine when the outcome has exactly two distinct values, or
the multiclass engine for three or more classes. This is the function CTA
nodes call at each split candidate.
}
\usage{
oda_fit(x, y, w = NULL,
        attr_type = c("auto","ordered","categorical","binary"),
        priors_on = TRUE, K_segments = NULL, degen = FALSE,
        miss_codes = NULL, missing_code = NULL,
        mcarlo = TRUE, mc_iter = 25000L, mc_target = 0.05,
        mc_stop = 99.9, mc_stopup = NA, mc_seed = NULL,
        loo = "off",
        boundary_mode = c("megaoda_halfopen","right_closed"),
        eval_order = c("mc_then_loo","loo_then_mc"),
        mindenom = 1L,
        direction = c("both","off","greater","less","ascending","descending"),
        direction_map = NULL)
}
\arguments{
  \item{x}{Attribute values (numeric, factor, character, or logical).}
  \item{y}{Class labels; must have 2 or 3+ distinct values.}
  \item{w}{Optional numeric case weights. Default: unit weights. These are
    economic or importance weights, distinct from prior-odds weighting which
    is controlled by \code{priors_on}.}
  \item{attr_type}{Attribute measurement type: \code{"auto"} (default),
    \code{"ordered"}, \code{"categorical"}, or \code{"binary"}.}
  \item{priors_on}{Logical; if \code{TRUE} (default), weight classes by the
    reciprocal of their sample frequency so that the objective maximizes mean
    PAC rather than overall PAC.}
  \item{K_segments}{Number of segments for multiclass ordered models.
    Default equals the number of classes \eqn{C}.}
  \item{degen}{Logical; if \code{FALSE} (default) require all \eqn{C}
    classes to appear in the predicted labels.}
  \item{miss_codes}{Numeric vector of values to treat as missing (excluded
    from analysis).}
  \item{missing_code}{Scalar alias for \code{miss_codes}; accepted for
    convenience.}
  \item{mcarlo}{Logical; run Monte Carlo Fisher-randomization p-value?
    Default \code{TRUE}.}
  \item{mc_iter}{Maximum Monte Carlo iterations. Default 25000.}
  \item{mc_target}{Significance threshold for STOP early stopping.
    Default 0.05.}
  \item{mc_stop}{Confidence level (percent) for lower-tail STOP.
    Default 99.9.}
  \item{mc_stopup}{Confidence level (percent) for upper-tail STOPUP.
    Default NA (disabled; matches MegaODA behavior).}
  \item{mc_seed}{Optional integer RNG seed for reproducibility.}
  \item{loo}{LOO mode. \code{"off"} (default): no LOO filter.
    \code{"on"}: synonym for \code{"pvalue"} when used with the multiclass
    engine.
    \code{"stable"}: binary only; accept when LOO ESS equals training ESS
    (|WESSL - WESS| <= 0.01 pp); split node reports
    \code{loo_status = "STABLE"}.
    \code{"pvalue"}: binary only; accept when LOO Fisher p is \emph{strictly
    less than} 0.05 (default threshold); split node reports
    \code{loo_status = "PVALUE"}.
    Numeric in (0, 1): binary only; accept when LOO Fisher p is \emph{strictly
    less than} the supplied value; must be a single finite value strictly in
    (0, 1); split node reports \code{loo_status = "PVALUE"}.
    Do not describe the p-value gate as "STABLE"  -  the two modes are distinct.}
  \item{boundary_mode}{Boundary convention for multiclass ordered rules.
    Default \code{"megaoda_halfopen"} matches MegaODA.exe behaviour.}
  \item{eval_order}{Controls whether Monte Carlo testing is run before LOO
    validation or whether eligible ordered-cut LOO stability is checked
    before Monte Carlo. The default \code{"mc_then_loo"} preserves
    standalone UniODA behaviour. CTA tree building uses \code{"loo_then_mc"}
    internally to reject LOO-unstable ordered-cut candidates before spending
    MC iterations.}
  \item{mindenom}{Minimum raw observation count required in each child node
    for a candidate cut to be evaluated. Default 1 (no enforcement).}
  \item{direction}{Directional hypothesis control.
    \code{"both"} (default, non-directional) or its synonym \code{"off"}.
    \code{"greater"} / \code{"less"}: MPE Chapter 2 binary ordered DIRECTIONAL
    (high / low attribute values predict class 1; binary class only  -  error on
    multiclass).
    \code{"ascending"} / \code{"descending"}: MPE Chapter 4 DIRECTIONAL for
    multiclass ordered (segment s maps to class s or C+1-s) or multiclass
    categorical with L == C (auto-creates identity or reverse mapping). Error on
    binary class (use \code{"greater"} / \code{"less"} for binary ordered).}
  \item{direction_map}{Named integer vector for categorical fixed-partition
    DIRECTIONAL (MPE Chapter 4). Names are attribute levels (character); values
    are predicted class labels. All attribute levels must be covered exactly once
    with at least two distinct target classes. When supplied, ODA evaluates only
    the specified mapping and skips the partition search. For binary class, values
    should be the original class labels (recoded to 0/1 internally). For
    multiclass, values should be 1..C class labels. Compatible with
    \code{direction = "both"} (default). Default \code{NULL}.}
}
\value{
A named list with components:
\describe{
  \item{\code{ok}}{Logical; \code{TRUE} if a valid model was found.}
  \item{\code{reason}}{Character reason string if \code{ok = FALSE}.}
  \item{\code{rule}}{The fitted rule (list; structure depends on
    \code{attr_type} and engine).}
  \item{\code{n_eff}}{Number of observations used (after missing removal).}
  \item{\code{ess}}{Effect Strength for Sensitivity (percent), scaled
    0--100.}
  \item{\code{pac}}{Percentage Accuracy in Classification (training).}
  \item{\code{p_mc}}{Monte Carlo p-value, or \code{NA} if
    \code{mcarlo = FALSE}.}
  \item{\code{loo}}{LOO results list, or \code{NULL} if \code{loo = "off"}.}
  \item{\code{engine}}{Character; \code{"binary"} or \code{"multiclass"}.}
  \item{\code{confusion}}{Confusion table. For the binary engine this is a
    list with integer counts \code{TP}, \code{TN}, \code{FP}, \code{FN}
    plus \code{sensitivity} and \code{specificity} as proportions [0,1].
    For the multiclass engine this is a numeric matrix of (possibly
    weighted) counts.}
}
}
\examples{
## Binary (C = 2)
x <- c(1,2,3,4,5,6,7,8)
y <- c(0L,0L,0L,0L,1L,1L,1L,1L)
fit <- oda_fit(x, y, mcarlo = FALSE)
fit$ok
fit$rule$cut_value

## Multiclass (C = 3)
x3 <- c(1,2,3,4,5,6,7,8,9)
y3 <- c(1L,1L,1L,2L,2L,2L,3L,3L,3L)
fit3 <- oda_fit(x3, y3, mcarlo = FALSE)
fit3$rule$cut_values
fit3$rule$seg_classes
}
