###############################################################################
# R/oda_fit.R
#
# oda_fit() — unified public entry point for ODA.
#
# Routes to oda_univariate_core() for binary-class problems (C=2)
# or oda_multiclass_unioda_core() for multiclass problems (C>=3).
# CTA nodes call oda_fit() at each split candidate.
#
# The rule object in the returned list has a consistent structure
# regardless of which engine was invoked:
#   $rule$type        "ordered_cut" | "binary_map" | "nominal_cut" |
#                     "multiclass_ordered" | "multiclass_nominal"
#   $rule$cut_value   (UniODA ordered only)
#   $rule$cut_values  (MultiODA ordered only)
#   $rule$seg_classes (MultiODA ordered only)
#   $rule$direction   (UniODA only)
#
# Returned list always contains:
#   ok, reason, rule, ess, pac, p_mc, loo, n_eff, attr_type
###############################################################################

#' Fit an ODA model (unified entry point for binary and multiclass problems)
#'
#' Dispatches to the binary-class engine when the outcome has exactly two
#' distinct values, or the multiclass engine when it has three or more.
#' This is the entry point CTA nodes should call.
#'
#' @param x Attribute values (numeric, factor, character, or logical).
#' @param y Class labels (any type coercible to integer; must have 2 or 3+ unique values).
#' @param w Optional case weights. Default: unit weights.
#' @param attr_type One of "auto", "ordered", "categorical", "binary".
#' @param priors_on Logical; if TRUE weight by inverse class frequency (default TRUE).
#' @param K_segments Number of segments for multiclass ordered models (default = C).
#' @param degen Allow degenerate multiclass solutions? Default FALSE.
#' @param miss_codes Additional values to treat as missing (scalar or vector).
#' @param mcarlo Run Monte Carlo p-value? Default TRUE.
#' @param mc_iter Maximum MC iterations. Default 25000.
#' @param mc_target Significance threshold. Default 0.05.
#' @param mc_stop Confidence for lower-tail early stop (percent). Default 99.9.
#' @param mc_stopup Confidence for upper-tail early stop (percent). Default 20.
#' @param mc_seed Optional RNG seed for reproducibility.
#' @param loo LOO mode. \code{"off"} (default): no LOO filter.
#'   \code{"on"}: synonym for \code{"pvalue"} when used with the multiclass engine.
#'   \code{"stable"}: binary only; accept when LOO ESS equals training ESS
#'   (|WESSL − WESS| ≤ 0.01 pp); split node reports \code{loo_status = "STABLE"}.
#'   \code{"pvalue"}: binary only; accept when LOO Fisher p is strictly less than
#'   0.05 (default threshold); split node reports \code{loo_status = "PVALUE"}.
#'   Numeric in (0, 1): binary only; accept when LOO Fisher p is strictly less than
#'   the supplied value; must be a single finite value strictly in (0, 1);
#'   split node reports \code{loo_status = "PVALUE"}.
#'   Do not describe the p-value gate as "STABLE" — the two modes are distinct.
#' @param direction Directional hypothesis control.
#'   \describe{
#'     \item{\code{"both"} (default)}{Non-directional; evaluates all directions.
#'       Backward-compatible synonym: \code{"off"}.}
#'     \item{\code{"greater"}, \code{"less"}}{MPE Chapter 2 binary ordered
#'       directional ODA. \code{"greater"}: high attribute values predict
#'       class 1 (MegaODA \code{DIRECTION < 0 1}). \code{"less"}: low values
#'       predict class 1 (MegaODA \code{DIRECTION > 0 1}). Binary/ordered
#'       attributes only; error on multiclass.}
#'     \item{\code{"ascending"}, \code{"descending"}}{MPE Chapter 4 ordered
#'       or categorical identity-map DIRECTIONAL. For multiclass ordered:
#'       constrains the segment-to-class assignment to be monotone ascending
#'       (segment s \eqn{\to} class s) or descending. For multiclass categorical
#'       with L == C: auto-creates identity or reverse \code{direction_map}.
#'       Error on binary class (use \code{"greater"}/\code{"less"} instead).}
#'   }
#' @param direction_map Named integer (or numeric) vector for fixed-partition
#'   categorical DIRECTIONAL (MPE Chapter 4). Names are attribute levels
#'   (as character); values are predicted class labels. All levels must be
#'   covered exactly once with at least two distinct target classes. When
#'   supplied, ODA evaluates only the specified mapping and skips the
#'   partition search. For binary class, values should be the original class
#'   labels (recoded automatically to 0/1 internally). For multiclass, values
#'   should be class labels 1..C. Compatible with \code{direction = "both"}
#'   (the default); do not combine with \code{"greater"} or \code{"less"}.
#' @param boundary_mode Boundary convention for multiclass ordered rules.
#' @return Named list with ok, rule, ess, pac, p_mc, loo, n_eff, attr_type, engine.
#' @export
oda_fit <- function(
    x,
    y,
    w            = NULL,
    attr_type    = c("auto","ordered","categorical","binary"),
    priors_on    = TRUE,
    K_segments   = NULL,
    degen        = FALSE,
    miss_codes   = NULL,
    missing_code = NULL,    # alias for miss_codes
    mcarlo       = TRUE,
    mc_iter      = 25000L,
    mc_target    = 0.05,
    mc_stop      = 99.9,
    mc_stopup    = 20,
    mc_seed      = NULL,
    loo          = "off",
    boundary_mode = c("megaoda_halfopen","right_closed"),
    eval_order   = c("mc_then_loo", "loo_then_mc"),
    mindenom     = 1L,
    direction    = c("both", "off", "greater", "less", "ascending", "descending"),
    direction_map = NULL
) {
  attr_type     <- match.arg(attr_type)
  boundary_mode <- match.arg(boundary_mode)
  eval_order    <- match.arg(eval_order)
  direction     <- match.arg(direction)
  if (direction == "both") direction <- "off"  # canonical synonym

  # Validate weights before any routing
  .validate_case_weights(w, length(x))

  # Resolve alias
  if (!is.null(missing_code))
    miss_codes <- unique(c(miss_codes, as.numeric(missing_code)))

  # Determine C from cleaned data
  clean <- oda_clean_xy(x, y, w, miss_codes)
  if (length(clean$y) == 0L)
    return(structure(
      list(ok = FALSE, reason = "all_missing", engine = NA_character_,
           priors_on = priors_on, miss_codes = miss_codes,
           has_weights = !is.null(w) && any(w != 1)),
      class = c("oda_fit_failed", "oda_fit")))

  C <- length(unique(as.integer(clean$y)))

  if (C < 2L)
    return(structure(
      list(ok = FALSE, reason = "pure_node", engine = NA_character_,
           priors_on = priors_on, miss_codes = miss_codes,
           has_weights = !is.null(w) && any(w != 1)),
      class = c("oda_fit_failed", "oda_fit")))

  if (C == 2L) {
    # ascending/descending are Chapter 4 multiclass; error for binary
    if (direction %in% c("ascending", "descending"))
      stop("direction = '", direction, "' is for multiclass ordered ODA ",
           "(MPE Chapter 4). For binary ordered ODA use direction = 'greater' ",
           "or 'less' (MPE Chapter 2).", call. = FALSE)

    # Binary engine requires y in {0, 1}. Recode from arbitrary {a, b} → {0, 1}
    # so the caller's label space (e.g. {1,2}) is transparent to unioda_core.
    bin_labels <- sort(unique(as.integer(clean$y)))   # e.g. c(1L, 2L)
    y_coded    <- ifelse(as.integer(y) == bin_labels[1L], 0L, 1L)
    # Propagate NAs (miss_codes handled inside clean; raw NAs stay NA)
    y_coded[is.na(y)] <- NA_integer_

    # Recode direction_map values from original label space to coded {0L, 1L}
    dm_coded <- NULL
    if (!is.null(direction_map)) {
      dm_vals  <- as.integer(direction_map)
      dm_coded <- stats::setNames(
        ifelse(dm_vals == bin_labels[1L], 0L, 1L),
        names(direction_map)
      )
    }

    # Numeric loo means "pvalue gate with this threshold" — do not convert to "off".
    # Validate numeric loo: must be a single finite value strictly in (0, 1).
    if (is.numeric(loo)) {
      if (length(loo) != 1L || is.na(loo) || !is.finite(loo) || loo <= 0 || loo >= 1)
        stop(
          "Numeric loo must be a single finite value strictly in (0, 1).",
          call. = FALSE
        )
    }
    # loo_alpha_val: the p-value threshold passed to oda_univariate_core().
    #   - numeric loo: threshold IS the supplied value
    #   - loo = "pvalue" string: default threshold 0.05
    #   - all other modes: 0.05 (unused by oda_univariate_core when loo = "off"/"stable")
    # Pass requires p < threshold (rejection is p >= threshold).
    loo_alpha_val <- if (is.numeric(loo)) as.double(loo) else 0.05
    loo_binary <- if (is.numeric(loo)) {
      "pvalue"
    } else {
      switch(as.character(loo),
        "off"    = "off",
        "on"     = "pvalue",
        "stable" = "stable",
        "pvalue" = "pvalue",
        "off"
      )
    }
    fit <- oda_univariate_core(
      x = x, y = y_coded, w = w,
      attr_type     = attr_type,
      priors_on     = priors_on,
      miss_codes    = miss_codes,
      loo           = loo_binary,
      loo_alpha     = loo_alpha_val,
      mcarlo        = isTRUE(mcarlo),
      mc_iter       = as.integer(mc_iter),
      mc_target     = mc_target,
      mc_stop       = mc_stop,
      mc_stopup     = mc_stopup,
      mc_seed       = mc_seed,
      eval_order    = eval_order,
      mindenom      = mindenom,
      direction     = direction,
      direction_map = dm_coded
    )

    # Remap rule and confusion back to original label space.
    # The binary engine returns direction "0->1" or "1->0"; remap those
    # sides to the caller's labels so CTA routing uses original values.
    if (isTRUE(fit$ok) && !is.null(fit$rule)) {
      fit$rule$label_0 <- bin_labels[1L]   # what coded 0 means in original space
      fit$rule$label_1 <- bin_labels[2L]   # what coded 1 means in original space
    }

    fit$engine      <- "binary"
    fit$priors_on   <- priors_on
    fit$miss_codes  <- miss_codes
    fit$has_weights <- !is.null(w) && any(w != 1)
    class(fit) <- c("oda_fit_binary", "oda_fit")
    return(fit)
  }

  # Multiclass engine (C >= 3)
  # "greater"/"less" are Chapter 2 binary-only; warn and ignore
  if (direction %in% c("greater", "less"))
    warning("direction = '", direction, "' (MPE Chapter 2) is for binary ordered ODA ",
            "only and will be ignored for multiclass problems. ",
            "For multiclass directional ODA use direction = 'ascending' or 'descending', ",
            "or supply a direction_map.", call. = FALSE)
  # "ascending"/"descending" are valid for multiclass — pass through
  direction_multi <- if (direction %in% c("ascending", "descending")) direction else "off"

  loo_multi <- if (loo == "off") "off" else "on"
  fit <- oda_multiclass_unioda_core(
    x = x, y = y, w = w,
    attr_type     = attr_type,
    priors_on     = priors_on,
    miss_codes    = miss_codes,
    K_segments    = K_segments,
    degen         = degen,
    mcarlo        = isTRUE(mcarlo),
    mc_iter       = as.integer(mc_iter),
    mc_target     = mc_target,
    mc_stop       = mc_stop,
    mc_stopup     = mc_stopup,
    mc_seed       = mc_seed,
    loo           = loo_multi,
    boundary_mode = boundary_mode,
    direction     = direction_multi,
    direction_map = direction_map
  )
  fit$engine        <- "multiclass"
  fit$priors_on     <- priors_on
  fit$miss_codes    <- miss_codes
  fit$has_weights   <- !is.null(w) && any(w != 1)
  fit$boundary_mode <- boundary_mode
  class(fit) <- c("oda_fit_multiclass", "oda_fit")
  return(fit)
}

# ---- cta_fit: public wrapper for CTA ----------------------------------------

cta_fit <- function(X, y, verbose = FALSE,
                    recursive        = FALSE,
                    min_n            = 30L,
                    max_depth        = 8L,
                    max_nodes        = 31L,
                    family_max_steps = 20L,
                    ...) {
  cls <- sort(unique(y[!is.na(y)]))
  if (length(cls) != 2L) {
    stop("cta_fit currently supports binary class variables only", call. = FALSE)
  }
  # recursive-only args: error if explicitly supplied with recursive = FALSE.
  if (!isTRUE(recursive)) {
    if (!missing(family_max_steps)) {
      stop("family_max_steps is only used when recursive = TRUE.", call. = FALSE)
    }
  }
  if (isTRUE(recursive)) {
    dots <- list(...)
    if ("mindenom" %in% names(dots)) {
      stop(paste0(
        "When recursive = TRUE, mindenom is selected by a per-node MDSA ",
        "family scan. Do not supply mindenom."
      ), call. = FALSE)
    }
    # Validate family_max_steps
    family_max_steps <- as.integer(family_max_steps)
    if (is.na(family_max_steps) || family_max_steps < 1L)
      stop("family_max_steps must be a positive integer.", call. = FALSE)
    return(.cta_ort_fit(
      X                = X,
      y                = y,
      w                = dots$w,
      mc_seed          = dots$mc_seed     %||% 42L,
      mc_iter          = dots$mc_iter     %||% 5000L,
      mc_stop          = dots$mc_stop     %||% 99.9,
      mc_stopup        = dots$mc_stopup   %||% 20,
      alpha_split      = dots$alpha_split %||% 0.05,
      prune_alpha      = dots$prune_alpha %||% 0.05,
      loo              = dots$loo         %||% "stable",
      min_n            = min_n,
      max_depth        = max_depth,
      max_nodes        = max_nodes,
      family_max_steps = family_max_steps,
      verbose          = verbose
    ))
  }
  oda_cta_fit(X = X, y = y, verbose = verbose, ...)
}

# ---- lort_fit: preferred explicit entry point for LORT ----------------------

#' Fit a Locally Optimal Recursive Tree (LORT)
#'
#' Preferred explicit entry point for the LORT workflow layer.  LORT is a
#' non-canonical workflow composition: at each recursive endpoint it runs a
#' full MDSA family scan (\code{\link{cta_descendant_family}}), selects the
#' min-D member, and recurses until no further structure is found or a compute
#' guard fires.  It uses canon CTA/MDSA components but is not itself a canon
#' CTA.exe behavior.
#'
#' \code{lort_fit()} is functionally equivalent to
#' \code{cta_fit(..., recursive = TRUE)}.  \code{cta_fit(..., recursive = TRUE)}
#' is retained as a legacy-compatible alias and will continue to work; prefer
#' \code{lort_fit()} for new code.  SORT and GORT are reserved and not
#' implemented.
#'
#' @param X Data frame or matrix of candidate predictor columns.
#' @param y Integer class variable vector.  Must have exactly two distinct values.
#' @param w Optional numeric case-weight vector.  Default \code{NULL} (unit weights).
#' @param mc_iter Integer; maximum Monte Carlo iterations per node.  Default \code{5000L}.
#' @param mc_seed Integer or \code{NULL}; RNG seed set once at LORT start.
#'   Child-node MDSA scans consume the stream in deterministic right-then-left
#'   traversal order without resetting the seed.  Default \code{42L}.
#' @param mc_stop Numeric; confidence bound for lower-tail early MC stopping
#'   (percent).  Default \code{99.9}.
#' @param mc_stopup Numeric; confidence bound for upper-tail early MC stopping
#'   (percent).  Default \code{20}.
#' @param alpha_split Numeric; node-level significance threshold.  Default \code{0.05}.
#' @param prune_alpha Numeric; pruning significance threshold.  Default \code{0.05}.
#' @param loo LOO gate mode per node: \code{"off"} (no gate), \code{"stable"}
#'   (MegaODA LOO STABLE; accept when |WESSL − WESS| ≤ 0.01 pp; default),
#'   \code{"pvalue"} (Fisher p strictly less than 0.05), or a single numeric
#'   in (0, 1) (Fisher p strictly less than the supplied threshold).
#' @param min_n Integer; minimum endpoint n to attempt recursion.  Endpoints
#'   smaller than \code{min_n} become terminal (stop reason \code{"min_n"}).
#'   Default \code{30L}.
#' @param max_depth Integer; safety cap on recursion depth.  Nodes at
#'   \code{depth >= max_depth} become terminal (stop reason
#'   \code{"max_depth"}).  Default \code{8L}.
#' @param max_nodes Integer; safety cap on total ORT nodes.  When node count
#'   exceeds \code{max_nodes} the current endpoint becomes terminal (stop
#'   reason \code{"max_nodes"}).  Default \code{31L}.
#' @param family_max_steps Integer; maximum MDSA family members evaluated at
#'   each recursive node.  Default \code{20L}.
#' @param verbose Logical; emit \code{[ORT]} progress messages.  Default
#'   \code{FALSE}.
#' @return A dual-tagged \code{cta_ort} / \code{cta_tree} object.
#'   \code{cta_ort}-aware methods (\code{predict.cta_ort},
#'   \code{print.cta_ort}, \code{summary.cta_ort}, \code{plot.cta_ort},
#'   \code{\link{cta_ort_node_table}}) operate on the full composite tree.
#'   \code{ort_settings$method} is always \code{"lort"}.
#' @seealso \code{\link{cta_fit}}, \code{\link{predict.cta_ort}},
#'   \code{\link{cta_ort_node_table}}, \code{\link{ort_plot_data}}
#' @examples
#' X <- data.frame(
#'   A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
#'   B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
#' )
#' y <- c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
#' fit <- lort_fit(X, y, mc_iter = 100L, mc_seed = 42L, loo = "off", min_n = 5L)
#' print(fit)
#' @export
lort_fit <- function(X, y, w = NULL,
                     mc_iter          = 5000L,
                     mc_seed          = 42L,
                     mc_stop          = 99.9,
                     mc_stopup        = 20,
                     alpha_split      = 0.05,
                     prune_alpha      = 0.05,
                     loo              = "stable",
                     min_n            = 30L,
                     max_depth        = 8L,
                     max_nodes        = 31L,
                     family_max_steps = 20L,
                     verbose          = FALSE) {
  cta_fit(
    X                = X,
    y                = y,
    w                = w,
    mc_iter          = mc_iter,
    mc_seed          = mc_seed,
    mc_stop          = mc_stop,
    mc_stopup        = mc_stopup,
    alpha_split      = alpha_split,
    prune_alpha      = prune_alpha,
    loo              = loo,
    min_n            = min_n,
    max_depth        = max_depth,
    max_nodes        = max_nodes,
    family_max_steps = family_max_steps,
    verbose          = verbose,
    recursive        = TRUE
  )
}
