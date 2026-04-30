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
#' @param loo LOO mode: "off", "on" (multiclass), "stable", or "pvalue" (binary).
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
    boundary_mode = c("megaoda_halfopen","right_closed")
) {
  attr_type     <- match.arg(attr_type)
  boundary_mode <- match.arg(boundary_mode)

  # Resolve alias
  if (!is.null(missing_code))
    miss_codes <- unique(c(miss_codes, as.numeric(missing_code)))

  # Determine C from cleaned data
  clean <- oda_clean_xy(x, y, w, miss_codes)
  if (length(clean$y) == 0L)
    return(list(ok = FALSE, reason = "all_missing", engine = NA_character_))

  C <- length(unique(as.integer(clean$y)))

  if (C < 2L)
    return(list(ok = FALSE, reason = "pure_node", engine = NA_character_))

  if (C == 2L) {
    # Binary engine requires y in {0, 1}. Recode from arbitrary {a, b} → {0, 1}
    # so the caller's label space (e.g. {1,2}) is transparent to unioda_core.
    bin_labels <- sort(unique(as.integer(clean$y)))   # e.g. c(1L, 2L)
    y_coded    <- ifelse(as.integer(y) == bin_labels[1L], 0L, 1L)
    # Propagate NAs (miss_codes handled inside clean; raw NAs stay NA)
    y_coded[is.na(y)] <- NA_integer_

    loo_binary <- switch(as.character(loo),
      "off"    = "off",
      "on"     = "pvalue",
      "stable" = "stable",
      "pvalue" = "pvalue",
      "off"
    )
    fit <- oda_univariate_core(
      x = x, y = y_coded, w = w,
      attr_type    = attr_type,
      priors_on    = priors_on,
      miss_codes   = miss_codes,
      loo          = loo_binary,
      mcarlo       = isTRUE(mcarlo),
      mc_iter      = as.integer(mc_iter),
      mc_target    = mc_target,
      mc_stop      = mc_stop,
      mc_stopup    = mc_stopup,
      mc_seed      = mc_seed
    )

    # Remap rule and confusion back to original label space.
    # The binary engine returns direction "0->1" or "1->0"; remap those
    # sides to the caller's labels so CTA routing uses original values.
    if (isTRUE(fit$ok) && !is.null(fit$rule)) {
      fit$rule$label_0 <- bin_labels[1L]   # what coded 0 means in original space
      fit$rule$label_1 <- bin_labels[2L]   # what coded 1 means in original space
    }

    fit$engine <- "binary"
    return(fit)
  }

  # Multiclass engine (C >= 3)
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
    boundary_mode = boundary_mode
  )
  fit$engine <- "multiclass"
  return(fit)
}
