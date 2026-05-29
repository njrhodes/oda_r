###############################################################################
# R/production_tools.R
#
# Production operator preflight and model-specific propensity tools.
#
# Functions:
#   oda_readiness_check()    -- preflight report: validates data before fitting
#   oda_clean_missing_codes() -- replace miss-code values with NA
#   oda_validate_group()     -- validate class/group variable
#   oda_validate_weights()   -- validate weight vector (report, not error)
#   oda_infer_attr_types()   -- infer attribute type per column
#   oda_propensity_weights() -- ODA rule strata propensity weights (binary only)
#   lort_propensity_weights() -- LORT terminal strata propensity weights
#
# Propensity doctrine:
#   - These are model-family-specific functions; no generic propensity API.
#   - SDA does not produce propensity weights.
#   - SORT/GORT propensity functions are future only.
###############################################################################

# ---------------------------------------------------------------------------
# oda_clean_missing_codes
# ---------------------------------------------------------------------------

#' Replace missing-code values with NA
#'
#' Replaces all values in \code{miss_codes} with \code{replacement} (default
#' \code{NA}).  Accepts a numeric vector or a data frame.  Does not modify
#' the class variable or weight vector --- pass those separately if needed.
#'
#' @param X Numeric vector or data frame of predictors.
#' @param miss_codes Numeric vector of values to treat as missing (e.g.
#'   \code{-9}).
#' @param replacement Replacement value (default \code{NA}).
#' @return Object of the same class and dimensions as \code{X} with
#'   \code{miss_codes} values replaced.
#' @export
oda_clean_missing_codes <- function(X, miss_codes, replacement = NA) {
  if (is.null(miss_codes) || length(miss_codes) == 0L)
    return(X)
  if (is.data.frame(X)) {
    for (col in names(X))
      X[[col]][X[[col]] %in% miss_codes] <- replacement
    return(X)
  }
  X[X %in% miss_codes] <- replacement
  X
}

# ---------------------------------------------------------------------------
# oda_validate_group
# ---------------------------------------------------------------------------

#' Validate a class / group variable
#'
#' Returns a structured report list rather than erroring.  Useful as a
#' preflight check before passing \code{y} to \code{oda_fit()} or
#' \code{cta_fit()}.
#'
#' @param y Integer (or coercible to integer) class vector.
#' @param binary_only Logical; if \code{TRUE}, flags more than 2 classes as
#'   an issue (useful for UniODA workflows).
#' @return Named list with: \code{ok} (logical), \code{n_classes} (integer),
#'   \code{class_levels} (integer vector), \code{class_counts} (named integer
#'   table), \code{issues} (character vector, empty if ok).
#' @export
oda_validate_group <- function(y, binary_only = FALSE) {
  issues <- character(0)

  if (is.null(y)) {
    return(list(ok = FALSE, n_classes = NA_integer_,
                class_levels = integer(0), class_counts = integer(0),
                issues = "y is NULL"))
  }
  if (!is.numeric(y) && !is.integer(y) && !is.factor(y))
    issues <- c(issues, "y must be numeric, integer, or factor")

  # coerce to integer for counting
  y_int <- tryCatch(as.integer(y), error = function(e) NULL)
  if (is.null(y_int)) {
    issues <- c(issues, "y cannot be coerced to integer")
    return(list(ok = FALSE, n_classes = NA_integer_,
                class_levels = integer(0), class_counts = integer(0),
                issues = issues))
  }

  non_missing <- y_int[!is.na(y_int)]
  if (length(non_missing) == 0L)
    issues <- c(issues, "y has no non-missing values")

  class_levels <- sort(unique(non_missing))
  n_classes    <- length(class_levels)
  class_counts <- as.integer(table(y_int))
  names(class_counts) <- names(table(y_int))

  if (n_classes == 0L)
    issues <- c(issues, "y has no unique class values")
  else if (n_classes == 1L)
    issues <- c(issues, "y has only one unique class (degenerate)")
  else if (binary_only && n_classes > 2L)
    issues <- c(issues, paste0("y has ", n_classes, " classes; binary_only = TRUE requires exactly 2"))

  if (anyNA(y_int))
    issues <- c(issues, paste0("y has ", sum(is.na(y_int)), " NA value(s)"))

  ok <- length(issues) == 0L
  list(ok          = ok,
       n_classes   = n_classes,
       class_levels = class_levels,
       class_counts = class_counts,
       issues       = issues)
}

# ---------------------------------------------------------------------------
# oda_validate_weights
# ---------------------------------------------------------------------------

#' Validate a case weight vector
#'
#' Returns a structured report rather than throwing an error.  \code{NULL}
#' weights are valid (interpreted as unit weights) and return
#' \code{ok = TRUE}.
#'
#' @param w Numeric weight vector or \code{NULL}.
#' @param n Expected length of \code{w}.
#' @return Named list with: \code{ok} (logical), \code{issues} (character
#'   vector, empty if ok), \code{n_weights} (integer or NA),
#'   \code{range} (numeric(2) or NULL).
#' @export
oda_validate_weights <- function(w, n) {
  if (is.null(w))
    return(list(ok = TRUE, issues = character(0),
                n_weights = NA_integer_, range = NULL))

  issues <- character(0)

  if (!is.numeric(w) && !is.integer(w))
    issues <- c(issues, "'w' must be numeric or integer")

  if (is.numeric(w) || is.integer(w)) {
    if (length(w) != n)
      issues <- c(issues,
        sprintf("'w' has length %d but n = %d", length(w), n))
    if (anyNA(w))
      issues <- c(issues, "'w' contains NA or NaN")
    if (any(!is.finite(w)))
      issues <- c(issues, "'w' contains Inf or -Inf")
    if (any(w <= 0, na.rm = TRUE))
      issues <- c(issues, "'w' contains zero or negative values")
  }

  ok    <- length(issues) == 0L
  rng   <- if (ok) range(w) else NULL
  list(ok        = ok,
       issues    = issues,
       n_weights = if (ok) as.integer(length(w)) else NA_integer_,
       range     = rng)
}

# ---------------------------------------------------------------------------
# oda_infer_attr_types
# ---------------------------------------------------------------------------

#' Infer attribute types from a predictor data frame
#'
#' Uses the same type-inference logic as \code{oda_fit()} (\dQuote{auto}
#' mode) to report the likely ODA attribute type for each column.
#'
#' @param X Data frame of predictors.
#' @param miss_codes Numeric vector of missing-code values to exclude when
#'   counting unique levels (default \code{NULL}).
#' @return Data frame with one row per column in \code{X}:
#'   \code{attribute} (character), \code{inferred_type} (one of
#'   \code{"ordered"}, \code{"categorical"}, \code{"binary"}),
#'   \code{n_unique} (integer, excluding miss_codes and NA),
#'   \code{n_missing} (integer, NA count),
#'   \code{n_miss_code} (integer, miss_code hit count).
#' @seealso \code{\link{oda_fit}}, \code{\link{oda_clean_missing_codes}}
#' @export
oda_infer_attr_types <- function(X, miss_codes = NULL) {
  if (!is.data.frame(X) && !is.matrix(X))
    stop("X must be a data.frame or matrix.", call. = FALSE)
  if (is.matrix(X)) X <- as.data.frame(X)

  nms  <- colnames(X)
  n    <- nrow(X)
  k    <- ncol(X)

  attr_type_v   <- character(k)
  n_unique_v    <- integer(k)
  n_missing_v   <- integer(k)
  n_miss_code_v <- integer(k)

  for (i in seq_len(k)) {
    col         <- X[[i]]
    n_missing   <- sum(is.na(col))
    n_miss_code <- if (!is.null(miss_codes)) sum(col %in% miss_codes, na.rm = TRUE) else 0L
    clean_vals  <- col[!is.na(col)]
    if (!is.null(miss_codes)) clean_vals <- clean_vals[!clean_vals %in% miss_codes]
    n_unique    <- length(unique(clean_vals))

    # Replicate oda_resolve_attr_type logic
    if (is.factor(col) || is.character(col) || is.logical(col)) {
      inferred <- if (n_unique <= 2L) "binary" else "categorical"
    } else if (n_unique <= 2L) {
      inferred <- "binary"
    } else {
      inferred <- "ordered"
    }

    attr_type_v[i]   <- inferred
    n_unique_v[i]    <- as.integer(n_unique)
    n_missing_v[i]   <- as.integer(n_missing)
    n_miss_code_v[i] <- as.integer(n_miss_code)
  }

  data.frame(
    attribute    = nms,
    inferred_type = attr_type_v,
    n_unique     = n_unique_v,
    n_missing    = n_missing_v,
    n_miss_code  = n_miss_code_v,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# oda_readiness_check
# ---------------------------------------------------------------------------

#' Preflight readiness check for ODA / CTA analysis
#'
#' Validates a predictor frame, class vector, and optional weight vector
#' before fitting.  Returns a structured report.  Does not modify inputs.
#'
#' Flags:
#' \itemize{
#'   \item Missing class/group variable.
#'   \item Non-binary group when \code{binary_only = TRUE}.
#'   \item Non-numeric weights, wrong-length weights, NA/Inf/zero weights.
#'   \item Missing-code patterns in predictors (if \code{miss_codes} supplied).
#'   \item Constant attributes (zero variance after miss-code removal).
#'   \item Insufficient class counts (< \code{min_class_n}).
#'   \item Attribute-type uncertainty (logical/factor columns).
#' }
#'
#' @param X Data frame of predictors.
#' @param y Integer class/group vector.
#' @param w Optional numeric weight vector.
#' @param miss_codes Numeric vector of missing-code values (default
#'   \code{NULL}).
#' @param binary_only Logical; flag > 2 classes as an issue (default
#'   \code{FALSE}).
#' @param min_class_n Minimum observations per class; flags if any class is
#'   below this threshold (default \code{5L}).
#' @return Named list with:
#'   \code{ok} (logical, TRUE if no issues),
#'   \code{issues} (character vector),
#'   \code{warnings} (character vector, non-fatal),
#'   \code{n_obs} (integer),
#'   \code{group_report} (from \code{oda_validate_group()}),
#'   \code{weight_report} (from \code{oda_validate_weights()}),
#'   \code{attr_types} (from \code{oda_infer_attr_types()}),
#'   \code{constant_attrs} (character vector of constant columns).
#' @seealso \code{\link{oda_validate_group}}, \code{\link{oda_validate_weights}},
#'   \code{\link{oda_infer_attr_types}}, \code{\link{oda_clean_missing_codes}}
#' @export
oda_readiness_check <- function(X, y, w = NULL, miss_codes = NULL,
                                binary_only = FALSE, min_class_n = 5L) {
  issues   <- character(0)
  warnings <- character(0)

  # --- dimensions ---
  if (!is.data.frame(X) && !is.matrix(X))
    stop("X must be a data.frame or matrix.", call. = FALSE)
  n_obs  <- nrow(X)
  n_attr <- ncol(X)

  if (length(y) != n_obs)
    issues <- c(issues, sprintf(
      "length(y) = %d does not match nrow(X) = %d", length(y), n_obs))

  # --- group validation ---
  grp <- oda_validate_group(y, binary_only = binary_only)
  if (!grp$ok)
    issues <- c(issues, grp$issues)
  else {
    # min class count check
    if (length(grp$class_counts) > 0L) {
      low_cls <- names(grp$class_counts)[grp$class_counts < as.integer(min_class_n)]
      if (length(low_cls) > 0L)
        warnings <- c(warnings, sprintf(
          "Class(es) %s have fewer than %d observations",
          paste(low_cls, collapse = ", "), min_class_n))
    }
  }

  # --- weight validation ---
  wt <- oda_validate_weights(w, n_obs)
  if (!wt$ok)
    issues <- c(issues, wt$issues)

  # --- attr type inference (after cleaning miss_codes for the report) ---
  X_clean <- if (!is.null(miss_codes)) oda_clean_missing_codes(X, miss_codes) else X
  atypes  <- tryCatch(
    oda_infer_attr_types(X_clean, miss_codes = NULL),
    error = function(e) NULL
  )

  # --- constant attributes ---
  constant_attrs <- character(0)
  if (!is.null(atypes)) {
    for (nm in atypes$attribute) {
      col_clean <- X_clean[[nm]]
      col_clean <- col_clean[!is.na(col_clean)]
      if (length(unique(col_clean)) <= 1L)
        constant_attrs <- c(constant_attrs, nm)
    }
  } else if (is.data.frame(X_clean)) {
    for (nm in names(X_clean)) {
      col_clean <- X_clean[[nm]]
      col_clean <- col_clean[!is.na(col_clean)]
      if (length(unique(col_clean)) <= 1L)
        constant_attrs <- c(constant_attrs, nm)
    }
  }
  if (length(constant_attrs) > 0L)
    warnings <- c(warnings, sprintf(
      "Constant attribute(s) (zero variance after miss-code removal): %s",
      paste(constant_attrs, collapse = ", ")))

  # --- factor/character columns (type uncertainty) ---
  if (!is.null(atypes)) {
    fc_cols <- character(0)
    for (nm in atypes$attribute) {
      col <- X[[nm]]
      if (is.factor(col) || is.character(col))
        fc_cols <- c(fc_cols, nm)
    }
    if (length(fc_cols) > 0L)
      warnings <- c(warnings, sprintf(
        "Factor/character attribute(s) will be treated as binary/categorical: %s",
        paste(fc_cols, collapse = ", ")))
  }

  ok <- length(issues) == 0L
  list(
    ok             = ok,
    issues         = issues,
    warnings       = warnings,
    n_obs          = as.integer(n_obs),
    n_attrs        = as.integer(n_attr),
    group_report   = grp,
    weight_report  = wt,
    attr_types     = atypes,
    constant_attrs = constant_attrs
  )
}

# ---------------------------------------------------------------------------
# oda_propensity_weights
# ---------------------------------------------------------------------------

#' ODA rule strata propensity weights
#'
#' Computes propensity weights from the two rule strata (left and right of
#' the ODA cutpoint) using stored training confusion counts.  Implements the
#' Yarnold/Linden stratum-weight formula:
#' \deqn{w = n_s \times \Pr(Z=z) / n_{z,s}}
#'
#' Currently implemented for \strong{binary (C=2) ODA fits only}.
#'
#' The fitted model must have been trained with the treatment/exposure/group
#' membership as the class variable (\code{y}), not a clinical outcome.
#' The user is responsible for this labeling decision.
#'
#' @param fit An \code{oda_fit} object with \code{ok == TRUE}.
#' @param adjusted Logical; if \code{TRUE} (default), applies a
#'   one-hypothetical-misclassification adjustment when a class is absent
#'   from a rule stratum.
#' @return Data frame with one row per (stratum, class) combination:
#'   \code{stratum_id} (1L = rule predicts class 0, 2L = rule predicts
#'   class 1), \code{predicted_class} (integer), \code{class} (character),
#'   \code{class_n} (integer), \code{stratum_n} (integer),
#'   \code{marginal_class_n} (integer), \code{marginal_total_n} (integer),
#'   \code{marginal_class_probability} (numeric),
#'   \code{propensity_weight} (numeric), \code{undefined_empirical}
#'   (logical), \code{adjusted} (logical),
#'   \code{adjusted_propensity_weight} (numeric),
#'   \code{model_family} (\code{"oda"}).
#' @seealso \code{\link{cta_propensity_weights}},
#'   \code{\link{lort_propensity_weights}}
#' @export
oda_propensity_weights <- function(fit, adjusted = TRUE) {
  if (!inherits(fit, "oda_fit"))
    stop("'fit' must be an oda_fit object.", call. = FALSE)

  empty_df <- function() {
    data.frame(
      stratum_id                  = integer(0),
      predicted_class             = integer(0),
      class                       = character(0),
      class_n                     = integer(0),
      stratum_n                   = integer(0),
      marginal_class_n            = integer(0),
      marginal_total_n            = integer(0),
      marginal_class_probability  = numeric(0),
      propensity_weight           = numeric(0),
      undefined_empirical         = logical(0),
      adjusted                    = logical(0),
      adjusted_propensity_weight  = numeric(0),
      model_family                = character(0),
      stringsAsFactors            = FALSE
    )
  }

  # Must be a successful binary fit
  if (!isTRUE(fit$ok)) return(empty_df())

  conf <- fit$confusion
  if (is.null(conf) || is.null(conf$TP))
    stop("'fit' does not contain a binary confusion matrix. ",
         "oda_propensity_weights() requires a binary (C=2) ODA fit.",
         call. = FALSE)

  # Raw integer counts from the confusion matrix
  TP <- as.integer(conf$TP)
  TN <- as.integer(conf$TN)
  FP <- as.integer(conf$FP)
  FN <- as.integer(conf$FN)
  N  <- TP + TN + FP + FN

  if (N == 0L) return(empty_df())

  # Class 0 and class 1 marginal counts
  n0 <- TN + FP   # total class-0 obs
  n1 <- TP + FN   # total class-1 obs

  # Stratum 1: rule predicts class 0 (left side)
  #   contains TN class-0 obs, FN class-1 obs
  n_s0 <- TN + FN

  # Stratum 2: rule predicts class 1 (right side)
  #   contains FP class-0 obs, TP class-1 obs
  n_s1 <- FP + TP

  # Build a table: 4 rows (2 strata x 2 classes)
  stratum_ids  <- c(1L, 1L, 2L, 2L)
  pred_cls     <- c(0L, 0L, 1L, 1L)
  class_labels <- c("0", "1", "0", "1")
  class_n_v    <- c(TN, FN, FP, TP)
  stratum_n_v  <- c(n_s0, n_s0, n_s1, n_s1)
  marg_cls_n_v <- c(n0, n1, n0, n1)

  marginal_total <- N
  marg_prob_v    <- marg_cls_n_v / marginal_total

  # Empirical weights (undefined when class_n == 0)
  undefined_emp <- class_n_v == 0L
  prop_wt       <- ifelse(
    undefined_emp,
    Inf,
    as.numeric(stratum_n_v) * marg_prob_v / as.numeric(class_n_v)
  )

  # Adjustment: one hypothetical obs added to each absent cell
  needs_adj <- isTRUE(adjusted) & undefined_emp

  adj_class_n    <- as.numeric(class_n_v)
  adj_class_n[needs_adj] <- adj_class_n[needs_adj] + 1.0

  # Adjusted stratum sizes: add 1 per absent class in each stratum
  adj_per_s0 <- sum(needs_adj[stratum_ids == 1L])
  adj_per_s1 <- sum(needs_adj[stratum_ids == 2L])
  adj_sn     <- ifelse(stratum_ids == 1L,
                       as.numeric(stratum_n_v) + adj_per_s0,
                       as.numeric(stratum_n_v) + adj_per_s1)

  # Adjusted marginals: per class, add n_adj_cells for that class
  adj_adds_0 <- sum(needs_adj[class_labels == "0"])
  adj_adds_1 <- sum(needs_adj[class_labels == "1"])
  adj_marg_n <- ifelse(class_labels == "0",
                       as.numeric(marg_cls_n_v) + adj_adds_0,
                       as.numeric(marg_cls_n_v) + adj_adds_1)
  adj_total  <- marginal_total + sum(needs_adj)
  adj_marg_prob <- adj_marg_n / adj_total

  adj_prop_wt <- adj_sn * adj_marg_prob / adj_class_n

  data.frame(
    stratum_id                  = stratum_ids,
    predicted_class             = pred_cls,
    class                       = class_labels,
    class_n                     = class_n_v,
    stratum_n                   = stratum_n_v,
    marginal_class_n            = marg_cls_n_v,
    marginal_total_n            = rep(as.integer(marginal_total), 4L),
    marginal_class_probability  = marg_prob_v,
    propensity_weight           = prop_wt,
    undefined_empirical         = undefined_emp,
    adjusted                    = needs_adj,
    adjusted_propensity_weight  = adj_prop_wt,
    model_family                = rep("oda", 4L),
    stringsAsFactors            = FALSE
  )
}

# ---------------------------------------------------------------------------
# lort_propensity_weights
# ---------------------------------------------------------------------------

#' LORT terminal strata propensity weights
#'
#' Computes propensity weights from the terminal strata of a fitted LORT
#' (Locally Optimal Recursive Tree) model.  Uses stored
#' \code{class_counts} per terminal node.  Implements the Yarnold/Linden
#' stratum-weight formula (same as \code{\link{cta_propensity_weights}}):
#' \deqn{w = n_s \times \Pr(Z=z) / n_{z,s}}
#'
#' The fitted model must have been trained with the treatment/exposure/group
#' membership as the class variable, not a clinical outcome.  The user is
#' responsible for this labeling decision.
#'
#' @param ort A \code{cta_ort} object produced by \code{\link{lort_fit}} or
#'   \code{cta_fit(recursive = TRUE)}.
#' @param target_class Integer target class for annotation (optional; if
#'   \code{NULL}, defaults to the numerically higher class in binary models).
#' @param adjusted Logical; if \code{TRUE} (default), applies a
#'   one-hypothetical-misclassification adjustment when a class is absent
#'   from a terminal stratum.
#' @return Data frame with one row per (stratum, class) combination.
#'   Columns: \code{stratum_id} (integer), \code{path} (character),
#'   \code{depth} (integer), \code{stratum_n} (integer),
#'   \code{terminal_class} (integer), \code{class} (character),
#'   \code{class_n} (integer), \code{target_class} (integer),
#'   \code{marginal_class_n} (integer), \code{marginal_total_n} (integer),
#'   \code{marginal_class_probability} (numeric),
#'   \code{propensity_weight} (numeric), \code{undefined_empirical}
#'   (logical), \code{adjusted} (logical),
#'   \code{adjusted_propensity_weight} (numeric), \code{model_family}
#'   (\code{"lort"}), \code{global_optimization} (\code{FALSE}),
#'   \code{sda_anchored} (\code{FALSE}).
#' @seealso \code{\link{cta_propensity_weights}},
#'   \code{\link{oda_propensity_weights}}, \code{\link{lort_fit}}
#' @export
lort_propensity_weights <- function(ort, target_class = NULL, adjusted = TRUE) {
  if (!inherits(ort, "cta_ort"))
    stop("'ort' must be a cta_ort object.", call. = FALSE)

  empty_df <- function() {
    data.frame(
      stratum_id                  = integer(0),
      path                        = character(0),
      depth                       = integer(0),
      stratum_n                   = integer(0),
      terminal_class              = integer(0),
      class                       = character(0),
      class_n                     = integer(0),
      target_class                = integer(0),
      marginal_class_n            = integer(0),
      marginal_total_n            = integer(0),
      marginal_class_probability  = numeric(0),
      propensity_weight           = numeric(0),
      undefined_empirical         = logical(0),
      adjusted                    = logical(0),
      adjusted_propensity_weight  = numeric(0),
      model_family                = character(0),
      global_optimization         = logical(0),
      sda_anchored                = logical(0),
      stringsAsFactors            = FALSE
    )
  }

  strata <- ort$strata
  if (is.null(strata) || nrow(strata) == 0L) return(empty_df())

  # --- Unnest class_counts into long form (one row per stratum x class) ---
  rows <- vector("list", nrow(strata))
  for (i in seq_len(nrow(strata))) {
    cc  <- strata$class_counts[[i]]    # named integer vector
    if (is.null(cc) || length(cc) == 0L) next
    cls_names <- names(cc)
    rows[[i]] <- data.frame(
      stratum_id     = strata$stratum_id[i],
      path           = strata$path[i],
      depth          = strata$depth[i],
      stratum_n      = as.integer(strata$n[i]),
      terminal_class = as.integer(strata$terminal_class[i]),
      class          = cls_names,
      class_n        = as.integer(unname(cc)),
      stringsAsFactors = FALSE
    )
  }
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(empty_df())
  ec <- do.call(rbind, rows)
  rownames(ec) <- NULL

  # --- Resolve target_class (annotation only) ---
  all_classes <- sort(unique(ec$class))
  if (is.null(target_class)) {
    if (length(all_classes) != 2L)
      stop("target_class must be specified for LORT models with ",
           length(all_classes), " classes (found: ",
           paste(all_classes, collapse = ", "), ")", call. = FALSE)
    target_class_int <- as.integer(max(as.integer(all_classes)))
  } else {
    target_class_int <- as.integer(target_class)
  }
  if (!as.character(target_class_int) %in% all_classes)
    stop("target_class ", target_class_int,
         " not found in LORT classes: ",
         paste(all_classes, collapse = ", "), call. = FALSE)

  # --- Marginal class counts aligned to ec rows ---
  marg_by_cls      <- tapply(ec$class_n, ec$class, sum)
  marginal_class_n <- as.integer(marg_by_cls[ec$class])
  marginal_total   <- as.integer(sum(ec$class_n))
  marg_prob        <- marginal_class_n / marginal_total

  # --- Empirical propensity weights ---
  undefined_emp <- ec$class_n == 0L
  prop_wt       <- ifelse(
    undefined_emp,
    Inf,
    as.numeric(ec$stratum_n) * marg_prob / as.numeric(ec$class_n)
  )

  # --- Adjustment: one hypothetical obs per absent class x stratum ---
  needs_adj <- isTRUE(adjusted) & undefined_emp

  adj_class_n <- as.numeric(ec$class_n)
  adj_class_n[needs_adj] <- adj_class_n[needs_adj] + 1.0

  adj_per_stratum   <- tapply(as.integer(needs_adj), ec$stratum_id, sum)
  adj_stratum_n     <- as.numeric(ec$stratum_n) +
                       as.numeric(adj_per_stratum[as.character(ec$stratum_id)])

  adj_additions  <- tapply(as.integer(needs_adj), ec$class, sum)
  adj_marg_cls_n <- marg_by_cls + adj_additions[names(marg_by_cls)]
  adj_marg_total <- sum(adj_marg_cls_n)

  adj_marg_n_row  <- adj_marg_cls_n[ec$class]
  adj_marg_prob   <- adj_marg_n_row / adj_marg_total
  adj_prop_wt     <- adj_stratum_n * adj_marg_prob / adj_class_n

  data.frame(
    stratum_id                  = ec$stratum_id,
    path                        = ec$path,
    depth                       = ec$depth,
    stratum_n                   = ec$stratum_n,
    terminal_class              = ec$terminal_class,
    class                       = ec$class,
    class_n                     = ec$class_n,
    target_class                = rep(target_class_int, nrow(ec)),
    marginal_class_n            = marginal_class_n,
    marginal_total_n            = rep(marginal_total, nrow(ec)),
    marginal_class_probability  = marg_prob,
    propensity_weight           = prop_wt,
    undefined_empirical         = undefined_emp,
    adjusted                    = needs_adj,
    adjusted_propensity_weight  = adj_prop_wt,
    model_family                = rep("lort", nrow(ec)),
    global_optimization         = rep(FALSE, nrow(ec)),
    sda_anchored                = rep(FALSE, nrow(ec)),
    stringsAsFactors            = FALSE
  )
}
