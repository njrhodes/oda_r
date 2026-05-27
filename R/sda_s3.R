###############################################################################
# R/sda_s3.R
#
# S3 methods for sda_fit objects.
#
# predict.sda_fit:  sequential selected-step application to newdata.
# print.sda_fit:    concise step-by-step summary.
# summary.sda_fit:  detailed object (returns sda_fit_summary).
###############################################################################

#' Predict from a fitted SDA model
#'
#' Applies the learned selected-step sequence to \code{newdata}. For each
#' observation, steps are applied in order; the first step whose rule
#' classifies the observation is authoritative. Observations not classified
#' by any step are returned as \code{NA} (\code{resolved = FALSE}).
#'
#' This is sequential selected-step application — it follows the learned SDA
#' structure, not a re-scan of X. It does not select a "first attribute" from
#' newdata; it replays \code{object$steps[[1]]}, \code{object$steps[[2]]}, …
#' in the order established at fit time.
#'
#' @param object A \code{sda_fit} object.
#' @param newdata Data frame or matrix. Must contain columns with names
#'   matching all selected attributes in \code{object$selected_attributes}.
#'   Extra columns are ignored (wide-newdata is supported).
#' @param type Output type. One of \code{"class"} (default), \code{"stage"},
#'   \code{"rule"}, or \code{"trace"}.
#'   \code{"propensity"} and \code{"weights"} are reserved for SDA-5 and
#'   error in SDA-1.
#' @param ... Unused.
#' @return
#' \describe{
#'   \item{\code{"class"}}{Integer vector of predicted class labels; \code{NA}
#'     for unresolved observations.}
#'   \item{\code{"stage"}}{Integer vector of step_id at which each observation
#'     was classified; \code{NA} for unresolved.}
#'   \item{\code{"rule"}}{Character vector of the selected attribute name at
#'     the classifying step; \code{NA} for unresolved.}
#'   \item{\code{"trace"}}{Data frame with one row per observation × step:
#'     \code{obs_id}, \code{step_id}, \code{attribute}, \code{classified},
#'     \code{class_pred}.}
#' }
#' @export
predict.sda_fit <- function(object, newdata, type = "class", ...) {
  type <- match.arg(type, c("class", "stage", "rule", "trace",
                             "propensity", "weights"))

  if (type %in% c("propensity", "weights"))
    stop(
      "predict.sda_fit(type = '", type, "') is not yet implemented in SDA-1. ",
      "Propensity and weight outputs are deferred to SDA-5.",
      call. = FALSE
    )

  steps <- object$steps
  if (length(steps) == 0L)
    stop("sda_fit object has no steps; cannot predict.", call. = FALSE)

  if (!is.data.frame(newdata) && !is.matrix(newdata))
    stop("newdata must be a data.frame or matrix.", call. = FALSE)

  # Name-based column lookup: require all selected attributes present.
  nd_names <- colnames(newdata)
  if (is.null(nd_names))
    stop("newdata must have column names for SDA prediction.", call. = FALSE)

  sel_attrs <- object$selected_attributes
  missing_cols <- setdiff(sel_attrs, nd_names)
  if (length(missing_cols) > 0L)
    stop(
      "newdata is missing required SDA split variable(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )

  n_new      <- nrow(newdata)
  out_class  <- rep(NA_integer_,   n_new)
  out_stage  <- rep(NA_integer_,   n_new)
  out_rule   <- rep(NA_character_, n_new)
  unresolved <- seq_len(n_new)     # indices not yet classified

  trace_rows <- if (type == "trace") vector("list", length(steps)) else NULL

  for (s in seq_along(steps)) {
    if (length(unresolved) == 0L) break

    step <- steps[[s]]
    attr <- step$attribute

    if (identical(step$mode, "novometric_min_d")) {
      # --- novometric: apply min-D CTA tree to unresolved rows
      best_tree <- step$model$members[[step$min_d_idx]]$tree
      nd_single <- newdata[unresolved, attr, drop = FALSE]
      pred_orig <- predict(best_tree, nd_single, missing_action = "na")
    } else {
      # --- unioda_max_ess: apply UniODA rule
      rule       <- step$model$rule
      x_step     <- newdata[unresolved, attr]
      pred_coded <- oda_rule_predict(x_step, rule)
      pred_orig  <- ifelse(pred_coded == 0L, rule$label_0, rule$label_1)
    }

    pred_orig <- as.integer(pred_orig)

    # classified = non-NA prediction
    classified_local  <- !is.na(pred_orig)
    classified_global <- unresolved[classified_local]

    out_class[classified_global] <- pred_orig[classified_local]
    out_stage[classified_global] <- s
    out_rule [classified_global] <- attr

    if (type == "trace") {
      trace_rows[[s]] <- data.frame(
        obs_id     = unresolved,
        step_id    = s,
        attribute  = attr,
        classified = classified_local,
        class_pred = pred_orig,
        stringsAsFactors = FALSE
      )
    }

    unresolved <- unresolved[!classified_local]
  }

  if (type == "class")  return(out_class)
  if (type == "stage")  return(out_stage)
  if (type == "rule")   return(out_rule)

  # type == "trace"
  do.call(rbind, trace_rows[!vapply(trace_rows, is.null, logical(1))])
}

#' Print an sda_fit object
#' @param x An \code{sda_fit} object.
#' @param ... Unused.
#' @export
print.sda_fit <- function(x, ...) {
  steps <- x$steps
  cat(sprintf("SDA fit  [mode: %s]\n", x$mode))
  cat(sprintf("  n_initial: %d  |  n_final_unresolved: %d\n",
              x$n_initial, x$n_final_unresolved))
  cat(sprintf("  steps: %d  |  stop_reason: %s\n",
              length(steps), x$stop_reason))
  if (length(steps) > 0L) {
    cat("  Selected sequence:\n")
    for (s in steps) {
      if (identical(s$mode, "novometric_min_d")) {
        cat(sprintf("    [%d] %-30s D=%.4f  ESS=%.2f%%  p=%.4f  n_in=%d  n_correct=%d\n",
                    s$step_id, s$attribute,
                    s$d   %||% NA_real_,
                    s$ess %||% NA_real_,
                    s$p_mc %||% NA_real_,
                    s$n_in, s$n_correct))
      } else {
        cat(sprintf("    [%d] %-30s ESS=%.2f%%  p=%.4f  n_in=%d  n_correct=%d\n",
                    s$step_id, s$attribute,
                    s$ess %||% NA_real_,
                    s$p_mc %||% NA_real_,
                    s$n_in, s$n_correct))
      }
    }
  }
  invisible(x)
}

#' Summarise an sda_fit object
#' @param object An \code{sda_fit} object.
#' @param ... Unused.
#' @export
summary.sda_fit <- function(object, ...) {
  steps <- object$steps
  step_df <- if (length(steps) > 0L) {
    data.frame(
      step_id   = vapply(steps, function(s) s$step_id,   integer(1)),
      attribute = vapply(steps, function(s) s$attribute, character(1)),
      n_in      = vapply(steps, function(s) s$n_in,      integer(1)),
      n_correct = vapply(steps, function(s) s$n_correct, integer(1)),
      ess       = vapply(steps, function(s) s$ess   %||% NA_real_, numeric(1)),
      p_mc      = vapply(steps, function(s) s$p_mc  %||% NA_real_, numeric(1)),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }

  out <- list(
    mode                = object$mode,
    n_initial           = object$n_initial,
    n_final_unresolved  = object$n_final_unresolved,
    stop_reason         = object$stop_reason,
    selected_attributes = object$selected_attributes,
    step_table          = step_df,
    settings            = object$settings
  )
  class(out) <- "sda_fit_summary"
  out
}

#' Print an sda_fit_summary object
#' @param x An \code{sda_fit_summary} object.
#' @param ... Unused.
#' @export
print.sda_fit_summary <- function(x, ...) {
  cat(sprintf("SDA summary  [mode: %s  |  stop: %s]\n",
              x$mode, x$stop_reason))
  cat(sprintf("  n_initial=%d  n_unresolved=%d  steps=%d\n",
              x$n_initial, x$n_final_unresolved, nrow(x$step_table)))
  if (nrow(x$step_table) > 0L) {
    cat("  Step table:\n")
    print(x$step_table, row.names = FALSE)
  }
  invisible(x)
}
