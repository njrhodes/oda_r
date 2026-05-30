###############################################################################
# R/sda_interop.R
#
# SDA-2: CTA / ORT interoperability accessors.
#
# These functions extract auditable outputs from sda_fit objects and bridge
# them into cta_fit() / cta_descendant_family() workflows.
#
# Canon reference: docs/SDA_AUTO_SDA_PLAN.md section 5.
###############################################################################

#' Return the selected attribute names from an SDA procedure result
#'
#' Returns the names of attributes selected across all SDA steps, in step
#' order. This is the constrained candidate set to pass to MDSA/CTA.
#'
#' @param fit An \code{sda_fit} object.
#' @return Character vector of selected attribute names (length = number of
#'   completed SDA steps). Empty character vector if no steps completed.
#' @export
sda_selected_attributes <- function(fit) {
  .sda_check_class(fit)
  fit$selected_attributes
}

#' Return a summary table of SDA steps
#'
#' One row per completed SDA step. Columns cover the key auditability fields
#' needed to review what was selected, why, and how the working sample changed.
#'
#' @param fit An \code{sda_fit} object.
#' @return Data frame with columns: \code{step_id}, \code{attribute},
#'   \code{n_in}, \code{n_correct}, \code{n_incorrect}, \code{ess}, \code{d},
#'   \code{p_mc}, \code{mindenom}.
#' @export
sda_step_table <- function(fit) {
  .sda_check_class(fit)
  steps <- fit$steps
  if (length(steps) == 0L)
    return(data.frame(
      step_id    = integer(0),
      attribute  = character(0),
      n_in       = integer(0),
      n_correct  = integer(0),
      n_incorrect = integer(0),
      ess        = numeric(0),
      d          = numeric(0),
      p_mc       = numeric(0),
      mindenom   = integer(0),
      stringsAsFactors = FALSE
    ))

  data.frame(
    step_id     = vapply(steps, function(s) s$step_id,     integer(1)),
    attribute   = vapply(steps, function(s) s$attribute,   character(1)),
    n_in        = vapply(steps, function(s) s$n_in,        integer(1)),
    n_correct   = vapply(steps, function(s) s$n_correct,   integer(1)),
    n_incorrect = vapply(steps, function(s) s$n_incorrect, integer(1)),
    ess         = vapply(steps, function(s) s$ess   %||% NA_real_, numeric(1)),
    d           = vapply(steps, function(s) s$d     %||% NA_real_, numeric(1)),
    p_mc        = vapply(steps, function(s) s$p_mc  %||% NA_real_, numeric(1)),
    mindenom    = vapply(steps, function(s) s$mindenom %||% NA_integer_, integer(1)),
    stringsAsFactors = FALSE
  )
}

#' Return the candidate table from one or all SDA steps
#'
#' The candidate table is the primary auditability record: one row per
#' candidate attribute evaluated at a step, showing ESS, p-value, eligibility,
#' and why a candidate was rejected or selected.
#'
#' @param fit An \code{sda_fit} object.
#' @param step Integer step index, or \code{NULL} (default) to return a list
#'   of candidate tables for all steps.
#' @return If \code{step} is an integer: the candidate table data frame for
#'   that step (with an added \code{step_id} column). If \code{step = NULL}:
#'   a named list of candidate table data frames, one per step.
#' @export
sda_candidate_table <- function(fit, step = NULL) {
  .sda_check_class(fit)
  steps <- fit$steps

  if (!is.null(step)) {
    s <- as.integer(step)
    if (is.na(s) || s < 1L || s > length(steps))
      stop(sprintf("step must be an integer in 1..%d (got %s).",
                   length(steps), step), call. = FALSE)
    ctab <- steps[[s]]$candidate_table
    ctab$step_id <- s
    return(ctab[, c("step_id", setdiff(names(ctab), "step_id"))])
  }

  # all steps
  out <- vector("list", length(steps))
  for (i in seq_along(steps)) {
    ctab <- steps[[i]]$candidate_table
    ctab$step_id <- i
    out[[i]] <- ctab[, c("step_id", setdiff(names(ctab), "step_id"))]
  }
  stats::setNames(out, paste0("step_", seq_along(steps)))
}

#' Subset a data frame to the SDA-selected candidate columns
#'
#' Returns \code{X} restricted to the columns identified by
#' \code{sda_selected_attributes(fit)}. Intended to produce the constrained
#' candidate frame for \code{\link{cta_fit}} or
#' \code{\link{cta_descendant_family}}.
#'
#' @param fit An \code{sda_fit} object.
#' @param X Data frame or matrix containing at least all selected attribute
#'   columns. Extra columns are dropped silently.
#' @return Data frame with columns matching \code{sda_selected_attributes(fit)},
#'   in SDA step order.
#' @export
as_cta_candidates <- function(fit, X) {
  .sda_check_class(fit)
  sel <- sda_selected_attributes(fit)
  if (length(sel) == 0L)
    stop("sda_fit has no selected attributes; cannot subset X.", call. = FALSE)
  if (!is.data.frame(X) && !is.matrix(X))
    stop("X must be a data.frame or matrix.", call. = FALSE)
  missing_cols <- setdiff(sel, colnames(X))
  if (length(missing_cols) > 0L)
    stop("X is missing SDA-selected column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  X[, sel, drop = FALSE]
}

#' Prepare X and y for CTA using SDA-selected attributes
#'
#' Returns a named list \code{list(X_cta, y_cta)} where \code{X_cta} contains
#' only the SDA-selected attribute columns and \code{y_cta} is the full
#' outcome vector (all observations, not just unresolved).
#'
#' This matches the Path B workflow from MPE Chapter 12: SDA identifies the
#' attribute subset; MDSA/CTA receives the full sample with a constrained
#' candidate frame. SDA resolution does not restrict which observations CTA
#' sees.
#'
#' @param fit An \code{sda_fit} object.
#' @param X Data frame of predictors (all observations).
#' @param y Integer class vector (all observations).
#' @return Named list with elements \code{X_cta} (data frame, selected columns
#'   only) and \code{y_cta} (integer vector, full length).
#' @export
sda_to_cta_data <- function(fit, X, y) {
  .sda_check_class(fit)
  if (nrow(X) != length(y))
    stop(sprintf("nrow(X) = %d != length(y) = %d.", nrow(X), length(y)),
         call. = FALSE)
  list(
    X_cta = as_cta_candidates(fit, X),
    y_cta = as.integer(y)
  )
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Check that an object is an sda_fit.
#' @keywords internal
.sda_check_class <- function(fit) {
  if (!inherits(fit, "sda_fit"))
    stop("'fit' must be an sda_fit object.", call. = FALSE)
  invisible(NULL)
}
