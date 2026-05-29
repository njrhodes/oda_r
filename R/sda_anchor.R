###############################################################################
# R/sda_anchor.R
#
# SDA anchor object: typed structural carrier for SDA selection history.
#
# The sda_anchor is an anchor for future SORT (Sequentially Optimal Recursive
# Tree) workflows. It preserves stage order, selected attributes, and metadata
# from an sda_fit so that downstream staged-CTA / SORT design can consume it
# without re-running SDA.
#
# sda_anchor IS NOT:
#   - a fitting object;
#   - a propensity-score estimator;
#   - an implementation of SORT or GORT.
#
# Canon reference: docs/SDA_ANCHOR_CONTRACT.md
###############################################################################

# ---------------------------------------------------------------------------
# Default task hook
# ---------------------------------------------------------------------------

.sda_anchor_task_hook <- function() {
  list(
    hook_type             = "sda_anchor",
    task_role             = c("candidate_set", "stage_order",
                              "sort_stage_seed", "reporting_only"),
    allowed_downstream    = c("sort_future", "cta_reporting", "agent_review"),
    future_downstream     = c("sort", "gort"),
    prohibited_downstream = c("propensity_weighting", "fraud_demo"),
    requires_human_review = TRUE,
    implementation_status = "anchor_only_no_sort",
    safety_notes          = c(
      "SDA anchors stages for future SORT / staged CTA",
      "SORT / staged CTA is not implemented",
      "GORT is not implemented",
      "LORT is not SDA-anchored",
      "SDA does not estimate propensity scores",
      "Human review is required",
      "Fraud / credit-card demos are out of scope until explicitly revived"
    )
  )
}

# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

#' Construct an \code{sda_anchor} object
#'
#' Low-level constructor.  Prefer \code{\link{as_sda_anchor}} when converting
#' from an \code{sda_fit}.  Use this constructor when building an explicit /
#' manual anchor from pre-specified fields (e.g. from a published attribute
#' ordering).
#'
#' An \code{sda_anchor} is a typed structural object that carries SDA
#' selection history for future SORT (staged CTA) workflows.  It is not a
#' fitting object and does not estimate propensity scores.
#'
#' @param anchor_type Character scalar: \code{"sda_fit"} for anchors derived
#'   from a fitted \code{sda_fit} object, or \code{"explicit"} for
#'   manually-declared anchors.
#' @param source_class Character vector: class of the source object, or
#'   \code{NULL}.
#' @param source_call Language object or \code{NULL}: the call used to produce
#'   the source object.
#' @param group_levels Integer vector of class/group levels, or \code{NULL}.
#' @param selected_attributes Non-empty character vector of selected attribute
#'   names in stage order.
#' @param candidate_universe Character vector of all attributes evaluated, or
#'   \code{NULL}.
#' @param stage_table Data frame with at least columns \code{stage_id} and
#'   \code{attribute}.
#' @param branch_candidate_map Named list for SORT branch-level candidates, or
#'   \code{NULL} (reserved for future SORT).
#' @param removal_history List of per-step removal records, or \code{NULL}.
#' @param weights_used Logical. \code{FALSE} unless weighted SDA is available.
#' @param weight_summary List or \code{NULL}.
#' @param loo_mode Character scalar or \code{NULL}: LOO mode string from
#'   SDA settings.
#' @param mc_iter Integer or \code{NULL}: MC iterations from SDA settings.
#' @param mc_seed Integer or \code{NULL}: RNG seed from SDA settings.
#' @param mindenom Integer or \code{NULL}: MINDENOM from SDA settings.
#' @param alpha Numeric or \code{NULL}: significance threshold from SDA
#'   settings.
#' @param stop_reason Character scalar or \code{NA}: SDA stop reason.
#' @param reproducibility_notes Character vector.
#' @param canon_notes Character vector.
#' @param task_hook List. Machine-readable metadata for future agent/pipeline
#'   consumers.  Defaults to the standard anchor task hook (see
#'   \code{?sda_anchor}).
#'
#' @return Object of class \code{c("sda_anchor", "list")}.
#'
#' @details
#' \strong{What an SDA anchor is not:}
#' \itemize{
#'   \item It is not a propensity-score estimator.  SDA produces stage order
#'     and selected attributes, not a propensity stratification.
#'   \item It is not an implementation of SORT or GORT.  Both remain future
#'     reserved workflows.
#'   \item Explicit / manual anchors are not SDA-derived and must be labeled
#'     \code{anchor_type = "explicit"}.
#' }
#'
#' \strong{Task hook:}
#' The default \code{task_hook} marks \code{implementation_status =
#' "anchor_only_no_sort"}, lists \code{prohibited_downstream =
#' c("propensity_weighting", "fraud_demo")}, and requires human review.
#'
#' @seealso \code{\link{as_sda_anchor}}, \code{\link{validate_sda_anchor}},
#'   \code{\link{sda_fit}}
#' @export
sda_anchor <- function(
    anchor_type,
    source_class          = NULL,
    source_call           = NULL,
    group_levels          = NULL,
    selected_attributes,
    candidate_universe    = NULL,
    stage_table,
    branch_candidate_map  = NULL,
    removal_history       = NULL,
    weights_used          = FALSE,
    weight_summary        = NULL,
    loo_mode              = NULL,
    mc_iter               = NULL,
    mc_seed               = NULL,
    mindenom              = NULL,
    alpha                 = NULL,
    stop_reason           = NA_character_,
    reproducibility_notes = character(0),
    canon_notes           = character(0),
    task_hook             = .sda_anchor_task_hook()
) {
  anchor <- list(
    anchor_type           = anchor_type,
    source_class          = source_class,
    source_call           = source_call,
    group_levels          = group_levels,
    selected_attributes   = selected_attributes,
    candidate_universe    = candidate_universe,
    stage_table           = stage_table,
    branch_candidate_map  = branch_candidate_map,
    removal_history       = removal_history,
    weights_used          = weights_used,
    weight_summary        = weight_summary,
    loo_mode              = loo_mode,
    mc_iter               = mc_iter,
    mc_seed               = mc_seed,
    mindenom              = mindenom,
    alpha                 = alpha,
    stop_reason           = stop_reason,
    reproducibility_notes = reproducibility_notes,
    canon_notes           = canon_notes,
    task_hook             = task_hook
  )
  structure(anchor, class = c("sda_anchor", "list"))
}

# ---------------------------------------------------------------------------
# Generic and methods: as_sda_anchor
# ---------------------------------------------------------------------------

#' Convert an object to an \code{sda_anchor}
#'
#' Generic converter.  Methods are provided for \code{sda_fit} and
#' \code{data.frame}.  Use \code{\link{sda_anchor}} for direct construction.
#'
#' @param x Object to convert.
#' @param ... Additional arguments passed to methods.
#' @return Object of class \code{c("sda_anchor", "list")}.
#' @seealso \code{\link{sda_anchor}}, \code{\link{validate_sda_anchor}}
#' @export
as_sda_anchor <- function(x, ...) UseMethod("as_sda_anchor")

#' @rdname as_sda_anchor
#' @param x An \code{sda_fit} object.
#' @export
as_sda_anchor.sda_fit <- function(x, ...) {
  # --- build stage_table from steps -----------------------------------
  steps <- x$steps
  if (length(steps) == 0L) {
    stage_table <- data.frame(
      stage_id   = integer(0),
      attribute  = character(0),
      n_in       = integer(0),
      n_correct  = integer(0),
      n_incorrect = integer(0),
      ess        = numeric(0),
      d          = numeric(0),
      p_mc       = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    stage_table <- data.frame(
      stage_id    = vapply(steps, function(s) s$step_id,     integer(1)),
      attribute   = vapply(steps, function(s) s$attribute,   character(1)),
      n_in        = vapply(steps, function(s) s$n_in,        integer(1)),
      n_correct   = vapply(steps, function(s) s$n_correct,   integer(1)),
      n_incorrect = vapply(steps, function(s) s$n_incorrect, integer(1)),
      ess         = vapply(steps, function(s) s$ess   %||% NA_real_, numeric(1)),
      d           = vapply(steps, function(s) s$d     %||% NA_real_, numeric(1)),
      p_mc        = vapply(steps, function(s) s$p_mc  %||% NA_real_, numeric(1)),
      stringsAsFactors = FALSE
    )
  }

  # --- removal history from per-step n_correct / n_incorrect ----------
  removal_history <- if (length(steps) > 0L) {
    lapply(steps, function(s) list(
      step_id     = s$step_id,
      attribute   = s$attribute,
      n_correct   = s$n_correct,
      n_incorrect = s$n_incorrect
    ))
  } else NULL

  # --- reproducibility notes ------------------------------------------
  repro <- character(0)
  if (!is.null(x$settings$mc_seed))
    repro <- c(repro, paste0("mc_seed = ", x$settings$mc_seed))
  if (!is.null(x$settings$mc_iter))
    repro <- c(repro, paste0("mc_iter = ", x$settings$mc_iter))
  if (!is.null(x$settings$loo))
    repro <- c(repro, paste0("loo = ", x$settings$loo))

  sda_anchor(
    anchor_type           = "sda_fit",
    source_class          = class(x),
    source_call           = x$call,
    group_levels          = x$class_levels,
    selected_attributes   = x$selected_attributes,
    candidate_universe    = x$candidate_names_initial,
    stage_table           = stage_table,
    branch_candidate_map  = NULL,
    removal_history       = removal_history,
    weights_used          = FALSE,
    weight_summary        = NULL,
    loo_mode              = x$settings$loo %||% NULL,
    mc_iter               = x$settings$mc_iter %||% NULL,
    mc_seed               = x$settings$mc_seed %||% NULL,
    mindenom              = x$settings$mindenom %||% NULL,
    alpha                 = x$settings$alpha %||% NULL,
    stop_reason           = x$stop_reason %||% NA_character_,
    reproducibility_notes = repro,
    canon_notes           = c(
      "SDA-derived anchor from sda_fit object",
      "SDA is not a propensity-score estimator",
      "This anchor is for future SORT / staged CTA workflows",
      "SORT is not implemented",
      "GORT is not implemented"
    )
  )
}

#' @rdname as_sda_anchor
#' @param x A data frame with at least columns \code{stage_id} (integer) and
#'   \code{attribute} (character).  Represents an explicit / manually-declared
#'   stage table.
#' @param selected_attributes Character vector of attribute names in stage
#'   order.  Must match \code{x$attribute} entries.
#' @param candidate_universe Character vector of all candidate attributes, or
#'   \code{NULL} (defaults to \code{selected_attributes}).
#' @param group_levels Integer vector, or \code{NULL}.
#' @param canon_notes Character vector describing the source.
#' @export
as_sda_anchor.data.frame <- function(
    x,
    selected_attributes,
    candidate_universe = NULL,
    group_levels       = NULL,
    canon_notes        = c(
      "Explicit / manual anchor - user-declared stage table",
      "Not derived from sda_fit",
      "This anchor is for future SORT / staged CTA workflows",
      "SORT is not implemented",
      "GORT is not implemented"
    ),
    ...
) {
  if (!all(c("stage_id", "attribute") %in% names(x)))
    stop("data.frame anchor requires columns 'stage_id' and 'attribute'.",
         call. = FALSE)
  if (missing(selected_attributes) || !is.character(selected_attributes) ||
      length(selected_attributes) == 0L)
    stop("'selected_attributes' must be a non-empty character vector.",
         call. = FALSE)
  sda_anchor(
    anchor_type           = "explicit",
    source_class          = "data.frame",
    source_call           = NULL,
    group_levels          = group_levels,
    selected_attributes   = selected_attributes,
    candidate_universe    = candidate_universe %||% selected_attributes,
    stage_table           = x,
    branch_candidate_map  = NULL,
    removal_history       = NULL,
    weights_used          = FALSE,
    weight_summary        = NULL,
    loo_mode              = NULL,
    mc_iter               = NULL,
    mc_seed               = NULL,
    mindenom              = NULL,
    alpha                 = NULL,
    stop_reason           = NA_character_,
    reproducibility_notes = c(
      "Explicit / manual anchor: user-declared",
      "No SDA run associated with this anchor"
    ),
    canon_notes           = canon_notes
  )
}

# ---------------------------------------------------------------------------
# Validator
# ---------------------------------------------------------------------------

#' Validate an \code{sda_anchor} object
#'
#' Checks that all required fields are present and well-formed.  Errors
#' clearly on any violation so that downstream SORT / staged-CTA code can rely
#' on the contract.
#'
#' @param anchor An object to validate.
#' @param strict Logical (default \code{TRUE}).  When \code{TRUE}, errors on
#'   any violation.  When \code{FALSE}, returns a character vector of issue
#'   messages (empty if valid).
#' @return \code{anchor} invisibly (on success).
#' @export
validate_sda_anchor <- function(anchor, strict = TRUE) {
  issues <- character(0)

  # list-like
  if (!is.list(anchor))
    issues <- c(issues, "anchor must be a list.")

  # anchor_type
  if (!is.character(anchor$anchor_type) || length(anchor$anchor_type) != 1L ||
      !anchor$anchor_type %in% c("sda_fit", "explicit"))
    issues <- c(issues,
      "anchor$anchor_type must be 'sda_fit' or 'explicit'.")

  # selected_attributes non-empty character
  if (!is.character(anchor$selected_attributes) ||
      length(anchor$selected_attributes) == 0L)
    issues <- c(issues,
      "anchor$selected_attributes must be a non-empty character vector.")

  # candidate_universe character
  if (!is.null(anchor$candidate_universe) &&
      !is.character(anchor$candidate_universe))
    issues <- c(issues,
      "anchor$candidate_universe must be a character vector or NULL.")

  # selected attributes within candidate universe
  if (is.character(anchor$selected_attributes) &&
      length(anchor$selected_attributes) > 0L &&
      is.character(anchor$candidate_universe) &&
      length(anchor$candidate_universe) > 0L) {
    outside <- setdiff(anchor$selected_attributes, anchor$candidate_universe)
    if (length(outside) > 0L)
      issues <- c(issues, paste0(
        "selected_attributes not in candidate_universe: ",
        paste(outside, collapse = ", ")))
  }

  # stage_table is data.frame
  if (!is.data.frame(anchor$stage_table))
    issues <- c(issues, "anchor$stage_table must be a data.frame.")

  # stage_table has required columns
  if (is.data.frame(anchor$stage_table)) {
    missing_cols <- setdiff(c("stage_id", "attribute"),
                            names(anchor$stage_table))
    if (length(missing_cols) > 0L)
      issues <- c(issues, paste0(
        "stage_table is missing column(s): ",
        paste(missing_cols, collapse = ", ")))

    # stage order preserved
    if ("attribute" %in% names(anchor$stage_table) &&
        nrow(anchor$stage_table) > 0L &&
        is.character(anchor$selected_attributes) &&
        length(anchor$selected_attributes) > 0L &&
        nrow(anchor$stage_table) == length(anchor$selected_attributes)) {
      if (!identical(anchor$stage_table$attribute,
                     anchor$selected_attributes))
        issues <- c(issues,
          "stage_table$attribute order does not match selected_attributes.")
    }
  }

  # weights_used is scalar logical
  if (!is.logical(anchor$weights_used) ||
      length(anchor$weights_used) != 1L ||
      is.na(anchor$weights_used))
    issues <- c(issues,
      "anchor$weights_used must be TRUE or FALSE (non-NA logical scalar).")

  # explicit anchor must not claim to be SDA-derived
  if (isTRUE(anchor$anchor_type == "explicit")) {
    sda_claim <- grepl("SDA-derived", anchor$canon_notes, fixed = TRUE)
    if (any(sda_claim))
      issues <- c(issues,
        "Explicit anchor must not claim 'SDA-derived' in canon_notes.")
  }

  # SDA anchor must not claim propensity weighting
  if (isTRUE(anchor$anchor_type == "sda_fit")) {
    prop_claim <- grepl("propensity weight", anchor$canon_notes,
                        ignore.case = TRUE)
    if (any(prop_claim))
      issues <- c(issues,
        "SDA anchor canon_notes must not claim propensity weighting.")
  }

  if (length(issues) > 0L) {
    if (strict)
      stop("sda_anchor validation failed:\n  ",
           paste(issues, collapse = "\n  "), call. = FALSE)
    return(issues)
  }
  invisible(anchor)
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print an \code{sda_anchor}
#'
#' Prints a concise summary: anchor type, number of stages, selected
#' attributes, implementation status.  Does not claim SORT or GORT are
#' implemented.
#'
#' @param x An \code{sda_anchor} object.
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @export
print.sda_anchor <- function(x, ...) {
  cat("SDA Anchor\n")
  cat("  anchor_type      :", x$anchor_type, "\n")
  cat("  stages           :", length(x$selected_attributes), "\n")
  cat("  selected_attrs   :", paste(x$selected_attributes, collapse = ", "), "\n")
  cat("  stop_reason      :", x$stop_reason %||% "(none)", "\n")
  cat("  weights_used     :", x$weights_used, "\n")
  status <- x$task_hook$implementation_status %||% "unknown"
  cat("  status           :", status, "\n")
  cat("  [SDA anchors future SORT / staged CTA; SORT is not implemented]\n")
  invisible(x)
}

#' Summarise an \code{sda_anchor}
#'
#' Returns a named list with the key structural fields needed to audit the
#' anchor or pass it to future SORT / staged-CTA pipelines.
#'
#' @param object An \code{sda_anchor} object.
#' @param ... Ignored.
#' @return Named list with fields: \code{anchor_type}, \code{n_stages},
#'   \code{selected_attributes}, \code{candidate_universe},
#'   \code{group_levels}, \code{stop_reason}, \code{weights_used},
#'   \code{loo_mode}, \code{mc_iter}, \code{mc_seed}, \code{mindenom},
#'   \code{alpha}, \code{stage_table}, \code{canon_notes},
#'   \code{implementation_status}, \code{safety_notes}.
#' @export
summary.sda_anchor <- function(object, ...) {
  out <- list(
    anchor_type           = object$anchor_type,
    n_stages              = length(object$selected_attributes),
    selected_attributes   = object$selected_attributes,
    candidate_universe    = object$candidate_universe,
    group_levels          = object$group_levels,
    stop_reason           = object$stop_reason,
    weights_used          = object$weights_used,
    loo_mode              = object$loo_mode,
    mc_iter               = object$mc_iter,
    mc_seed               = object$mc_seed,
    mindenom              = object$mindenom,
    alpha                 = object$alpha,
    stage_table           = object$stage_table,
    canon_notes           = object$canon_notes,
    implementation_status = object$task_hook$implementation_status,
    safety_notes          = object$task_hook$safety_notes
  )
  structure(out, class = "sda_anchor_summary")
}

#' @export
print.sda_anchor_summary <- function(x, ...) {
  cat("SDA Anchor Summary\n")
  cat("  type             :", x$anchor_type, "\n")
  cat("  stages           :", x$n_stages, "\n")
  cat("  attributes       :", paste(x$selected_attributes, collapse = ", "), "\n")
  if (!is.null(x$candidate_universe))
    cat("  universe         :", length(x$candidate_universe), "attributes\n")
  cat("  weights_used     :", x$weights_used, "\n")
  cat("  stop_reason      :", x$stop_reason %||% "(none)", "\n")
  cat("  status           :", x$implementation_status %||% "unknown", "\n")
  if (!is.null(x$stage_table) && nrow(x$stage_table) > 0L) {
    cat("  Stage table      :\n")
    print(x$stage_table, row.names = FALSE)
  }
  invisible(x)
}
