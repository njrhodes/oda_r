###############################################################################
# R/sda_core.R
#
# sda_fit() — Sequential Discriminant Analysis (SDA).
#
# SDA-1 scope:
#   mode = "unioda_max_ess" only.
#   mode = "novometric_min_d" reserved; errors with not-yet-implemented.
#   Case weights not supported; weights != NULL errors.
#   Returns class c("sda_fit", "odacore_sda").
#
# Canon reference: docs/ORT_SELECTION_METHODS.md §3, docs/SDA_AUTO_SDA_PLAN.md.
###############################################################################

#' Fit a Sequential Discriminant Analysis (SDA) model
#'
#' Executes staged attribute-set identification on binary class data.
#' Evaluates candidate attributes sequentially, selecting the best eligible
#' attribute at each step, removing correctly classified observations, and
#' repeating on the unresolved sample until a stopping condition is met.
#'
#' @param X Data frame of candidate attribute columns.
#' @param y Integer class vector. Must have exactly two distinct values.
#' @param mode SDA mode. \code{"novometric_min_d"} (MPE-canon; not yet
#'   implemented in SDA-1) or \code{"unioda_max_ess"} (legacy/iterative
#'   UniODA; default for SDA-1). Must be declared explicitly; do not mix modes.
#' @param attr_types Named character vector of attribute types
#'   (\code{"ordered"}, \code{"categorical"}, \code{"binary"}), or \code{NULL}
#'   for auto-detection. Names must match column names of \code{X}.
#' @param weights Case weights. Must be \code{NULL} in SDA-1; weighted SDA is
#'   not yet implemented and will error if non-\code{NULL}.
#' @param mindenom Integer MINDENOM (novometric mode only; ignored with warning
#'   in unioda_max_ess mode).
#' @param mc_iter Maximum Monte Carlo iterations per attribute fit. Default 5000L.
#' @param mc_seed RNG seed set once before the SDA run. Default 42L.
#' @param mc_stop Lower-tail early-stop confidence (percent). Default 99.9.
#' @param mc_stopup Upper-tail early-stop confidence (percent). Default 20.
#' @param alpha Significance threshold for p-value gate. Default 0.05.
#' @param loo LOO mode passed to \code{oda_fit()}. Default \code{"off"}.
#' @param max_steps Maximum number of SDA steps (safety cap). Default \code{NULL}
#'   (no cap beyond candidate exhaustion).
#' @param min_n Minimum working-sample size. If unresolved n drops below this,
#'   stop with \code{"min_n"}. Default \code{NULL}.
#' @param min_class_n Minimum per-class count. Stop with \code{"min_class_n"}
#'   if either class falls below this. Default \code{NULL}.
#' @param remove_correct Logical. If \code{TRUE} (canonical SDA), remove
#'   correctly classified observations after each step. If \code{FALSE},
#'   diagnostic dry-run: step logic executes but working sample is not modified.
#'   Default \code{TRUE}.
#' @param collinearity How to handle duplicate candidate columns:
#'   \code{"skip"} (silent), \code{"warn"}, or \code{"allow"}. Default
#'   \code{"skip"}.
#' @param verbose Logical. Emit \code{[SDA]} progress messages. Default
#'   \code{FALSE}.
#' @return Object of class \code{c("sda_fit", "odacore_sda")}.
#' @export
sda_fit <- function(
    X,
    y,
    mode           = c("novometric_min_d", "unioda_max_ess"),
    attr_types     = NULL,
    weights        = NULL,
    mindenom       = NULL,
    mc_iter        = 5000L,
    mc_seed        = 42L,
    mc_stop        = 99.9,
    mc_stopup      = 20,
    alpha          = 0.05,
    loo            = "off",
    max_steps      = NULL,
    min_n          = NULL,
    min_class_n    = NULL,
    remove_correct = TRUE,
    collinearity   = c("skip", "warn", "allow"),
    verbose        = FALSE
) {
  mode        <- match.arg(mode)
  collinearity <- match.arg(collinearity)

  # --- mode gate -------------------------------------------------------
  if (mode == "novometric_min_d")
    stop(
      "sda_fit(mode = 'novometric_min_d') is not yet implemented. ",
      "Use mode = 'unioda_max_ess' for the current SDA-1 implementation.",
      call. = FALSE
    )

  # --- weights gate ----------------------------------------------------
  if (!is.null(weights))
    stop(
      "Weighted SDA is not implemented yet; use weights = NULL. ",
      "Weighted SDA requires explicit weighted removal, WESS labeling, ",
      "and weighted candidate-table semantics.",
      call. = FALSE
    )

  # --- mindenom warning (unioda_max_ess ignores it) --------------------
  if (!is.null(mindenom))
    warning(
      "mindenom is ignored when mode = 'unioda_max_ess' (no MINDENOM gate ",
      "in iterative UniODA mode). Supplied value discarded.",
      call. = FALSE
    )

  # --- basic validation ------------------------------------------------
  val <- .sda_validate_candidate_frame(X, y, min_class_n = min_class_n)
  if (!val$ok)
    stop("sda_fit validation failed: ", val$reason, call. = FALSE)
  if (length(val$warnings) > 0L && verbose)
    for (w in val$warnings) message("[SDA] warning: ", w)

  # --- collinearity ----------------------------------------------------
  cand_names   <- val$candidate_names
  excl_names   <- character(0)
  excl_reasons <- character(0)
  dupes        <- .sda_find_collinear(X[, cand_names, drop = FALSE])
  if (length(dupes) > 0L) {
    if (collinearity == "skip") {
      cand_names   <- setdiff(cand_names, dupes)
      excl_names   <- dupes
      excl_reasons <- rep("collinear", length(dupes))
    } else if (collinearity == "warn") {
      warning("sda_fit: duplicate/collinear columns excluded: ",
              paste(dupes, collapse = ", "), call. = FALSE)
      cand_names   <- setdiff(cand_names, dupes)
      excl_names   <- dupes
      excl_reasons <- rep("collinear", length(dupes))
    }
    # "allow" — keep as-is
  }

  if (length(cand_names) == 0L)
    stop("sda_fit: no candidate attributes remain after validation.", call. = FALSE)

  # --- resolve attr_types ----------------------------------------------
  resolved_types <- .sda_resolve_attr_types(X, cand_names, attr_types)

  # --- initialize state ------------------------------------------------
  n_total     <- nrow(X)
  class_levels <- sort(unique(as.integer(y[!is.na(y)])))
  active_rows  <- seq_len(n_total)          # global indices of unresolved obs
  candidates   <- cand_names
  steps        <- list()
  stop_reason  <- NA_character_
  diagnostics  <- list(warnings = val$warnings,
                       excluded = data.frame(name   = excl_names,
                                             reason = excl_reasons,
                                             stringsAsFactors = FALSE))

  # --- set RNG seed once -----------------------------------------------
  if (!is.null(mc_seed)) set.seed(as.integer(mc_seed))

  settings <- list(
    mode        = mode,
    mc_iter     = as.integer(mc_iter),
    mc_seed     = mc_seed,
    mc_stop     = mc_stop,
    mc_stopup   = mc_stopup,
    alpha       = alpha,
    loo         = loo,
    max_steps   = max_steps,
    min_n       = min_n,
    min_class_n = min_class_n,
    remove_correct = remove_correct,
    collinearity   = collinearity
  )

  step_cap <- if (!is.null(max_steps)) as.integer(max_steps) else .Machine$integer.max

  # --- main SDA loop ---------------------------------------------------
  for (step_id in seq_len(step_cap)) {

    n_active <- length(active_rows)

    # stopping: no candidates left
    if (length(candidates) == 0L) {
      stop_reason <- "no_candidates"; break
    }

    # stopping: min_n
    if (!is.null(min_n) && n_active < as.integer(min_n)) {
      stop_reason <- "min_n"; break
    }

    # stopping: min_class_n / class composition
    y_active      <- y[active_rows]
    class_tab     <- tabulate(match(as.integer(y_active), class_levels))
    if (any(class_tab == 0L)) {
      stop_reason <- "class_resolved"; break
    }
    if (!is.null(min_class_n) && any(class_tab < as.integer(min_class_n))) {
      stop_reason <- "min_class_n"; break
    }

    if (verbose)
      message(sprintf("[SDA] step %d: n=%d, candidates=%d (%s)",
                      step_id, n_active, length(candidates),
                      paste(candidates, collapse = ", ")))

    # run one unioda_max_ess step
    step_result <- .sda_step_unioda_max_ess(
      X            = X,
      y            = y,
      candidates   = candidates,
      active_rows  = active_rows,
      class_levels = class_levels,
      attr_types   = resolved_types,
      settings     = settings
    )

    # stopping: no eligible candidate at this step
    if (is.null(step_result$winner_attr)) {
      stop_reason <- "p_gate"; break
    }

    winner_attr <- step_result$winner_attr
    winner_fit  <- step_result$winner_fit

    if (verbose)
      message(sprintf("[SDA]   selected: %s (ESS=%.4f, p=%.4f)",
                      winner_attr,
                      winner_fit$ess %||% NA_real_,
                      winner_fit$p_mc %||% NA_real_))

    # apply winner rule to active rows → compute removal
    x_active   <- X[active_rows, winner_attr]
    pred_coded <- oda_rule_predict(x_active, winner_fit$rule)
    # remap 0/1 to original class labels
    pred_orig  <- ifelse(pred_coded == 0L,
                         winner_fit$rule$label_0,
                         winner_fit$rule$label_1)

    removal <- .sda_remove_correctly_classified(
      y_true    = y[active_rows],
      y_pred    = pred_orig,
      active_rows = active_rows
    )

    n_correct   <- length(removal$correct_indices_global)
    n_incorrect <- length(removal$incorrect_indices_global)

    # build rule_summary
    rule_summary <- list(
      type      = winner_fit$rule$type,
      cut_value = winner_fit$rule$cut_value %||% NULL,
      direction = winner_fit$rule$direction %||% NULL,
      label_0   = winner_fit$rule$label_0,
      label_1   = winner_fit$rule$label_1
    )

    # record step
    steps[[step_id]] <- list(
      step_id                   = step_id,
      attribute                 = winner_attr,
      mode                      = mode,
      n_in                      = n_active,
      class_counts_in           = class_tab,
      model                     = winner_fit,
      rule_summary              = rule_summary,
      ess                       = winner_fit$ess,
      d                         = NA_real_,
      p_mc                      = winner_fit$p_mc %||% NA_real_,
      ge_count                  = winner_fit$mc_info$ge_count %||% NA_integer_,
      iter_used                 = winner_fit$mc_info$iter_used %||% NA_integer_,
      mindenom                  = NA_integer_,
      n_correct                 = n_correct,
      n_incorrect               = n_incorrect,
      correct_indices_global    = removal$correct_indices_global,
      incorrect_indices_global  = removal$incorrect_indices_global,
      selected                  = TRUE,
      reason                    = "max_ess",
      candidate_table           = step_result$candidate_table
    )

    # update working state
    if (isTRUE(remove_correct)) {
      active_rows <- removal$incorrect_indices_global
    }
    candidates <- setdiff(candidates, winner_attr)

    # stopping: all resolved
    if (length(active_rows) == 0L) {
      stop_reason <- "all_resolved"; break
    }

    # stopping: one class fully resolved
    y_remain   <- y[active_rows]
    remain_tab <- tabulate(match(as.integer(y_remain), class_levels))
    if (any(remain_tab == 0L)) {
      stop_reason <- "class_resolved"; break
    }
  }

  if (is.na(stop_reason)) stop_reason <- "max_steps"
  if (!isTRUE(remove_correct)) stop_reason <- "dry_run"

  selected_attributes <- vapply(steps, function(s) s$attribute, character(1))

  structure(
    list(
      call                    = match.call(),
      mode                    = mode,
      outcome_name            = NA_character_,
      class_levels            = class_levels,
      n_initial               = n_total,
      n_final_unresolved      = length(active_rows),
      candidate_names_initial = cand_names,
      selected_attributes     = selected_attributes,
      steps                   = steps,
      unresolved_indices      = active_rows,
      resolved_indices        = setdiff(seq_len(n_total), active_rows),
      settings                = settings,
      stop_reason             = stop_reason,
      diagnostics             = diagnostics
    ),
    class = c("sda_fit", "odacore_sda")
  )
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Run one unioda_max_ess SDA step: evaluate all candidates, select max ESS.
#' @keywords internal
.sda_step_unioda_max_ess <- function(X, y, candidates, active_rows,
                                     class_levels, attr_types, settings) {
  y_sub <- y[active_rows]

  results <- lapply(candidates, function(attr) {
    x_sub <- X[active_rows, attr]

    # raw class counts before oda cleaning
    class_counts <- tabulate(match(as.integer(y_sub), class_levels))

    at <- attr_types[[attr]] %||% "auto"

    fit <- tryCatch(
      oda_fit(
        x         = x_sub,
        y         = y_sub,
        w         = NULL,
        attr_type = at,
        mcarlo    = TRUE,
        mc_iter   = settings$mc_iter,
        mc_stop   = settings$mc_stop,
        mc_stopup = settings$mc_stopup,
        mc_target = settings$alpha,
        loo       = settings$loo
        # mc_seed intentionally absent: RNG stream flows from top-level set.seed()
      ),
      error = function(e) list(ok = FALSE, reason = conditionMessage(e))
    )

    eligible          <- FALSE
    ineligible_reason <- NA_character_
    ess               <- NA_real_
    p_mc              <- NA_real_
    n_attr            <- as.integer(sum(!is.na(x_sub)))

    if (!isTRUE(fit$ok)) {
      r <- fit$reason %||% "no_model"
      ineligible_reason <- switch(r,
        pure_node   = "pure_node",
        all_missing = "all_missing",
        "no_model"
      )
    } else {
      ess  <- fit$ess %||% NA_real_
      p_mc <- fit$p_mc %||% NA_real_
      if (is.na(p_mc) || p_mc > settings$alpha) {
        ineligible_reason <- "p_gate"
      } else {
        eligible <- TRUE
      }
    }

    list(
      attr              = attr,
      fit               = if (isTRUE(fit$ok)) fit else NULL,
      ess               = ess,
      p_mc              = p_mc,
      n                 = n_attr,
      class_counts      = class_counts,
      eligible          = eligible,
      ineligible_reason = ineligible_reason
    )
  })

  # winner: max ESS among eligible (first-max = enumeration order tie-break)
  eligible_idx <- which(vapply(results, function(r) r$eligible, logical(1)))

  winner_attr <- NULL
  winner_fit  <- NULL
  winner_pos  <- NA_integer_

  if (length(eligible_idx) > 0L) {
    ess_vals   <- vapply(results[eligible_idx],
                         function(r) r$ess %||% -Inf, numeric(1))
    local_best <- eligible_idx[which.max(ess_vals)]
    winner_attr <- results[[local_best]]$attr
    winner_fit  <- results[[local_best]]$fit
    winner_pos  <- local_best
  }

  # build candidate table
  n_cands  <- length(results)
  statuses <- character(n_cands)
  for (i in seq_len(n_cands)) {
    r <- results[[i]]
    statuses[i] <- if (!is.na(winner_pos) && i == winner_pos) "selected"
                   else if (r$eligible) "eligible"
                   else "ineligible"
  }

  ctab <- data.frame(
    attribute         = vapply(results, function(r) r$attr,              character(1)),
    status            = statuses,
    n                 = vapply(results, function(r) r$n,                 integer(1)),
    ess               = vapply(results, function(r) r$ess  %||% NA_real_, numeric(1)),
    d                 = NA_real_,
    p_mc              = vapply(results, function(r) r$p_mc %||% NA_real_, numeric(1)),
    mindenom          = NA_integer_,
    min_ep_n          = NA_integer_,
    eligible          = vapply(results, function(r) r$eligible,          logical(1)),
    ineligible_reason = vapply(results, function(r) r$ineligible_reason %||% NA_character_,
                               character(1)),
    selected          = vapply(seq_len(n_cands),
                               function(i) !is.na(winner_pos) && i == winner_pos,
                               logical(1)),
    stringsAsFactors  = FALSE
  )
  # list column for class_counts
  ctab$class_counts <- lapply(results, function(r) r$class_counts)

  list(
    winner_attr     = winner_attr,
    winner_fit      = winner_fit,
    candidate_table = ctab
  )
}

#' Remove correctly classified observations and return global index vectors.
#' @keywords internal
.sda_remove_correctly_classified <- function(y_true, y_pred, active_rows) {
  # y_true: true class labels for active_rows (length = length(active_rows))
  # y_pred: predicted class labels in original label space (same length)
  # active_rows: global row indices
  correct_mask <- !is.na(y_pred) & (as.integer(y_pred) == as.integer(y_true))
  list(
    correct_indices_global   = active_rows[ correct_mask],
    incorrect_indices_global = active_rows[!correct_mask]
  )
}

#' Validate the candidate frame for sda_fit.
#' @keywords internal
.sda_validate_candidate_frame <- function(X, y, min_class_n = NULL) {
  warnings <- character(0)

  if (!is.data.frame(X) && !is.matrix(X))
    return(list(ok = FALSE, reason = "X must be a data.frame or matrix",
                warnings = warnings, candidate_names = character(0)))

  if (nrow(X) != length(y))
    return(list(ok = FALSE,
                reason = sprintf("nrow(X) = %d != length(y) = %d",
                                 nrow(X), length(y)),
                warnings = warnings, candidate_names = character(0)))

  if (ncol(X) == 0L)
    return(list(ok = FALSE, reason = "X has no columns",
                warnings = warnings, candidate_names = character(0)))

  cls <- sort(unique(as.integer(y[!is.na(y)])))
  if (length(cls) != 2L)
    return(list(ok = FALSE,
                reason = sprintf("y must have exactly 2 distinct values (found %d)",
                                 length(cls)),
                warnings = warnings, candidate_names = character(0)))

  if (!is.null(min_class_n)) {
    ctab <- tabulate(match(as.integer(y), cls))
    if (any(ctab < as.integer(min_class_n)))
      return(list(ok = FALSE,
                  reason = sprintf("initial class counts below min_class_n = %d",
                                   as.integer(min_class_n)),
                  warnings = warnings, candidate_names = character(0)))
  }

  # zero-variance columns: warn, keep (caller decides exclusion)
  cnames <- if (!is.null(colnames(X))) colnames(X) else paste0("V", seq_len(ncol(X)))
  for (nm in cnames) {
    col <- X[, nm]
    uv  <- unique(col[!is.na(col)])
    if (length(uv) <= 1L)
      warnings <- c(warnings, paste0("column '", nm, "' has zero variance"))
  }

  list(ok = TRUE, reason = NULL, warnings = warnings, candidate_names = cnames)
}

#' Identify duplicate (collinear) columns in X.
#' Returns character vector of column names to drop (keeps first occurrence).
#' @keywords internal
.sda_find_collinear <- function(X) {
  nms  <- colnames(X)
  seen <- list()
  drop <- character(0)
  for (nm in nms) {
    key <- paste(X[, nm], collapse = "|~|")
    if (key %in% names(seen)) {
      drop <- c(drop, nm)
    } else {
      seen[[key]] <- nm
    }
  }
  drop
}

#' Resolve attr_types for SDA candidates (fill in "auto" for unspecified).
#' @keywords internal
.sda_resolve_attr_types <- function(X, cand_names, attr_types) {
  out <- stats::setNames(rep("auto", length(cand_names)), cand_names)
  if (!is.null(attr_types)) {
    matches <- intersect(names(attr_types), cand_names)
    out[matches] <- attr_types[matches]
  }
  out
}
