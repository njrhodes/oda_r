###############################################################################
# R/sda_core.R
#
# sda_fit()  -  Structural Decomposition Analysis (SDA).
#
# SDA-1 scope:  mode = "unioda_max_ess"   -  iterative UniODA.
# SDA-4B scope: mode = "novometric_min_d"  -  per-attribute MDSA via
#               cta_descendant_family(); MINDENOM gate + p gate + min-D selection.
#   Case weights not supported; weights != NULL errors.
#   Returns class c("sda_fit", "odacore_sda").
#
# Canon reference: docs/SDA_AUTO_SDA_PLAN.md section SDA-4 Novometric Acceptance Contract.
###############################################################################

#' Run a Structural Decomposition Analysis (SDA) procedure
#'
#' Executes staged attribute-set identification on binary class data.
#' Traverses the attribute space by class, selecting the best eligible
#' attribute at each step, removing correctly classified observations, and
#' repeating on the unresolved sample until a stopping condition is met.
#' The result identifies which attributes to pass to downstream CTA or MDSA.
#'
#' @param X Data frame of candidate attribute columns.
#' @param y Integer class vector. Must have exactly two distinct values.
#' @param mode SDA mode. \code{"novometric_min_d"} (MPE-canon; per-attribute
#'   MDSA via \code{cta_descendant_family()}; requires \code{mindenom}) or
#'   \code{"unioda_max_ess"} (iterative UniODA; default for SDA-1). Must be
#'   declared explicitly; do not mix modes.
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

  # --- weights gate ----------------------------------------------------
  if (!is.null(weights))
    stop(
      "Weighted SDA is not implemented yet; use weights = NULL. ",
      "Weighted SDA requires explicit weighted removal, WESS labeling, ",
      "and weighted candidate-table semantics.",
      call. = FALSE
    )

  # --- novometric: mindenom is mandatory -------------------------------
  if (mode == "novometric_min_d") {
    if (is.null(mindenom))
      stop(
        "sda_fit(mode = 'novometric_min_d') requires explicit mindenom.\n",
        "MINDENOM must be tied to statistical power (MPE):\n",
        "  moderate effect: N \u2248 32 per class per category \u2192 mindenom = 64\n",
        "  strong effect:   N \u2248 12 per class per category \u2192 mindenom = 24\n",
        "Do not pass mindenom = 1 and claim power requirements are satisfied.",
        call. = FALSE
      )
    mindenom <- as.integer(mindenom)
  }

  # --- unioda_max_ess: mindenom ignored (warn if supplied) -------------
  if (mode == "unioda_max_ess" && !is.null(mindenom))
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
    # "allow"  -  keep as-is
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
    collinearity   = collinearity,
    mindenom       = if (mode == "novometric_min_d") as.integer(mindenom) else NA_integer_
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

    # run one SDA step (dispatch on mode)
    if (mode == "novometric_min_d") {
      step_result <- .sda_step_novometric_min_d(
        X            = X,
        y            = y,
        candidates   = candidates,
        active_rows  = active_rows,
        class_levels = class_levels,
        settings     = settings
      )
    } else {
      step_result <- .sda_step_unioda_max_ess(
        X            = X,
        y            = y,
        candidates   = candidates,
        active_rows  = active_rows,
        class_levels = class_levels,
        attr_types   = resolved_types,
        settings     = settings
      )
    }

    # stopping: no eligible candidate at this step
    if (is.null(step_result$winner_attr)) {
      stop_reason <- step_result$stop_reason_hint %||% "p_gate"
      break
    }

    winner_attr <- step_result$winner_attr

    # ---- novometric path ------------------------------------------------
    if (mode == "novometric_min_d") {
      wr          <- step_result$winner_result   # list(family, min_d_idx, ess, d, p_mc, ...)
      best_family <- wr$family
      best_midx   <- wr$min_d_idx
      best_tree   <- best_family$members[[best_midx]]$tree

      if (verbose)
        message(sprintf("[SDA]   selected: %s (D=%.4f, ESS=%.4f%%, p=%.4f)",
                        winner_attr,
                        wr$d %||% NA_real_,
                        wr$ess %||% NA_real_,
                        wr$p_mc %||% NA_real_))

      # apply best CTA tree to active rows
      X_active_single  <- X[active_rows, winner_attr, drop = FALSE]
      pred_orig        <- predict(best_tree, X_active_single, missing_action = "na")

      removal <- .sda_remove_correctly_classified(
        y_true      = y[active_rows],
        y_pred      = pred_orig,
        active_rows = active_rows
      )

      n_correct   <- length(removal$correct_indices_global)
      n_incorrect <- length(removal$incorrect_indices_global)

      rule_summary <- list(
        type      = "novometric_min_d",
        mindenom  = settings$mindenom,
        min_d_idx = best_midx,
        strata    = best_family$members[[best_midx]]$strata,
        ess       = wr$ess,
        d         = wr$d
      )

      steps[[step_id]] <- list(
        step_id                   = step_id,
        attribute                 = winner_attr,
        mode                      = mode,
        n_in                      = n_active,
        class_counts_in           = class_tab,
        model                     = best_family,   # full cta_family object
        min_d_idx                 = best_midx,
        rule_summary              = rule_summary,
        ess                       = wr$ess,
        d                         = wr$d,
        p_mc                      = wr$p_mc,
        ge_count                  = wr$ge_count  %||% NA_integer_,
        iter_used                 = wr$iter_used %||% NA_integer_,
        mindenom                  = settings$mindenom,
        n_correct                 = n_correct,
        n_incorrect               = n_incorrect,
        correct_indices_global    = removal$correct_indices_global,
        incorrect_indices_global  = removal$incorrect_indices_global,
        selected                  = TRUE,
        reason                    = "min_d",
        candidate_table           = step_result$candidate_table
      )

    } else {
      # ---- unioda_max_ess path --------------------------------------------
      winner_fit  <- step_result$winner_fit

      if (verbose)
        message(sprintf("[SDA]   selected: %s (ESS=%.4f, p=%.4f)",
                        winner_attr,
                        winner_fit$ess %||% NA_real_,
                        winner_fit$p_mc %||% NA_real_))

      # apply winner rule to active rows -> compute removal
      x_active   <- X[active_rows, winner_attr]
      pred_coded <- oda_rule_predict(x_active, winner_fit$rule)
      pred_orig  <- ifelse(pred_coded == 0L,
                           winner_fit$rule$label_0,
                           winner_fit$rule$label_1)

      removal <- .sda_remove_correctly_classified(
        y_true      = y[active_rows],
        y_pred      = pred_orig,
        active_rows = active_rows
      )

      n_correct   <- length(removal$correct_indices_global)
      n_incorrect <- length(removal$incorrect_indices_global)

      rule_summary <- list(
        type      = winner_fit$rule$type,
        cut_value = winner_fit$rule$cut_value %||% NULL,
        direction = winner_fit$rule$direction %||% NULL,
        label_0   = winner_fit$rule$label_0,
        label_1   = winner_fit$rule$label_1
      )

      steps[[step_id]] <- list(
        step_id                   = step_id,
        attribute                 = winner_attr,
        mode                      = mode,
        n_in                      = n_active,
        class_counts_in           = class_tab,
        model                     = winner_fit,
        min_d_idx                 = NA_integer_,
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
    }

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

# ---------------------------------------------------------------------------
# novometric_min_d step engine
# ---------------------------------------------------------------------------

#' Run one novometric_min_d SDA step via per-attribute MDSA.
#'
#' For each candidate attribute, runs \code{cta_descendant_family()} on the
#' active working sample starting at \code{settings$mindenom}. Applies the
#' MINDENOM gate (Axiom 1) then the p gate, and selects the eligible candidate
#' with minimum D, using the tie-breaking hierarchy from the SDA-4A contract.
#'
#' @return Named list: \code{winner_attr}, \code{winner_result},
#'   \code{candidate_table}, \code{stop_reason_hint}.
#' @keywords internal
.sda_step_novometric_min_d <- function(X, y, candidates, active_rows,
                                        class_levels, settings) {
  y_sub    <- y[active_rows]
  mindenom <- as.integer(settings$mindenom)
  alpha    <- settings$alpha

  results <- lapply(candidates, function(attr) {
    x_sub        <- X[active_rows, attr, drop = FALSE]
    n_attr       <- as.integer(sum(!is.na(x_sub[, 1L])))
    class_counts <- tabulate(match(as.integer(y_sub), class_levels))

    X_single <- x_sub    # single-column data.frame; column name = attr

    family <- tryCatch(
      cta_descendant_family(
        X              = X_single,
        y              = y_sub,
        w              = NULL,
        start_mindenom = mindenom,
        alpha_split    = alpha,
        mc_iter        = settings$mc_iter,
        mc_stop        = settings$mc_stop,
        mc_stopup      = settings$mc_stopup,
        loo            = settings$loo
        # mc_seed intentionally absent: RNG flows from top-level set.seed()
      ),
      error = function(e) list(.error = conditionMessage(e))
    )

    # error in family construction
    if (!is.null(family$.error)) {
      return(list(
        attr = attr, family = NULL, eligible = FALSE,
        ineligible_reason = "no_model",
        n = n_attr, class_counts = class_counts,
        ess = NA_real_, d = NA_real_, p_mc = NA_real_,
        strata = NA_integer_, min_terminal_denom = NA_integer_,
        min_d_idx = NA_integer_
      ))
    }

    # MINDENOM / Axiom 1 gate: no feasible member in the family?
    midx <- family$min_d_idx
    if (is.na(midx)) {
      return(list(
        attr = attr, family = family, eligible = FALSE,
        ineligible_reason = "axiom1",
        n = n_attr, class_counts = class_counts,
        ess = NA_real_, d = NA_real_, p_mc = NA_real_,
        strata = NA_integer_, min_terminal_denom = NA_integer_,
        min_d_idx = NA_integer_
      ))
    }

    # Extract min-D family member stats
    member <- family$members[[midx]]
    ess    <- member$overall_ess %||% NA_real_
    d      <- member$d           %||% NA_real_
    strata <- member$strata      %||% NA_integer_
    min_td <- member$min_terminal_denom %||% NA_integer_

    # p_mc, ge_count, iter_used from root split node of the min-D tree
    p_mc      <- NA_real_
    ge_count  <- NA_integer_
    iter_used <- NA_integer_
    tr   <- member$tree
    if (!is.null(tr) && !isTRUE(tr$no_tree) && !is.null(tr$nodes)) {
      root_nd   <- tr$nodes[[tr$root_id]]
      if (!is.null(root_nd)) {
        p_mc      <- root_nd$p_mc      %||% NA_real_
        ge_count  <- root_nd$ge_count  %||% NA_integer_
        iter_used <- root_nd$iter_used %||% NA_integer_
      }
    }

    # p gate
    if (is.na(p_mc) || p_mc > alpha) {
      return(list(
        attr = attr, family = family, eligible = FALSE,
        ineligible_reason = "p_gate",
        n = n_attr, class_counts = class_counts,
        ess = ess, d = d, p_mc = p_mc,
        ge_count = ge_count, iter_used = iter_used,
        strata = as.integer(strata),
        min_terminal_denom = as.integer(min_td),
        min_d_idx = midx
      ))
    }

    list(
      attr = attr, family = family, eligible = TRUE,
      ineligible_reason = NA_character_,
      n = n_attr, class_counts = class_counts,
      ess = ess, d = d, p_mc = p_mc,
      ge_count = ge_count, iter_used = iter_used,
      strata = as.integer(strata),
      min_terminal_denom = as.integer(min_td),
      min_d_idx = midx
    )
  })

  # ---- winner selection: min-D with tie-breaking --------------------------
  eligible_idx <- which(vapply(results, function(r) r$eligible, logical(1)))

  winner_attr    <- NULL
  winner_result  <- NULL
  winner_pos     <- NA_integer_

  if (length(eligible_idx) > 0L) {
    d_eligible <- vapply(results[eligible_idx],
                         function(r) r$d %||% Inf, numeric(1))

    # tick-level tie detection (1e9 precision)
    d_ticks  <- vapply(d_eligible, tick, integer(1))
    min_tick <- min(d_ticks)
    tied_loc <- which(d_ticks == min_tick)   # positions within eligible_idx

    if (length(tied_loc) == 1L) {
      winner_pos <- eligible_idx[tied_loc]
    } else {
      # Tie-break 2: fewer strata
      tied_glob   <- eligible_idx[tied_loc]
      strata_tied <- vapply(results[tied_glob],
                            function(r) r$strata %||% .Machine$integer.max,
                            integer(1))
      min_s   <- min(strata_tied)
      pass_s  <- which(strata_tied == min_s)

      if (length(pass_s) == 1L) {
        winner_pos <- tied_glob[pass_s]
      } else {
        # Tie-break 3: larger min terminal denominator
        pass_s_glob <- tied_glob[pass_s]
        mtd_tied    <- vapply(results[pass_s_glob],
                              function(r) r$min_terminal_denom %||% 0L,
                              integer(1))
        max_mtd  <- max(mtd_tied)
        pass_mtd <- which(mtd_tied == max_mtd)

        # Tie-break 4: first in column order of X (pass_mtd already index-ordered)
        winner_pos <- pass_s_glob[pass_mtd[1L]]
      }
    }

    winner_attr   <- results[[winner_pos]]$attr
    winner_result <- results[[winner_pos]]
  }

  # ---- tied_objective / selected_by_tie_break flags -----------------------
  n_cands                   <- length(results)
  tied_objective_flags      <- rep(FALSE, n_cands)
  selected_by_tie_break     <- rep(FALSE, n_cands)

  if (!is.na(winner_pos) && length(eligible_idx) > 0L) {
    win_d_tick <- tick(results[[winner_pos]]$d %||% Inf)
    for (i in eligible_idx) {
      if (tick(results[[i]]$d %||% Inf) == win_d_tick)
        tied_objective_flags[i] <- TRUE
    }
    if (sum(tied_objective_flags) > 1L)
      selected_by_tie_break[winner_pos] <- TRUE
  }

  # ---- stop_reason_hint when no winner ------------------------------------
  stop_hint <- NULL
  if (is.null(winner_attr)) {
    any_p_gate <- any(vapply(results,
      function(r) identical(r$ineligible_reason, "p_gate"), logical(1)))
    stop_hint <- if (any_p_gate) "p_gate" else "axiom1_violated"
  }

  # ---- candidate table ----------------------------------------------------
  statuses <- character(n_cands)
  for (i in seq_len(n_cands)) {
    r <- results[[i]]
    statuses[i] <- if (!is.na(winner_pos) && i == winner_pos) "selected"
                   else if (r$eligible) "eligible"
                   else "ineligible"
  }

  ctab <- data.frame(
    attribute             = vapply(results, function(r) r$attr,                       character(1)),
    status                = statuses,
    eligible              = vapply(results, function(r) r$eligible,                   logical(1)),
    ineligible_reason     = vapply(results, function(r) r$ineligible_reason %||% NA_character_, character(1)),
    n                     = vapply(results, function(r) r$n,                          integer(1)),
    mindenom              = rep(mindenom, n_cands),
    min_terminal_denom    = vapply(results, function(r) r$min_terminal_denom %||% NA_integer_, integer(1)),
    ess                   = vapply(results, function(r) r$ess  %||% NA_real_,         numeric(1)),
    d                     = vapply(results, function(r) r$d    %||% NA_real_,         numeric(1)),
    p_mc                  = vapply(results, function(r) r$p_mc %||% NA_real_,         numeric(1)),
    ge_count              = vapply(results, function(r) r$ge_count  %||% NA_integer_, integer(1)),
    iter_used             = vapply(results, function(r) r$iter_used %||% NA_integer_, integer(1)),
    strata                = vapply(results, function(r) r$strata %||% NA_integer_,    integer(1)),
    selected              = vapply(seq_len(n_cands),
                                   function(i) !is.na(winner_pos) && i == winner_pos, logical(1)),
    tied_objective        = tied_objective_flags,
    selected_by_tie_break = selected_by_tie_break,
    stringsAsFactors      = FALSE
  )
  ctab$class_counts <- lapply(results, function(r) r$class_counts)

  list(
    winner_attr      = winner_attr,
    winner_result    = winner_result,
    candidate_table  = ctab,
    stop_reason_hint = stop_hint
  )
}
