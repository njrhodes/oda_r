###############################################################################
# R/sda_auto_plan.R
#
# auto_sda_plan() -- dry-run planning and validation layer for SDA.
#
# SDA-3 scope:
#   dry_run = TRUE only; no fitting.
#   No weights. No propensity weights. No causal/timing inference.
#   Exact-duplicate detection only (collinearity_threshold = 1.0).
#   Returns class c("auto_sda_plan", "odacore_plan").
#
# Canon reference: docs/SDA_AUTO_SDA_PLAN.md Sec. 6.
###############################################################################

#' Dry-run planning and validation layer for SDA
#'
#' Validates and constructs a candidate set for \code{\link{sda_fit}} without
#' fitting. Returns an auditable plan object that records which columns were
#' accepted, which were excluded and why, and what settings would be passed to
#' \code{sda_fit()}.
#'
#' \strong{Agent principle:} \code{auto_sda_plan()} proposes and validates.
#' It does not silently decide causal validity, temporal ordering, exposure
#' roles, or outcome roles.  If temporal or causal structure is required,
#' declare it via \code{role_map}, \code{time_map}, or \code{stage_map}.
#'
#' @param data A data frame.
#' @param outcome Character scalar: name of the binary outcome column in
#'   \code{data}.
#' @param candidates Character vector of candidate column names, or \code{NULL}
#'   (default) to use all non-outcome columns after exclusions.
#' @param exclude Character vector of column names to force-exclude from
#'   candidates regardless of other checks.
#' @param role_map Named list mapping column names to declared roles:
#'   \code{"assignment_mechanism"}, \code{"outcome"}, \code{"id"},
#'   or \code{"leakage"}.  Columns declared as \code{"id"} or \code{"leakage"}
#'   are excluded from the candidate set.
#' @param time_map Named numeric/integer vector mapping column names to time
#'   indices. Columns with \code{time_map} value greater than the outcome's
#'   time index generate a warning flagging potential leakage. Temporal
#'   validity is a scientific judgment; auto_sda_plan flags but does not
#'   auto-exclude on this basis.
#' @param stage_map Named integer vector mapping column names to stage
#'   assignments. Stored for downstream use; not used for exclusions.
#' @param attr_types Named character vector mapping column names to declared
#'   attribute types (\code{"ordered"}, \code{"categorical"}, \code{"binary"}).
#'   Overrides type inference for named columns.
#' @param collinearity_threshold Numeric threshold for collinearity detection.
#'   Default \code{1.0} detects exact-duplicate columns only.
#' @param min_n Passed through to \code{proposed_call} for \code{sda_fit()}.
#' @param min_class_n Passed through to \code{proposed_call}.
#' @param mode SDA mode: \code{"unioda_max_ess"} (legacy/iterative UniODA) or
#'   \code{"novometric_min_d"} (MPE-canon; per-attribute MDSA via
#'   \code{cta_descendant_family()}; requires \code{mindenom}).
#' @param dry_run Logical. Must be \code{TRUE} (default). Fitting is not
#'   performed in SDA-3; \code{dry_run = FALSE} errors.
#' @return Object of class \code{c("auto_sda_plan", "odacore_plan")}.
#' @export
auto_sda_plan <- function(
    data,
    outcome,
    candidates             = NULL,
    exclude                = NULL,
    role_map               = NULL,
    time_map               = NULL,
    stage_map              = NULL,
    attr_types             = NULL,
    collinearity_threshold = 1.0,
    min_n                  = NULL,
    min_class_n            = NULL,
    mode                   = c("unioda_max_ess", "novometric_min_d"),
    dry_run                = TRUE
) {
  mode <- match.arg(mode)

  # --- dry_run gate ---------------------------------------------------------
  if (!isTRUE(dry_run))
    stop("auto_sda_plan() currently supports dry_run = TRUE only.",
         call. = FALSE)

  # --- validate data --------------------------------------------------------
  if (!is.data.frame(data) && !is.matrix(data))
    stop("'data' must be a data.frame.", call. = FALSE)
  if (!is.data.frame(data)) data <- as.data.frame(data)

  all_names <- colnames(data)
  if (is.null(all_names))
    stop("'data' must have column names.", call. = FALSE)

  n_rows <- nrow(data)
  n_cols <- ncol(data)

  # --- validate outcome -----------------------------------------------------
  outcome <- as.character(outcome)[1L]
  if (!(outcome %in% all_names))
    stop(sprintf("outcome '%s' not found in data.", outcome), call. = FALSE)

  # --- class counts ---------------------------------------------------------
  y_raw   <- as.integer(data[[outcome]])
  cls     <- sort(unique(y_raw[!is.na(y_raw)]))
  if (length(cls) != 2L)
    warning(sprintf(
      "outcome '%s' has %d distinct values; SDA expects exactly 2.",
      outcome, length(cls)), call. = FALSE)
  class_counts          <- tabulate(match(y_raw, cls))
  names(class_counts)   <- as.character(cls)

  # --- mutable exclusion state ----------------------------------------------
  excl_names   <- character(0)
  excl_reasons <- character(0)
  plan_warnings <- character(0)

  add_excl <- function(nms, reason) {
    excl_names   <<- c(excl_names,   nms)
    excl_reasons <<- c(excl_reasons, rep(reason, length(nms)))
  }

  # --- build initial pool ---------------------------------------------------
  if (is.null(candidates)) {
    pool <- setdiff(all_names, outcome)
    add_excl(outcome, "is_outcome")
  } else {
    candidates <- as.character(candidates)
    unknown_c  <- setdiff(candidates, all_names)
    if (length(unknown_c) > 0L)
      stop(sprintf("candidates not found in data: %s",
                   paste(unknown_c, collapse = ", ")), call. = FALSE)
    if (outcome %in% candidates) {
      plan_warnings <- c(plan_warnings,
                         paste0("outcome '", outcome,
                                "' was listed in candidates and has been removed"))
    }
    pool <- setdiff(candidates, outcome)
    add_excl(outcome, "is_outcome")
  }

  # --- force exclude --------------------------------------------------------
  if (!is.null(exclude)) {
    exclude     <- as.character(exclude)
    unknown_e   <- setdiff(exclude, all_names)
    if (length(unknown_e) > 0L)
      plan_warnings <- c(plan_warnings, paste0(
        "exclude names not found in data: ",
        paste(unknown_e, collapse = ", ")))
    to_excl <- intersect(exclude, pool)
    if (length(to_excl) > 0L) {
      add_excl(to_excl, "force_excluded")
      pool <- setdiff(pool, to_excl)
    }
  }

  # --- role_map exclusions --------------------------------------------------
  if (!is.null(role_map)) {
    unknown_r <- setdiff(names(role_map), all_names)
    if (length(unknown_r) > 0L)
      plan_warnings <- c(plan_warnings, paste0(
        "role_map names not found in data: ",
        paste(unknown_r, collapse = ", ")))
    for (nm in names(role_map)) {
      if (!(nm %in% pool)) next
      role <- role_map[[nm]]
      if (identical(role, "id")) {
        add_excl(nm, "role_id")
        pool <- setdiff(pool, nm)
      } else if (identical(role, "leakage")) {
        add_excl(nm, "role_leakage")
        pool <- setdiff(pool, nm)
      }
    }
  }

  # --- time_map checks (warn only; no auto-exclusion) -----------------------
  if (!is.null(time_map)) {
    unknown_t <- setdiff(names(time_map), all_names)
    if (length(unknown_t) > 0L)
      plan_warnings <- c(plan_warnings, paste0(
        "time_map names not found in data: ",
        paste(unknown_t, collapse = ", ")))
    outcome_time <- if (outcome %in% names(time_map)) time_map[[outcome]] else NULL
    if (!is.null(outcome_time)) {
      for (nm in pool) {
        nm_time <- if (nm %in% names(time_map)) time_map[[nm]] else NULL
        if (!is.null(nm_time) && nm_time > outcome_time)
          plan_warnings <- c(plan_warnings, sprintf(
            "'%s' time index (%s) > outcome time (%s) -- potential leakage; review before including",
            nm, nm_time, outcome_time))
      }
    }
  }

  # --- stage_map checks -----------------------------------------------------
  if (!is.null(stage_map)) {
    unknown_s <- setdiff(names(stage_map), all_names)
    if (length(unknown_s) > 0L)
      plan_warnings <- c(plan_warnings, paste0(
        "stage_map names not found in data: ",
        paste(unknown_s, collapse = ", ")))
  }

  # --- invalid column types -------------------------------------------------
  invalid <- Filter(function(nm) {
    col <- data[[nm]]
    is.list(col) || (is.vector(col) && typeof(col) == "raw")
  }, pool)
  if (length(invalid) > 0L) {
    add_excl(invalid, "invalid_type")
    pool <- setdiff(pool, invalid)
  }

  # --- all-missing columns --------------------------------------------------
  all_miss <- Filter(function(nm) all(is.na(data[[nm]])), pool)
  if (length(all_miss) > 0L) {
    add_excl(all_miss, "all_missing")
    pool <- setdiff(pool, all_miss)
  }

  # --- zero-variance (constant) columns -------------------------------------
  zero_var <- Filter(function(nm) {
    uv <- unique(data[[nm]][!is.na(data[[nm]])])
    length(uv) <= 1L
  }, pool)
  if (length(zero_var) > 0L) {
    add_excl(zero_var, "zero_variance")
    pool <- setdiff(pool, zero_var)
  }

  # --- exact-duplicate detection (collinearity_threshold >= 1.0) -----------
  if (collinearity_threshold >= 1.0 && length(pool) > 1L) {
    seen  <- list()
    dupes <- character(0)
    for (nm in pool) {
      key <- paste(data[[nm]], collapse = "|~|")
      if (key %in% names(seen)) {
        dupes <- c(dupes, nm)
      } else {
        seen[[key]] <- nm
      }
    }
    if (length(dupes) > 0L) {
      plan_warnings <- c(plan_warnings, paste0(
        "exact duplicate columns detected and excluded: ",
        paste(dupes, collapse = ", ")))
      add_excl(dupes, "collinear")
      pool <- setdiff(pool, dupes)
    }
  }

  # --- infer attr_types -----------------------------------------------------
  resolved_types <- .auto_sda_infer_types(data, pool, attr_types)

  # --- build proposed_call --------------------------------------------------
  proposed_call <- list(
    mode           = mode,
    attr_types     = resolved_types,
    mc_iter        = 5000L,
    mc_seed        = 42L,
    mc_stop        = 99.9,
    mc_stopup      = 20,
    alpha          = 0.05,
    loo            = "off",
    min_n          = min_n,
    min_class_n    = min_class_n,
    remove_correct = TRUE,
    collinearity   = "skip"
  )

  # --- exclusion reasons data.frame ----------------------------------------
  excl_df <- data.frame(
    name   = excl_names,
    reason = excl_reasons,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      call              = match.call(),
      outcome           = outcome,
      mode              = mode,
      dry_run           = dry_run,
      candidate_names   = pool,
      excluded_names    = excl_names,
      exclusion_reasons = excl_df,
      attr_types        = resolved_types,
      role_map          = role_map,
      time_map          = time_map,
      stage_map         = stage_map,
      warnings          = plan_warnings,
      proposed_call     = proposed_call,
      settings          = list(
        collinearity_threshold = collinearity_threshold,
        min_n                  = min_n,
        min_class_n            = min_class_n
      ),
      n_rows        = n_rows,
      n_cols        = n_cols,
      class_counts  = class_counts
    ),
    class = c("auto_sda_plan", "odacore_plan")
  )
}

#' Print an auto_sda_plan object
#' @param x An \code{auto_sda_plan} object.
#' @param ... Unused.
#' @export
print.auto_sda_plan <- function(x, ...) {
  cat(sprintf("auto_sda_plan  [mode: %s  |  dry_run: %s]\n",
              x$mode, x$dry_run))
  cat(sprintf("  outcome: '%s'  |  n=%d  |  class_counts: %s\n",
              x$outcome, x$n_rows,
              paste(names(x$class_counts), x$class_counts, sep = "=",
                    collapse = ", ")))
  cat(sprintf("  candidates accepted: %d  |  excluded: %d\n",
              length(x$candidate_names), nrow(x$exclusion_reasons)))
  if (length(x$candidate_names) > 0L)
    cat("  Candidates:", paste(x$candidate_names, collapse = ", "), "\n")
  if (nrow(x$exclusion_reasons) > 0L) {
    cat("  Excluded:\n")
    for (i in seq_len(nrow(x$exclusion_reasons)))
      cat(sprintf("    %-30s  [%s]\n",
                  x$exclusion_reasons$name[i],
                  x$exclusion_reasons$reason[i]))
  }
  if (length(x$warnings) > 0L) {
    cat("  Warnings:\n")
    for (w in x$warnings) cat("    *", w, "\n")
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Infer attribute types for SDA candidate columns.
#' User-declared types in attr_types take precedence.
#' @keywords internal
.auto_sda_infer_types <- function(data, pool, attr_types) {
  out <- stats::setNames(character(length(pool)), pool)
  for (nm in pool) {
    if (!is.null(attr_types) && nm %in% names(attr_types)) {
      out[nm] <- attr_types[[nm]]
      next
    }
    col <- data[[nm]]
    uv  <- unique(col[!is.na(col)])
    if (is.logical(col)) {
      out[nm] <- "binary"
    } else if (is.factor(col) || is.character(col)) {
      out[nm] <- if (length(uv) == 2L) "binary" else "categorical"
    } else if (is.numeric(col) || is.integer(col)) {
      out[nm] <- if (length(uv) == 2L) "binary" else "ordered"
    } else {
      out[nm] <- "auto"
    }
  }
  out
}
