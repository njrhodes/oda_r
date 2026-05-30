###############################################################################
# R/balance.R  -  Covariate balance analysis (univariate ODA + SMD companion)
#
# Public API:
#   oda_balance_table()      -  univariate ODA balance diagnostics
#   smd_balance_table()      -  conventional SMD companion (no p-values)
#   oda_balance_plot_data()  -  renderer-ready plot data (no fitting)
#
# Canon:
#   Linden & Yarnold (2016). Using machine learning to assess covariate
#   balance in matching studies. J Eval Clin Pract, 22(6), 861-867.
#
# Scope:
#   group   = binary class/study-arm/exposure membership (plays y in ODA)
#   X       = observed baseline covariates (play x in ODA)
#   outcome = explicitly out of scope; never pass the outcome as group or
#             as a covariate in X
#
# Implementation invariants:
#   - oda_balance_table() calls oda_fit() per covariate; it does not
#     duplicate ODA search or MC p-value logic.
#   - smd_balance_table() is descriptive arithmetic only; no ODA/CTA.
#   - oda_balance_plot_data() consumes oda_balance_table() output and an
#     optional smd_balance_table(); it does NOT call oda_fit() and does NOT
#     accept group or X arguments.
#   - SMD is a conventional companion; it is not the odacore balance objective.
#   - Deferred to v3B2: cta_balance_table(), cta_balance_plot_data().
###############################################################################

# ---- Internal helpers ------------------------------------------------------- #

# Human-readable summary of one ODA rule for balance reporting.
# attr_nm is passed explicitly because oda_fit() does not set fit$attr_name.
.oda_balance_rule_summary <- function(fit, attr_nm = "?") {
  if (is.null(fit) || !isTRUE(fit$ok) || is.null(fit$rule))
    return(NA_character_)
  rule <- fit$rule
  tryCatch({
    t <- rule$type
    if (identical(t, "ordered_cut")) {
      dir <- if (identical(rule$direction, "0->1")) "<=" else ">"
      sprintf("%s %s %g", attr_nm, dir, rule$cut_value)
    } else if (t %in% c("binary_map", "nominal_cut")) {
      lv <- paste(rule$left_levels  %||% character(0), collapse = ",")
      rv <- paste(rule$right_levels %||% character(0), collapse = ",")
      sprintf("%s: {%s} vs {%s}", attr_nm, lv, rv)
    } else {
      attr_nm
    }
  }, error = function(e) attr_nm)
}

# Sidak correction: 1 - (1 - p)^k.  NA-safe.
.sidak_p <- function(p, k) {
  if (is.na(p) || !is.finite(p) || k < 1L) return(NA_real_)
  min(1 - (1 - p)^k, 1)
}

# Bonferroni correction: p * k, capped at 1.  NA-safe.
.bonferroni_p <- function(p, k) {
  if (is.na(p) || !is.finite(p) || k < 1L) return(NA_real_)
  min(p * k, 1)
}

# ---- oda_balance_table ------------------------------------------------------ #

#' Univariate ODA covariate balance diagnostics
#'
#' Fits a univariate ODA model for each covariate in \code{X} with
#' \code{group} as the class variable.  Returns one row per covariate
#' summarizing ODA-based balance diagnostics: rule, sensitivity, specificity,
#' Mean PAC, ESS/WESS, and permutation p-value with Sidak and Bonferroni
#' multiplicity corrections.
#'
#' Balance asks whether group membership (treatment, exposure, or study arm)
#' can be predicted from observed baseline covariates.  When no covariate
#' predicts group membership above chance, the groups are considered balanced
#' on those covariates under the declared analytic constraints.
#'
#' \strong{group vs. outcome:} \code{group} is the binary class variable in
#' every ODA call.  The scientific outcome of interest is strictly out of
#' scope; do not pass the outcome as \code{group} or as a column of \code{X}.
#'
#' \strong{SMD:} conventional standardized mean difference is a companion
#' diagnostic, not the odacore balance objective.  Use
#' \code{\link{smd_balance_table}} for the conventional companion table.
#'
#' @param group Integer (or coercible) binary group indicator.  Must have
#'   exactly two distinct non-missing values.  Plays the role of the class
#'   variable (\code{y}) in ODA.
#' @param X Data frame of baseline covariate columns.  Plays the role of
#'   attributes (\code{x}) in ODA.  Must have \code{nrow(X) == length(group)}.
#' @param w Optional numeric case-weight vector (length \code{length(group)}).
#'   When supplied, ODA fits use case weights and \code{ess_label = "WESS"}.
#' @param covariate_types Optional named character vector mapping column names
#'   to ODA attribute types (\code{"auto"}, \code{"ordered"},
#'   \code{"categorical"}, \code{"binary"}).  Unmapped columns use
#'   \code{"auto"}.
#' @param loo LOO gate mode passed to each \code{\link{oda_fit}} call.
#'   Default \code{"off"}.  See \code{\link{oda_fit}} for options.
#' @param mcarlo Logical; run Monte Carlo permutation p-value?  Default
#'   \code{TRUE}.  Set \code{FALSE} to skip p-value computation.
#' @param mc_iter Integer; maximum MC iterations per covariate.  Default
#'   \code{1000L}.  Higher values (e.g., \code{25000L}) are recommended for
#'   final published analyses but increase runtime.
#' @param alpha Numeric significance threshold for the \code{significant}
#'   flag.  Default \code{0.05}.
#' @param adjust Character; which p-value drives the primary \code{significant}
#'   column: \code{"none"} (raw p, default), \code{"sidak"}, or
#'   \code{"bonferroni"}.  All three significance flags are always returned
#'   regardless of this choice.
#' @param ... Additional arguments forwarded to each \code{\link{oda_fit}}
#'   call (e.g., \code{mc_seed}, \code{priors_on}, \code{mc_stop}).
#' @return A list of class \code{"oda_balance_table"} with elements:
#' \describe{
#'   \item{\code{rows}}{Data frame; one row per covariate.  Key columns:
#'     \code{attribute}, \code{attr_type}, \code{n_total},
#'     \code{n_group_0}, \code{n_group_1}, \code{sensitivity},
#'     \code{specificity}, \code{mean_pac}, \code{ess}, \code{wess},
#'     \code{ess_display} (operative measure), \code{p_mc},
#'     \code{p_sidak}, \code{p_bonferroni}, \code{significant_raw},
#'     \code{significant_sidak}, \code{significant_bonferroni},
#'     \code{significant} (driven by \code{adjust}), \code{rule_type},
#'     \code{rule_summary}, \code{loo_status}, \code{ess_loo},
#'     \code{has_weights}, \code{fit_ok}, \code{fit_reason}.}
#'   \item{\code{meta}}{List of metadata: \code{n_covariates},
#'     \code{n_obs}, \code{has_weights}, \code{ess_label},
#'     \code{alpha}, \code{adjust}, \code{k_valid} (number of covariates
#'     with valid p_mc used for multiplicity correction), \code{loo_mode},
#'     \code{mcarlo}, \code{mc_iter}.}
#' }
#' @references
#' Linden A, Yarnold PR (2016). Using machine learning to assess covariate
#' balance in matching studies. \emph{Journal of Evaluation in Clinical
#' Practice}, \strong{22}(6), 861-867.
#' @seealso \code{\link{smd_balance_table}}, \code{\link{oda_balance_plot_data}},
#'   \code{\link{oda_fit}}
#' @examples
#' set.seed(1)
#' n <- 60
#' group <- c(rep(0L, 30), rep(1L, 30))
#' X <- data.frame(
#'   age    = c(rnorm(30, 45, 8), rnorm(30, 55, 8)),   # imbalanced
#'   score  = rnorm(60, 50, 10)                         # balanced
#' )
#' bt <- oda_balance_table(group, X, mcarlo = TRUE, mc_iter = 500L)
#' bt$rows[, c("attribute", "ess_display", "p_mc", "significant_raw")]
#' @export
oda_balance_table <- function(group,
                               X,
                               w               = NULL,
                               covariate_types = NULL,
                               loo             = "off",
                               mcarlo          = TRUE,
                               mc_iter         = 1000L,
                               alpha           = 0.05,
                               adjust          = c("none", "sidak", "bonferroni"),
                               ...) {
  adjust  <- match.arg(adjust)
  group   <- as.integer(group)
  X       <- as.data.frame(X)
  n       <- nrow(X)
  n_cov   <- ncol(X)
  col_nms <- names(X)

  if (length(group) != n)
    stop("group must have the same length as nrow(X).", call. = FALSE)

  grp_vals <- sort(unique(group[!is.na(group)]))
  if (length(grp_vals) != 2L)
    stop("group must be a binary variable with exactly 2 distinct non-missing values.",
         call. = FALSE)

  if (!is.null(w)) .validate_case_weights(w, n)

  has_wts <- !is.null(w) && any(w != 1, na.rm = TRUE)

  # ---- Per-covariate ODA fits -----------------------------------------------
  row_list <- vector("list", n_cov)

  for (j in seq_len(n_cov)) {
    col_nm <- col_nms[j]

    at <- if (!is.null(covariate_types) && col_nm %in% names(covariate_types))
      covariate_types[[col_nm]] else "auto"

    fit <- tryCatch(
      oda_fit(x         = X[[j]],
              y         = group,
              w         = w,
              attr_type = at,
              mcarlo    = isTRUE(mcarlo),
              mc_iter   = as.integer(mc_iter),
              loo       = loo,
              ...),
      error = function(e) list(ok = FALSE, reason = conditionMessage(e))
    )

    ok <- isTRUE(fit$ok)

    row_list[[j]] <- list(
      attribute   = col_nm,
      attr_type   = if (ok) fit$attr_type %||% NA_character_ else NA_character_,
      n_total     = if (ok) as.integer(fit$n_eff %||% NA_integer_) else NA_integer_,
      n_group_0   = if (ok && !is.null(fit$confusion$n0))
                      as.integer(fit$confusion$n0) else NA_integer_,
      n_group_1   = if (ok && !is.null(fit$confusion$n1))
                      as.integer(fit$confusion$n1) else NA_integer_,
      sensitivity = if (ok) fit$confusion$sensitivity %||% NA_real_ else NA_real_,
      specificity = if (ok) fit$confusion$specificity %||% NA_real_ else NA_real_,
      mean_pac    = if (ok) fit$confusion$mean_pac    %||% NA_real_ else NA_real_,
      ess         = if (ok && !isTRUE(fit$has_weights)) fit$ess %||% NA_real_ else NA_real_,
      wess        = if (ok && isTRUE(fit$has_weights))  fit$ess %||% NA_real_ else NA_real_,
      ess_display = if (ok) fit$ess %||% NA_real_ else NA_real_,
      p_mc        = if (ok) fit$p_mc %||% NA_real_ else NA_real_,
      p_sidak         = NA_real_,   # filled after loop
      p_bonferroni    = NA_real_,   # filled after loop
      significant_raw         = FALSE,
      significant_sidak       = FALSE,
      significant_bonferroni  = FALSE,
      significant             = FALSE,
      rule_type    = if (ok && !is.null(fit$rule))
                       fit$rule$type %||% NA_character_ else NA_character_,
      rule_summary = .oda_balance_rule_summary(fit, attr_nm = col_nm),
      loo_status   = if (ok) fit$loo$status  %||% NA_character_ else NA_character_,
      ess_loo      = if (ok) fit$loo$ess_loo %||% NA_real_      else NA_real_,
      has_weights  = has_wts,
      fit_ok       = ok,
      fit_reason   = if (!ok) as.character(fit$reason %||% NA_character_) else NA_character_
    )
  }

  tbl <- do.call(rbind, lapply(row_list, as.data.frame, stringsAsFactors = FALSE))

  # ---- Multiplicity corrections -----------------------------------------------
  # k = number of covariates with a valid (non-NA) p_mc
  valid_p <- !is.na(tbl$p_mc)
  k       <- sum(valid_p)

  if (k > 0L) {
    tbl$p_sidak[valid_p] <-
      vapply(tbl$p_mc[valid_p], .sidak_p,      double(1L), k = k)
    tbl$p_bonferroni[valid_p] <-
      vapply(tbl$p_mc[valid_p], .bonferroni_p, double(1L), k = k)
  }

  tbl$significant_raw        <- !is.na(tbl$p_mc)         & tbl$p_mc         < alpha
  tbl$significant_sidak      <- !is.na(tbl$p_sidak)      & tbl$p_sidak      < alpha
  tbl$significant_bonferroni <- !is.na(tbl$p_bonferroni) & tbl$p_bonferroni < alpha
  tbl$significant <- switch(adjust,
    "none"        = tbl$significant_raw,
    "sidak"       = tbl$significant_sidak,
    "bonferroni"  = tbl$significant_bonferroni
  )

  out <- list(
    rows = tbl,
    meta = list(
      n_covariates = n_cov,
      n_obs        = n,
      has_weights  = has_wts,
      ess_label    = if (has_wts) "WESS" else "ESS",
      alpha        = alpha,
      adjust       = adjust,
      k_valid      = k,
      loo_mode     = as.character(loo),
      mcarlo       = isTRUE(mcarlo),
      mc_iter      = as.integer(mc_iter)
    )
  )
  class(out) <- "oda_balance_table"
  out
}

# ---- smd_balance_table ------------------------------------------------------ #

#' Conventional SMD companion table for covariate balance
#'
#' Computes standardized mean differences (SMD) between two groups for each
#' covariate in \code{X}.  Returns one row per covariate with group means,
#' standard deviations, raw and absolute SMD, and conventional balance
#' thresholds.
#'
#' \strong{SMD is a conventional companion diagnostic, not the odacore
#' balance objective.}  The primary odacore balance assessment uses
#' \code{\link{oda_balance_table}}.  This function is intended for comparison
#' with non-ODA balance reports.
#'
#' No p-values are computed.  SMD is a descriptive statistic.  For a variable
#' with zero within-group variance in both groups, \code{smd} is \code{NA}.
#'
#' @param group Integer (or coercible) binary group indicator.  Must have
#'   exactly two distinct non-missing values.
#' @param X Data frame of baseline covariate columns.
#' @param w Optional numeric case-weight vector.  When supplied, weighted
#'   group means (\code{wmean_0}, \code{wmean_1}) and weighted SMD
#'   (\code{wsmd}) are added.  Weighted SMD uses the unweighted pooled SD as
#'   the standardizer (Rubin simplification).
#' @return A \code{data.frame} of class \code{c("smd_balance_table",
#'   "data.frame")} with one row per covariate and columns:
#'   \code{attribute}, \code{n_group_0}, \code{n_group_1},
#'   \code{mean_0}, \code{sd_0}, \code{mean_1}, \code{sd_1},
#'   \code{smd}, \code{abs_smd},
#'   \code{balanced_020} (\code{abs_smd < 0.20}),
#'   \code{balanced_010} (\code{abs_smd < 0.10}),
#'   \code{wmean_0}, \code{wmean_1}, \code{wsmd}, \code{wabs_smd},
#'   \code{wbalanced_020}, \code{wbalanced_010}
#'   (weighted variants; \code{NA} when \code{w = NULL}).
#' @seealso \code{\link{oda_balance_table}}, \code{\link{oda_balance_plot_data}}
#' @examples
#' group <- c(rep(0L, 30), rep(1L, 30))
#' X     <- data.frame(age = c(rep(45, 30), rep(55, 30)),
#'                     score = rnorm(60, 50, 10))
#' smd_balance_table(group, X)
#' @export
smd_balance_table <- function(group, X, w = NULL) {
  group <- as.integer(group)
  X     <- as.data.frame(X)
  n     <- nrow(X)

  if (length(group) != n)
    stop("group must have the same length as nrow(X).", call. = FALSE)

  grp_vals <- sort(unique(group[!is.na(group)]))
  if (length(grp_vals) != 2L)
    stop("group must be a binary variable with exactly 2 distinct non-missing values.",
         call. = FALSE)

  if (!is.null(w)) .validate_case_weights(w, n)

  g0   <- grp_vals[1L]; g1 <- grp_vals[2L]
  idx0 <- which(!is.na(group) & group == g0)
  idx1 <- which(!is.na(group) & group == g1)

  col_nms  <- names(X)
  row_list <- vector("list", ncol(X))

  for (j in seq_along(col_nms)) {
    x0 <- suppressWarnings(as.numeric(X[[j]][idx0]))
    x1 <- suppressWarnings(as.numeric(X[[j]][idx1]))

    n0_obs <- sum(!is.na(x0)); n1_obs <- sum(!is.na(x1))
    m0 <- if (n0_obs >= 1L) mean(x0, na.rm = TRUE) else NA_real_
    m1 <- if (n1_obs >= 1L) mean(x1, na.rm = TRUE) else NA_real_
    s0 <- if (n0_obs >= 2L) sd(x0,   na.rm = TRUE) else NA_real_
    s1 <- if (n1_obs >= 2L) sd(x1,   na.rm = TRUE) else NA_real_

    # Pooled SD (Hedges/Cohen formula)
    pooled_sd <- NA_real_
    if (!is.na(s0) && !is.na(s1) && n0_obs >= 2L && n1_obs >= 2L) {
      pvar     <- ((n0_obs - 1L) * s0^2 + (n1_obs - 1L) * s1^2) /
                  (n0_obs + n1_obs - 2L)
      pooled_sd <- sqrt(pvar)
    }

    smd_val <- if (!is.na(m0) && !is.na(m1) && !is.na(pooled_sd) && pooled_sd > 0)
      (m1 - m0) / pooled_sd else NA_real_
    abs_smd_val <- if (!is.na(smd_val)) abs(smd_val) else NA_real_

    # Weighted means and SMD
    wm0 <- NA_real_; wm1 <- NA_real_
    wsmd_val <- NA_real_; wabs_smd_val <- NA_real_

    if (!is.null(w)) {
      w0 <- w[idx0]; w1 <- w[idx1]
      k0w <- !is.na(x0) & !is.na(w0) & w0 > 0
      k1w <- !is.na(x1) & !is.na(w1) & w1 > 0
      if (sum(k0w) >= 1L)
        wm0 <- sum(w0[k0w] * x0[k0w]) / sum(w0[k0w])
      if (sum(k1w) >= 1L)
        wm1 <- sum(w1[k1w] * x1[k1w]) / sum(w1[k1w])
      # Weighted SMD uses unweighted pooled SD as standardizer (Rubin simplification)
      if (!is.na(wm0) && !is.na(wm1) && !is.na(pooled_sd) && pooled_sd > 0) {
        wsmd_val     <- (wm1 - wm0) / pooled_sd
        wabs_smd_val <- abs(wsmd_val)
      }
    }

    row_list[[j]] <- data.frame(
      attribute     = col_nms[j],
      n_group_0     = n0_obs,
      n_group_1     = n1_obs,
      mean_0        = m0,
      sd_0          = s0,
      mean_1        = m1,
      sd_1          = s1,
      smd           = smd_val,
      abs_smd       = abs_smd_val,
      balanced_020  = if (!is.na(abs_smd_val)) abs_smd_val < 0.20 else NA,
      balanced_010  = if (!is.na(abs_smd_val)) abs_smd_val < 0.10 else NA,
      wmean_0       = wm0,
      wmean_1       = wm1,
      wsmd          = wsmd_val,
      wabs_smd      = wabs_smd_val,
      wbalanced_020 = if (!is.na(wabs_smd_val)) wabs_smd_val < 0.20 else NA,
      wbalanced_010 = if (!is.na(wabs_smd_val)) wabs_smd_val < 0.10 else NA,
      stringsAsFactors = FALSE
    )
  }

  tbl <- do.call(rbind, row_list)
  class(tbl) <- c("smd_balance_table", "data.frame")
  tbl
}

# ---- oda_balance_plot_data -------------------------------------------------- #

#' Renderer-ready plot data for univariate ODA covariate balance
#'
#' Transforms an \code{\link{oda_balance_table}} result (and optionally an
#' \code{\link{smd_balance_table}} result) into a renderer-independent data
#' structure suitable for Graphics v3 plotting.
#'
#' \strong{This function does not fit any ODA models and does not accept
#' \code{group} or \code{X} arguments.}  It is a pure transformation of
#' pre-computed balance tables.
#'
#' @param balance_table An \code{"oda_balance_table"} object from
#'   \code{\link{oda_balance_table}}.
#' @param smd_table Optional \code{"smd_balance_table"} from
#'   \code{\link{smd_balance_table}}.  When supplied, SMD columns are joined
#'   by \code{attribute} name and included in \code{rows}.
#' @param p_col Character; which p-value column to use for the \code{p_plot}
#'   and \code{significant} columns in the output rows.
#'   One of \code{"p_mc"} (default), \code{"p_sidak"},
#'   \code{"p_bonferroni"}.
#' @param rank_by Character; how to rank covariates for display order.
#'   \code{"abs_ess"} (default): descending ESS/WESS (most imbalanced first).
#'   \code{"p"}: ascending p (most significant first).
#'   \code{"abs_smd"}: descending absolute SMD (requires \code{smd_table}).
#' @return A list of class \code{"oda_balance_plot_data"} with elements:
#' \describe{
#'   \item{\code{rows}}{Data frame; one row per covariate, sorted by
#'     \code{rank_by}.  Columns: \code{attribute}, \code{attr_type},
#'     \code{ess_display}, \code{ess_display_bar} (clipped to [0, 100]),
#'     \code{p_plot} (selected p column), \code{significant},
#'     \code{significance_label} (\code{"*"} or \code{""}),
#'     \code{rule_summary}, \code{abs_smd}, \code{wsmd_available},
#'     \code{abs_smd_display} (weighted if active), \code{fit_ok},
#'     \code{rank}.}
#'   \item{\code{has_weights}}{Logical.}
#'   \item{\code{ess_label}}{Character; \code{"WESS"} or \code{"ESS"}.}
#'   \item{\code{p_col_used}}{Character; selected p column name.}
#'   \item{\code{alpha}}{Numeric; significance threshold from metadata.}
#'   \item{\code{n_covariates}}{Integer.}
#'   \item{\code{n_significant}}{Integer; covariates significant on
#'     \code{p_col_used}.}
#'   \item{\code{rank_by}}{Character.}
#' }
#' @seealso \code{\link{oda_balance_table}}, \code{\link{smd_balance_table}}
#' @examples
#' set.seed(1)
#' group <- c(rep(0L, 30), rep(1L, 30))
#' X     <- data.frame(age   = c(rnorm(30, 45, 8), rnorm(30, 55, 8)),
#'                     score = rnorm(60, 50, 10))
#' bt  <- oda_balance_table(group, X, mcarlo = TRUE, mc_iter = 500L)
#' smd <- smd_balance_table(group, X)
#' pd  <- oda_balance_plot_data(bt, smd_table = smd)
#' pd$rows[, c("attribute", "ess_display", "p_plot", "significant", "abs_smd")]
#' @export
oda_balance_plot_data <- function(balance_table,
                                   smd_table = NULL,
                                   p_col     = c("p_mc", "p_sidak", "p_bonferroni"),
                                   rank_by   = c("abs_ess", "p", "abs_smd")) {
  stopifnot(inherits(balance_table, "oda_balance_table"))
  p_col   <- match.arg(p_col)
  rank_by <- match.arg(rank_by)

  tbl  <- balance_table$rows
  meta <- balance_table$meta
  nl   <- nrow(tbl)

  if (!p_col %in% names(tbl))
    stop("p_col '", p_col, "' not found in balance_table$rows.", call. = FALSE)

  # ---- SMD join -------------------------------------------------------------
  abs_smd_raw    <- rep(NA_real_, nl)
  abs_smd_disp   <- rep(NA_real_, nl)
  wsmd_avail     <- rep(FALSE,    nl)

  if (!is.null(smd_table)) {
    if (!inherits(smd_table, c("smd_balance_table", "data.frame")))
      stop("smd_table must be a smd_balance_table or data.frame.", call. = FALSE)
    smd_idx <- match(tbl$attribute, smd_table$attribute)
    valid   <- !is.na(smd_idx)
    if (any(valid)) {
      if ("abs_smd" %in% names(smd_table))
        abs_smd_raw[valid] <- smd_table$abs_smd[smd_idx[valid]]
      has_wsmd <- "wabs_smd" %in% names(smd_table) &&
                  any(!is.na(smd_table$wabs_smd))
      wsmd_avail[valid] <- has_wsmd
      if (isTRUE(meta$has_weights) && has_wsmd) {
        abs_smd_disp[valid] <- smd_table$wabs_smd[smd_idx[valid]]
      } else {
        abs_smd_disp[valid] <- abs_smd_raw[valid]
      }
    }
  }

  # ---- Significance ---------------------------------------------------------
  p_plot     <- tbl[[p_col]]
  sig_col    <- switch(p_col,
    p_mc         = "significant_raw",
    p_sidak      = "significant_sidak",
    p_bonferroni = "significant_bonferroni"
  )
  significant <- if (sig_col %in% names(tbl)) tbl[[sig_col]]
                 else !is.na(p_plot) & p_plot < (meta$alpha %||% 0.05)

  sig_lbl <- ifelse(!is.na(significant) & significant, "*", "")

  # ---- Rows data frame ------------------------------------------------------
  ess_disp <- tbl$ess_display
  rows_df  <- data.frame(
    attribute          = tbl$attribute,
    attr_type          = tbl$attr_type,
    ess_display        = ess_disp,
    ess_display_bar    = pmax(0, pmin(100, ifelse(is.na(ess_disp), NA_real_, ess_disp))),
    p_plot             = p_plot,
    significant        = significant,
    significance_label = sig_lbl,
    rule_summary       = tbl$rule_summary,
    abs_smd            = abs_smd_raw,
    wsmd_available     = wsmd_avail,
    abs_smd_display    = abs_smd_disp,
    fit_ok             = tbl$fit_ok,
    stringsAsFactors   = FALSE
  )

  # ---- Ranking --------------------------------------------------------------
  rank_key <- switch(rank_by,
    abs_ess = { v <- rows_df$ess_display; -ifelse(is.na(v), -Inf, v) },
    p       = { v <- rows_df$p_plot;      ifelse(is.na(v),  Inf,  v) },
    abs_smd = { v <- rows_df$abs_smd_display; -ifelse(is.na(v), -Inf, v) }
  )
  ord     <- order(rank_key, na.last = TRUE)
  rows_df <- rows_df[ord, ]
  rows_df$rank <- seq_len(nl)
  rownames(rows_df) <- NULL

  list(
    rows          = rows_df,
    has_weights   = isTRUE(meta$has_weights),
    ess_label     = meta$ess_label %||% "ESS",
    p_col_used    = p_col,
    alpha         = meta$alpha %||% 0.05,
    n_covariates  = nl,
    n_significant = sum(!is.na(significant) & significant),
    rank_by       = rank_by
  ) |> structure(class = "oda_balance_plot_data")
}

# ---- cta_balance_table ------------------------------------------------------ #

#' Multivariate CTA covariate balance diagnostics
#'
#' Fits a single \code{\link{cta_fit}} model with \code{group} as the class
#' variable and all columns of \code{X} as candidate predictors.  Returns a
#' structured summary of the CTA balance result.
#'
#' A \code{status = "no_tree"} result means no combination of baseline
#' covariates in \code{X} predicted group membership at the declared
#' significance level, LOO constraint, and minimum endpoint denominator.
#' This is \strong{favorable evidence of multivariable covariate balance}
#' under the declared analytic constraints.  It must not be interpreted as
#' a model failure; in balance analysis, inability to discriminate groups is
#' the goal.
#'
#' \strong{group vs. outcome:} \code{group} is the binary class variable.
#' The scientific outcome is strictly out of scope.
#'
#' \strong{Implementation constraint:} this function calls \code{\link{cta_fit}}
#' once; it does not reimplement ENUMERATE or node-growth logic.
#'
#' @param group Integer (or coercible) binary group indicator.  Must have
#'   exactly two distinct non-missing values.
#' @param X Data frame of baseline covariate columns.
#' @param w Optional numeric case-weight vector.  When supplied, CTA uses case
#'   weights and \code{has_weights = TRUE} in the result.
#' @param mindenom Integer minimum endpoint denominator passed to
#'   \code{\link{cta_fit}}.  Default \code{1L}.
#' @param alpha Numeric significance threshold stored in the result and used
#'   in the \code{no_tree_message} of \code{\link{cta_balance_plot_data}}.
#'   Default \code{0.05}.  Does not override \code{alpha_split}; pass
#'   \code{alpha_split} via \code{...} to change the CTA node-level threshold.
#' @param loo LOO gate mode passed to \code{\link{cta_fit}}.  Default
#'   \code{"off"}.
#' @param mc_iter Integer MC iterations per CTA node.  Default \code{5000L}.
#' @param mc_seed Integer RNG seed; \code{NULL} for unseeded.
#' @param ... Additional arguments forwarded to \code{\link{cta_fit}} (e.g.,
#'   \code{alpha_split}, \code{prune_alpha}, \code{priors_on}).
#' @return A list of class \code{"cta_balance_table"} with fields:
#' \describe{
#'   \item{\code{status}}{Character: \code{"valid_tree"}, \code{"stump"},
#'     \code{"no_tree"}, or \code{"fit_error"}.}
#'   \item{\code{balance_interpretation}}{Character: \code{"discriminating"} or
#'     \code{"no_discriminating_combinations"} (when \code{no_tree});
#'     \code{NA} on fit error.}
#'   \item{\code{root_attribute}}{Character; root split variable name;
#'     \code{NA} when \code{no_tree}.}
#'   \item{\code{n_endpoints}}{Integer; number of terminal endpoints;
#'     \code{NA} when \code{no_tree}.}
#'   \item{\code{overall_ess}}{Numeric; full-tree ESS (\%) when weights not
#'     active; \code{NA} otherwise.}
#'   \item{\code{overall_wess}}{Numeric; full-tree WESS (\%) when weights
#'     active; \code{NA} otherwise.}
#'   \item{\code{ess_display}}{Numeric; operative measure (\code{overall_wess}
#'     when weights active, else \code{overall_ess}); \code{NA} for no_tree.}
#'   \item{\code{d_stat}}{Numeric; parsimony-adjusted D statistic;
#'     \code{NA} for no_tree.}
#'   \item{\code{mindenom}}{Integer; MINDENOM used.}
#'   \item{\code{alpha}}{Numeric; significance threshold stored for downstream
#'     use.}
#'   \item{\code{has_weights}}{Logical; whether case weights were active.}
#'   \item{\code{tree}}{The raw \code{cta_tree} object; \code{NULL} on fit
#'     error.}
#'   \item{\code{endpoint_table}}{Data frame from
#'     \code{\link{cta_endpoint_table}}; zero-row for no_tree.}
#'   \item{\code{node_table}}{Data frame from \code{\link{cta_node_table}}.}
#'   \item{\code{fit_error}}{Logical; \code{TRUE} when \code{cta_fit} threw.}
#'   \item{\code{fit_reason}}{Character; error message when \code{fit_error};
#'     \code{NA} otherwise.}
#' }
#' @references
#' Linden A, Yarnold PR (2016). Using machine learning to assess covariate
#' balance in matching studies. \emph{Journal of Evaluation in Clinical
#' Practice}, \strong{22}(6), 861-867.
#' @seealso \code{\link{cta_balance_plot_data}}, \code{\link{oda_balance_table}},
#'   \code{\link{cta_fit}}
#' @examples
#' X <- data.frame(
#'   A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
#'   B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
#' )
#' group <- c(rep(0L, 40), rep(1L, 20))
#' ct <- cta_balance_table(group, X, mindenom = 5L,
#'                          mc_iter = 200L, mc_seed = 42L)
#' ct$status
#' ct$balance_interpretation
#' @export
cta_balance_table <- function(group,
                               X,
                               w        = NULL,
                               mindenom = 1L,
                               alpha    = 0.05,
                               loo      = "off",
                               mc_iter  = 5000L,
                               mc_seed  = NULL,
                               ...) {
  group <- as.integer(group)
  X     <- as.data.frame(X)
  n     <- nrow(X)

  if (length(group) != n)
    stop("group must have the same length as nrow(X).", call. = FALSE)

  grp_vals <- sort(unique(group[!is.na(group)]))
  if (length(grp_vals) != 2L)
    stop("group must be a binary variable with exactly 2 distinct non-missing values.",
         call. = FALSE)

  if (!is.null(w)) .validate_case_weights(w, n)

  has_wts <- !is.null(w) && any(w != 1, na.rm = TRUE)

  # ---- Fit CTA ---------------------------------------------------------------
  fit_args <- list(X        = X,
                   y        = group,
                   mindenom = as.integer(mindenom),
                   loo      = loo,
                   mc_iter  = as.integer(mc_iter))
  if (!is.null(w))       fit_args$w       <- w
  if (!is.null(mc_seed)) fit_args$mc_seed <- as.integer(mc_seed)
  fit_args <- modifyList(fit_args, list(...))

  tree <- tryCatch(
    do.call(cta_fit, fit_args),
    error = function(e) list(.cta_fit_error = TRUE, reason = conditionMessage(e))
  )

  # Handle fit failure
  if (isTRUE(tree$.cta_fit_error)) {
    empty_ep <- cta_endpoint_table(
      structure(list(no_tree = TRUE, nodes = list(), training_confusion = NULL,
                     has_weights = has_wts), class = "cta_tree"))
    out <- list(
      status                 = "fit_error",
      balance_interpretation = NA_character_,
      root_attribute         = NA_character_,
      n_endpoints            = NA_integer_,
      overall_ess            = NA_real_,
      overall_wess           = NA_real_,
      ess_display            = NA_real_,
      d_stat                 = NA_real_,
      mindenom               = as.integer(mindenom),
      alpha                  = alpha,
      has_weights            = has_wts,
      tree                   = NULL,
      endpoint_table         = empty_ep,
      node_table             = data.frame(),
      fit_error              = TRUE,
      fit_reason             = as.character(tree$reason)
    )
    class(out) <- "cta_balance_table"
    return(out)
  }

  # ---- Classify status -------------------------------------------------------
  no_tree  <- isTRUE(tree$no_tree)
  n_strata <- if (no_tree) NA_integer_ else cta_strata(tree)

  status <- if (no_tree) {
    "no_tree"
  } else if (!is.na(n_strata) && n_strata == 2L) {
    "stump"
  } else {
    "valid_tree"
  }

  balance_interp <- if (no_tree) "no_discriminating_combinations" else "discriminating"

  # ---- Root attribute (node_id == 1 has parent_id == 0L) ---------------------
  root_attr <- NA_character_
  if (!no_tree && !is.null(tree$nodes) && length(tree$nodes) >= 1L) {
    nd1 <- tree$nodes[[1L]]
    if (!is.null(nd1) && !isTRUE(nd1$leaf))
      root_attr <- nd1$attribute %||% NA_character_
  }

  # ---- ESS -------------------------------------------------------------------
  ess_val      <- if (!no_tree) tree$overall_ess %||% NA_real_ else NA_real_
  overall_ess  <- if (!has_wts) ess_val else NA_real_
  overall_wess <- if (has_wts)  ess_val else NA_real_

  # ---- Tables ----------------------------------------------------------------
  ep_tbl <- cta_endpoint_table(tree)
  nd_tbl <- cta_node_table(tree)

  out <- list(
    status                 = status,
    balance_interpretation = balance_interp,
    root_attribute         = root_attr,
    n_endpoints            = n_strata,
    overall_ess            = overall_ess,
    overall_wess           = overall_wess,
    ess_display            = ess_val,
    d_stat                 = cta_d_stat(tree),
    mindenom               = as.integer(mindenom),
    alpha                  = alpha,
    has_weights            = has_wts,
    tree                   = tree,
    endpoint_table         = ep_tbl,
    node_table             = nd_tbl,
    fit_error              = FALSE,
    fit_reason             = NA_character_
  )
  class(out) <- "cta_balance_table"
  out
}

# ---- cta_balance_plot_data -------------------------------------------------- #

#' Renderer-ready plot data for CTA covariate balance
#'
#' Transforms a \code{\link{cta_balance_table}} result into a
#' renderer-independent data structure suitable for Graphics v3 plotting.
#' For \code{no_tree} results, populates \code{no_tree_message} with the
#' favorable-balance interpretation.
#'
#' \strong{This function does not fit any CTA models.}  It is a pure
#' transformation of the pre-computed \code{cta_balance_table} result.
#'
#' @param cta_balance A \code{"cta_balance_table"} object from
#'   \code{\link{cta_balance_table}}.
#' @param target_class Integer; target class for endpoint coloring in the
#'   embedded tree diagram.  Default \code{1L} (group 1 = treated).
#' @param digits Integer; decimal digits passed to \code{\link{cta_plot_data}}.
#'   Default \code{1L}.
#' @return A list of class \code{"cta_balance_plot_data"} with elements:
#' \describe{
#'   \item{\code{status}}{Character; \code{"valid_tree"}, \code{"stump"},
#'     \code{"no_tree"}, or \code{"fit_error"}.}
#'   \item{\code{balance_interpretation}}{Character.}
#'   \item{\code{no_tree_message}}{Character; human-readable no-tree
#'     annotation for renderers; \code{NA} when status is not
#'     \code{"no_tree"}.}
#'   \item{\code{cta_pd}}{List from \code{\link{cta_plot_data}} when a valid
#'     tree or stump was found; \code{NULL} for no_tree or fit_error.}
#'   \item{\code{ess_display}}{Numeric; full-tree ESS/WESS (\%);
#'     \code{NA} for no_tree.}
#'   \item{\code{d_stat}}{Numeric; \code{NA} for no_tree.}
#'   \item{\code{has_weights}}{Logical.}
#'   \item{\code{ess_label}}{Character; \code{"WESS"} or \code{"ESS"}.}
#' }
#' @seealso \code{\link{cta_balance_table}}, \code{\link{cta_plot_data}}
#' @examples
#' X <- data.frame(
#'   A = c(rep(0L, 20), rep(1L, 20), rep(1L, 20)),
#'   B = c(rep(0L, 20), rep(0L, 20), rep(1L, 20))
#' )
#' group <- c(rep(0L, 40), rep(1L, 20))
#' ct  <- cta_balance_table(group, X, mindenom = 5L,
#'                           mc_iter = 200L, mc_seed = 42L)
#' cpd <- cta_balance_plot_data(ct)
#' cpd$status
#' @export
cta_balance_plot_data <- function(cta_balance, target_class = 1L, digits = 1L) {
  stopifnot(inherits(cta_balance, "cta_balance_table"))

  status  <- cta_balance$status
  no_tree <- status %in% c("no_tree", "fit_error")

  # no_tree_message ----------------------------------------------------------------
  no_tree_msg <- NA_character_
  if (identical(status, "no_tree")) {
    no_tree_msg <- sprintf(
      paste0("No combination of covariates predicted group membership\n",
             "(MINDENOM = %d, alpha = %g).\n",
             "This is favorable evidence of multivariable balance\n",
             "under the declared constraints."),
      cta_balance$mindenom,
      cta_balance$alpha %||% 0.05
    )
  } else if (identical(status, "fit_error")) {
    no_tree_msg <- paste0("CTA fitting error: ",
                          cta_balance$fit_reason %||% "unknown")
  }

  # cta_pd: populate only when a real tree exists --------------------------------
  cta_pd <- NULL
  if (!no_tree && !is.null(cta_balance$tree) &&
      inherits(cta_balance$tree, "cta_tree")) {
    cta_pd <- tryCatch(
      cta_plot_data(cta_balance$tree,
                    target_class = as.integer(target_class),
                    digits       = as.integer(digits)),
      error = function(e) NULL
    )
  }

  out <- list(
    status                 = status,
    balance_interpretation = cta_balance$balance_interpretation,
    no_tree_message        = no_tree_msg,
    cta_pd                 = cta_pd,
    ess_display            = cta_balance$ess_display,
    d_stat                 = cta_balance$d_stat,
    has_weights            = cta_balance$has_weights,
    ess_label              = if (isTRUE(cta_balance$has_weights)) "WESS" else "ESS"
  )
  class(out) <- "cta_balance_plot_data"
  out
}
