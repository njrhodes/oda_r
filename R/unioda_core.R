###############################################################################
# R/unioda_core.R
#
# UniODA Core - univariate binary-class ODA engine.
#
# Spec (MegaODA-faithful):
#   PRIMARY   = MAXSENS  (if priors_on)  / MEANSENS (if priors_off)
#   SECONDARY = SAMPLEREP
#   TERTIARY  = FIRST IDENTIFIED  (enumeration order)
#   LOO       = true refit-per-fold (no global rule reuse)
#   MC        = Fisher randomization with Clopper-Pearson STOP/STOPUP
#
# Exported entry points:
#   oda_univariate_core()   - fit one attribute
#   oda_loo_for_rule()      - LOO given a fixed rule (UniODA)
#   oda_mc_p_value()        - MC p-value (called internally; also exported)
#   oda_confusion_binary()  - confusion table
#   oda_rule_predict()      - apply a UniODA rule to new x
#   oda_mean_pac()          - (sens + spec) / 2
#   oda_ess_from_meanpac()  - ESS from mean PAC
###############################################################################

# ---- Basic performance metrics -------------------------------------------- #

#' Mean PAC from sensitivity and specificity
oda_mean_pac <- function(sens, spec) (sens + spec) / 2

#' Effect Strength for Sensitivity (ESS) from mean PAC
#' ESS(%) = 100 * (mean_pac - chance) / (1 - chance)
oda_ess_from_meanpac <- function(mean_pac, chance) {
  if (is.na(mean_pac) || is.na(chance) || chance >= 1) return(NA_real_)
  100 * (mean_pac - chance) / (1 - chance)
}

# ---- Missing-value handling ------------------------------------------------ #

#' Drop observations missing on x, y, or w for this attribute only.
oda_clean_xy <- function(x, y, w = NULL, miss_codes = NULL) {
  stopifnot(length(x) == length(y))
  n <- length(y)
  if (is.null(w)) w <- rep(1, n)
  stopifnot(length(w) == n)

  miss_mask <- is.na(x) | is.na(y) | is.na(w)
  if (!is.null(miss_codes)) {
    miss_mask <- miss_mask |
      (x %in% miss_codes) | (y %in% miss_codes) | (w %in% miss_codes)
  }

  keep <- !miss_mask
  list(x = x[keep], y = y[keep], w = as.numeric(w[keep]), idx = which(keep))
}

# ---- Attribute-type resolver ---------------------------------------------- #

#' Resolve attribute type from data if "auto".
oda_resolve_attr_type <- function(x,
    attr_type = c("auto","ordered","categorical","binary")) {
  attr_type <- match.arg(attr_type)
  if (attr_type != "auto") return(attr_type)

  if (is.factor(x) || is.character(x) || is.logical(x)) {
    k <- length(unique(x[!is.na(x)]))
    return(if (k == 2L) "binary" else "categorical")
  }
  ux <- unique(x[!is.na(x)])
  if (length(ux) == 2L) return("binary")
  "ordered"
}

# ---- Priors (binary class) ------------------------------------------------ #

#' Apply prior-odds weighting: weight class-0 and class-1 observations so each
#' class contributes equally to the objective (PRIORS ON).
oda_apply_priors <- function(y, w, priors_on = TRUE) {
  y <- as.integer(y)
  if (!priors_on) return(as.numeric(w))
  n0 <- sum(w[y == 0])
  n1 <- sum(w[y == 1])
  if (n0 == 0 || n1 == 0) return(as.numeric(w))   # pure node; caller handles
  w_adj        <- as.numeric(w)
  w_adj[y == 0] <- w_adj[y == 0] / n0
  w_adj[y == 1] <- w_adj[y == 1] / n1
  w_adj
}

# ---- Rule representation & prediction ------------------------------------- #

#' Assign each observation to the left (0) or right (1) side of a rule.
oda_rule_side <- function(x, rule) {
  if (rule$type == "ordered_cut") {
    as.integer(!(x <= rule$cut_value))    # x <= cut -> left(0), else right(1)
  } else if (rule$type %in% c("nominal_cut", "binary_map")) {
    as.integer(!(x %in% rule$left_levels))
  } else {
    stop("oda_rule_side: unknown rule$type '", rule$type, "'")
  }
}

#' Predict class label (0 or 1) for each observation given a UniODA rule.
oda_rule_predict <- function(x, rule) {
  side <- oda_rule_side(x, rule)          # 0 = left/Yes, 1 = right/No
  if (rule$direction == "0->1") {
    side                                   # left->0, right->1
  } else if (rule$direction == "1->0") {
    1L - side                              # left->1, right->0
  } else {
    stop("oda_rule_predict: unknown direction '", rule$direction, "'")
  }
}

# ---- Binary confusion table ----------------------------------------------- #

#' Compute weighted binary confusion table and sensitivity/specificity/meanPAC.
oda_confusion_binary <- function(y, y_pred, w = NULL) {
  y      <- as.integer(y)
  y_pred <- as.integer(y_pred)
  if (is.null(w)) w <- rep(1, length(y))
  w <- as.numeric(w)
  stopifnot(length(w) == length(y))

  TP <- sum(w[y_pred == 1L & y == 1L])
  TN <- sum(w[y_pred == 0L & y == 0L])
  FP <- sum(w[y_pred == 1L & y == 0L])
  FN <- sum(w[y_pred == 0L & y == 1L])

  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  mp   <- oda_mean_pac(sens, spec)

  list(TP = TP, TN = TN, FP = FP, FN = FN,
       sensitivity = sens, specificity = spec, mean_pac = mp)
}

# ---- Ordered-attribute block builder -------------------------------------- #

#' Collapse sorted x into unique-value blocks; return cumulative class counts.
oda_make_blocks_ordered <- function(x, y, w) {
  ord      <- order(x)
  x_ord    <- x[ord]; y_ord <- y[ord]; w_ord <- w[ord]
  block_id <- cumsum(c(TRUE, diff(x_ord) != 0))
  m        <- max(block_id)

  z <- numeric(m); o <- numeric(m); x_rep <- numeric(m); n_v <- integer(m)
  for (j in seq_len(m)) {
    idx   <- which(block_id == j)
    z[j]  <- sum(w_ord[idx][y_ord[idx] == 0])
    o[j]  <- sum(w_ord[idx][y_ord[idx] == 1])
    x_rep[j] <- x_ord[idx[1L]]
    n_v[j] <- length(idx)
  }
  list(ord = ord, z = z, o = o, x_rep = x_rep, n_v = n_v)
}

# ---- Primary / secondary tie-breaking ------------------------------------- #
#
# Spec (strictly book-faithful):
#   1. Keep all candidates at maximum ESS  (within tol_ess)
#   2. PRIMARY  heuristic further filters
#   3. SECONDARY = SAMPLEREP  (sample representativeness)
#   4. FIRST IDENTIFIED in enumeration order
#
# SAMPLEREP: maximise min(N_predicted_0, N_predicted_1)
#   === minimise |predicted_freq_0 - observed_freq_0|  for binary class.
#   We proxy via max( min(N0_pred, N1_pred) ) which is equivalent when total N
#   is fixed.  For direct L1 distance see multioda_core which has raw counts.
#' Apply primary/secondary tie-breaking to candidates
#' @keywords internal

oda_apply_primary_secondary <- function(cand_df, primary, secondary, y, w,
                                         preds_list) {
  w <- as.numeric(w)

  # 1. ESS filter
  if (!"ess" %in% names(cand_df))
    stop("cand_df must have column 'ess'.")
  max_ess <- max(cand_df$ess, na.rm = TRUE)
  tol_ess <- 1e-10
  keep    <- which(cand_df$ess >= max_ess - tol_ess)
  cand    <- cand_df[keep, , drop = FALSE]

  # 2. SAMPLEREP side counts (needed for secondary even when primary finishes)
  if (!"min_side_n" %in% names(cand)) {
    min_side <- numeric(nrow(cand))
    for (i in seq_len(nrow(cand))) {
      pred         <- preds_list[[ cand$row_id[i] ]]
      N0           <- sum(w[pred == 0L])
      N1           <- sum(w[pred == 1L])
      min_side[i]  <- min(N0, N1)
    }
    cand$min_side_n <- min_side
  }

  # 3. Primary filter
  primary <- match.arg(primary,
                       c("maxsens","maxspec","meansens","balanced","samplerep"))
  if (primary == "maxsens") {
    m    <- max(cand$sens_1,   na.rm = TRUE)
    cand <- cand[cand$sens_1   >= m - 1e-12, , drop = FALSE]
  } else if (primary == "maxspec") {
    m    <- max(cand$spec_0,   na.rm = TRUE)
    cand <- cand[cand$spec_0   >= m - 1e-12, , drop = FALSE]
  } else if (primary == "meansens") {
    m    <- max(cand$mean_pac, na.rm = TRUE)
    cand <- cand[cand$mean_pac >= m - 1e-12, , drop = FALSE]
  } else if (primary == "balanced") {
    m    <- min(cand$balance,  na.rm = TRUE)
    cand <- cand[cand$balance  <= m + 1e-12, , drop = FALSE]
  } else if (primary == "samplerep") {
    m    <- max(cand$min_side_n, na.rm = TRUE)
    cand <- cand[cand$min_side_n >= m - 1e-12, , drop = FALSE]
  }

  if (nrow(cand) <= 1L) return(cand[1L, , drop = FALSE])

  # 4. SAMPLEREP secondary (unless primary already was samplerep)
  secondary <- match.arg(secondary, c("samplerep","none"))
  if (secondary == "samplerep") {
    m    <- max(cand$min_side_n, na.rm = TRUE)
    cand <- cand[cand$min_side_n >= m - 1e-12, , drop = FALSE]
  }

  # 5. First identified (row_id is enumeration order)
  cand[which.min(cand$row_id), , drop = FALSE]
}

# ---- Monte Carlo p-value (Fisher randomization) --------------------------- #

# --------------------------------------------------------------------------- #
# Fast MC permutation ESS helpers (internal)                                  #
# --------------------------------------------------------------------------- #

# Precompute x/w structure for fast per-permutation ESS scoring.
# x_clean and w_clean must already have miss codes removed and be the same
# length as y passed to oda_mc_p_value() (which is pre-cleaned by the caller).
# attr_type: resolved string ("ordered", "binary", or "categorical").
# mindenom: MINDENOM for admissible cuts. Use 1L here to match the default
#   used when oda_univariate_core() is called from within oda_mc_p_value().
# Returns a precomp list or NULL when the fast path is not applicable.
.oda_mc_precomp <- function(x_clean, w_clean, attr_type, mindenom = 1L) {
  if (attr_type == "categorical") return(NULL)
  ord      <- order(x_clean)
  x_s      <- x_clean[ord]
  w_s      <- w_clean[ord]
  block_id <- as.integer(cumsum(c(TRUE, diff(x_s) != 0)))
  m        <- block_id[length(block_id)]
  if (m < 2L) return(NULL)
  tw_b  <- as.numeric(rowsum(w_s, block_id, reorder = FALSE))
  TW    <- cumsum(tw_b)
  n_b   <- as.integer(tabulate(block_id, nbins = m))
  N_v   <- cumsum(n_b)
  Nm    <- N_v[m]
  adm_j <- seq_len(m - 1L)
  adm_j <- adm_j[N_v[adm_j] >= mindenom & (Nm - N_v[adm_j]) >= mindenom]
  if (length(adm_j) == 0L) return(NULL)
  list(ord = ord, block_id = block_id, w_s = w_s,
       TW = TW, total_w = TW[m], m = m, adm_j = adm_j)
}

# Per-permutation max ESS for chance_model = "class" (chance = 0.5).
# y_star: permuted integer y (length == length(precomp$ord)).
# direction: "off" (non-directional), "0->1", or "1->0".
# Priors are irrelevant here: sens = class-1-right / total-class-1 and
# spec = class-0-left / total-class-0 are identical whether weights are
# priors-scaled or not (the class totals cancel in both numerator and
# denominator). Returns the max ESS scalar >= 0.
.oda_fast_ess_perm <- function(y_star, precomp, direction) {
  y_s <- y_star[precomp$ord]
  c1  <- as.numeric(rowsum(precomp$w_s * (y_s == 1L), precomp$block_id,
                           reorder = FALSE))
  O   <- cumsum(c1)
  W1  <- O[precomp$m]
  W0  <- precomp$total_w - W1
  if (W1 <= 0 || W0 <= 0) return(0)
  # ESS_01 at each admissible cut (direction "0->1"; ESS_10 = -ESS_01).
  ess_01 <- ((W1 - O[precomp$adm_j]) / W1 +
              (precomp$TW[precomp$adm_j] - O[precomp$adm_j]) / W0 - 1) * 100
  if (direction == "off") max(abs(ess_01))
  else if (direction == "0->1") max(ess_01)
  else max(-ess_01)   # "1->0"
}

#' Monte Carlo Fisher-randomization p-value with Clopper-Pearson early stopping.
#'
#' @param x,y,w  Data for the current attribute (already cleaned).
#' @param attr_type  "ordered", "categorical", or "binary".
#' @param priors_on  Logical.
#' @param primary,secondary  Tie-break heuristic strings.
#' @param miss_codes  Optional numeric vector of additional missing codes.
#' @param chance_model  "class" (1/2) or "attribute" (1/k_attr).
#' @param mc_iter  Maximum iterations.
#' @param mc_target  Significance threshold (e.g. 0.05).
#' @param mc_stop  Confidence level for lower-tail stop (e.g. 99.9).
#' @param mc_stopup  Confidence level for upper-tail stop (e.g. 20 -> 0.20). Default NA (disabled).
#' @param mc_adjust  Kept for API compatibility; not used.
#' @param seed  Optional RNG seed.
#' @param ess_obs  Observed ESS (must be supplied).
#' @param direction  Directional constraint forwarded from oda_univariate_core():
#'   "both" (canonical non-directional default), "off" (synonym for "both"),
#'   "greater", or "less". Each permutation refit uses the same constraint.
#' @param direction_map  Named integer vector for categorical fixed-partition
#'   DIRECTIONAL. When supplied, each permutation evaluates the SAME fixed
#'   mapping on permuted y labels. Default NULL.
#' @return List with p_mc, ge_count, iter_used, ess_obs.
oda_mc_p_value <- function(
    x, y,
    w            = NULL,
    attr_type,
    priors_on,
    primary,
    secondary,
    miss_codes   = NULL,
    chance_model = c("class","attribute"),
    mc_iter      = 25000L,
    mc_target    = 0.05,
    mc_stop      = 99.9,
    mc_stopup    = NA,
    mc_adjust    = FALSE,
    seed         = NULL,
    ess_obs      = NULL,
    direction    = c("both", "off", "greater", "less"),
    direction_map = NULL
) {
  chance_model <- match.arg(chance_model)
  direction    <- match.arg(direction)
  if (direction == "both") direction <- "off"  # canonical synonym
  y <- as.integer(y)
  n <- length(y)
  if (is.null(w)) w <- rep(1, n) else w <- as.numeric(w)
  if (!is.null(seed)) set.seed(seed)

  if (is.null(ess_obs) || !is.finite(ess_obs))
    stop("oda_mc_p_value: ess_obs must be supplied and finite.")

  # Convert stop confidence levels to probability scale.
  # STOP  (lower-tail): early accept when upper CP bound < mc_target at conf_stop  confidence.
  # STOPUP (upper-tail): early reject when lower CP bound > mc_target at conf_stopup confidence.
  conf_stop   <- if (!is.na(mc_stop)   && mc_stop   > 1) mc_stop   / 100 else mc_stop
  conf_stopup <- if (!is.na(mc_stopup) && mc_stopup > 1) mc_stopup / 100 else mc_stopup

  ge_count  <- 0L
  iter_used <- 0L
  min_check   <- 50L
  check_every <- 50L

  # Precompute for fast ordered/binary ESS scoring (chance_model = "class" only).
  # x is pre-cleaned by the caller (oda_univariate_core), so no miss-code filter needed.
  .fast_pc <- if (chance_model == "class") {
    .oda_mc_precomp(x, w, oda_resolve_attr_type(x, attr_type), mindenom = 1L)
  } else NULL

  for (b in seq_len(mc_iter)) {
    iter_used <- b
    y_star    <- sample(y, size = n, replace = FALSE)

    if (!is.null(.fast_pc)) {
      ess_b <- .oda_fast_ess_perm(y_star, .fast_pc, direction)
    } else {
      fit_b <- oda_univariate_core(
        x = x, y = y_star, w = w,
        attr_type     = attr_type,
        priors_on     = priors_on,
        primary       = primary,
        secondary     = secondary,
        miss_codes    = miss_codes,
        loo           = "off",
        mcarlo        = FALSE,
        mc_iter       = 0L,
        mc_target     = mc_target,
        mc_stop       = mc_stop,
        mc_stopup     = mc_stopup,
        mc_adjust     = mc_adjust,
        mc_seed       = NULL,
        chance_model  = chance_model,
        direction     = direction,
        direction_map = direction_map
      )
      ess_b <- if (!isTRUE(fit_b$ok)) 0 else fit_b$ess
    }

    if (ess_b >= ess_obs - 1e-12) ge_count <- ge_count + 1L

    # Clopper-Pearson early stopping
    if (b >= min_check && (b %% check_every == 0L)) {
      if (ge_count == 0L) {
        lower <- 0
        upper <- if (!is.na(conf_stop))   stats::qbeta(conf_stop,   1,           b)           else NA_real_
      } else if (ge_count == b) {
        upper <- 1
        lower <- if (!is.na(conf_stopup)) stats::qbeta(1 - conf_stopup, b,       1)           else NA_real_
      } else {
        upper <- if (!is.na(conf_stop))   stats::qbeta(conf_stop,   ge_count + 1, b - ge_count)      else NA_real_
        lower <- if (!is.na(conf_stopup)) stats::qbeta(1 - conf_stopup, ge_count, b - ge_count + 1)  else NA_real_
      }
      if (!is.na(mc_target) && !is.na(upper) && upper < mc_target) break   # STOP:   conf_stop   sure p < CUTOFF
      if (!is.na(mc_target) && !is.na(lower) && lower > mc_target) break   # STOPUP: conf_stopup sure p > CUTOFF
    }
  }

  list(p_mc      = if (iter_used <= 0L) NA_real_
                   else if (ge_count == 0L) 0.0
                   else ge_count / iter_used,
       ge_count  = ge_count,
       iter_used = iter_used,
       ess_obs   = ess_obs)
}

# ---- Algebraic ordered-cut LOO helper ------------------------------------ #

# Fast count-table LOO for ordered_cut rules.
#
# Instead of calling oda_univariate_core() n times (once per fold), this
# function builds a sorted unique-value count table and for each unique
# (x value, class) bin simulates removing one observation algebraically:
# decrement the bin count, rescan all cuts in ascending order using the same
# priors-adjusted ESS formula as training, and predict the held-out x.  All
# observations sharing the same (x value, class) bin receive the same LOO
# prediction, so the inner loop runs over at most 2k bins rather than n folds.
#
# Conservative fallback: returns NULL when raw weights are non-uniform within
# any bin (non-canonical case), leaving the caller to use full refits.
#
# Tie-breaking: strict > throughout preserves FIRST IDENTIFIED (ascending cut
# scan) as the tertiary rule, matching oda_univariate_core().
#
# Returns: integer vector y_pred_loo of length n (same 0L/1L encoding as
# oda_loo_for_rule), or NULL if the algebraic path is not applicable.
oda_loo_ordered_cut_counts <- function(x, y, w, priors_on, rule) {

  if (!identical(rule$type, "ordered_cut")) return(NULL)

  y <- as.integer(y)
  n <- length(y)

  # Extended to handle non-uniform case weights via per-bin class weight sums.
  # The ESS formula (cum0/tc0 + right1/tc1 - 1)*100 is equivalent to
  # (specificity + sensitivity - 1)*100 when cum0 and right1 are class weight
  # sums and tc0, tc1 are class weight totals.  For binary class, priors_on=TRUE
  # and FALSE produce the same sensitivity/specificity ratios (the priors scaling
  # cancels), so the same formula works for both settings.
  # For uniform weights, the per-obs loop gives the same result as the former
  # per-(bin, class) loop because every obs in the same (val, class) cell has
  # the same weight and thus the same adjusted table.

  obs  <- !is.na(x)

  # Class weight totals over observed rows.
  tc0_w <- sum(w[obs & y == 0L])
  tc1_w <- sum(w[obs & y == 1L])
  if (tc0_w <= 0 || tc1_w <= 0) return(NULL)   # pure node; no split possible

  vals <- sort(unique(x[obs]))
  k    <- length(vals)
  if (k < 2L) return(NULL)                      # no candidate cut exists

  # Precompute bin index for each observation.
  # bin_id[i] is NA when x[i] is NA (obs[i]=FALSE); those rows are skipped below.
  bin_id <- match(x, vals)

  # Build per-bin class weight sums: z0v[j], z1v[j] = sum of w at vals[j]
  # for class 0 and class 1 respectively.
  z0v <- numeric(k)
  z1v <- numeric(k)
  for (vi in seq_len(k)) {
    m      <- obs & (bin_id == vi)
    z0v[vi] <- sum(w[m & y == 0L])
    z1v[vi] <- sum(w[m & y == 1L])
  }

  # Scan helper: find best cut/direction on a (possibly adjusted) bin table.
  # ESS formula (works for both integer counts and float weight sums):
  #   ESS("0->1") = (cum0/tc0a + right1/tc1a - 1) * 100
  #   ESS("1->0") = (right0/tc0a + cum1/tc1a  - 1) * 100
  # Strict > keeps first-identified cut/direction when ESS ties.
  # SAMPLEREP secondary is not implemented here; first-identified is the
  # tiebreaker.  This matches the uniform-weight algebraic path (A1-A3).
  .find_best <- function(c0a, c1a, tc0a, tc1a, vls) {
    ka     <- length(vls)
    best_e <- -Inf
    best_c <- NA_real_
    best_d <- NA_character_
    cum0   <- 0L
    cum1   <- 0L
    for (i in seq_len(ka - 1L)) {
      cum0 <- cum0 + c0a[i]
      cum1 <- cum1 + c1a[i]
      r0   <- tc0a - cum0
      r1   <- tc1a - cum1
      if ((cum0 + cum1) == 0L || (r0 + r1) == 0L) next
      cut  <- (vls[i] + vls[i + 1L]) / 2.0
      e01  <- (cum0 / tc0a + r1   / tc1a - 1.0) * 100.0  # "0->1"
      e10  <- (r0   / tc0a + cum1 / tc1a - 1.0) * 100.0  # "1->0"
      if (e01 > best_e) { best_e <- e01; best_c <- cut; best_d <- "0->1" }
      if (e10 > best_e) { best_e <- e10; best_c <- cut; best_d <- "1->0" }
    }
    if (is.na(best_c)) NULL else list(cut = best_c, dir = best_d)
  }

  eps        <- 1e-12
  y_pred_loo <- integer(n)

  # Per-observation LOO: subtract each obs's weight from its bin/class cell and
  # refit the best cut on the adjusted n-1 weight-sum table.
  # For uniform weights every obs in the same (val, class) cell has the same wi,
  # producing the same adjusted table -- identical to the former per-(bin, class)
  # loop behavior that proved A1-A3.
  for (i in seq_len(n)) {
    if (!obs[i]) { y_pred_loo[i] <- NA_integer_; next }

    vi  <- bin_id[i]; wi <- w[i]; ci <- y[i]; val <- x[i]

    z0a <- z0v; z1a <- z1v
    if (ci == 0L) z0a[vi] <- z0a[vi] - wi else z1a[vi] <- z1a[vi] - wi
    tc0a <- tc0_w - wi * (ci == 0L)
    tc1a <- tc1_w - wi * (ci == 1L)

    keep_v <- (z0a + z1a) > eps
    vls_a  <- vals[keep_v]
    z0a_k  <- z0a[keep_v]
    z1a_k  <- z1a[keep_v]

    # Degenerate fold: one class eliminated or fewer than 2 values remain.
    # Fall back to the training rule for this observation.
    if (tc0a <= eps || tc1a <= eps || length(vls_a) < 2L) {
      y_pred_loo[i] <- oda_rule_predict(val, rule)
    } else {
      best <- .find_best(z0a_k, z1a_k, tc0a, tc1a, vls_a)
      if (is.null(best)) {
        y_pred_loo[i] <- oda_rule_predict(val, rule)
      } else if (identical(best$dir, "0->1")) {
        y_pred_loo[i] <- if (val <= best$cut) 0L else 1L
      } else {
        y_pred_loo[i] <- if (val <= best$cut) 1L else 0L
      }
    }
  }

  # NA-x observations: match oda_rule_predict(NA, rule) = NA_integer_.
  if (any(!obs)) y_pred_loo[!obs] <- NA_integer_

  y_pred_loo
}

# ---- Algebraic delete-one LOO for binary_map rules ------------------------ #
#
# MPE Chapter 2 (pp. 32-33) canonical LOO contract:
#   For each observation i: hold out i, obtain the ODA model from all
#   observations except i, classify i with that model, store the result.
#
# This function implements that contract for binary_map rules using a fast
# algebraic approach: build a weighted 2x2 table (x-side x class), subtract
# the held-out observation's weight from its cell, and recompute the optimal
# binary direction from the adjusted table.  Results match explicit delete-one
# oda_univariate_core(mcarlo=FALSE) refits for all tested configurations,
# including non-uniform case weights.
#
# Tie-breaking when ESS("0->1") == ESS("1->0") within tol=1e-10:
#   1. Primary MAXSENS: higher class-1 sensitivity wins (c11a/tc1a).
#   2. SAMPLEREP is symmetric for binary (both directions predict the same
#      side-weight totals); SAMPLEREP always ties here.
#   3. First-identified: "0->1" (matches oda_univariate_core enumeration order).
#
# Weighted binary LOO is canonical.  Weighted categorical LOO is out of scope
# and is separately guarded elsewhere.
#
# Returns integer vector y_pred_loo of length n, or NULL if not applicable
# (rule type mismatch or pure-node training data).
oda_loo_binary_map_counts <- function(x, y, w, priors_on, rule) {

  if (!identical(rule$type, "binary_map")) return(NULL)

  y <- as.integer(y)
  n <- length(y)
  if (is.null(w)) w <- rep(1.0, n) else w <- as.numeric(w)

  obs    <- !is.na(x)
  x_chr  <- as.character(x)
  left   <- obs & (x_chr %in% as.character(rule$left_levels))

  tc0_w  <- sum(w[obs & y == 0L])
  tc1_w  <- sum(w[obs & y == 1L])
  if (tc0_w <= 0 || tc1_w <= 0) return(NULL)   # pure node; no split possible

  c00_w  <- sum(w[left         & y == 0L])   # left side, class 0
  c01_w  <- sum(w[left         & y == 1L])   # left side, class 1
  c10_w  <- sum(w[obs & !left  & y == 0L])   # right side, class 0
  c11_w  <- sum(w[obs & !left  & y == 1L])   # right side, class 1

  # Choose direction from adjusted 2x2 table.
  # ESS("0->1") = c00a/tc0a + c11a/tc1a - 1
  # ESS("1->0") = c10a/tc0a + c01a/tc1a - 1
  .best_dir <- function(c00a, c10a, tc0a, c01a, c11a, tc1a) {
    if (tc0a <= 0 || tc1a <= 0) return(rule$direction)   # degenerate fold
    e01  <- c00a / tc0a + c11a / tc1a
    e10  <- c10a / tc0a + c01a / tc1a
    tol  <- 1e-10
    if (e01 > e10 + tol)                    return("0->1")
    if (e10 > e01 + tol)                    return("1->0")
    # ESS tie: primary MAXSENS (class-1 sensitivity)
    s01  <- c11a / tc1a
    s10  <- c01a / tc1a
    tol2 <- 1e-12
    if (s01 > s10 + tol2)                   return("0->1")
    if (s10 > s01 + tol2)                   return("1->0")
    # SAMPLEREP symmetric for binary; first-identified wins
    "0->1"
  }

  y_pred_loo <- integer(n)

  for (i in seq_len(n)) {
    if (!obs[i]) { y_pred_loo[i] <- NA_integer_; next }

    wi  <- w[i]; yi <- y[i]; li <- left[i]

    c00a <- if ( li && yi == 0L) c00_w - wi else c00_w
    c01a <- if ( li && yi == 1L) c01_w - wi else c01_w
    c10a <- if (!li && yi == 0L) c10_w - wi else c10_w
    c11a <- if (!li && yi == 1L) c11_w - wi else c11_w
    tc0a <- if (yi == 0L) tc0_w - wi else tc0_w
    tc1a <- if (yi == 1L) tc1_w - wi else tc1_w

    dir              <- .best_dir(c00a, c10a, tc0a, c01a, c11a, tc1a)
    y_pred_loo[i]    <- if (dir == "0->1") {
      if (li) 0L else 1L
    } else {
      if (li) 1L else 0L
    }
  }

  if (any(!obs)) y_pred_loo[!obs] <- NA_integer_

  y_pred_loo
}

# ---- LOO for a fixed UniODA rule ------------------------------------------ #

#' True leave-one-out cross-validation for a UniODA rule.
#' Each fold refits the full ODA model on n-1 observations then predicts the
#' held-out case.  LOO p-value uses Fisher's exact test.
#'
#' @param x,y  Predictors and class labels.
#' @param rule  A UniODA rule list (from oda_univariate_core()$rule).
#' @param w  Optional case weights.
#' @param chance_model  "class" or "attribute".
#' @param k_attr  Number of attribute levels (for ESS baseline).
#' @param attr_type  Attribute type string.
#' @param allow_weighted_categorical_loo  Logical (default FALSE per spec).
#' @param ...  Additional args forwarded to oda_univariate_core for each fold.
#' @return List with \code{allowed}, \code{confusion}, \code{ess_loo},
#'   \code{p_value}, and \code{alternative} (the \code{fisher.test} alternative
#'   used: \code{"two.sided"} when \code{direction = "off"}, \code{"greater"}
#'   when a directional hypothesis is declared).
oda_loo_for_rule <- function(
    x, y,
    rule,
    w            = NULL,
    chance_model = c("class","attribute"),
    k_attr       = NULL,
    attr_type    = c("ordered","categorical","binary"),
    allow_weighted_categorical_loo = FALSE,
    priors_on  = TRUE,
    primary    = "maxsens",
    secondary  = "samplerep",
    miss_codes = NULL,
    mc_iter    = 25000L,
    mc_target  = 0.05,
    mc_stop    = 99.9,
    mc_stopup  = NA_real_,
    mc_adjust  = FALSE,
    mc_seed    = NULL,
    direction  = c("both", "off", "greater", "less"),
    direction_map = NULL
) {
  chance_model <- match.arg(chance_model)
  attr_type    <- match.arg(attr_type)
  direction    <- match.arg(direction)
  if (direction == "both") direction <- "off"  # canonical synonym

  y <- as.integer(y)
  n <- length(y)
  if (is.null(w)) w <- rep(1, n) else w <- as.numeric(w)

  # Weighted categorical LOO is forbidden
  if (attr_type == "categorical" &&
      !allow_weighted_categorical_loo &&
      any(w != w[1L])) {
    return(list(allowed = FALSE,
                reason  = "weighted_categorical_loo_not_supported",
                confusion = NULL, ess_loo = NA_real_, p_value = NA_real_))
  }

  # Engineering guard: categorical LOO is not supported when every observed
  # attribute level is unique.  In that case every LOO fold's held-out level
  # is absent from training, so the absent-level fallback would drive all
  # predictions  -  not a meaningful cross-validation.
  # NOTE: this is NOT the MPE canon "superfluous LOO" case.  The MPE condition
  # (political affiliation, rater agreement) requires a declared directional
  # diagonal/agreement hypothesis (k_attr == C + DIRECTIONAL), which is a
  # separate deferred feature in issue #6.
  if (attr_type == "categorical") {
    x_obs <- x[!is.na(x)]
    if (!is.null(miss_codes)) x_obs <- x_obs[!(x_obs %in% miss_codes)]
    if (length(x_obs) > 0L && length(unique(x_obs)) == length(x_obs)) {
      return(list(allowed   = FALSE,
                  reason    = "loo_not_supported_all_unique_categories",
                  confusion = NULL,
                  ess_loo   = NA_real_,
                  p_value   = NA_real_))
    }
  }

  y_pred_loo <- integer(n)

  # Algebraic paths (fast, exact, no per-fold oda_univariate_core calls):
  #   binary_map: oda_loo_binary_map_counts() -- weighted 2x2 table adjustment.
  #   ordered_cut (uniform w): oda_loo_ordered_cut_counts() -- bin count adjustment.
  # Both return NULL when not applicable; fall through to explicit per-fold refits.
  #
  # All other rule types: full per-row refit via oda_univariate_core().

  .loo_done <- FALSE

  if (identical(rule$type, "binary_map")) {
    .alg <- oda_loo_binary_map_counts(x, y, w, priors_on, rule)
    if (!is.null(.alg)) {
      y_pred_loo <- .alg
      .loo_done  <- TRUE
    }
  }

  if (!.loo_done && identical(rule$type, "ordered_cut")) {
    .alg <- oda_loo_ordered_cut_counts(x, y, w, priors_on, rule)
    if (!is.null(.alg)) {
      y_pred_loo <- .alg
      .loo_done  <- TRUE
    }
  }

  if (!.loo_done) {
    for (i in seq_len(n)) {
      keep <- seq_len(n) != i

      # Fixed categorical rule (direction_map supplied): predictions are
      # determined entirely by x, not training data.  Apply training rule.
      if (identical(rule$type, "nominal_cut") && !is.null(direction_map)) {
        if (!is.na(x[i])) {
          y_pred_loo[i] <- oda_rule_predict(x[i], rule)
          next
        }
      }

      fit_i <- oda_univariate_core(
        x = x[keep], y = y[keep], w = w[keep],
        attr_type     = attr_type,
        priors_on     = priors_on,
        primary       = primary,
        secondary     = secondary,
        miss_codes    = miss_codes,
        loo           = "off",
        mcarlo        = FALSE,
        mc_iter       = mc_iter,
        mc_target     = mc_target,
        mc_stop       = mc_stop,
        mc_stopup     = mc_stopup,
        mc_adjust     = mc_adjust,
        mc_seed       = if (is.null(mc_seed)) NULL else mc_seed + i,
        chance_model  = chance_model,
        direction     = direction,
        direction_map = direction_map
      )

      if (!isTRUE(fit_i$ok)) {
        return(list(allowed = FALSE,
                    reason  = paste0("loo_fit_failed_case_", i,
                                     "_reason_", fit_i$reason),
                    confusion = NULL, ess_loo = NA_real_, p_value = NA_real_))
      }
      # Fixture-derived MegaODA compatibility: when the held-out observation's category was
      # absent from the n-1 LOO training fold (singleton or rare category
      # dropped out), MegaODA assigns it to the left-side class rather than
      # defaulting to the right side via !(x %in% left_levels).
      # Fixture evidence: MPE Chapter 4 Bray-Curtis Table 4.2, category E
      # (one S30 observation); absent-level override recovers LOO ESS = 43.5.
      # Applies to nominal_cut rules only; oda_rule_side() is unchanged.
      if (fit_i$rule$type == "nominal_cut" && !is.na(x[i])) {
        known  <- c(fit_i$rule$left_levels, fit_i$rule$right_levels)
        absent <- !(as.character(x[i]) %in% known)
      } else {
        absent <- FALSE
      }
      if (absent) {
        y_pred_loo[i] <- if (fit_i$rule$direction == "0->1") 0L else 1L
      } else {
        y_pred_loo[i] <- oda_rule_predict(x[i], fit_i$rule)
      }
    }
  }

  # ESS for LOO  -  use the same priors+case-weight scheme as training so that
  # ess_loo can be compared to ess_obj with all.equal().  When every fold uses
  # the training rule (admissible ordered/binary case), the weighted confusion
  # here is identical to the training confusion, giving ess_loo == ess_obj
  # exactly (within floating-point noise), matching MegaODA's WESSL = WESS.
  w_loo_adj <- oda_apply_priors(y, w, priors_on)
  conf_loo  <- oda_confusion_binary(y, y_pred_loo, w_loo_adj)

  chance_class <- 0.5
  if (is.null(k_attr)) k_attr <- max(2L, length(unique(x[!is.na(x)])))
  chance_attr  <- 1 / k_attr
  chance       <- if (chance_model == "class") chance_class else chance_attr
  ess_loo      <- oda_ess_from_meanpac(conf_loo$mean_pac, chance)

  # Fisher's exact test on raw (unit-weight) counts  -  requires integer inputs.
  conf_loo_r <- oda_confusion_binary(y, y_pred_loo)   # unit weights
  tab <- matrix(c(conf_loo_r$TP, conf_loo_r$FP,
                  conf_loo_r$FN, conf_loo_r$TN),
                nrow = 2, byrow = TRUE)
  fisher_alt <- if (direction == "off") "two.sided" else "greater"
  p_fisher <- tryCatch(
    stats::fisher.test(tab, alternative = fisher_alt)$p.value,
    error = function(e) NA_real_
  )

  list(allowed = TRUE, reason = NULL,
       confusion = conf_loo, ess_loo = ess_loo,
       p_value = p_fisher, alternative = fisher_alt)
}

# ---- Main UniODA entry point ---------------------------------------------- #

#' Fit a univariate binary-class ODA model.
#'
#' @param x  Attribute values (numeric, factor, character, or logical).
#' @param y  Binary class labels (0/1 integer-coercible).
#' @param w  Optional case weights (default: unit weights).
#' @param attr_type  One of "auto","ordered","categorical","binary".
#' @param priors_on  Logical; if TRUE, weight by inverse class frequency.
#' @param primary  Primary tie-break heuristic (NULL = default by priors).
#' @param secondary  Secondary tie-break heuristic (NULL = "samplerep").
#' @param miss_codes  Additional values to treat as missing.
#' @param loo  "off", "stable" (must match training ESS), or "pvalue".
#' @param loo_alpha  Alpha for LOO pvalue filter.
#' @param mcarlo  Logical; run Monte Carlo p-value.
#' @param mc_iter  Max MC iterations.
#' @param mc_target,mc_stop,mc_stopup  MC stopping parameters.
#' @param mc_adjust  Legacy parameter (unused).
#' @param mc_seed  RNG seed for MC.
#' @param chance_model  "class" (1/C) or "attribute" (1/k_attr).
#' @param direction Directional hypothesis (MPE Chapter 2 scope for ordered/binary):
#'   \code{"both"} (default, non-directional; evaluates both "0->1" and "1->0") or
#'   its backward-compatible synonym \code{"off"}, \code{"greater"} (Chapter 2
#'   greater-than direction: x > cut predicts class 1; MegaODA Appendix A
#'   \code{DIRECTION < 0 1}; internal "0->1"), or \code{"less"} (Chapter 2
#'   less-than direction: x <= cut predicts class 1; MegaODA Appendix A
#'   \code{DIRECTION > 0 1}; internal "1->0").
#'   Ordered and binary attributes only. Categorical returns ok=FALSE with
#'   reason "direction_not_supported_for_categorical".
#'   MPE Chapter 4 categorical/table DIRECTIONAL is Phase 6C (not yet implemented).
#' @return Named list with ok, rule, confusion, ess, pac, p_mc, loo, tie_block.
oda_univariate_core <- function(
    x,
    y,
    w            = NULL,
    attr_type    = c("auto","ordered","categorical","binary"),
    priors_on    = TRUE,
    primary      = NULL,
    secondary    = NULL,
    miss_codes   = NULL,
    missing_code = NULL,   # alias for miss_codes (scalar or vector)
    loo        = c("off","stable","pvalue"),
    loo_alpha  = 0.05,
    mcarlo     = TRUE,
    mc_iter    = 25000L,
    mc_target  = 0.05,
    mc_stop    = 99.9,
    mc_stopup  = NA_real_,
    mc_adjust  = FALSE,
    mc_seed    = NULL,
    chance_model = c("class","attribute"),
    eval_order   = c("mc_then_loo", "loo_then_mc"),
    mindenom     = 1L,
    direction    = c("both", "off", "greater", "less"),
    direction_map = NULL
) {
  chance_model <- match.arg(chance_model)
  loo          <- match.arg(loo)
  eval_order   <- match.arg(eval_order)
  mindenom     <- max(1L, as.integer(mindenom))
  direction    <- match.arg(direction)
  if (direction == "both") direction <- "off"  # canonical synonym

  # Resolve missing_code alias -> miss_codes
  if (!is.null(missing_code)) {
    miss_codes <- unique(c(miss_codes, as.numeric(missing_code)))
  }

  # 1. Clean
  clean <- oda_clean_xy(x, y, w, miss_codes)
  x     <- clean$x
  y     <- as.integer(clean$y)
  w_raw <- clean$w
  n_eff <- length(y)
  if (n_eff == 0L) return(list(ok = FALSE, reason = "all_missing", type = "leaf"))

  # 2. Priors
  w <- oda_apply_priors(y, w_raw, priors_on)

  # 3. Attribute type
  attr_type <- oda_resolve_attr_type(x, attr_type)

  # 4. Both classes present?
  n0 <- sum(w[y == 0L]); n1 <- sum(w[y == 1L])
  if (n0 == 0 || n1 == 0) return(list(ok = FALSE, reason = "pure_node", type = "leaf"))

  # 4a. Categorical DIRECTION not supported in Phase 6A.
  # MPE Chapter 4 TABLE/DIRECTIONAL semantics are deferred to Phase 6C.
  if (direction != "off" && attr_type == "categorical")
    return(list(ok = FALSE,
                reason = "direction_not_supported_for_categorical",
                type   = "leaf"))

  # 5. Chance baselines
  C_class     <- 2L
  chance_class <- 1 / C_class
  k_attr      <- NA_integer_

  # 6. Enumerate admissible rules -------------------------------------------

  # Placeholder for rule built from direction_map (categorical fixed-partition).
  # Set in the categorical branch; used in step 9 to skip reconstruction.
  rule_from_dmap <- NULL

  # Accumulate candidate fields as parallel vectors; build cand_df once after
  # enumeration. Avoids repeated data.frame() + deparse + make.names overhead
  # in the inner candidate loop (the hot path for large ordered attributes).
  cand_j         <- integer(0)
  cand_direction <- character(0)
  cand_cut_value <- numeric(0)
  cand_pac       <- numeric(0)
  cand_ess_class <- numeric(0)
  cand_ess_attr  <- numeric(0)
  cand_ess       <- numeric(0)
  cand_sens_1    <- numeric(0)
  cand_spec_0    <- numeric(0)
  cand_mean_pac  <- numeric(0)
  cand_balance   <- numeric(0)
  preds_list     <- list()   # parallel to cand_* vectors; indexed by row_id

  add_candidate <- function(rule, j_index = NA_integer_) {
    y_pred    <- oda_rule_predict(x, rule)
    conf      <- oda_confusion_binary(y, y_pred, w)
    mp        <- conf$mean_pac
    pac       <- mp * 100
    chance_a  <- 1 / k_attr
    ess_class <- oda_ess_from_meanpac(mp, chance_class)
    ess_attr  <- oda_ess_from_meanpac(mp, chance_a)
    ess_obj   <- if (chance_model == "class") ess_class else ess_attr

    preds_list[[length(preds_list) + 1L]] <<- y_pred
    cand_j         <<- c(cand_j,        j_index)
    cand_direction <<- c(cand_direction, rule$direction)
    cand_cut_value <<- c(cand_cut_value,
                          if (!is.null(rule$cut_value)) rule$cut_value else NA_real_)
    cand_pac       <<- c(cand_pac,       pac)
    cand_ess_class <<- c(cand_ess_class, ess_class)
    cand_ess_attr  <<- c(cand_ess_attr,  ess_attr)
    cand_ess       <<- c(cand_ess,       ess_obj)
    cand_sens_1    <<- c(cand_sens_1,    conf$sensitivity)
    cand_spec_0    <<- c(cand_spec_0,    conf$specificity)
    cand_mean_pac  <<- c(cand_mean_pac,  mp)
    cand_balance   <<- c(cand_balance,   abs(conf$sensitivity - conf$specificity))
  }

  # --- binary attribute ---
  if (attr_type == "binary") {
    levs <- sort(unique(as.character(x)))
    if (length(levs) != 2L)
      return(list(ok = FALSE, reason = "not_binary", type = "leaf"))
    k_attr  <- 2L
    n_left  <- sum(x == levs[1L])
    n_right <- sum(x == levs[2L])
    if (n_left >= mindenom && n_right >= mindenom) {
      rules  <- list(
        list(type="binary_map", direction="0->1", left_levels=levs[1],right_levels=levs[2]),
        list(type="binary_map", direction="1->0", left_levels=levs[1],right_levels=levs[2])
      )
      for (r in rules) add_candidate(r, j_index = NA_integer_)
    }
  }

  # --- categorical attribute ---
  if (attr_type == "categorical") {
    x_fac <- factor(x)
    levs   <- levels(x_fac)
    k_attr <- length(levs)
    if (k_attr < 2L) return(list(ok = FALSE, reason = "no_levels", type = "leaf"))

    if (!is.null(direction_map)) {
      # Fixed-partition directional: evaluate only the specified mapping.
      dm_names <- as.character(names(direction_map))
      dm_vals  <- as.integer(direction_map)
      if (!setequal(dm_names, levs))
        return(list(ok = FALSE,
                    reason = "direction_map_levels_mismatch",
                    type   = "leaf"))
      if (!all(dm_vals %in% c(0L, 1L)))
        return(list(ok = FALSE,
                    reason = "direction_map_values_not_binary",
                    type   = "leaf"))
      left_levels  <- dm_names[dm_vals == 0L]
      right_levels <- dm_names[dm_vals == 1L]
      if (length(left_levels) == 0L || length(right_levels) == 0L)
        return(list(ok = FALSE,
                    reason = "direction_map_one_sided",
                    type   = "leaf"))
      # Build and register the single fixed rule
      rule_from_dmap <- list(type        = "nominal_cut",
                             direction   = "0->1",
                             left_levels = left_levels,
                             right_levels = right_levels)
      add_candidate(rule_from_dmap, j_index = NA_integer_)
    } else {
      # Nondirectional: exhaustive ordered partition search
      # Order by weighted class-1 rate
      p1 <- sapply(levs, function(l) {
        wl <- w[x_fac == l]; yl <- y[x_fac == l]
        if (sum(wl) == 0) NA_real_ else sum(wl * yl) / sum(wl)
      })
      levs_ord <- levs[order(p1, na.last = TRUE)]

      for (j in seq_len(k_attr - 1L)) {
        n_left  <- sum(x_fac %in% levs_ord[seq_len(j)])
        n_right <- sum(x_fac %in% levs_ord[(j + 1L):k_attr])
        if (n_left < mindenom || n_right < mindenom) next
        left  <- levs_ord[seq_len(j)]
        right <- levs_ord[(j + 1L):k_attr]
        for (dir in c("0->1","1->0")) {
          rule <- list(type="nominal_cut", direction=dir,
                       left_levels=left, right_levels=right)
          add_candidate(rule, j_index = j)
        }
      }
    }
  }

  # --- ordered attribute ---
  if (attr_type == "ordered") {
    blocks <- oda_make_blocks_ordered(x, y, w)
    z      <- blocks$z; o <- blocks$o; x_rep <- blocks$x_rep; n_v <- blocks$n_v
    m      <- length(z)
    if (m < 2L) return(list(ok = FALSE, reason = "no_blocks", type = "leaf"))
    k_attr <- m

    Z  <- cumsum(z); O  <- cumsum(o)
    Zm <- sum(z);    Om <- sum(o)
    N_v <- cumsum(n_v); Nm <- sum(n_v)

    admissible <- function(j, dir) {
      if (dir == "0->1") return(Z[j] > 0 && (Om - O[j]) > 0)
      if (dir == "1->0") return(O[j] > 0 && (Zm - Z[j]) > 0)
      FALSE
    }

    for (j in seq_len(m - 1L)) {
      if (N_v[j] < mindenom || (Nm - N_v[j]) < mindenom) next
      for (dir in c("0->1","1->0")) {
        if (!admissible(j, dir)) next
        cut_v <- (x_rep[j] + x_rep[j + 1L]) / 2
        rule  <- list(type="ordered_cut", direction=dir, cut_value=cut_v)
        add_candidate(rule, j_index = j)
      }
    }
  }

  if (length(preds_list) == 0L)
    return(list(ok = FALSE, reason = "no_admissible_cut", type = "leaf"))

  cand_df <- data.frame(
    row_id    = seq_along(preds_list),
    j         = cand_j,
    direction = cand_direction,
    cut_value = cand_cut_value,
    pac       = cand_pac,
    ess_class = cand_ess_class,
    ess_attr  = cand_ess_attr,
    ess       = cand_ess,
    sens_1    = cand_sens_1,
    spec_0    = cand_spec_0,
    mean_pac  = cand_mean_pac,
    balance   = cand_balance,
    stringsAsFactors = FALSE
  )

  # 6b. Direction filter (ordered and binary attributes; categorical already
  # returned early above). Restricts the candidate set to the declared
  # direction before tie-breaking, so the same constraint applies to
  # training, MC permutation refits, and LOO fold refits.
  if (direction != "off") {
    target_dir <- if (direction == "greater") "0->1" else "1->0"
    cand_df    <- cand_df[cand_df$direction == target_dir, , drop = FALSE]
    if (nrow(cand_df) == 0L)
      return(list(ok = FALSE,
                  reason = "no_admissible_cut_for_direction",
                  type   = "leaf"))
  }

  # 7. Tie-breaking defaults
  if (is.null(primary))   primary   <- if (priors_on) "maxsens" else "meansens"
  if (is.null(secondary)) secondary <- "samplerep"

  # 8. Select best candidate
  cand_best <- oda_apply_primary_secondary(cand_df, primary, secondary,
                                           y = y, w = w, preds_list = preds_list)
  best_row  <- cand_best[1L, ]
  best_pred <- preds_list[[ best_row$row_id ]]

  # 9. Reconstruct winning rule
  if (attr_type == "binary") {
    levs      <- sort(unique(as.character(x)))
    rule_best <- list(type="binary_map", direction=best_row$direction,
                      left_levels=levs[1], right_levels=levs[2])
  } else if (attr_type == "categorical") {
    if (!is.null(rule_from_dmap)) {
      # Fixed-partition rule already fully constructed above
      rule_best <- rule_from_dmap
    } else {
      x_fac     <- factor(x)
      levs      <- levels(x_fac)
      p1        <- sapply(levs, function(l) {
        wl <- w[x_fac == l]; yl <- y[x_fac == l]
        if (sum(wl) == 0) NA_real_ else sum(wl * yl) / sum(wl)
      })
      levs_ord  <- levs[order(p1, na.last = TRUE)]
      j         <- best_row$j
      left_raw  <- levs_ord[seq_len(j)]
      right_raw <- levs_ord[(j + 1L):length(levs_ord)]
      # MegaODA convention: left always = side predicting 0
      if (best_row$direction == "0->1") {
        left_lev <- left_raw;  right_lev <- right_raw
      } else {
        left_lev <- right_raw; right_lev <- left_raw
      }
      rule_best <- list(type="nominal_cut", direction=best_row$direction,
                        left_levels=left_lev, right_levels=right_lev)
    }
  } else {
    rule_best <- list(type="ordered_cut", direction=best_row$direction,
                      cut_value=best_row$cut_value)
  }

  # 10. Recompute confusion for winner
  # confusion: raw counts (w_raw = unit weights); confusion_wt: priors-weighted
  conf_best_raw <- oda_confusion_binary(y, best_pred, w_raw)
  conf_best     <- oda_confusion_binary(y, best_pred, w)   # weighted (used for ESS)
  mp_best    <- conf_best$mean_pac
  chance_a   <- 1 / k_attr
  ess_class  <- oda_ess_from_meanpac(mp_best, chance_class)
  ess_attr   <- oda_ess_from_meanpac(mp_best, chance_a)
  ess_obj    <- if (chance_model == "class") ess_class else ess_attr
  pac        <- mp_best * 100

  # 10a. Early LOO stability gate (loo_then_mc mode only).
  # For ordered_cut rules with uniform weights, the algebraic count-table LOO
  # can determine stability in O(k^2) time  -  far cheaper than MC. If it proves
  # instability, reject before burning MC iterations.
  # Conditions: eval_order=="loo_then_mc", loo=="stable", mcarlo==TRUE,
  #             ordered_cut rule, algebraic helper applicable (uniform w).
  if (eval_order == "loo_then_mc" && loo == "stable" &&
      isTRUE(mcarlo) && identical(rule_best$type, "ordered_cut")) {
    .pre_alg <- oda_loo_ordered_cut_counts(x, y, w_raw, priors_on, rule_best)
    if (!is.null(.pre_alg)) {
      .w_pre      <- oda_apply_priors(y, w_raw, priors_on)
      .conf_pre   <- oda_confusion_binary(y, .pre_alg, .w_pre)
      .chance_pre <- if (chance_model == "class") chance_class else chance_a
      .ess_pre    <- oda_ess_from_meanpac(.conf_pre$mean_pac, .chance_pre)
      if (!isTRUE(all.equal(ess_obj, .ess_pre, tolerance = 1e-12)))
        return(list(ok = FALSE, reason = "loo_not_stable", type = "leaf"))
    }
  }

  # 11. Monte Carlo p-value
  mc_res <- NULL
  p_mc   <- NA_real_
  if (isTRUE(mcarlo)) {
    mc_res <- oda_mc_p_value(
      x = x, y = y, w = w_raw,  # raw weights; MC loop re-applies priors per permutation
      attr_type     = attr_type,
      priors_on     = priors_on,
      primary       = primary,
      secondary     = secondary,
      miss_codes    = miss_codes,
      chance_model  = chance_model,
      mc_iter       = mc_iter,
      mc_target     = mc_target,
      mc_stop       = mc_stop,
      mc_stopup     = mc_stopup,
      mc_adjust     = mc_adjust,
      seed          = mc_seed,
      ess_obs       = ess_obj,
      direction     = direction,
      direction_map = direction_map
    )
    p_mc <- mc_res$p_mc
  }

  # 12. LOO
  loo_out <- NULL
  if (loo != "off") {
    loo_out <- oda_loo_for_rule(
      x = x, y = y, rule = rule_best, w = w_raw,
      chance_model  = chance_model,
      k_attr        = k_attr,
      attr_type     = attr_type,
      priors_on     = priors_on,
      primary       = primary,
      secondary     = secondary,
      miss_codes    = miss_codes,
      mc_iter       = mc_iter,
      mc_target     = mc_target,
      mc_stop       = mc_stop,
      mc_stopup     = mc_stopup,
      mc_adjust     = mc_adjust,
      mc_seed       = if (is.null(mc_seed)) NULL else mc_seed + 1L,
      direction     = direction,
      direction_map = direction_map
    )

    if (isTRUE(loo_out$allowed)) {
      if (loo == "stable") {
        stable <- isTRUE(all.equal(ess_obj, loo_out$ess_loo, tolerance = 1e-12))
        if (!stable) return(list(ok=FALSE, reason="loo_not_stable", type="leaf"))
      }
      if (loo == "pvalue") {
        if (is.na(loo_out$p_value) || loo_out$p_value >= loo_alpha)
          return(list(ok=FALSE, reason="loo_p_not_significant", type="leaf"))
      }
    }
  }

  list(
    ok           = TRUE,
    reason       = NULL,
    type         = rule_best$type,
    rule         = rule_best,
    attr_type    = attr_type,
    k_attr       = k_attr,
    n_eff        = n_eff,
    weights      = w,
    # confusion: raw integer counts (TP/TN/FP/FN) + sensitivity/specificity as proportions
    confusion    = conf_best_raw,
    confusion_wt = conf_best,   # priors-weighted version (used internally for ESS)
    ess_class = ess_class,
    ess_attr  = ess_attr,
    ess       = ess_obj,
    pac       = pac,
    p_mc      = p_mc,
    mc_info   = mc_res,
    loo       = loo_out,
    tie_block = cand_df
  )
}
