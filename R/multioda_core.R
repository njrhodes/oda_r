###############################################################################
# R/multioda_core.R
#
# MultiODA Core - univariate multiclass ODA engine.
#
# Spec (MegaODA-faithful):
#   PRIMARY   = MAXSENS  (overall PAC in priors-weighted objective space)
#   SECONDARY = SAMPLEREP  (L1 distance between predicted and observed
#                           class relative frequencies)
#   TERTIARY  = FIRST IDENTIFIED  (strict enumeration order via tick())
#   LOO       = true refit-per-fold  ("refit") OR fixed global grid ("fixed")
#   MC        = Fisher randomization against number-correct statistic
#
# Key invariants:
#   - PRIORS ON:  objective = overall PAC on priors-weighted counts
#                 -> equivalent to maximising mean PAC on raw counts
#   - DEGEN OFF (default): all C classes must appear in predicted labels
#   - DEGEN ON:  degenerate solutions allowed; PRIORS forced OFF
#
# Depends on: R/utils.R, R/unioda_core.R
###############################################################################

ODA_CRITERIA <- c("maxsens","meansens","samplerep","balanced","distance","random")

# ---- Weighting policy ------------------------------------------------------ #
#' Enforce weighting policy for multiclass ODA
#' @keywords internal

oda_enforce_weighting_policy <- function(
    attr_type_res, priors_on, loo, w_case, reason_prefix = ""
) {
  priors_on_eff      <- isTRUE(priors_on)
  loo_eff            <- loo
  notes              <- character(0)
  weighted_requested <- priors_on_eff || any(as.numeric(w_case) != 1)
  allow_weighted_loo <- (attr_type_res == "ordered")

  if (loo_eff == "on" && !allow_weighted_loo && weighted_requested) {
    notes <- c(notes, paste0(reason_prefix,
        "categorical LOO requested with weighting; must be disallowed"))
  }
  list(priors_on_eff      = priors_on_eff,
       loo_eff            = loo_eff,
       allow_weighted_loo = allow_weighted_loo,
       weighted_requested = weighted_requested,
       notes              = notes)
}

# ---- Priors (multiclass) --------------------------------------------------- #

#' Class-normalise weights so each class sums to 1 (PRIORS ON).
oda_apply_priors_multiclass <- function(y, w, priors_on = TRUE) {
  y <- as.integer(y); w <- as.numeric(w)
  if (!isTRUE(priors_on)) return(w)
  classes <- sort(unique(y))
  w_adj   <- w
  for (cl in classes) {
    idx <- which(y == cl)
    if (!length(idx)) next
    tot <- sum(w_adj[idx])
    if (tot > 0) w_adj[idx] <- w_adj[idx] / tot
  }
  w_adj
}

# ---- Multiclass confusion table -------------------------------------------- #

#' Compute weighted multiclass confusion matrix, PAC and PV by class.
oda_confusion_multiclass <- function(y, y_pred, w = NULL) {
  y <- as.integer(y); y_pred <- as.integer(y_pred)
  n <- length(y)
  if (is.null(w)) w <- rep(1, n)
  w <- as.numeric(w)

  classes <- sort(unique(c(y, y_pred)))
  C       <- length(classes)

  conf <- matrix(0, nrow = C, ncol = C,
                 dimnames = list(actual = classes, predicted = classes))
  for (i in seq_len(n)) {
    ai <- match(y[i],      classes)
    pi <- match(y_pred[i], classes)
    conf[ai, pi] <- conf[ai, pi] + w[i]
  }

  total       <- sum(conf)
  correct     <- sum(diag(conf))
  overall_acc <- if (total > 0) correct / total else NA_real_

  row_tot      <- rowSums(conf)
  pac_by_class <- ifelse(row_tot > 0, diag(conf) / row_tot, NA_real_)
  mean_pac     <- mean(pac_by_class, na.rm = TRUE)

  col_tot      <- colSums(conf)
  pv_by_class  <- ifelse(col_tot > 0, diag(conf) / col_tot, NA_real_)
  mean_pv      <- mean(pv_by_class, na.rm = TRUE)

  list(classes      = classes,
       confusion    = conf,
       total        = total,
       correct      = correct,
       overall_acc  = overall_acc,
       pac_by_class = pac_by_class,
       mean_pac     = mean_pac,
       pv_by_class  = pv_by_class,
       mean_pv      = mean_pv)
}

#' ESS from mean metric (PAC or PV) for C-class problem.
oda_ess_from_mean <- function(mean_metric, C) {
  chance <- 1 / C
  if (!is.finite(mean_metric) || chance >= 1) return(NA_real_)
  100 * (mean_metric - chance) / (1 - chance)
}

#' Non-degeneracy check: predicted labels must cover all C classes.
oda_is_nondegenerate_labels <- function(pred_labels, C) {
  length(unique(as.integer(pred_labels))) == C
}

# ---- Ordered block builders ------------------------------------------------ #

#' Build per-class count matrix over unique-x blocks (objective or raw space).
oda_make_blocks_ordered_multiclass <- function(x, y, w) {
  ord      <- order(x)
  x_ord    <- x[ord]; y_ord <- as.integer(y[ord]); w_ord <- as.numeric(w[ord])
  block_id <- cumsum(c(TRUE, diff(x_ord) != 0))
  m        <- max(block_id)
  classes  <- sort(unique(y_ord))
  C        <- length(classes)

  counts <- matrix(0, nrow = m, ncol = C)
  x_rep  <- numeric(m)
  for (j in seq_len(m)) {
    idx    <- which(block_id == j)
    x_rep[j] <- x_ord[idx[1L]]
    for (c_idx in seq_len(C)) {
      cl <- classes[c_idx]
      counts[j, c_idx] <- sum(w_ord[idx][y_ord[idx] == cl])
    }
  }
  list(ord = ord, x_rep = x_rep, counts = counts, classes = classes)
}

#' Build global block grid on full x (for LOO fixed-grid mode).
oda_make_global_blocks_ordered <- function(x) {
  ord      <- order(x)
  x_ord    <- x[ord]
  block_id <- cumsum(c(TRUE, diff(x_ord) != 0))
  m        <- max(block_id)
  first_idx <- which(c(TRUE, diff(x_ord) != 0))
  x_rep    <- x_ord[first_idx]

  block_of_ord <- as.integer(block_id)
  block_of_obs <- integer(length(x))
  block_of_obs[ord] <- block_of_ord
  list(ord = ord, x_rep = x_rep, block_of_obs = block_of_obs, m = m)
}

# ---- Candidate metrics ----------------------------------------------------- #
#' Compute PAC, SAMPLEREP, and other metrics for a partition
#' @keywords internal

oda_metrics_candidate <- function(NP_obj, NP_raw, NA_raw, x_rep,
                                   cuts_idx, priors_on_eff = TRUE) {
  TP              <- diag(NP_obj)
  NA_actual_obj   <- rowSums(NP_obj)
  pac_by_class    <- ifelse(NA_actual_obj > 0, TP / NA_actual_obj, NA_real_)
  mean_pac        <- mean(pac_by_class, na.rm = TRUE)
  den_obj         <- sum(NP_obj)
  overall_pac     <- if (den_obj > 0) sum(TP) / den_obj else -Inf

  # SAMPLEREP: always in raw-count space (guidebook: compare predicted to observed
  # sample relative frequencies  -  these are raw counts, independent of weighting).
  NP_raw_pred     <- colSums(NP_raw)
  den_pred_raw    <- sum(NP_raw_pred)
  den_act_raw     <- sum(NA_raw)
  srep_dist_raw   <- if (den_pred_raw <= 0 || den_act_raw <= 0) Inf else
    sum(abs((NP_raw_pred / den_pred_raw) - (NA_raw / den_act_raw)))

  # Keep obj-space version for diagnostics only; never used for selection
  NP_obj_pred     <- colSums(NP_obj)
  den_pred_obj    <- sum(NP_obj_pred)
  den_act_obj     <- sum(NA_actual_obj)
  srep_dist_obj   <- if (den_pred_obj <= 0 || den_act_obj <= 0) Inf else
    sum(abs((NP_obj_pred / den_pred_obj) - (NA_actual_obj / den_act_obj)))

  # ALWAYS use raw-space SREP for selection regardless of priors setting
  srep_dist <- srep_dist_raw

  pac_ok   <- pac_by_class[is.finite(pac_by_class)]
  bal_diff <- if (length(pac_ok) >= 2L) (max(pac_ok) - min(pac_ok)) else Inf

  dist <- Inf
  if (length(cuts_idx)) {
    d <- numeric(length(cuts_idx))
    for (k in seq_along(cuts_idx)) {
      j      <- cuts_idx[k]
      left_d  <- if (j > 1L)                 x_rep[j] - x_rep[j - 1L] else Inf
      right_d <- if (j < length(x_rep))      x_rep[j + 1L] - x_rep[j] else Inf
      d[k]   <- min(left_d, right_d)
    }
    dist <- max(d)
  }

  list(overall_pac    = overall_pac,
       mean_pac       = mean_pac,
       pac_by_class   = pac_by_class,
       srep_dist      = srep_dist,
       srep_dist_raw  = srep_dist_raw,
       srep_dist_obj  = srep_dist_obj,
       bal_diff       = bal_diff,
       neg_srep       = -srep_dist,
       neg_bal        = -bal_diff,
       neg_dist       = -dist)
}

# ---- Ordered multiclass partition selector --------------------------------- #

#' Select the best K-segment ordered partition by MegaODA spec:
#'   PRIMARY -> SECONDARY -> FIRST IDENTIFIED (enum order via tick()).
#'
#' @param x_rep  Representative x value per unique block.
#' @param counts_obj  m x C count matrix in objective (priors-weighted) space.
#' @param counts_raw  m x C count matrix in raw (case-weight) space.
#' @param K  Number of segments (cuts = K-1).
#' @param priors_on_eff  Logical.
#' @param degen  Allow degenerate solutions?
#' @param primary,secondary  Heuristic strings (NULL = spec defaults).
#' @param cut_value_mode  "midpoint", "lower", or "upper".
#' @param debug_return_ties  Return all primary-tied candidates for diagnostics.
#' @param debug_max_ties  Cap on number of ties stored.
#' @param direction  Directional constraint (MPE Chapter 4 ordered DIRECTIONAL).
#'   "ascending" forces segment s to map to class s; "descending" forces class
#'   C+1-s. Default "off" (nondirectional; all assignments evaluated).
#' @return List with ok, cuts_idx, cut_values, seg_cls_idx, primary_obj,
#'   secondary_obj, best_enum_id, ties, classes.
oda_best_ordered_multiclass_partition <- function(
    x_rep, counts_obj, counts_raw, K,
    priors_on_eff   = TRUE,
    degen           = FALSE,
    primary         = NULL,
    secondary       = NULL,
    cut_value_mode  = c("midpoint","lower","upper"),
    debug_return_ties = FALSE,
    debug_max_ties  = 200L,
    direction       = "off"
) {
  cut_value_mode <- match.arg(cut_value_mode)

  if (is.null(primary))   primary   <- if (isTRUE(priors_on_eff)) "maxsens" else "meansens"
  if (is.null(secondary)) secondary <- "samplerep"
  if (!(primary   %in% ODA_CRITERIA)) stop("Invalid PRIMARY: ",   primary)
  if (!(secondary %in% ODA_CRITERIA)) stop("Invalid SECONDARY: ", secondary)

  m       <- nrow(counts_obj)
  C       <- ncol(counts_obj)
  classes <- seq_len(C)

  K <- max(2L, min(as.integer(K), m))
  if (!isTRUE(degen) && K < C)
    return(list(ok = FALSE, reason = "no_solution_degen_off_K_lt_C"))

  pref_obj <- apply(counts_obj, 2, cumsum)
  pref_raw <- apply(counts_raw, 2, cumsum)
  seg_sum  <- function(pref, i, j)
    if (i == 1L) pref[j, ] else pref[j, ] - pref[i - 1L, ]

  NA_raw <- colSums(counts_raw)

  best_primary   <- -Inf
  best_secondary <- -Inf
  best_enum_id   <- 2147483647L
  best_cuts      <- NULL
  best_seg       <- NULL

  need_ties <- isTRUE(debug_return_ties) || identical(secondary, "samplerep")
  ties       <- if (need_ties) list() else NULL

  enum_id  <- 0L
  current  <- integer(K)

  cut_grid <- if (K == 2L) {
    matrix(seq_len(m - 1L), ncol = 1L)
  } else {
    t(utils::combn(seq_len(m - 1L), K - 1L))
  }

  cuts_to_values <- function(ci)
    vapply(ci, function(j) switch(cut_value_mode,
      midpoint = (x_rep[j] + x_rep[j + 1L]) / 2,
      lower    = x_rep[j],
      upper    = x_rep[j + 1L]), numeric(1))

  # TWO-LEVEL SELECTION (MegaODA-faithful):
  # Level 2 (inner): within a given cut position, pick the best seg assignment using
  #   primary -> SECONDARY (SAMPLEREP) -> first assignment enumerated
  # Level 1 (outer): across cut positions, compare level-2 winners using
  #   primary -> FIRST CUT POSITION ENUMERATED  (NOT secondary/SAMPLEREP)
  #
  # Evidence: iris V1 cuts 6.15 vs 6.25 both have identical primary and only one
  # valid assignment each -> no assignment tie -> first-identified-cut picks 6.15.
  # DAT5 cuts (1.5,2.5) with assignments (1,2,3) vs (2,1,3) tied on primary ->
  # SAMPLEREP (assignment-level secondary) picks (2,1,3).

  for (r in seq_len(nrow(cut_grid))) {
    cuts_idx <- sort(as.integer(cut_grid[r, ]))
    starts   <- c(1L, cuts_idx + 1L)
    ends     <- c(cuts_idx, m)

    cand    <- vector("list", K)
    seg_obj <- vector("list", K)
    seg_raw <- vector("list", K)

    for (s in seq_len(K)) {
      v_obj      <- seg_sum(pref_obj, starts[s], ends[s])
      v_raw      <- seg_sum(pref_raw, starts[s], ends[s])
      mx         <- max(v_obj)
      cand[[s]]  <- sort(which(abs(v_obj - mx) <= 1e-12))
      seg_obj[[s]] <- v_obj
      seg_raw[[s]] <- v_raw
    }

    # Level-2 state: best assignment found at THIS cut position
    cut_best_primary   <- -Inf
    cut_best_secondary <- -Inf
    cut_best_enum_id   <- 2147483647L
    cut_best_seg       <- NULL
    cut_best_NP_raw    <- NULL
    cut_best_NP_obj    <- NULL
    cut_best_met       <- NULL

    recurse <- function(s) {
      if (s > K) {
        enum_id <<- enum_id + 1L
        if (!isTRUE(degen) && length(unique(current)) < C) return()

        NP_obj <- matrix(0, C, C); NP_raw <- matrix(0, C, C)
        for (t in seq_len(K)) {
          cls <- current[t]
          NP_obj[, cls] <- NP_obj[, cls] + seg_obj[[t]]
          NP_raw[, cls] <- NP_raw[, cls] + seg_raw[[t]]
        }

        met <- oda_metrics_candidate(NP_obj, NP_raw, NA_raw,
                                     x_rep, cuts_idx, priors_on_eff)

        sp <- switch(primary,
          maxsens   = met$overall_pac,
          meansens  = met$mean_pac,
          samplerep = met$neg_srep,
          balanced  = met$neg_bal,
          distance  = met$neg_dist,
          random    = stats::runif(1))

        ss <- switch(secondary,
          maxsens   = met$overall_pac,
          meansens  = met$mean_pac,
          samplerep = met$neg_srep,
          balanced  = met$neg_bal,
          distance  = met$neg_dist,
          random    = stats::runif(1))

        sp_t <- tick(sp); ss_t <- tick(ss)
        cp_t <- tick(cut_best_primary); cs_t <- tick(cut_best_secondary)

        # Level-2 comparison: primary -> secondary (SAMPLEREP) -> first enum
        cut_better <- FALSE
        if (!is.na(sp_t) && !is.na(cp_t)) {
          if      (sp_t > cp_t)  cut_better <- TRUE
          else if (sp_t == cp_t && !is.na(ss_t) && !is.na(cs_t)) {
            if      (ss_t > cs_t)  cut_better <- TRUE
            else if (ss_t == cs_t && enum_id < cut_best_enum_id) cut_better <- TRUE
          }
        }

        if (cut_better) {
          cut_best_primary   <<- sp
          cut_best_secondary <<- ss
          cut_best_enum_id   <<- enum_id
          cut_best_seg       <<- current
          cut_best_NP_raw    <<- NP_raw
          cut_best_NP_obj    <<- NP_obj
          cut_best_met       <<- met
        }

        # Store in global ties list for diagnostics (keyed on primary score)
        if (need_ties && tick(sp) == tick(best_primary) &&
            length(ties) < debug_max_ties) {
          ties[[length(ties) + 1L]] <<- list(
            enum_id        = enum_id,
            primary        = sp,
            secondary      = ss,
            primary_tick   = sp_t,
            secondary_tick = ss_t,
            overall_pac    = met$overall_pac,
            mean_pac       = met$mean_pac,
            srep_dist      = met$srep_dist,
            srep_dist_raw  = met$srep_dist_raw,
            srep_dist_obj  = met$srep_dist_obj,
            bal_diff       = met$bal_diff,
            pac_by_class   = met$pac_by_class,
            cuts_idx       = cuts_idx,
            cut_values     = cuts_to_values(cuts_idx),
            seg_cls_idx    = current,
            pred_raw_counts = colSums(NP_raw),
            act_raw_counts  = NA_raw
          )
        }
        return()
      }
      if (direction %in% c("ascending", "descending")) {
        # Constrain segment-to-class assignment (MPE Chapter 4 DIRECTIONAL)
        cl <- if (direction == "ascending") s else C + 1L - s
        current[s] <<- cl
        recurse(s + 1L)
      } else {
        for (cl in cand[[s]]) { current[s] <<- cl; recurse(s + 1L) }
      }
    }

    recurse(1L)

    # Level-1 comparison: compare this cut's winner to global best
    # Use: primary -> FIRST CUT POSITION (no secondary/SAMPLEREP across cuts)
    # NOTE: plain <- (not <<-) because this code is in the for-loop body
    # (not inside recurse), so we are already in the function's local environment.
    if (!is.null(cut_best_seg)) {
      cp_t <- tick(cut_best_primary); bp_t <- tick(best_primary)
      if (!is.na(cp_t) && !is.na(bp_t)) {
        if (cp_t > bp_t) {
          # Strictly better primary: this cut wins globally
          best_primary   <- cut_best_primary
          best_secondary <- cut_best_secondary
          best_enum_id   <- cut_best_enum_id
          best_cuts      <- cuts_idx
          best_seg       <- cut_best_seg
        }
        # If cp_t == bp_t: first-identified cut wins -> do nothing (keep earlier cut)
        # If cp_t < bp_t: this cut is worse -> do nothing
      }
    }
  }

  if (is.null(best_cuts)) return(list(ok = FALSE, reason = "no_solution"))

  list(ok            = TRUE,
       cuts_idx      = best_cuts,
       cut_values    = cuts_to_values(best_cuts),
       seg_cls_idx   = best_seg,
       primary_obj   = best_primary,
       secondary_obj = best_secondary,
       best_enum_id  = best_enum_id,
       ties          = ties,
       classes       = classes)
}

# ---- Multiclass prediction ------------------------------------------------- #

#' Predict class labels given a multiclass ODA rule.
oda_rule_predict_multiclass <- function(
    x, rule, boundary = c("megaoda_halfopen","right_closed")
) {
  boundary <- match.arg(boundary)

  if (rule$type == "multiclass_ordered") {
    cuts    <- rule$cut_values
    seg_cls <- rule$seg_classes
    K       <- length(seg_cls)
    if (K == 1L) return(rep(seg_cls[1L], length(x)))

    idx <- if (boundary == "megaoda_halfopen")
      findInterval(x, vec = cuts, rightmost.closed = FALSE, left.open = FALSE) + 1L
    else
      findInterval(x, vec = cuts, rightmost.closed = TRUE,  left.open = TRUE)  + 1L

    idx[idx < 1L] <- 1L
    idx[idx > K]  <- K
    return(seg_cls[idx])
  }

  if (rule$type == "multiclass_nominal") {
    levs <- rule$levels
    map  <- rule$level_class
    xf   <- factor(x, levels = levs)
    li   <- as.integer(xf)
    li[is.na(li)] <- 0L
    out  <- integer(length(li))
    for (i in seq_along(li))
      out[i] <- if (li[i] == 0L) rule$default_class else map[li[i]]
    return(out)
  }

  stop("oda_rule_predict_multiclass: unsupported rule$type = '", rule$type, "'")
}

# ---- Monte Carlo p-value (multiclass) ------------------------------------- #
#' Monte Carlo p-value for multiclass ODA
#' @keywords internal

oda_mc_p_value_multiclass <- function(
    x, y, w, attr_type, priors_on, degen, K_segments,
    mc_iter   = 25000L, mc_target = 0.05,
    mc_stop   = 99.9,   mc_stopup = NA,
    mc_adjust = FALSE,  seed = NULL, observed_mean_pac,
    direction = "off", direction_map = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  n           <- length(y)
  ge          <- 0L; used <- 0L
  conf_level  <- if (!is.na(mc_stop)   && mc_stop   > 1) mc_stop   / 100 else mc_stop
  stopup_lv   <- if (!is.na(mc_stopup) && mc_stopup > 1) mc_stopup / 100 else mc_stopup
  min_check   <- 50L; check_every <- 50L

  for (b in seq_len(mc_iter)) {
    used   <- b
    y_star <- sample(y, size = n, replace = FALSE)

    fit_b <- oda_multiclass_unioda_core(
      x = x, y = y_star, w = w,
      attr_type  = attr_type, priors_on = priors_on, degen = degen,
      miss_codes = NULL,       K_segments = K_segments,
      mcarlo = FALSE,          loo = "off",
      mc_iter = 0L, mc_target = mc_target, mc_stop = mc_stop,
      mc_stopup = mc_stopup,   mc_adjust = mc_adjust, mc_seed = NULL,
      direction = direction,   direction_map = direction_map
    )

    # Compare on the training-objective scale (mean PAC, percent units) so that
    # priors-weighted and unbalanced designs are assessed consistently.
    # Raw correct count is NOT used: it diverges from mean PAC when class sizes
    # are unbalanced (issue #7).
    mp_b <- if (isTRUE(fit_b$ok)) fit_b$mean_pac else 0
    if (mp_b >= observed_mean_pac - 1e-7) ge <- ge + 1L

    if (!is.na(conf_level) && used >= min_check && (used %% check_every == 0L)) {
      if (ge == 0L) {
        lower <- 0; upper <- stats::qbeta(conf_level, 1, used)
      } else if (ge == used) {
        lower <- stats::qbeta(1 - conf_level, used, 1); upper <- 1
      } else {
        upper <- stats::qbeta(conf_level, ge + 1, used - ge)
        lower <- stats::qbeta(1 - conf_level, ge, used - ge + 1)
      }
      if (!is.na(mc_target) && upper < mc_target) break
      if (!is.na(stopup_lv) && lower > stopup_lv) break
    }
  }

  p_mc <- if (used <= 0L) NA_real_ else if (ge == 0L) 0.0 else ge / used
  list(p_mc = p_mc, ge_count = ge, iter_used = used)
}

# ---- LOO: categorical multiclass (unweighted only) ------------------------ #
#' Leave-one-out for categorical multiclass ODA
#' @keywords internal

oda_loo_multiclass <- function(
    x, y, w0, attr_type, priors_on_eff, weighted_requested,
    degen, K_segments, miss_codes = NULL
) {
  # Block only on actual non-unit case weights. priors_on_eff is an objective
  # weighting device; it does not invalidate the LOO math (each fold refit
  # uses priors_on = FALSE). Blocking on priors_on_eff wrongly suppresses
  # categorical multiclass LOO for all unweighted fits. (MULTICAT-LOO1 fix)
  if (any(as.numeric(w0) != 1))
    return(list(allowed = FALSE,
                reason  = "categorical_weighted_loo_not_allowed"))

  # Engineering guard: same rationale as oda_loo_for_rule().
  # See comment there for distinction from MPE canon "superfluous LOO".
  if (attr_type == "categorical") {
    x_obs <- x[!is.na(x)]
    if (!is.null(miss_codes)) x_obs <- x_obs[!(x_obs %in% miss_codes)]
    if (length(x_obs) > 0L && length(unique(x_obs)) == length(x_obs))
      return(list(allowed = FALSE,
                  reason  = "loo_not_supported_all_unique_categories"))
  }

  n      <- length(y); y <- as.integer(y)
  y_pred <- integer(n)

  for (i in seq_len(n)) {
    keep   <- seq_len(n) != i
    fit_i  <- oda_multiclass_unioda_core(
      x = x[keep], y = y[keep], w = rep(1, sum(keep)),
      attr_type = attr_type, priors_on = FALSE, degen = degen,
      miss_codes = miss_codes, mcarlo = FALSE, loo = "off",
      K_segments = K_segments
    )
    if (!isTRUE(fit_i$ok))
      return(list(allowed = FALSE,
                  reason  = paste0("loo_fit_failed_case_", i)))
    y_pred[i] <- oda_rule_predict_multiclass(x[i], fit_i$rule)
  }

  conf_raw <- oda_confusion_multiclass(y, y_pred, w = rep(1, n))
  list(allowed = TRUE, reason = NULL,
       confusion = conf_raw, y_pred = y_pred)
}

# ---- LOO: ordered multiclass (refit or fixed grid) ------------------------ #

#' Leave-one-out cross-validation for ordered multiclass ODA.
#'
#' @param x,y  Attribute and class vectors.
#' @param w0   Raw case weights.
#' @param priors_on_eff  Logical.
#' @param degen  Logical.
#' @param K_segments  Number of segments.
#' @param miss_codes  Optional missing codes.
#' @param cut_value_mode  "midpoint","lower","upper".
#' @param grid_mode  "refit" (true per-fold rebuild) or "fixed" (global grid).
#' @param boundary_mode  "megaoda_halfopen" or "right_closed".
#' @param loo_use_samplerep  Include samplerep in fold selection.
#' @param loo_return_folds  Return per-fold rules and debug info.
#' @param loo_priors_mode  "fold" (renorm each fold) or "global" (global wts).
#' @return List with allowed, confusion_raw, confusion_weighted, y_pred, and
#'   optional fold_rule, fold_debug, fold_best_enum_id.
oda_loo_multiclass_ordered <- function(
    x, y, w0,
    priors_on_eff,
    degen,
    K_segments,
    miss_codes        = NULL,
    cut_value_mode    = c("midpoint","lower","upper"),
    grid_mode         = c("fixed","refit"),
    boundary_mode     = c("megaoda_halfopen","right_closed"),
    loo_use_samplerep = FALSE,
    loo_return_folds  = FALSE,
    loo_priors_mode   = c("fold","global")
) {
  cut_value_mode  <- match.arg(cut_value_mode)
  grid_mode       <- match.arg(grid_mode)
  boundary_mode   <- match.arg(boundary_mode)
  loo_priors_mode <- match.arg(loo_priors_mode)

  y    <- as.integer(y); w0 <- as.numeric(w0)
  n    <- length(y);      C  <- length(sort(unique(y)))
  need_ties <- isTRUE(loo_return_folds) || isTRUE(loo_use_samplerep)

  fold_rule        <- if (loo_return_folds) vector("list", n) else NULL
  fold_debug       <- if (loo_return_folds) vector("list", n) else NULL
  fold_best_enum_id <- if (loo_return_folds) integer(n) else NULL

  w_obj_global <- NULL
  if (isTRUE(priors_on_eff) && loo_priors_mode == "global")
    w_obj_global <- oda_apply_priors_multiclass(y, w0, priors_on = TRUE)

  y_pred <- integer(n)

  # ---- refit mode: rebuild blocks per fold --------------------------------
  if (grid_mode == "refit") {
    for (i in seq_len(n)) {
      keep <- seq_len(n) != i

      w_fold_obj <- if (isTRUE(priors_on_eff)) {
        if (loo_priors_mode == "global") w_obj_global[keep]
        else oda_apply_priors_multiclass(y[keep], w0[keep], priors_on = TRUE)
      } else w0[keep]

      blocks_obj  <- oda_make_blocks_ordered_multiclass(x[keep], y[keep], w_fold_obj)
      blocks_raw  <- oda_make_blocks_ordered_multiclass(x[keep], y[keep], w0[keep])
      x_rep       <- blocks_obj$x_rep
      counts_obj  <- blocks_obj$counts
      counts_raw  <- blocks_raw$counts
      classes_fold <- blocks_obj$classes
      C_fold      <- length(classes_fold)
      m           <- nrow(counts_obj)

      if (m < 2L)
        return(list(allowed = FALSE,
                    reason = paste0("loo_no_blocks_case_", i)))

      K <- as.integer(K_segments)
      if (K > m) K <- m
      if (K < 2L) K <- 2L
      if (!isTRUE(degen) && K < C_fold)
        return(list(allowed = FALSE,
                    reason = paste0("loo_no_solution_case_", i, "_K_lt_C")))

      sel <- oda_best_ordered_multiclass_partition(
        x_rep = x_rep, counts_obj = counts_obj, counts_raw = counts_raw,
        K = K, priors_on_eff = priors_on_eff, degen = degen,
        cut_value_mode = cut_value_mode,
        primary = NULL, secondary = NULL,
        debug_return_ties = need_ties, debug_max_ties = 200L
      )

      if (!isTRUE(sel$ok))
        return(list(allowed = FALSE,
                    reason = paste0("loo_no_solution_case_", i, "_",
                                    sel$reason %||% "no_solution")))

      if (!is.null(fold_debug))
        fold_debug[[i]] <- if (!is.null(sel$ties) && length(sel$ties))
          lapply(sel$ties, function(tr) {
            tr$seg_classes <- classes_fold[as.integer(tr$seg_cls_idx)]; tr
          }) else NULL

      if (!is.null(fold_best_enum_id))
        fold_best_enum_id[i] <- sel$best_enum_id %||% NA_integer_

      rule_i <- list(
        type        = "multiclass_ordered",
        cut_values  = as.numeric(sel$cut_values),
        seg_classes = as.integer(classes_fold[sel$seg_cls_idx]),
        K           = length(sel$seg_cls_idx),
        boundary    = boundary_mode
      )

      if (!is.null(fold_rule)) fold_rule[[i]] <- rule_i
      y_pred[i] <- as.integer(
        oda_rule_predict_multiclass(x[i], rule_i, boundary = boundary_mode))
    }

    conf_raw <- oda_confusion_multiclass(y, y_pred, w = rep(1, n))
    conf_w   <- oda_confusion_multiclass(
      y, y_pred,
      w = if (isTRUE(priors_on_eff)) oda_apply_priors_multiclass(y, w0, TRUE) else w0)

    return(list(allowed = TRUE, reason = NULL,
                confusion_raw = conf_raw, confusion_weighted = conf_w,
                y_pred = y_pred,
                fold_rule = fold_rule, fold_debug = fold_debug,
                fold_best_enum_id = fold_best_enum_id))
  }

  # ---- fixed global grid --------------------------------------------------
  gb          <- oda_make_global_blocks_ordered(x)
  x_rep_g     <- gb$x_rep
  block_of_obs <- gb$block_of_obs
  m           <- gb$m

  for (i in seq_len(n)) {
    keep        <- seq_len(n) != i
    counts_obj  <- matrix(0, nrow = m, ncol = C)
    counts_raw  <- matrix(0, nrow = m, ncol = C)

    w_obj_full  <- numeric(n)
    if (isTRUE(priors_on_eff)) {
      w_obj_full[keep] <- if (loo_priors_mode == "global") w_obj_global[keep]
        else oda_apply_priors_multiclass(y[keep], w0[keep], priors_on = TRUE)
    } else {
      w_obj_full[keep] <- w0[keep]
    }

    for (j in which(keep)) {
      b <- block_of_obs[j]; cl <- y[j]
      counts_obj[b, cl] <- counts_obj[b, cl] + w_obj_full[j]
      counts_raw[b, cl] <- counts_raw[b, cl] + w0[j]
    }

    if (length(which(rowSums(counts_raw) > 0)) < 2L)
      return(list(allowed = FALSE,
                  reason = paste0("loo_no_blocks_case_", i)))

    K <- as.integer(K_segments)
    if (K > m) K <- m; if (K < 2L) K <- 2L
    if (!isTRUE(degen) && K < C)
      return(list(allowed = FALSE,
                  reason = paste0("loo_no_solution_case_", i, "_K_lt_C")))

    sel <- oda_best_ordered_multiclass_partition(
      x_rep = x_rep_g, counts_obj = counts_obj, counts_raw = counts_raw,
      K = K, priors_on_eff = priors_on_eff, degen = degen,
      cut_value_mode = cut_value_mode,
      primary = NULL, secondary = NULL,
      debug_return_ties = need_ties, debug_max_ties = 200L
    )

    if (!isTRUE(sel$ok))
      return(list(allowed = FALSE,
                  reason = paste0("loo_no_solution_case_", i, "_",
                                  sel$reason %||% "no_solution")))

    if (!is.null(fold_debug))
      fold_debug[[i]] <- if (!is.null(sel$ties) && length(sel$ties))
        lapply(sel$ties, function(tr) {
          tr$seg_classes <- as.integer(tr$seg_cls_idx); tr
        }) else NULL

    if (!is.null(fold_best_enum_id))
      fold_best_enum_id[i] <- sel$best_enum_id %||% NA_integer_

    rule_i <- list(
      type        = "multiclass_ordered",
      cut_values  = as.numeric(sel$cut_values),
      seg_classes = as.integer(sel$seg_cls_idx),
      K           = length(sel$seg_cls_idx),
      boundary    = boundary_mode
    )

    if (!is.null(fold_rule)) fold_rule[[i]] <- rule_i
    y_pred[i] <- as.integer(
      oda_rule_predict_multiclass(x[i], rule_i, boundary = boundary_mode))
  }

  conf_raw <- oda_confusion_multiclass(y, y_pred, w = rep(1, n))
  conf_w   <- oda_confusion_multiclass(
    y, y_pred,
    w = if (isTRUE(priors_on_eff)) oda_apply_priors_multiclass(y, w0, TRUE) else w0)

  list(allowed = TRUE, reason = NULL,
       confusion_raw = conf_raw, confusion_weighted = conf_w,
       y_pred = y_pred,
       fold_rule = fold_rule, fold_debug = fold_debug,
       fold_best_enum_id = fold_best_enum_id)
}

# ---- Main MultiODA entry point -------------------------------------------- #

#' Fit a univariate multiclass ODA model.
#'
#' @param x  Attribute values.
#' @param y  Integer class labels (will be re-coded to 1..C internally).
#' @param w  Optional case weights.
#' @param attr_type  "auto","ordered","categorical","binary".
#' @param priors_on  Logical; inverse class-frequency weighting.
#' @param miss_codes  Additional missing-value codes (canonical name).
#' @param missing_code  Alias for \code{miss_codes} (scalar or vector).
#' @param K_segments  Number of segments (default = C).
#' @param degen  Allow degenerate solutions (forces priors_on = FALSE).
#' @param mcarlo  Run Monte Carlo p-value?
#' @param mc_iter,mc_target,mc_stop,mc_stopup,mc_adjust,mc_seed  MC params.
#' @param loo  "off" or "on".
#' @param boundary_mode  "megaoda_halfopen" or "right_closed".
#' @param loo_opts  Named list of LOO detail options. Keys:
#'   \code{cut_value_mode} ("midpoint","lower","upper"),
#'   \code{grid_mode} ("fixed","refit"),
#'   \code{use_samplerep} (logical),
#'   \code{return_folds} (logical),
#'   \code{priors_mode} ("fold","global").
#' @return Named list with ok, rule, confusion, pac, mean_pac, ess_pac,
#'   p_mc, loo, and metadata.
oda_multiclass_unioda_core <- function(
    x,
    y,
    w                = NULL,
    attr_type        = c("auto","ordered","categorical","binary"),
    priors_on        = TRUE,
    miss_codes       = NULL,
    missing_code     = NULL,   # alias for miss_codes (scalar or vector)
    K_segments       = NULL,
    degen            = FALSE,
    mcarlo           = TRUE,
    mc_iter          = 25000L,
    mc_target        = 0.05,
    mc_stop          = 99.9,
    mc_stopup        = NA,
    mc_adjust        = FALSE,
    mc_seed          = NULL,
    loo              = c("off","on"),
    boundary_mode    = c("megaoda_halfopen","right_closed"),
    loo_opts         = list(),
    direction        = "off",
    direction_map    = NULL
) {
  attr_type     <- match.arg(attr_type)
  loo           <- match.arg(loo)
  boundary_mode <- match.arg(boundary_mode)

  if (!is.null(missing_code))
    miss_codes <- unique(c(miss_codes, as.numeric(missing_code)))

  .loo_defaults <- list(
    cut_value_mode = "midpoint",
    grid_mode      = "fixed",
    use_samplerep  = FALSE,
    return_folds   = FALSE,
    priors_mode    = "fold"
  )
  loo_opts <- modifyList(.loo_defaults, loo_opts)
  stopifnot(loo_opts$cut_value_mode %in% c("midpoint","lower","upper"))
  stopifnot(loo_opts$grid_mode      %in% c("fixed","refit"))
  stopifnot(loo_opts$priors_mode    %in% c("fold","global"))

  if (isTRUE(degen) && isTRUE(priors_on)) {
    warning("DEGEN=TRUE forces PRIORS=FALSE.")
    priors_on <- FALSE
  }

  if (is.null(w)) w <- rep(1, length(y))
  clean  <- oda_clean_xy(x, y, w, miss_codes)
  x      <- clean$x; y_raw <- clean$y; w_case <- as.numeric(clean$w)
  n_eff  <- length(y_raw)
  if (n_eff == 0L)
    return(list(ok = FALSE, reason = "all_missing", type = "leaf"))

  y <- as.integer(y_raw)

  # Stable 1..C encoding
  u <- sort(unique(y))
  if (!all(u == seq_along(u))) {
    map <- setNames(seq_along(u), as.character(u))
    y   <- as.integer(map[as.character(y)])
    u   <- seq_along(u)
    # Recode direction_map values through same map
    if (!is.null(direction_map)) {
      direction_map <- stats::setNames(
        as.integer(map[as.character(direction_map)]),
        names(direction_map)
      )
    }
  }
  C <- length(u)
  if (C < 3L)
    return(list(ok = FALSE, reason = "not_multiclass", type = "leaf"))

  attr_type_res <- oda_resolve_attr_type(x, attr_type)
  if (attr_type_res == "binary") attr_type_res <- "categorical"

  pol              <- oda_enforce_weighting_policy(attr_type_res, priors_on, loo, w_case)
  priors_on_eff    <- pol$priors_on_eff
  loo_eff          <- pol$loo_eff
  weighted_req     <- pol$weighted_requested
  policy_notes     <- pol$notes

  w_obj  <- oda_apply_priors_multiclass(y, w_case, priors_on = priors_on_eff)
  w_eval <- if (isTRUE(priors_on_eff)) w_obj else w_case

  # ---- ORDERED ----
  if (attr_type_res == "ordered") {
    blocks_obj <- oda_make_blocks_ordered_multiclass(x, y, w_obj)
    blocks_raw <- oda_make_blocks_ordered_multiclass(x, y, w_case)
    x_rep      <- blocks_obj$x_rep
    counts_obj <- blocks_obj$counts
    counts_raw <- blocks_raw$counts
    m          <- nrow(counts_obj)
    if (m < 2L)
      return(list(ok = FALSE, reason = "no_blocks", type = "leaf"))

    K <- as.integer(K_segments %||% C)
    if (K > m) K <- m; if (K < 2L) K <- 2L
    if (!isTRUE(degen) && K < C)
      return(list(ok = FALSE, reason = "no_solution_degen_off_K_lt_C", type = "leaf"))

    sel <- oda_best_ordered_multiclass_partition(
      x_rep = x_rep, counts_obj = counts_obj, counts_raw = counts_raw,
      K = K, priors_on_eff = priors_on_eff, degen = degen,
      primary = NULL, secondary = NULL,
      cut_value_mode = loo_opts$cut_value_mode,
      debug_return_ties = isTRUE(loo_opts$return_folds), debug_max_ties = 200L,
      direction = direction
    )

    if (!isTRUE(sel$ok))
      return(list(ok = FALSE, reason = sel$reason %||% "no_solution", type = "leaf"))

    seg_classes <- blocks_obj$classes[sel$seg_cls_idx]
    if (!isTRUE(degen) && !oda_is_nondegenerate_labels(seg_classes, C))
      return(list(ok = FALSE, reason = "degenerate_solution_disallowed", type = "leaf"))

    rule <- list(type = "multiclass_ordered",
                 cut_values  = sel$cut_values,
                 seg_classes = seg_classes,
                 K = K, boundary = boundary_mode)

    y_pred    <- oda_rule_predict_multiclass(x, rule, boundary = boundary_mode)
    conf_raw  <- oda_confusion_multiclass(y, y_pred, w = rep(1, n_eff))   # raw counts
    conf_wt   <- oda_confusion_multiclass(y, y_pred, w_eval)               # weighted

    out <- list(
      ok = TRUE, reason = NULL, type = rule$type, rule = rule,
      attr_type = attr_type_res, n_eff = n_eff,
      weights_case = w_case, weights_obj = w_obj,
      priors_on_eff = priors_on_eff, policy_notes = policy_notes,
      classes = conf_raw$classes,
      # confusion: raw integer counts (rows = actual, cols = predicted)
      confusion      = conf_raw$confusion,
      confusion_wt   = conf_wt$confusion,
      correct        = conf_raw$correct,             # raw count
      accuracy       = conf_raw$overall_acc,          # raw overall PAC, proportion
      pac            = conf_raw$overall_acc * 100,    # raw overall PAC, percent
      pac_by_class   = conf_raw$pac_by_class * 100,  # raw per-class PAC, percent
      mean_pac       = conf_raw$mean_pac * 100,       # raw mean PAC, percent
      pv_by_class    = conf_raw$pv_by_class * 100,
      mean_pv        = conf_raw$mean_pv * 100,
      ess     = oda_ess_from_mean(conf_raw$mean_pac, C),   # public ESS slot (Patch 1)
      ess_pac = oda_ess_from_mean(conf_raw$mean_pac, C),   # compat alias; retire in Patch 5
      ess_pv  = oda_ess_from_mean(conf_raw$mean_pv,  C),
      p_mc = NA_real_, mc_info = NULL, loo = NULL
    )

    if (isTRUE(mcarlo)) {
      mc     <- oda_mc_p_value_multiclass(
        x = x, y = y, w = w_eval, attr_type = attr_type_res,
        priors_on = priors_on_eff, degen = degen, K_segments = K,
        mc_iter = mc_iter, mc_target = mc_target, mc_stop = mc_stop,
        mc_stopup = mc_stopup, mc_adjust = mc_adjust, seed = mc_seed,
        observed_mean_pac = out$mean_pac,
        direction = direction, direction_map = NULL)
      out$p_mc    <- mc$p_mc
      out$mc_info <- mc
    }

    if (loo_eff == "on") {
      out$loo <- oda_loo_multiclass_ordered(
        x = x, y = y, w0 = w_case,
        priors_on_eff = priors_on_eff, degen = degen,
        K_segments = K, miss_codes = miss_codes,
        cut_value_mode = loo_opts$cut_value_mode, grid_mode = loo_opts$grid_mode,
        boundary_mode = boundary_mode,
        loo_use_samplerep = loo_opts$use_samplerep,
        loo_return_folds = loo_opts$return_folds,
        loo_priors_mode = loo_opts$priors_mode)
    }

    return(out)
  }

  # ---- CATEGORICAL ----
  if (attr_type_res == "categorical") {
    xf   <- factor(x); levs <- levels(xf); L <- length(levs)

    if (!isTRUE(degen) && C > L)
      return(list(ok = FALSE,
                  reason = "error68_classes_exceed_attribute_values",
                  type = "leaf"))

    counts_obj <- matrix(0, nrow = L, ncol = C)
    counts_raw <- matrix(0, nrow = L, ncol = C)
    for (i in seq_along(xf)) {
      li <- as.integer(xf[i]); ci <- y[i]
      counts_obj[li, ci] <- counts_obj[li, ci] + w_obj[i]
      counts_raw[li, ci] <- counts_raw[li, ci] + w_case[i]
    }

    # --- MPE Chapter 4 DIRECTIONAL: auto-create direction_map when L == C ---
    if (direction %in% c("ascending", "descending") && is.null(direction_map)) {
      if (L != C)
        stop("direction = '", direction, "' for categorical multiclass requires ",
             "L == C (one attribute level per class) or an explicit direction_map. ",
             "L = ", L, ", C = ", C, ". Supply direction_map for a fixed-partition ",
             "DIRECTIONAL hypothesis.", call. = FALSE)
      direction_map <- stats::setNames(
        if (direction == "ascending") seq_len(C) else rev(seq_len(C)),
        levs
      )
    }

    # --- Determine level_class: fixed map or optimal search ---
    if (!is.null(direction_map)) {
      # Fixed-partition directional: validate and extract level_class directly
      dm_names <- as.character(names(direction_map))
      dm_vals  <- as.integer(direction_map)
      if (!setequal(dm_names, levs))
        return(list(ok = FALSE, reason = "direction_map_levels_mismatch", type = "leaf"))
      if (!all(dm_vals %in% seq_len(C)))
        return(list(ok = FALSE, reason = "direction_map_values_not_class_labels", type = "leaf"))
      level_class <- dm_vals[match(levs, dm_names)]
      if (!isTRUE(degen) && length(unique(level_class)) < C)
        return(list(ok = FALSE, reason = "degenerate_solution_disallowed", type = "leaf"))
    } else {
      # Nondirectional: exhaustive or greedy search
      NA_raw     <- colSums(counts_raw)
      best_primary <- -Inf; best_overall <- -Inf; best_mean <- -Inf
      best_map   <- NULL;   best_srep    <- Inf

      use_exhaustive <- (L <= 9L)
      if (use_exhaustive) {
        all_maps <- as.matrix(expand.grid(rep(list(seq_len(C)), L)))
        for (r in seq_len(nrow(all_maps))) {
          map_r <- as.integer(all_maps[r, ])
          if (!isTRUE(degen) && length(unique(map_r)) < C) next

          NP_obj <- matrix(0, C, C); NP_raw <- matrix(0, C, C)
          for (l in seq_len(L)) {
            pc <- map_r[l]
            NP_obj[, pc] <- NP_obj[, pc] + counts_obj[l, ]
            NP_raw[, pc] <- NP_raw[, pc] + counts_raw[l, ]
          }

          TP_obj     <- diag(NP_obj); NA_o <- rowSums(NP_obj)
          mp_obj     <- mean(ifelse(NA_o > 0, TP_obj / NA_o, NA_real_), na.rm = TRUE)
          tot_obj    <- sum(NP_obj)
          ov_obj     <- if (tot_obj > 0) sum(TP_obj) / tot_obj else NA_real_
          if (!is.finite(mp_obj)) mp_obj <- -Inf
          if (!is.finite(ov_obj)) ov_obj <- -Inf

          prim <- if (isTRUE(priors_on_eff)) ov_obj else mp_obj
          sec  <- if (isTRUE(priors_on_eff)) mp_obj else ov_obj

          NP_raw_pred <- colSums(NP_raw)
          dp <- sum(NP_raw_pred); da <- sum(NA_raw)
          sr <- if (dp <= 0 || da <= 0) Inf else
            sum(abs((NP_raw_pred / dp) - (NA_raw / da)))

          better <- FALSE
          if      (prim > best_primary + 1e-12) better <- TRUE
          else if (abs(prim - best_primary) <= 1e-12) {
            cur_sec <- if (isTRUE(priors_on_eff)) best_mean else best_overall
            if      (sec > cur_sec + 1e-12) better <- TRUE
            else if (abs(sec - cur_sec) <= 1e-12 && sr < best_srep - 1e-12) better <- TRUE
          }

          if (better) {
            best_primary <- prim; best_overall <- ov_obj
            best_mean    <- mp_obj; best_map <- map_r; best_srep <- sr
          }
        }
      }

      level_class <- if (!is.null(best_map)) best_map else {
        lc <- apply(counts_obj, 1, which.max)
        if (!isTRUE(degen) && length(unique(lc)) < C && L >= C) {
          present <- sort(unique(lc)); missing <- setdiff(seq_len(C), present)
          for (mc_cl in missing) {
            losses <- vapply(seq_len(L), function(l) {
              cur <- lc[l]; if (cur == mc_cl) return(Inf)
              counts_obj[l, cur] - counts_obj[l, mc_cl]
            }, numeric(1))
            l_star <- which.min(losses)
            if (length(l_star) && is.finite(losses[l_star])) lc[l_star] <- mc_cl
          }
        }
        if (!isTRUE(degen) && length(unique(lc)) < C)
          return(list(ok=FALSE, reason="degenerate_solution_disallowed", type="leaf"))
        lc
      }
    }

    overall_class <- which.max(colSums(counts_obj))
    rule <- list(type = "multiclass_nominal",
                 levels = levs,
                 level_class   = as.integer(level_class),
                 default_class = as.integer(overall_class))

    y_pred    <- oda_rule_predict_multiclass(x, rule)
    conf_raw  <- oda_confusion_multiclass(y, y_pred, w = rep(1, n_eff))
    conf_wt   <- oda_confusion_multiclass(y, y_pred, w_eval)

    out <- list(
      ok = TRUE, reason = NULL, type = rule$type, rule = rule,
      attr_type = attr_type_res, n_eff = n_eff,
      weights_case = w_case, weights_obj = w_obj,
      priors_on_eff = priors_on_eff, policy_notes = policy_notes,
      classes = conf_raw$classes,
      confusion      = conf_raw$confusion,
      confusion_wt   = conf_wt$confusion,
      correct        = conf_raw$correct,
      accuracy       = conf_raw$overall_acc,
      pac            = conf_raw$overall_acc * 100,
      pac_by_class   = conf_raw$pac_by_class * 100,
      mean_pac       = conf_raw$mean_pac * 100,
      pv_by_class    = conf_raw$pv_by_class * 100,
      mean_pv        = conf_raw$mean_pv * 100,
      ess     = oda_ess_from_mean(conf_raw$mean_pac, C),   # public ESS slot (Patch 1)
      ess_pac = oda_ess_from_mean(conf_raw$mean_pac, C),   # compat alias; retire in Patch 5
      ess_pv  = oda_ess_from_mean(conf_raw$mean_pv,  C),
      p_mc = NA_real_, mc_info = NULL, loo = NULL
    )

    if (isTRUE(mcarlo)) {
      mc     <- oda_mc_p_value_multiclass(
        x = x, y = y, w = w_eval, attr_type = attr_type_res,
        priors_on = priors_on_eff, degen = degen, K_segments = 1L,
        mc_iter = mc_iter, mc_target = mc_target, mc_stop = mc_stop,
        mc_stopup = mc_stopup, mc_adjust = mc_adjust, seed = mc_seed,
        observed_mean_pac = out$mean_pac,
        direction = direction, direction_map = direction_map)
      out$p_mc    <- mc$p_mc
      out$mc_info <- mc
    }

    if (loo_eff == "on")
      out$loo <- oda_loo_multiclass(
        x = x, y = y, w0 = w_case,
        attr_type = attr_type_res, priors_on_eff = priors_on_eff,
        weighted_requested = weighted_req, degen = degen,
        K_segments = 1L, miss_codes = miss_codes)

    return(out)
  }

  list(ok = FALSE, reason = paste0("unsupported_attr_type:", attr_type_res),
       type = "leaf")
}
