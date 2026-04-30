###############################################################################
# R/unioda_core.R
#
# UniODA Core – univariate binary-class ODA engine.
#
# Spec (MegaODA-faithful):
#   PRIMARY   = MAXSENS  (if priors_on)  / MEANSENS (if priors_off)
#   SECONDARY = SAMPLEREP
#   TERTIARY  = FIRST IDENTIFIED  (enumeration order)
#   LOO       = true refit-per-fold (no global rule reuse)
#   MC        = Fisher randomization with Clopper-Pearson STOP/STOPUP
#
# Exported entry points:
#   oda_univariate_core()   – fit one attribute
#   oda_loo_for_rule()      – LOO given a fixed rule (UniODA)
#   oda_mc_p_value()        – MC p-value (called internally; also exported)
#   oda_confusion_binary()  – confusion table
#   oda_rule_predict()      – apply a UniODA rule to new x
#   oda_mean_pac()          – (sens + spec) / 2
#   oda_ess_from_meanpac()  – ESS from mean PAC
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

  z <- numeric(m); o <- numeric(m); x_rep <- numeric(m)
  for (j in seq_len(m)) {
    idx   <- which(block_id == j)
    z[j]  <- sum(w_ord[idx][y_ord[idx] == 0])
    o[j]  <- sum(w_ord[idx][y_ord[idx] == 1])
    x_rep[j] <- x_ord[idx[1L]]
  }
  list(ord = ord, z = z, o = o, x_rep = x_rep)
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
#   ≡ minimise |predicted_freq_0 - observed_freq_0|  for binary class.
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
#' @param mc_stopup  Confidence level for upper-tail stop (e.g. 20 → 0.20).
#' @param mc_adjust  Kept for API compatibility; not used.
#' @param seed  Optional RNG seed.
#' @param ess_obs  Observed ESS (must be supplied).
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
    mc_stopup    = 20,
    mc_adjust    = FALSE,
    seed         = NULL,
    ess_obs      = NULL
) {
  chance_model <- match.arg(chance_model)
  y <- as.integer(y)
  n <- length(y)
  if (is.null(w)) w <- rep(1, n) else w <- as.numeric(w)
  if (!is.null(seed)) set.seed(seed)

  if (is.null(ess_obs) || !is.finite(ess_obs))
    stop("oda_mc_p_value: ess_obs must be supplied and finite.")

  # Convert stop confidence to probability scale.
  # STOPUP uses the same CUTOFF threshold as STOP (mc_target), not mc_stopup.
  # The canonical MegaODA command is: MC ITER N CUTOFF p STOP c
  # There is no separate STOPUP threshold — both directions compare against CUTOFF.
  conf_level   <- if (!is.na(mc_stop) && mc_stop > 1) mc_stop / 100 else mc_stop

  ge_count  <- 0L
  iter_used <- 0L
  min_check   <- 50L
  check_every <- 50L

  for (b in seq_len(mc_iter)) {
    iter_used <- b
    y_star    <- sample(y, size = n, replace = FALSE)

    fit_b <- oda_univariate_core(
      x = x, y = y_star, w = w,
      attr_type    = attr_type,
      priors_on    = priors_on,
      primary      = primary,
      secondary    = secondary,
      miss_codes   = miss_codes,
      loo          = "off",
      mcarlo       = FALSE,
      mc_iter      = 0L,
      mc_target    = mc_target,
      mc_stop      = mc_stop,
      mc_stopup    = mc_stopup,
      mc_adjust    = mc_adjust,
      mc_seed      = NULL,
      chance_model = chance_model
    )

    ess_b <- if (!isTRUE(fit_b$ok)) 0 else fit_b$ess
    if (ess_b >= ess_obs - 1e-12) ge_count <- ge_count + 1L

    # Clopper-Pearson early stopping
    if (!is.na(conf_level) && b >= min_check && (b %% check_every == 0L)) {
      if (ge_count == 0L) {
        lower <- 0; upper <- stats::qbeta(conf_level, 1, b)
      } else if (ge_count == b) {
        lower <- stats::qbeta(1 - conf_level, b, 1); upper <- 1
      } else {
        upper <- stats::qbeta(conf_level, ge_count + 1, b - ge_count)
        lower <- stats::qbeta(1 - conf_level, ge_count, b - ge_count + 1)
      }
      if (!is.na(mc_target) && upper < mc_target) break   # STOP:   99.9% sure p < CUTOFF
      if (!is.na(mc_target) && lower > mc_target) break   # STOPUP: 99.9% sure p > CUTOFF
    }
  }

  list(p_mc      = (ge_count + 1) / (iter_used + 1),
       ge_count  = ge_count,
       iter_used = iter_used,
       ess_obs   = ess_obs)
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
#' @return List with allowed, confusion, ess_loo, p_value.
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
    mc_seed    = NULL
) {
  chance_model <- match.arg(chance_model)
  attr_type    <- match.arg(attr_type)

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

  y_pred_loo <- integer(n)

  for (i in seq_len(n)) {
    keep <- seq_len(n) != i

    # MegaODA-faithful LOO STABLE: when the training rule is still valid on the
    # n-1 fold, apply it directly to predict the held-out case instead of
    # refitting.  This produces WESSL = WESS (LOO ESS = training ESS) for stable
    # nodes, matching every STABLE node in MegaODA output.
    #
    # For ordered cuts: admissible iff both classes still appear on their
    # correctly-predicted side of the cut in the n-1 data.  If inadmissible
    # (one class is entirely removed from one side), fall through to a full
    # refit — the resulting different rule may yield a different prediction for
    # obs i, signalling genuine instability.
    #
    # For binary maps: only one split is possible; the training rule is always
    # the optimal rule as long as both binary values remain in the n-1 data.
    # Apply the training rule directly (avoids direction-flip artefacts from
    # case-weighted LOO refits).

    if (identical(rule$type, "ordered_cut")) {
      x_k   <- x[keep]
      y_k   <- y[keep]
      w_k   <- w[keep]
      w_k_a <- oda_apply_priors(y_k, w_k, priors_on)

      left  <- !is.na(x_k) & x_k <= rule$cut_value
      right <- !is.na(x_k) & x_k >  rule$cut_value

      admissible <- if (rule$direction == "0->1") {
        sum(w_k_a[left  & y_k == 0L]) > 0 &&
        sum(w_k_a[right & y_k == 1L]) > 0
      } else {                              # "1->0"
        sum(w_k_a[left  & y_k == 1L]) > 0 &&
        sum(w_k_a[right & y_k == 0L]) > 0
      }

      if (admissible) {
        y_pred_loo[i] <- oda_rule_predict(x[i], rule)
        next
      }
      # Fall through to full refit when training cut is inadmissible on n-1 obs

    } else if (identical(rule$type, "binary_map")) {
      # Binary attribute: only one possible split.  Apply the training rule
      # directly as long as obs i has a non-missing value.
      if (!is.na(x[i])) {
        y_pred_loo[i] <- oda_rule_predict(x[i], rule)
        next
      }
    }

    fit_i <- oda_univariate_core(
      x = x[keep], y = y[keep], w = w[keep],
      attr_type    = attr_type,
      priors_on    = priors_on,
      primary      = primary,
      secondary    = secondary,
      miss_codes   = miss_codes,
      loo          = "off",
      mcarlo       = FALSE,
      mc_iter      = mc_iter,
      mc_target    = mc_target,
      mc_stop      = mc_stop,
      mc_stopup    = mc_stopup,
      mc_adjust    = mc_adjust,
      mc_seed      = if (is.null(mc_seed)) NULL else mc_seed + i,
      chance_model = chance_model
    )

    if (!isTRUE(fit_i$ok)) {
      return(list(allowed = FALSE,
                  reason  = paste0("loo_fit_failed_case_", i,
                                   "_reason_", fit_i$reason),
                  confusion = NULL, ess_loo = NA_real_, p_value = NA_real_))
    }
    y_pred_loo[i] <- oda_rule_predict(x[i], fit_i$rule)
  }

  # ESS for LOO — use the same priors+case-weight scheme as training so that
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

  # Fisher's exact test on raw (unit-weight) counts — requires integer inputs.
  conf_loo_r <- oda_confusion_binary(y, y_pred_loo)   # unit weights
  tab <- matrix(c(conf_loo_r$TP, conf_loo_r$FP,
                  conf_loo_r$FN, conf_loo_r$TN),
                nrow = 2, byrow = TRUE)
  p_fisher <- tryCatch(
    stats::fisher.test(tab, alternative = "greater")$p.value,
    error = function(e) NA_real_
  )

  list(allowed = TRUE, reason = NULL,
       confusion = conf_loo, ess_loo = ess_loo, p_value = p_fisher)
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
    chance_model = c("class","attribute")
) {
  chance_model <- match.arg(chance_model)
  loo          <- match.arg(loo)

  # Resolve missing_code alias → miss_codes
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

  # 5. Chance baselines
  C_class     <- 2L
  chance_class <- 1 / C_class
  k_attr      <- NA_integer_

  # 6. Enumerate admissible rules -------------------------------------------

  cand_rows  <- list()
  preds_list <- list()   # index matches cand_rows for SAMPLEREP

  add_candidate <- function(rule, j_index = NA_integer_) {
    y_pred    <- oda_rule_predict(x, rule)
    conf      <- oda_confusion_binary(y, y_pred, w)
    mp        <- conf$mean_pac
    pac       <- mp * 100
    chance_a  <- 1 / k_attr
    ess_class <- oda_ess_from_meanpac(mp, chance_class)
    ess_attr  <- oda_ess_from_meanpac(mp, chance_a)
    ess_obj   <- if (chance_model == "class") ess_class else ess_attr

    row <- data.frame(
      row_id    = length(preds_list) + 1L,
      j         = j_index,
      direction = rule$direction,
      cut_value = if (!is.null(rule$cut_value)) rule$cut_value else NA_real_,
      pac       = pac,
      ess_class = ess_class,
      ess_attr  = ess_attr,
      ess       = ess_obj,
      sens_1    = conf$sensitivity,
      spec_0    = conf$specificity,
      mean_pac  = mp,
      balance   = abs(conf$sensitivity - conf$specificity),
      stringsAsFactors = FALSE
    )
    cand_rows[[length(cand_rows) + 1L]] <<- row
    preds_list[[length(preds_list) + 1L]] <<- y_pred
  }

  # --- binary attribute ---
  if (attr_type == "binary") {
    levs <- sort(unique(as.character(x)))
    if (length(levs) != 2L)
      return(list(ok = FALSE, reason = "not_binary", type = "leaf"))
    k_attr <- 2L
    rules  <- list(
      list(type="binary_map", direction="0->1", left_levels=levs[1],right_levels=levs[2]),
      list(type="binary_map", direction="1->0", left_levels=levs[1],right_levels=levs[2])
    )
    for (r in rules) add_candidate(r, j_index = NA_integer_)
  }

  # --- categorical attribute ---
  if (attr_type == "categorical") {
    x_fac <- factor(x)
    levs   <- levels(x_fac)
    k_attr <- length(levs)
    if (k_attr < 2L) return(list(ok = FALSE, reason = "no_levels", type = "leaf"))

    # Order by weighted class-1 rate
    p1 <- sapply(levs, function(l) {
      wl <- w[x_fac == l]; yl <- y[x_fac == l]
      if (sum(wl) == 0) NA_real_ else sum(wl * yl) / sum(wl)
    })
    levs_ord <- levs[order(p1, na.last = TRUE)]

    for (j in seq_len(k_attr - 1L)) {
      left  <- levs_ord[seq_len(j)]
      right <- levs_ord[(j + 1L):k_attr]
      for (dir in c("0->1","1->0")) {
        rule <- list(type="nominal_cut", direction=dir,
                     left_levels=left, right_levels=right)
        add_candidate(rule, j_index = j)
      }
    }
  }

  # --- ordered attribute ---
  if (attr_type == "ordered") {
    blocks <- oda_make_blocks_ordered(x, y, w)
    z      <- blocks$z; o <- blocks$o; x_rep <- blocks$x_rep
    m      <- length(z)
    if (m < 2L) return(list(ok = FALSE, reason = "no_blocks", type = "leaf"))
    k_attr <- m

    Z  <- cumsum(z); O  <- cumsum(o)
    Zm <- sum(z);    Om <- sum(o)

    admissible <- function(j, dir) {
      if (dir == "0->1") return(Z[j] > 0 && (Om - O[j]) > 0)
      if (dir == "1->0") return(O[j] > 0 && (Zm - Z[j]) > 0)
      FALSE
    }

    for (j in seq_len(m - 1L)) {
      for (dir in c("0->1","1->0")) {
        if (!admissible(j, dir)) next
        cut_v <- (x_rep[j] + x_rep[j + 1L]) / 2
        rule  <- list(type="ordered_cut", direction=dir, cut_value=cut_v)
        add_candidate(rule, j_index = j)
      }
    }
  }

  if (length(cand_rows) == 0L)
    return(list(ok = FALSE, reason = "no_admissible_cut", type = "leaf"))

  cand_df <- do.call(rbind, cand_rows)

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

  # 11. Monte Carlo p-value
  mc_res <- NULL
  p_mc   <- NA_real_
  if (isTRUE(mcarlo)) {
    mc_res <- oda_mc_p_value(
      x = x, y = y, w = w_raw,  # raw weights; MC loop re-applies priors per permutation
      attr_type    = attr_type,
      priors_on    = priors_on,
      primary      = primary,
      secondary    = secondary,
      miss_codes   = miss_codes,
      chance_model = chance_model,
      mc_iter      = mc_iter,
      mc_target    = mc_target,
      mc_stop      = mc_stop,
      mc_stopup    = mc_stopup,
      mc_adjust    = mc_adjust,
      seed         = mc_seed,
      ess_obs      = ess_obj
    )
    p_mc <- mc_res$p_mc
  }

  # 12. LOO
  loo_out <- NULL
  if (loo != "off") {
    loo_out <- oda_loo_for_rule(
      x = x, y = y, rule = rule_best, w = w_raw,
      chance_model = chance_model,
      k_attr       = k_attr,
      attr_type    = attr_type,
      priors_on    = priors_on,
      primary      = primary,
      secondary    = secondary,
      miss_codes   = miss_codes,
      mc_iter      = mc_iter,
      mc_target    = mc_target,
      mc_stop      = mc_stop,
      mc_stopup    = mc_stopup,
      mc_adjust    = mc_adjust,
      mc_seed      = if (is.null(mc_seed)) NULL else mc_seed + 1L
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
