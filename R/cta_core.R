###############################################################################
# R/cta_core.R - Classification Tree Analysis (CTA)
#
# Implements EO-CTA (Enumerated Optimal CTA), corresponding to MegaODA's
# ENUMERATE command.
#
# Algorithm (Yarnold & Soltysik 2016, Chapters 10-11):
#
#   ENUMERATE: Evaluate every valid (root A x left-child B x right-child C)
#     combination.  Each candidate split must pass MC significance AND LOO
#     STABLE locally.  Below the top-3 nodes, grow greedily (HO-CTA).
#     Compute full-tree WESS over all n observations.
#     Retain the combination with the highest full-tree WESS.
#
#   HO-CTA (nodes 4+):  At each node, take the attribute with the HIGHEST
#     local ESS that passes MC significance AND LOO STABLE.  Recurse.
#
# Implementation note on B-subtree caching:
#   B's subtree (children of node 2) does not depend on C.  It is grown
#   once per (A, B) pair; then reused for every C option.  C's subtree
#   node IDs are offset past B's last used node ID to avoid collisions.
###############################################################################

# ---- Internal diagnostic helpers ------------------------------------------ #

# .row_hash(idx)
# Lightweight fingerprint of a row-index set for diagnostic cache audits.
# Not cryptographic; designed to detect exact repeats in ENUMERATE sub-nodes.
# idx: integer vector of row positions (order-independent).
.row_hash <- function(idx) {
  idx <- sort.int(as.integer(idx))
  n   <- length(idx)
  paste0(n, ".", sum(idx), ".", (sum(idx * seq_len(n)) %% 999983L))
}

# .cta_predictions_degenerate(preds, required_classes)
#
# Returns TRUE when a prediction vector does not contain at least
# `required_classes` distinct non-NA classes.  Used in two places:
#   1. Expanded ENUMERATE post-prune gate: candidate is skipped.
#   2. Root-only stump defensive guard.
.cta_predictions_degenerate <- function(preds, required_classes = 2L) {
  vals <- unique(preds[!is.na(preds)])
  length(vals) < required_classes
}

# ---- CTA-specific ordered-cut selection ------------------------------------ #

# .cta_ordered_scan()
#
# CTA-faithful ordered-cut selection for binary class.
#
# The rule differs from generic oda_univariate_core():
#   "0->1" direction: select the RIGHTMOST cut where the class-1
#     priors-adjusted PAC on the right branch (sens_wa) exceeds 0.5.
#     Since sens_wa_01 = 1 - C1a[j] is non-increasing, this is the last j
#     before sens_wa drops to <= 0.5 - equivalently, maximise spec_wa subject
#     to sens_wa > 0.5.
#   "1->0" direction: select the LEFTMOST cut where the class-1
#     priors-adjusted PAC on the LEFT branch (sens_wa_10 = C1a[j]) exceeds 0.5.
#     Since C1a[j] is non-decreasing, this is the first j where sens_wa_10
#     crosses above 0.5 - equivalently, maximise spec_wa_10 subject to
#     sens_wa_10 > 0.5.
#   The two rules are symmetric: each picks the "outermost" cut that maintains
#   class-1 majority on the prediction side.
#   The direction yielding the higher ESS is returned (ties -> "0->1").
#
# Parameters
#   x, y, w    raw data vectors (same length); y must be in {0, 1}.
#   priors_on  if TRUE apply prior-odds weighting.
#   miss_codes numeric miss-code vector (may be NULL).
#   mindenom   minimum raw obs count per child node.
#
# Returns list(rule, ess, sens_wa, spec_wa) or NULL if no eligible cut exists.
.cta_ordered_scan <- function(x, y, w, priors_on, miss_codes, mindenom) {
  # 1. Clean missing
  miss_mask <- is.na(x) | is.na(y)
  if (!is.null(miss_codes)) miss_mask <- miss_mask | (x %in% miss_codes)
  keep <- !miss_mask
  x <- x[keep]; y <- as.integer(y[keep]); w <- as.numeric(w[keep])
  n <- length(x)
  if (n < 2L * as.integer(mindenom)) return(NULL)
  n0 <- sum(y == 0L); n1 <- sum(y == 1L)
  if (n0 == 0L || n1 == 0L) return(NULL)

  # 2. Priors-adjusted weights
  W0 <- sum(w[y == 0L]); W1 <- sum(w[y == 1L])
  if (W0 <= 0 || W1 <= 0) return(NULL)
  w_adj <- w
  if (isTRUE(priors_on)) {
    w_adj[y == 0L] <- w_adj[y == 0L] / W0
    w_adj[y == 1L] <- w_adj[y == 1L] / W1
  }

  # 3. Unique-value blocks (sorted left to right)
  uvals <- sort(unique(x))
  m     <- length(uvals)
  if (m < 2L) return(NULL)

  z0a <- numeric(m); z1a <- numeric(m); nv <- integer(m)
  for (b in seq_len(m)) {
    idx    <- x == uvals[b]
    z0a[b] <- sum(w_adj[idx][y[idx] == 0L])
    z1a[b] <- sum(w_adj[idx][y[idx] == 1L])
    nv[b]  <- sum(idx)
  }
  C0a <- cumsum(z0a); C1a <- cumsum(z1a); Nv <- cumsum(nv); Nm <- sum(nv)

  # 4. Direction "0->1": sens_wa = 1 - C1a[j] (non-increasing).
  #    Condition: sens_wa > 0.5 <=> C1a[j] < 0.5.
  #    Rule: RIGHTMOST eligible j -> overwrite on every valid j.
  bj01 <- NA_integer_; bes01 <- NA_real_; bse01 <- NA_real_; bsp01 <- NA_real_
  for (j in seq_len(m - 1L)) {
    nL <- Nv[j]; nR <- Nm - nL
    if (nL < mindenom || nR < mindenom) next
    sens <- 1 - C1a[j]
    if (sens <= 0.5) next
    spec <- C0a[j]                        # non-decreasing: rightmost = max
    bj01 <- j; bes01 <- (spec + sens - 1) * 100; bse01 <- sens; bsp01 <- spec
  }

  # 5. Direction "1->0": sens_wa_10 = C1a[j] (non-decreasing).
  #    Condition: sens_wa_10 > 0.5 <=> C1a[j] > 0.5.
  #    Rule: LEFTMOST eligible j -> break on first valid j.
  bj10 <- NA_integer_; bes10 <- NA_real_; bse10 <- NA_real_; bsp10 <- NA_real_
  for (j in seq_len(m - 1L)) {
    nL <- Nv[j]; nR <- Nm - nL
    if (nL < mindenom || nR < mindenom) next
    sens <- C1a[j]
    if (sens <= 0.5) next
    spec <- 1 - C0a[j]                   # non-increasing: leftmost = max
    bj10 <- j; bes10 <- (spec + sens - 1) * 100; bse10 <- sens; bsp10 <- spec
    break                                 # leftmost = first found
  }

  # 6. Build result for each direction, pick direction with higher ESS (ties -> "0->1")
  .make_r <- function(dir, bj, bes, bse, bsp)
    list(rule = list(type = "ordered_cut", direction = dir,
                     cut_value = (uvals[bj] + uvals[bj + 1L]) / 2),
         ess = bes, sens_wa = bse, spec_wa = bsp)

  r01 <- if (!is.na(bj01)) .make_r("0->1", bj01, bes01, bse01, bsp01) else NULL
  r10 <- if (!is.na(bj10)) .make_r("1->0", bj10, bes10, bse10, bsp10) else NULL

  if (is.null(r01) && is.null(r10)) return(NULL)
  if (is.null(r01)) return(r10)
  if (is.null(r10)) return(r01)
  if (r01$ess >= r10$ess) r01 else r10
}

# .cta_mc_ordered()
#
# CTA-specific Monte Carlo Fisher-randomization p-value for ordered attributes.
# Mirrors oda_mc_p_value() but applies .cta_ordered_scan() to each permutation
# so permuted fits respect the same CTA cut-selection rule as the observed fit.
#
# Parameters
#   x, y_coded   raw vectors (with miss_codes intact); y must be in {0,1}.
#   w            case weights.
#   obs_ess      observed ESS from .cta_ordered_scan() - required.
#   priors_on, miss_codes, mindenom  forwarded to .cta_ordered_scan().
#   mc_iter, mc_target, mc_stop, mc_stopup  Fisher randomization params.
#   mc_seed      optional RNG seed.
#
# Returns list(p_mc, ge_count, iter_used).
.cta_mc_ordered <- function(x, y_coded, w, obs_ess,
                             priors_on, miss_codes, mindenom,
                             mc_iter, mc_target, mc_stop, mc_stopup, mc_seed,
                             diag_env   = NULL,
                             row_hash   = NULL,
                             attr_name  = NA_character_,
                             pos_id     = NA,
                             n_obs      = NA_integer_) {
  if (!is.null(mc_seed)) set.seed(mc_seed)
  y_coded   <- as.integer(y_coded)
  n         <- length(y_coded)
  mc_iter   <- as.integer(mc_iter)

  conf_stop   <- if (!is.na(mc_stop)   && mc_stop   > 1) mc_stop   / 100 else mc_stop
  conf_stopup <- if (!is.na(mc_stopup) && mc_stopup > 1) mc_stopup / 100 else mc_stopup

  ge_count         <- 0L
  iter_used        <- 0L
  min_check        <- 50L
  check_every      <- 50L
  stop_reason_code <- "max_iter"
  t0_mc            <- if (!is.null(diag_env)) proc.time()[["elapsed"]] else NULL

  for (b in seq_len(mc_iter)) {
    iter_used <- b
    y_star    <- sample(y_coded, size = n, replace = FALSE)
    scan_b    <- .cta_ordered_scan(x, y_star, w, priors_on, miss_codes, mindenom)
    ess_b     <- if (is.null(scan_b) || is.null(scan_b$ess)) 0 else scan_b$ess %||% 0
    if (ess_b >= obs_ess - 1e-12) ge_count <- ge_count + 1L

    if (b >= min_check && (b %% check_every == 0L)) {
      if (ge_count == 0L) {
        upper <- if (!is.na(conf_stop))   stats::qbeta(conf_stop,   1L,          b)          else NA_real_
        lower <- 0
      } else if (ge_count == b) {
        upper <- 1
        lower <- if (!is.na(conf_stopup)) stats::qbeta(1 - conf_stopup, b,       1L)         else NA_real_
      } else {
        upper <- if (!is.na(conf_stop))   stats::qbeta(conf_stop,   ge_count + 1L, b - ge_count)     else NA_real_
        lower <- if (!is.na(conf_stopup)) stats::qbeta(1 - conf_stopup, ge_count, b - ge_count + 1L) else NA_real_
      }
      if (!is.na(mc_target) && !is.na(upper) && upper < mc_target) {
        stop_reason_code <- "stop_significant"; break
      }
      if (!is.na(mc_target) && !is.na(lower) && lower > mc_target) {
        stop_reason_code <- "stop_nonsignificant"; break
      }
    }
  }

  if (!is.null(diag_env)) {
    elapsed_mc <- proc.time()[["elapsed"]] - t0_mc
    p_pseudo   <- (ge_count + 1L) / (iter_used + 1L)
    p_raw      <- if (iter_used <= 0L) NA_real_
                  else if (ge_count == 0L) 0.0
                  else ge_count / iter_used
    diag_env$mc_log[[length(diag_env$mc_log) + 1L]] <- list(
      path             = "cta_specific",
      attr_name        = attr_name,
      pos_id           = pos_id,
      row_hash         = row_hash,
      n_obs            = n_obs,
      obs_ess          = obs_ess,
      mc_target        = mc_target,
      ge_count         = ge_count,
      iter_used        = iter_used,
      stop_reason      = stop_reason_code,
      p_mc_current     = p_raw,
      p_mc_raw         = p_raw,
      p_mc_pseudocount = p_pseudo,
      signif_current   = !is.na(p_raw)    && p_raw    < mc_target,
      signif_raw       = !is.na(p_raw)    && p_raw    < mc_target,
      elapsed_sec      = elapsed_mc
    )
  }

  list(p_mc      = if (iter_used <= 0L) NA_real_
                   else if (ge_count == 0L) 0.0
                   else ge_count / iter_used,
       ge_count  = ge_count,
       iter_used = iter_used)
}

# Normalize a binary_map rule on a numeric {0,1} attribute to ordered_cut.
.cta_norm_rule <- function(rule, x_j, miss_codes = NULL) {
  if (!identical(rule$type, "binary_map")) return(rule)
  x_num <- as.numeric(x_j[!is.na(x_j)])
  if (!is.null(miss_codes)) x_num <- x_num[!(x_num %in% miss_codes)]
  uvals <- sort(unique(x_num))
  if (length(uvals) != 2L || !isTRUE(all.equal(uvals, c(0, 1)))) return(rule)
  low_left <- 0 %in% as.numeric(rule$left_levels %||% integer(0))
  list(type = "ordered_cut", cut_value = 0.5,
       direction = if (low_left) "0->1" else "1->0")
}

# Compute weighted ESS from predictions.  Returns ESS in [0, 100].
.wess_classes <- function(y_actual, y_pred, w) {
  classes <- sort(unique(as.integer(y_actual)))
  C_local <- length(classes)
  if (C_local < 2L) return(0)
  pac_sum <- 0
  for (cl in classes) {
    m  <- y_actual == cl
    wt <- sum(w[m])
    if (wt <= 0) next
    pac_sum <- pac_sum + sum(w[m & y_pred == cl]) / wt
  }
  (pac_sum / C_local - 1.0 / C_local) / (1.0 - 1.0 / C_local) * 100
}


#' Fit a Classification Tree Analysis (CTA) model using ENUMERATE
#'
#' Builds a classification tree via EO-CTA (Enumerated Optimal CTA) as
#' described in Yarnold & Soltysik (2016).  All valid combinations of root,
#' left child, and right child are evaluated; below the top three nodes the
#' tree is grown greedily (HO-CTA).  The combination with the highest
#' full-tree WESS is returned.
#'
#' @param X Data frame or matrix of attribute columns.
#' @param y Class variable vector.
#' @param w Optional numeric case weights (MegaODA WEIGHT).
#' @param priors_on Logical; use prior-odds weighting at each node.
#' @param miss_codes Numeric vector of missing-value codes.
#' @param alpha_split p threshold to split (MegaODA MC CUTOFF).
#' @param mindenom Minimum raw observation count to attempt a split.
#' @param prune_alpha Don't add nodes with p >= prune_alpha (MegaODA PRUNE).
#' @param max_depth Maximum depth (root = depth 1).
#' @param ess_min Minimum ESS to split.
#' @param mc_iter Maximum MC iterations per fit.
#' @param mc_target MC significance threshold for early stopping.
#' @param mc_stop Confidence level (pct) for STOP bound.
#' @param mc_stopup Confidence level (pct) for STOPUP bound.
#' @param mc_seed Integer base seed; NULL = unseeded.
#' @param loo LOO mode per node: \code{"off"} (default), \code{"stable"}
#'   (MegaODA LOO STABLE; accept when |WESSL − WESS| ≤ 0.01 pp; reports
#'   \code{loo_status = "STABLE"}), \code{"pvalue"} (Fisher p strictly less
#'   than 0.05; reports \code{loo_status = "PVALUE"}), or a single numeric in
#'   (0, 1) (Fisher p strictly less than the supplied threshold; reports
#'   \code{loo_status = "PVALUE"}).
#' @param attr_names Optional character vector of attribute names.
#' @param K_segments Segments for multiclass ordered splits.
#' @return An object of class \code{cta_tree}.
#' @export
oda_cta_fit <- function(
    X,
    y,
    w           = NULL,
    priors_on   = TRUE,
    miss_codes  = NULL,
    alpha_split = 0.05,
    mindenom    = 5L,
    prune_alpha = 1.0,
    max_depth   = 10L,
    ess_min     = 0,
    mc_iter     = 25000L,
    mc_target   = 0.05,
    mc_stop     = 99.9,
    mc_stopup   = 20,
    mc_seed     = NULL,
    loo         = "off",
    attr_names  = NULL,
    K_segments  = NULL,
    verbose     = FALSE,
    diag_env    = NULL
) {
  # ---- Validate and prep ----------------------------------------------------
  if (!is.null(diag_env)) {
    if (!exists("mc_log",   envir = diag_env, inherits = FALSE)) diag_env$mc_log   <- list()
    if (!exists("loo_log",  envir = diag_env, inherits = FALSE)) diag_env$loo_log  <- list()
    if (!exists("fast_log", envir = diag_env, inherits = FALSE)) diag_env$fast_log <- list()
    diag_env$fit_start <- proc.time()[["elapsed"]]
  }
  if (!is.data.frame(X)) X <- as.data.frame(X)
  n <- nrow(X)
  stopifnot(length(y) == n)
  y <- as.integer(y)
  .validate_case_weights(w, n)
  if (is.null(w)) w <- rep(1.0, n) else w <- as.numeric(w)
  if (is.null(attr_names))
    attr_names <- colnames(X) %||% paste0("V", seq_len(ncol(X)))
  n_attrs <- ncol(X)
  C            <- length(unique(y))
  class_levels <- sort(unique(as.integer(y)))   # fit-level class universe
  loo_arg <- if (is.numeric(loo)) "pvalue" else as.character(loo)

  .vmsg <- if (isTRUE(verbose)) function(...) message(...) else function(...) invisible(NULL)
  .vmsg("[CTA] fit: n=", n, " attrs=", n_attrs,
        " mindenom=", mindenom, " loo=", loo_arg,
        " mc_iter=", mc_iter)

  # ---- Inner helpers (closures over X, y, w, parameters) -------------------

  .get_ess <- function(fit) {
    if (!is.null(fit$ess)     && is.finite(fit$ess))     return(fit$ess)
    if (!is.null(fit$ess_pac) && is.finite(fit$ess_pac)) return(fit$ess_pac)
    NA_real_
  }

  .majority_class <- function(idx) {
    ctab <- tapply(w[idx], y[idx], sum)
    as.integer(names(which.max(ctab)))
  }

  # Leaf node record.
  # pred_class: if provided (not NULL), use it as majority_class (= class
  #   assigned by the parent ODA rule for the branch leading here).
  #   If NULL, fall back to local weighted majority (root-only case).
  .leaf_nd <- function(nid, parent_id, depth, idx, pred_class = NULL) {
    maj <- if (!is.null(pred_class)) as.integer(pred_class)
           else .majority_class(idx)
    yy <- as.integer(y[idx])
    ww <- w[idx]
    class_counts_raw <- setNames(
      as.integer(table(factor(yy, levels = class_levels))),
      as.character(class_levels)
    )
    class_counts_weighted <- setNames(
      vapply(class_levels, function(cl) sum(ww[yy == cl]), numeric(1L)),
      as.character(class_levels)
    )
    list(
      node_id = nid, parent_id = parent_id, depth = depth,
      n_obs = length(idx), n_weighted = sum(ww),
      leaf = TRUE, majority_class = maj,
      class_counts_raw = class_counts_raw,
      class_counts_weighted = class_counts_weighted,
      attribute = NA_character_, attr_col = NA_integer_,
      attr_type = NA_character_, rule = NULL,
      ess = NA_real_, ess_weighted = NA_real_, p_mc = NA_real_,
      loo_status = NA_character_, loo_ess = NA_real_, loo_p = NA_real_,
      confusion = NULL, split_labels = integer(0), child_ids = integer(0)
    )
  }

  # Apply candidate's rule to subset idx.
  # Returns list(rule, engine, y_pred, y_pred_raw):
  #   y_pred     = NA for miss_coded obs (used for child routing)
  #   y_pred_raw = includes miss_coded obs (used for node confusion)
  .apply_cand <- function(cand, idx) {
    engine <- cand$fit$engine %||% "binary"
    x_j    <- X[[cand$j]][idx]
    rule   <- if (engine != "multiclass")
      .cta_norm_rule(cand$fit$rule, x_j, miss_codes)
    else
      cand$fit$rule
    y_pred_raw <- if (engine == "multiclass")
      as.integer(oda_rule_predict_multiclass(
        x_j, rule, boundary = rule$boundary %||% "megaoda_halfopen"))
    else
      as.integer(oda_rule_predict(x_j, rule))
    y_pred <- y_pred_raw
    # Remap y_pred from rule-side 0/1 to actual class labels for binary rules.
    # y_pred_raw stays as 0/1 (for .conf_binary which expects 0/1 convention).
    if (engine != "multiclass") {
      cls <- sort(unique(as.integer(y[idx])))
      if (length(cls) == 2L) y_pred <- cls[y_pred_raw + 1L]
    }
    miss_mask <- is.na(x_j)
    if (!is.null(miss_codes)) miss_mask <- miss_mask | (x_j %in% miss_codes)
    y_pred[miss_mask] <- NA_integer_
    list(rule = rule, engine = engine,
         y_pred = y_pred, y_pred_raw = y_pred_raw)
  }

  # 2x2 binary confusion matrix
  .conf_binary <- function(y_node, y_pred_raw) {
    conf <- oda_confusion_binary(y_node, y_pred_raw)
    matrix(c(conf$TN, conf$FP, conf$FN, conf$TP), nrow = 2L, ncol = 2L,
           dimnames = list(actual = c("0","1"), predicted = c("0","1")))
  }

  # Final full-tree training confusion matrix from actual and predicted labels.
  # y_actual and y_pred must already be filtered to classified observations
  # (NA rows excluded by the caller).
  # levels defaults to class_levels (the fit-level class universe) so that
  # the stored matrix always spans all fitted classes - preventing silent
  # class dropping in edge cases where path-local scoring excludes a class.
  .make_training_conf <- function(y_actual, y_pred, levels = class_levels) {
    y_actual <- as.integer(y_actual)
    y_pred   <- as.integer(y_pred)
    levels   <- as.integer(levels)
    tbl <- table(
      actual    = factor(y_actual, levels = levels),
      predicted = factor(y_pred,   levels = levels)
    )
    matrix(as.integer(tbl), nrow = length(levels), ncol = length(levels),
           dimnames = list(actual    = as.character(levels),
                           predicted = as.character(levels)))
  }

  # LOO metadata from a fit result.
  # status vocabulary:
  #   "OFF"    — LOO mode is "off" or LOO did not run
  #   "STABLE" — candidate passed the LOO STABLE gate (|WESSL - WESS| <= 0.01 pp)
  #   "PVALUE" — candidate passed the LOO p-value gate
  # Uses enclosing loo_arg (closure over oda_cta_fit scope).
  .loo_info <- function(fit) {
    lr <- fit$loo
    status <- if (is.null(lr) || !isTRUE(lr$allowed)) {
      "OFF"
    } else if (identical(loo_arg, "stable")) {
      "STABLE"
    } else {
      "PVALUE"
    }
    list(
      status  = status,
      ess_loo = if (!is.null(lr)) lr$ess_loo %||% NA_real_ else NA_real_,
      p_value = if (!is.null(lr)) lr$p_value %||% NA_real_ else NA_real_
    )
  }

  # CTA-specific candidate evaluation for weighted ordered binary predictors.
  #
  # Canon (docs/CTA_ORDERED_CUT_AUDIT.md, MPE.pdf):
  #   Full WESS: .cta_ordered_scan() - rightmost cut with class-1 right-branch
  #              priors-adjusted PAC (sens_wa) > 0.5.
  #   MC p-value: .cta_mc_ordered() - Fisher randomization on CTA scan.
  #   LOO STABLE: generic ODA per-fold refit (oda_loo_for_rule()).
  #               Canon: WESSL must equal WESS (|WESS - WESSL| <= 0.01 pp).
  #               Signif T alone is insufficient; WESSL = WESS is required.
  #
  # Trigger: non-uniform case weights AND binary class AND >2 unique attr values.
  # Binary attributes (<=2 unique values) use the generic ODA path unchanged.
  # Returns a fit-like list or NULL if rejected (Signif F, no eligible cut, or
  # LOO UNSTABLE). When this path is triggered, callers must not fall back to
  # generic ODA - the candidate is either accepted here or rejected entirely.
  .cta_full_fit_ordered <- function(j, x_j, y_n, w_n, pos_id,
                                     row_hash = NULL, n_obs_parent = NA_integer_) {
    # Guard: uniform weights -> generic ODA applies
    if (all(w_n == w_n[1L])) return(NULL)
    # Guard: binary class (C = 2) only
    y_uniq <- sort(unique(as.integer(y_n[!is.na(y_n)])))
    if (length(y_uniq) != 2L) return(NULL)
    # Guard: attribute must have >2 distinct non-missing values
    x_v <- x_j[!is.na(x_j)]
    if (!is.null(miss_codes)) x_v <- x_v[!(x_v %in% miss_codes)]
    if (length(unique(x_v)) <= 2L) return(NULL)

    # Recode y to {0, 1} using sorted label order
    y_01 <- ifelse(as.integer(y_n) == y_uniq[1L], 0L, 1L)

    # Full-model WESS: CTA-specific ordered scan
    cta_scan <- .cta_ordered_scan(x_j, y_01, w_n,
                                  priors_on  = priors_on,
                                  miss_codes = miss_codes,
                                  mindenom   = mindenom)
    if (is.null(cta_scan) || cta_scan$ess <= 0) return(NULL)
    obs_ess <- cta_scan$ess

    # MC significance: CTA-specific permutation test
    seed_j <- if (is.null(mc_seed)) NULL else
      (as.integer(mc_seed) + as.integer(pos_id) * 100L + j) %% .Machine$integer.max
    mc_res <- .cta_mc_ordered(
      x          = x_j,
      y_coded    = y_01,
      w          = w_n,
      obs_ess    = obs_ess,
      priors_on  = priors_on,
      miss_codes = miss_codes,
      mindenom   = mindenom,
      mc_iter    = as.integer(mc_iter),
      mc_target  = mc_target,
      mc_stop    = mc_stop,
      mc_stopup  = mc_stopup,
      mc_seed    = seed_j,
      diag_env   = diag_env,
      row_hash   = row_hash,
      attr_name  = attr_names[j],
      pos_id     = pos_id,
      n_obs      = n_obs_parent
    )
    p_mc <- mc_res$p_mc
    .vmsg("  [CTA-scan] ", attr_names[j],
          " cut=", round(cta_scan$rule$cut_value, 4),
          " WESS=", round(obs_ess, 4), "%  p=", round(p_mc, 4))
    if (is.na(p_mc) || p_mc >= alpha_split) {
      .vmsg("  -> rejected (CTA path): Signif F  p=", round(p_mc %||% NA_real_, 4))
      return(NULL)
    }

    # LOO STABLE gate (MPE.pdf canon): WESSL must equal WESS.
    # LOO uses generic ODA per-fold refit - true refit-per-fold LOO.
    loo_result <- NULL
    if (!identical(loo_arg, "off")) {
      t0_loo <- if (!is.null(diag_env)) proc.time()[["elapsed"]] else NULL
      loo_result <- tryCatch(
        oda_loo_for_rule(
          x          = x_j,
          y          = y_01,
          w          = w_n,
          rule       = cta_scan$rule,
          attr_type  = "ordered",
          priors_on  = priors_on,
          miss_codes = miss_codes
        ),
        error = function(e) NULL
      )
      if (!is.null(diag_env)) {
        loo_elapsed <- proc.time()[["elapsed"]] - t0_loo
        wessl_d     <- if (!is.null(loo_result)) loo_result$ess_loo %||% NA_real_ else NA_real_
        delta_d     <- abs(obs_ess - (wessl_d %||% Inf))
        loo_stable_flag <- if (identical(loo_arg, "stable"))
          !is.na(wessl_d) && delta_d <= 0.01
        else NA
        diag_env$loo_log[[length(diag_env$loo_log) + 1L]] <- list(
          path        = "cta_specific",
          attr_name   = attr_names[j],
          pos_id      = pos_id,
          row_hash    = row_hash,
          n_obs       = n_obs_parent,
          loo_mode    = if (is.null(loo_result)) "failed" else "n_fold_or_algebraic",
          ess_loo     = wessl_d,
          stable      = loo_stable_flag,
          elapsed_sec = loo_elapsed
        )
      }
      if (identical(loo_arg, "stable")) {
        wessl <- if (!is.null(loo_result)) loo_result$ess_loo %||% NA_real_
                 else NA_real_
        delta <- abs(obs_ess - (wessl %||% Inf))
        .vmsg("  [CTA-scan] WESSL=", round(wessl %||% NA_real_, 4),
              "%  |delta|=", round(delta, 4), " pp")
        if (is.na(wessl) || delta > 0.01) {
          .vmsg("  -> rejected (CTA path): LOO UNSTABLE",
                "  WESS=", round(obs_ess, 4), "%",
                "  WESSL=", round(wessl %||% NA_real_, 4), "%")
          return(NULL)
        }
      } else if (identical(loo_arg, "pvalue")) {
        # Numeric loo passes loo_arg = "pvalue"; threshold is the original numeric
        # value when is.numeric(loo), otherwise the default 0.05.
        loo_pv     <- if (!is.null(loo_result)) loo_result$p_value %||% NA_real_
                      else NA_real_
        loo_thresh <- if (is.numeric(loo)) as.double(loo) else 0.05
        .vmsg("  [CTA-scan] LOO p=", round(loo_pv %||% NA_real_, 4),
              "  threshold=", loo_thresh)
        if (is.na(loo_pv) || loo_pv >= loo_thresh) {
          .vmsg("  -> rejected (CTA path): LOO P_FAIL  p=",
                round(loo_pv %||% NA_real_, 4),
                "  threshold=", loo_thresh)
          return(NULL)
        }
      }
    }

    # Candidate accepted
    .vmsg("  -> accepted (CTA path): ", attr_names[j],
          "  WESS=", round(obs_ess, 2), "%  p=", round(p_mc, 4))
    list(
      ok        = TRUE,
      rule      = cta_scan$rule,
      ess       = obs_ess,
      ess_pac   = obs_ess,
      p_mc      = p_mc,
      ge_count  = mc_res$ge_count,
      iter_used = mc_res$iter_used,
      attr_type = "ordered",
      engine    = "binary",
      confusion = NULL,
      loo       = if (!is.null(loo_result))
        list(allowed = TRUE,
             ess_loo = loo_result$ess_loo %||% obs_ess,
             p_value = loo_result$p_value %||% NA_real_)
      else
        NULL
    )
  }

  # Split node record
  .split_nd <- function(nid, parent_id, depth, idx, cand, appl) {
    sl  <- sort(unique(appl$y_pred[!is.na(appl$y_pred)]))
    conf <- if (appl$engine == "multiclass") cand$fit$confusion
            else .conf_binary(y[idx], appl$y_pred_raw)
    li  <- .loo_info(cand$fit)
    list(
      node_id = nid, parent_id = parent_id, depth = depth,
      n_obs = length(idx), n_weighted = sum(w[idx]),
      leaf = FALSE, majority_class = .majority_class(idx),
      attribute = cand$name, attr_col = cand$j,
      attr_type = cand$fit$attr_type %||% "ordered",
      rule = appl$rule, ess = cand$ess, ess_weighted = cand$ess,
      p_mc      = cand$p_mc,
      ge_count  = cand$ge_count  %||% NA_integer_,
      iter_used = cand$iter_used %||% NA_integer_,
      loo_status = li$status,
      loo_ess = li$ess_loo, loo_p = li$p_value,
      confusion = conf,
      split_labels = as.integer(sl),
      child_ids = integer(0)
    )
  }

  # Fast screen: mcarlo=FALSE for all attributes at idx.
  # Returns list of list(j, ess), sorted by ess descending.
  .fast_screen <- function(idx, pos_label = "?") {
    .vmsg("[CTA node ", pos_label, "] fast screen: ", length(idx), " obs, ", n_attrs, " attrs")
    if (length(idx) < mindenom) return(list())
    y_n <- y[idx]; w_n <- w[idx]
    if (length(unique(y_n)) < 2L) return(list())
    rh_fs <- if (!is.null(diag_env)) .row_hash(idx) else NULL
    result <- list()
    for (j in seq_len(n_attrs)) {
      x_j <- X[[j]][idx]
      x_v <- x_j[!is.na(x_j)]
      if (!is.null(miss_codes)) x_v <- x_v[!(x_v %in% miss_codes)]
      if (length(x_v) < 2L || length(unique(x_v)) < 2L) next
      t0_fs <- if (!is.null(diag_env)) proc.time()[["elapsed"]] else NULL
      fit_f <- tryCatch(
        oda_fit(x = x_j, y = y_n, w = w_n, priors_on = priors_on,
                miss_codes = miss_codes, K_segments = K_segments,
                mcarlo = FALSE, loo = "off", mindenom = mindenom),
        error = function(e) list(ok = FALSE))
      if (!is.null(diag_env)) {
        elapsed_fs <- proc.time()[["elapsed"]] - t0_fs
        ess_fs     <- if (isTRUE(fit_f$ok)) .get_ess(fit_f) else NA_real_
        passed_fs  <- isTRUE(fit_f$ok) && !is.na(ess_fs) && ess_fs >= ess_min
        diag_env$fast_log[[length(diag_env$fast_log) + 1L]] <- list(
          pos_id      = pos_label,
          attr_name   = attr_names[j],
          row_hash    = rh_fs,
          n_obs       = length(idx),
          ess_screen  = ess_fs,
          passed      = passed_fs,
          elapsed_sec = elapsed_fs
        )
      }
      if (!isTRUE(fit_f$ok)) next
      ess_f <- .get_ess(fit_f)
      if (is.na(ess_f) || ess_f < ess_min) next
      result[[length(result) + 1L]] <- list(j = j, ess = ess_f)
    }
    if (length(result) > 1L) {
      ord    <- order(vapply(result, `[[`, 0.0, "ess"), decreasing = TRUE)
      result <- result[ord]
    }
    result
  }

  # Full oda_fit (mcarlo=TRUE, loo=loo_arg) for attribute j at idx.
  # Returns a validated candidate list(j,name,fit,ess,p_mc) or NULL if rejected.
  .full_fit_one <- function(j, idx, pos_id) {
    y_n <- y[idx]; w_n <- w[idx]; x_j <- X[[j]][idx]
    # ---- CTA-specific ordered-cut path (weighted binary ordered predictors) ----
    if (any(w_n != w_n[1L]) &&
        length(unique(as.integer(y_n[!is.na(y_n)]))) == 2L) {
      x_v <- x_j[!is.na(x_j)]
      if (!is.null(miss_codes)) x_v <- x_v[!(x_v %in% miss_codes)]
      if (length(unique(x_v)) > 2L) {
        .vmsg("[CTA node ", pos_id, "] CTA-scan path: ", attr_names[j])
        rh      <- .row_hash(idx)
        cta_fit <- .cta_full_fit_ordered(j, x_j, y_n, w_n, pos_id,
                                          row_hash = rh, n_obs_parent = length(idx))
        if (!is.null(cta_fit))
          return(list(j = j, name = attr_names[j], fit = cta_fit,
                      ess = cta_fit$ess, p_mc = cta_fit$p_mc,
                      ge_count  = cta_fit$ge_count  %||% NA_integer_,
                      iter_used = cta_fit$iter_used %||% NA_integer_))
        return(NULL)  # triggered and rejected; do not fall back to generic ODA
      }
    }
    # ---- end CTA path ---------------------------------------------------------
    seed_j <- if (is.null(mc_seed)) NULL else
      (as.integer(mc_seed) + as.integer(pos_id) * 100L + j) %% .Machine$integer.max
    t0 <- proc.time()[["elapsed"]]
    .vmsg("[CTA node ", pos_id, "] MC+LOO: ", attr_names[j])
    fit <- tryCatch(
      oda_fit(x = x_j, y = y_n, w = w_n, priors_on = priors_on,
              miss_codes = miss_codes, K_segments = K_segments,
              mcarlo = TRUE, mc_iter = as.integer(mc_iter),
              mc_target = mc_target, mc_stop = mc_stop,
              mc_stopup = mc_stopup, mc_seed = seed_j,
              loo = loo_arg, eval_order = "loo_then_mc",
              mindenom = mindenom),
      error = function(e) list(ok = FALSE))
    elapsed <- round(proc.time()[["elapsed"]] - t0, 2)
    rh_gen  <- if (!is.null(diag_env)) .row_hash(idx) else NULL
    # Diagnostic helper: log this generic-ODA call regardless of outcome.
    .log_gen <- function(p_mc_val, ess_val, accepted, reject_reason) {
      if (is.null(diag_env)) return(invisible(NULL))
      diag_env$mc_log[[length(diag_env$mc_log) + 1L]] <- list(
        path             = "generic_oda",
        attr_name        = attr_names[j],
        pos_id           = pos_id,
        row_hash         = rh_gen,
        n_obs            = length(idx),
        obs_ess          = ess_val,
        mc_target        = mc_target,
        ge_count         = NA_integer_,
        iter_used        = NA_integer_,
        stop_reason      = NA_character_,
        p_mc_current     = p_mc_val,
        p_mc_raw         = p_mc_val,
        p_mc_pseudocount = NA_real_,
        signif_current   = !is.na(p_mc_val) && p_mc_val < alpha_split,
        signif_raw       = NA,
        elapsed_sec      = elapsed,
        accepted         = accepted,
        reject_reason    = reject_reason
      )
    }
    if (!isTRUE(fit$ok)) {
      .vmsg("  -> rejected: ", fit$reason %||% "unknown", " (", elapsed, "s)")
      .log_gen(NA_real_, NA_real_, FALSE, "fit_failed")
      return(NULL)
    }
    p_mc <- fit$p_mc
    if (is.na(p_mc) || p_mc >= alpha_split) {
      .vmsg("  -> rejected: p=", round(p_mc %||% NA_real_, 4), " (", elapsed, "s)")
      .log_gen(p_mc, NA_real_, FALSE, "signif_fail")
      return(NULL)
    }
    ess <- .get_ess(fit)
    if (is.na(ess) || ess < ess_min) {
      .vmsg("  -> rejected: ESS=", round(ess %||% NA_real_, 2), "% (", elapsed, "s)")
      .log_gen(p_mc, ess, FALSE, "ess_low")
      return(NULL)
    }
    if (identical(loo_arg, "stable")) {
      if (is.null(fit$loo) || !isTRUE(fit$loo$allowed)) {
        .vmsg("  -> rejected: LOO unstable (", elapsed, "s)")
        .log_gen(p_mc, ess, FALSE, "loo_unstable")
        return(NULL)
      }
    } else if (is.numeric(loo)) {
      loo_pv <- fit$loo$p_value %||% NA_real_
      if (is.na(loo_pv) || loo_pv >= loo) {
        .vmsg("  -> rejected: LOO p=", round(loo_pv, 4), " (", elapsed, "s)")
        .log_gen(p_mc, ess, FALSE, "loo_p_fail")
        return(NULL)
      }
    }
    .vmsg("  -> accepted: ESS=", round(ess, 2), "% p=", round(p_mc, 4),
          " (", elapsed, "s)")
    .log_gen(p_mc, ess, TRUE, NA_character_)
    list(j = j, name = attr_names[j], fit = fit, ess = ess, p_mc = p_mc,
         ge_count  = fit$mc_info$ge_count  %||% NA_integer_,
         iter_used = fit$mc_info$iter_used %||% NA_integer_)
  }

  # Get ALL valid candidates at a node position (fast screen then ALL full MC+LOO).
  # pos_id: unique position for seed derivation (1=root,2=left,3=right).
  .all_cands <- function(idx, pos_id) {
    fast <- .fast_screen(idx, pos_label = pos_id)
    if (length(fast) == 0L) return(list())
    valid <- list()
    for (fp in fast) {
      cand <- .full_fit_one(fp$j, idx, pos_id)
      if (!is.null(cand)) valid[[length(valid) + 1L]] <- cand
    }
    if (length(valid) > 1L) {
      ord   <- order(vapply(valid, `[[`, 0.0, "ess"), decreasing = TRUE)
      valid <- valid[ord]
    }
    valid
  }

  # HO-CTA growth: greedy (highest fast-ESS passing full MC+LOO) at each node.
  # Fills env$nodes (list indexed by integer node_id), increments env$counter.
  # Returns the root node_id of the grown subtree.
  # pred_class: class assigned by the parent ODA rule for the branch leading here.
  .ho_grow <- function(idx, depth, parent_id, env, pred_class = NULL) {
    env$counter <- env$counter + 1L
    nid <- env$counter

    nd     <- .leaf_nd(nid, parent_id, depth, idx, pred_class = pred_class)
    y_node <- y[idx]
    n_obs  <- length(idx)

    if (length(unique(y_node)) < 2L || n_obs < mindenom || depth >= max_depth) {
      env$nodes[[nid]] <- nd
      return(nid)
    }

    fast <- .fast_screen(idx, pos_label = nid)
    if (length(fast) == 0L) {
      env$nodes[[nid]] <- nd
      return(nid)
    }

    # HO-CTA greedy: accept the first (highest fast-ESS) candidate passing full MC+LOO
    best <- NULL
    for (fp in fast) {
      cand <- .full_fit_one(fp$j, idx, pos_id = nid)
      if (!is.null(cand)) { best <- cand; break }
    }

    if (is.null(best)) {
      env$nodes[[nid]] <- nd
      return(nid)
    }

    appl <- .apply_cand(best, idx)
    sl   <- sort(unique(appl$y_pred[!is.na(appl$y_pred)]))
    nd   <- .split_nd(nid, parent_id, depth, idx, best, appl)

    child_ids <- integer(length(sl))
    for (k in seq_along(sl)) {
      child_idx <- idx[!is.na(appl$y_pred) & appl$y_pred == sl[k]]
      if (length(child_idx) == 0L) next
      cid           <- .ho_grow(child_idx, depth + 1L, nid, env,
                                pred_class = as.integer(sl[k]))
      child_ids[k]  <- cid
    }
    nd$split_labels  <- as.integer(sl)
    nd$child_ids     <- child_ids[child_ids > 0L]
    env$nodes[[nid]] <- nd
    return(nid)
  }

  # Predict all n observations through a nodes_list (root at root_id).
  # Missing obs are routed to the node's majority class (majority fallback).
  .predict_all <- function(nodes_list, root_id) {
    pred_one <- function(i) {
      nid <- root_id
      repeat {
        nd <- nodes_list[[nid]]
        if (is.null(nd) || isTRUE(nd$leaf) || length(nd$child_ids) == 0L)
          return(nd$majority_class %||% NA_integer_)
        j     <- nd$attr_col
        x_val <- X[[j]][i]
        miss_here <- is.na(x_val)
        if (!is.null(miss_codes) && !miss_here) miss_here <- x_val %in% miss_codes
        if (miss_here) {
          return(nd$majority_class)
        }
        rule  <- nd$rule
        if (!is.null(rule$type) &&
            rule$type %in% c("multiclass_ordered","multiclass_nominal")) {
          # Multiclass: y_hat is actual class label; look up in split_labels.
          y_hat <- as.integer(oda_rule_predict_multiclass(
            x_val, rule, boundary = rule$boundary %||% "megaoda_halfopen"))
          sl  <- nd$split_labels
          ic  <- which(sl == y_hat)
          if (length(ic) == 0L) return(nd$majority_class)
        } else {
          # Binary: y_hat is rule-side 0/1; route by position (0->child 1, 1->child 2).
          y_hat <- as.integer(oda_rule_predict(x_val, rule))
          if (!(y_hat %in% 0:1)) return(nd$majority_class)
          ic <- y_hat + 1L
        }
        nid <- nd$child_ids[ic[1L]]
        if (is.na(nid) || nid < 1L) return(nd$majority_class)
      }
    }
    vapply(seq_len(n), pred_one, integer(1L))
  }

  # Sidak-Bonferroni + maximum-accuracy post-growth pruning.
  #
  # Algorithm (MPE Chapters 10-11; MODEL1.TXT canon):
  #   Repeat until stable:
  #     1. Count reachable split nodes K in current tree.
  #     2. Compute per-comparison Sidak threshold: alpha_k = 1-(1-prune_alpha)^(1/K).
  #     3. Flag every split node whose p_mc >= alpha_k (Sidak-unreliable).
  #     4. For each flagged node: tentatively collapse it to a leaf; score the
  #        resulting tree with .predict_all() / .wess_classes().
  #     5. Apply the collapse that yields the highest WESS >= current WESS
  #        (equal-WESS tie: keep first-found in BFS order -> shallowest/leftmost).
  #     6. If no collapse improves or maintains WESS: stop.
  #
  # Captures: y, w, prune_alpha, .predict_all, .wess_classes from enclosing scope.
  # Reusable: operates only on a node list; does not depend on enumeration structure.
  # Returns: a (possibly modified) copy of nodes_list.
  .prune_tree <- function(nodes_list, root_id = 1L) {
    if (prune_alpha >= 1.0) return(nodes_list)

    # BFS traversal - returns reachable node IDs in breadth-first (depth) order.
    .bfs_ids <- function(nlist, rid) {
      visited <- integer(0)
      queue   <- rid
      while (length(queue) > 0L) {
        nid   <- queue[1L]; queue <- queue[-1L]
        if (nid %in% visited) next
        visited <- c(visited, nid)
        nd <- nlist[[nid]]
        if (!is.null(nd) && !isTRUE(nd$leaf) && length(nd$child_ids) > 0L)
          queue <- c(queue, nd$child_ids[nd$child_ids > 0L])
      }
      visited
    }

    current <- nodes_list

    repeat {
      all_ids   <- .bfs_ids(current, root_id)
      split_ids <- all_ids[vapply(all_ids, function(nid) {
        nd <- current[[nid]]
        !is.null(nd) && !isTRUE(nd$leaf) && length(nd$child_ids) > 0L
      }, logical(1L))]

      K <- length(split_ids)
      if (K == 0L) break

      alpha_k <- 1.0 - (1.0 - prune_alpha)^(1.0 / K)

      # Sidak-flagged split nodes (BFS order = shallowest first for tie-breaking)
      flagged <- split_ids[vapply(split_ids, function(nid) {
        p <- current[[nid]]$p_mc
        !is.na(p) && p >= alpha_k
      }, logical(1L))]

      if (length(flagged) == 0L) break

      current_wess <- tryCatch({
        .wess_classes(y, .predict_all(current, root_id), w)
      }, error = function(e) -Inf)

      # Find the best-WESS eligible collapse (must not decrease WESS).
      # Initialise threshold just below current so equal-WESS pruning is accepted.
      best_nid    <- NA_integer_
      best_wess_p <- current_wess - 1e-8

      for (nid in flagged) {
        trial <- current
        trial[[nid]]$leaf      <- TRUE
        trial[[nid]]$child_ids <- integer(0)
        tw <- tryCatch({
          .wess_classes(y, .predict_all(trial, root_id), w)
        }, error = function(e) -Inf)
        if (tw > best_wess_p) {
          best_wess_p <- tw
          best_nid    <- nid
        }
      }

      if (is.na(best_nid)) break   # no eligible collapse - stop

      current[[best_nid]]$leaf      <- TRUE
      current[[best_nid]]$child_ids <- integer(0)
      .vmsg("[PRUNE] node ", best_nid, ": WESS ",
            round(current_wess, 2), "% -> ", round(best_wess_p, 2), "%")
    }

    current
  }

  # ---- ENUMERATE: top-3-node enumeration ------------------------------------

  root_cands <- .all_cands(seq_len(n), pos_id = 1L)

  # No valid root -> single-leaf tree (no_tree = TRUE)
  if (length(root_cands) == 0L) {
    .vmsg("[CTA] no valid root found, returning leaf node")
    nodes_out <- list(.leaf_nd(1L, 0L, 1L, seq_len(n)))
    return(structure(
      list(nodes = nodes_out, root_id = 1L, n_nodes = 1L, no_tree = TRUE,
           n = n, C = C, attr_names = attr_names, n_attrs = n_attrs,
           miss_codes = miss_codes, priors_on = priors_on,
           alpha_split = alpha_split, mindenom = mindenom,
           prune_alpha = prune_alpha, max_depth = max_depth,
           ess_min = ess_min, loo = loo,
           # Final tree objective on the ESS/WESS scale: WESS when weights are
           # active, ESS otherwise.  NA for no-tree fits.  Retained for the
           # cta_d_stat() contract; do not rename without updating accessors.
           overall_ess       = NA_real_,
           has_weights       = any(w != 1),
           training_confusion = NULL),
      class = "cta_tree"))
  }

  # ---- ENUMERATE: full AxBxC top-3-node enumeration ---------------------------
  # Canonical algorithm: evaluate every valid combination of:
  #   A = root split (all valid root candidates, sorted by root ESS desc)
  #   B = left-child option  (all valid split candidates + explicit leaf sentinel)
  #   C = right-child option (all valid split candidates + explicit leaf sentinel)
  #
  # Leaf sentinel is included even when split candidates exist: a locally
  # significant child split may still reduce full-tree WESS after pruning,
  # missingness, or interaction with the other branch.
  #
  # (B=leaf, C=leaf) is excluded from this expanded loop.  That configuration
  # is handled by the root-only stump phase below under path-local scoring.
  # Including it here under majority-fallback .predict_all() would duplicate
  # it with different missingness behaviour and risk conflicting WESS values
  # on datasets with missing root-attribute observations (Myeloma V17).
  #
  # B-subtree caching: B's sub-tree is grown once per (A, B) pair and reused
  # across all C options for that (A, B).  C's sub-tree is re-grown per
  # (A, B, C) triple; re-indexing optimisation deferred.
  #
  # Tie-breaking (FIRST IDENTIFIED): A sorted by root ESS desc; within each A,
  # B split candidates sorted by left-branch ESS desc then leaf last; within
  # each (A, B), C split candidates sorted by right-branch ESS desc then leaf last.

  best_wess       <- -Inf
  best_nodes      <- NULL
  best_n_nodes    <- 0L
  best_confusion  <- NULL

  .vmsg("[ENUMERATE] ", length(root_cands), " root candidates (A\u00d7B\u00d7C)")

  for (i in seq_along(root_cands)) {
    A_cand <- root_cands[[i]]
    .vmsg("[ENUMERATE ", i, "/", length(root_cands), "] root=", A_cand$name,
          " ESS=", round(A_cand$ess, 2), "%")

    appl_A  <- .apply_cand(A_cand, seq_len(n))
    valid_A <- !is.na(appl_A$y_pred)
    sl_A    <- sort(unique(appl_A$y_pred[valid_A]))
    if (length(sl_A) < 2L) next

    left_idx  <- which(valid_A & appl_A$y_pred == sl_A[1L])
    right_idx <- which(valid_A & appl_A$y_pred == sl_A[2L])

    # All valid split candidates for B (left) and C (right) positions.
    B_cands <- .all_cands(left_idx,  pos_id = 2L)
    C_cands <- .all_cands(right_idx, pos_id = 3L)

    # B/C options: split candidates (ESS desc) then leaf sentinel (NULL = leaf)
    B_options <- c(B_cands, list(NULL))
    C_options <- c(C_cands, list(NULL))

    .vmsg("  B: ", length(B_cands), " split(s) + leaf")
    .vmsg("  C: ", length(C_cands), " split(s) + leaf")

    for (bi in seq_along(B_options)) {
      B_opt     <- B_options[[bi]]
      B_is_leaf <- is.null(B_opt)

      # ---- Build B sub-tree (cached for all C at this (A, B)) -----------------
      if (B_is_leaf) {
        nd2   <- .leaf_nd(2L, 1L, 2L, left_idx, pred_class = as.integer(sl_A[1L]))
        B_max <- 2L
      } else {
        B_cand <- B_opt
        appl_B <- .apply_cand(B_cand, left_idx)
        sl_B   <- sort(unique(appl_B$y_pred[!is.na(appl_B$y_pred)]))
        if (length(sl_B) < 2L) next   # degenerate B split - skip

        nd2 <- .split_nd(2L, 1L, 2L, left_idx, B_cand, appl_B)

        # Grow HO-CTA below B's children (depth 3+); node 2 is already placed.
        B_sub_env         <- new.env(parent = emptyenv())
        B_sub_env$nodes   <- vector("list", n + 4L)
        B_sub_env$counter <- 2L   # children start at 3

        B_left_idx  <- left_idx[!is.na(appl_B$y_pred) & appl_B$y_pred == sl_B[1L]]
        B_right_idx <- left_idx[!is.na(appl_B$y_pred) & appl_B$y_pred == sl_B[2L]]

        cid_BL <- .ho_grow(B_left_idx,  3L, 2L, B_sub_env,
                           pred_class = as.integer(sl_B[1L]))
        cid_BR <- .ho_grow(B_right_idx, 3L, 2L, B_sub_env,
                           pred_class = as.integer(sl_B[2L]))
        B_max  <- B_sub_env$counter
        nd2$child_ids <- c(cid_BL, cid_BR)
      }

      # ---- Inner C loop -------------------------------------------------------
      for (ci in seq_along(C_options)) {
        C_opt     <- C_options[[ci]]
        C_is_leaf <- is.null(C_opt)

        # (B=leaf, C=leaf) handled by root-only stump phase - skip here
        if (B_is_leaf && C_is_leaf) next

        c_root_id <- B_max + 1L

        if (C_is_leaf) {
          nd_C  <- .leaf_nd(c_root_id, 1L, 2L, right_idx,
                            pred_class = as.integer(sl_A[2L]))
          C_max <- c_root_id
        } else {
          C_cand <- C_opt
          appl_C <- .apply_cand(C_cand, right_idx)
          sl_C   <- sort(unique(appl_C$y_pred[!is.na(appl_C$y_pred)]))
          if (length(sl_C) < 2L) next   # degenerate C split - skip

          nd_C <- .split_nd(c_root_id, 1L, 2L, right_idx, C_cand, appl_C)

          # Grow HO-CTA below C's children (depth 3+); c_root_id is already placed.
          C_sub_env         <- new.env(parent = emptyenv())
          C_sub_env$nodes   <- vector("list", n + B_max + 4L)
          C_sub_env$counter <- c_root_id   # children start at c_root_id + 1

          C_left_idx  <- right_idx[!is.na(appl_C$y_pred) & appl_C$y_pred == sl_C[1L]]
          C_right_idx <- right_idx[!is.na(appl_C$y_pred) & appl_C$y_pred == sl_C[2L]]

          cid_CL <- .ho_grow(C_left_idx,  3L, c_root_id, C_sub_env,
                             pred_class = as.integer(sl_C[1L]))
          cid_CR <- .ho_grow(C_right_idx, 3L, c_root_id, C_sub_env,
                             pred_class = as.integer(sl_C[2L]))
          C_max <- C_sub_env$counter
          nd_C$child_ids <- c(cid_CL, cid_CR)
        }

        B_label <- if (B_is_leaf) "leaf" else B_opt$name
        C_label <- if (C_is_leaf) "leaf" else C_opt$name
        .vmsg("  [B=", B_label, " C=", C_label, "]")

        # ---- Assemble full candidate node list --------------------------------
        cand_nodes <- vector("list", C_max)

        nd1 <- .split_nd(1L, 0L, 1L, seq_len(n), A_cand, appl_A)
        nd1$split_labels <- as.integer(sl_A)
        nd1$child_ids    <- c(2L, c_root_id)
        cand_nodes[[1L]] <- nd1
        cand_nodes[[2L]] <- nd2

        if (!B_is_leaf) {
          for (bid in seq(3L, B_max)) {
            nb <- B_sub_env$nodes[[bid]]
            if (!is.null(nb)) cand_nodes[[bid]] <- nb
          }
        }

        cand_nodes[[c_root_id]] <- nd_C

        if (!C_is_leaf) {
          for (cid_v in seq(c_root_id + 1L, C_max)) {
            nc <- C_sub_env$nodes[[cid_v]]
            if (!is.null(nc)) cand_nodes[[cid_v]] <- nc
          }
        }

        # ---- Prune, score, compare -------------------------------------------
        pruned_nodes <- .prune_tree(cand_nodes)

        preds     <- NULL
        wess_cand <- tryCatch({
          preds <- .predict_all(pruned_nodes, root_id = 1L)
          .wess_classes(y, preds, w)
        }, error = function(e) -Inf)

        .vmsg("    -> WESS=", round(wess_cand, 2), "%")

        # Reject degenerate trees: pruning may collapse a class-1 branch to a
        # class-0 majority node, leaving all terminals predicting the same class.
        # Such a tree is invalid regardless of its WESS value.
        if (!is.null(preds) && .cta_predictions_degenerate(preds)) {
          .vmsg("    -> degenerate (all predictions class ",
                unique(preds[!is.na(preds)]), ") - skipped")
          next
        }

        if (wess_cand > best_wess) {
          best_wess    <- wess_cand
          best_nodes   <- pruned_nodes
          best_n_nodes <- sum(!vapply(pruned_nodes, is.null, logical(1L)))
          # Capture final-tree confusion at the moment of selection.
          # preds is from majority-fallback .predict_all(); filter classified.
          if (!is.null(preds)) {
            ok_exp         <- !is.na(preds)
            best_confusion <- .make_training_conf(y[ok_exp], preds[ok_exp])
          }
        }
      }  # end C loop
    }  # end B loop
  }  # end A loop

  # ---- ENUMERATE root-only (stump) candidate phase ----------------------------
  # CTA.exe canonical: MODEL1.TXT Trees 5-7 show CTA.exe explicitly evaluating
  # each root candidate as a stump (root split + two leaves), scored path-locally
  # (observations whose root attribute is missing are excluded, not majority-routed).
  # This phase runs after the expanded phase; root-only candidates compete against
  # expanded candidates for the overall best WESS.
  #
  # For MINDENOM=30: V17 stump scores 186 obs -> WESS=16.51%;
  #                  V14 stump scores 255 obs -> WESS=14.06%.
  #                  V17 wins (16.51% > 14.06% and > any majority-fallback expanded score).
  # For MINDENOM=1:  expanded V14->V15 WESS=27.69% > all stumps -> expanded wins.
  .vmsg("[ENUMERATE root-only] ", length(root_cands), " stump candidates (path-local)")
  for (i in seq_along(root_cands)) {
    A_cand  <- root_cands[[i]]
    appl_A  <- .apply_cand(A_cand, seq_len(n))
    valid_A <- !is.na(appl_A$y_pred)
    sl_A    <- sort(unique(appl_A$y_pred[valid_A]))
    if (length(sl_A) < 2L) next

    # appl_A$y_pred already carries NA_integer_ for obs whose root attribute is
    # missing (value in miss_codes).  Scoring over ok=!is.na(stump_preds) is
    # therefore path-local: missing-root obs are excluded from WESS computation.
    left_idx  <- which(valid_A & appl_A$y_pred == sl_A[1L])
    right_idx <- which(valid_A & appl_A$y_pred == sl_A[2L])
    if (length(left_idx) < mindenom || length(right_idx) < mindenom) {
      .vmsg("  [root-only] ", A_cand$name,
            " skipped: child sizes ", length(left_idx), "/", length(right_idx),
            " < MINDENOM=", mindenom)
      next
    }

    stump_preds <- appl_A$y_pred
    ok          <- !is.na(stump_preds)
    # Canonical guard: stump predictions must cover both classes.
    # sl_A already has >= 2 levels (checked above), so this should never fire
    # in practice; the check is retained as a defensive invariant.
    if (.cta_predictions_degenerate(stump_preds)) {
      .vmsg("  [root-only] ", A_cand$name,
            " degenerate stump (all predictions same class) - skipped")
      next
    }
    wess_stump  <- if (sum(ok) < 2L) -Inf else
                     .wess_classes(y[ok], stump_preds[ok], w[ok])

    .vmsg("  [root-only] ", A_cand$name,
          " n_classified=", sum(ok),
          " WESS=", round(wess_stump, 2), "%")

    if (wess_stump > best_wess) {

      stump_nodes       <- vector("list", 3L)
      nd1               <- .split_nd(1L, 0L, 1L, seq_len(n), A_cand, appl_A)
      nd1$split_labels  <- as.integer(sl_A)
      nd1$child_ids     <- c(2L, 3L)
      stump_nodes[[1L]] <- nd1
      stump_nodes[[2L]] <- .leaf_nd(2L, 1L, 2L, left_idx,
                                    pred_class = as.integer(sl_A[1L]))
      stump_nodes[[3L]] <- .leaf_nd(3L, 1L, 2L, right_idx,
                                    pred_class = as.integer(sl_A[2L]))

      best_wess    <- wess_stump
      best_nodes   <- stump_nodes
      best_n_nodes <- 3L
      # Capture final-tree confusion at the moment of selection.
      # stump_preds[ok] is path-local: missing-root obs already excluded.
      best_confusion <- .make_training_conf(y[ok], stump_preds[ok])
    }
  }  # end root-only phase

  # Fallback: enumeration produced no valid tree
  no_tree <- is.null(best_nodes)
  if (no_tree) {
    .vmsg("[CTA] no valid tree from enumeration, returning leaf node")
    best_nodes   <- list(.leaf_nd(1L, 0L, 1L, seq_len(n)))
    best_n_nodes <- 1L
  } else {
    root_attr  <- best_nodes[[1L]]$attribute %||% "?"
    n_split    <- sum(!vapply(best_nodes,
                              function(nd) is.null(nd) || isTRUE(nd$leaf),
                              logical(1L)))
    .vmsg("[CTA] selected: root=", root_attr,
          " WESS=", round(best_wess, 2), "%",
          " split_nodes=", n_split)
  }

  structure(
    list(
      nodes       = best_nodes,
      root_id     = 1L,
      n_nodes     = best_n_nodes,
      no_tree     = no_tree,
      n           = n,
      C           = C,
      attr_names  = attr_names,
      n_attrs     = n_attrs,
      miss_codes  = miss_codes,
      priors_on   = priors_on,
      alpha_split = alpha_split,
      mindenom    = mindenom,
      prune_alpha = prune_alpha,
      max_depth   = max_depth,
      ess_min     = ess_min,
      loo         = loo,
      # Final tree objective on the ESS/WESS scale: WESS when weights are
      # active, ESS otherwise.  NA for no-tree fits.  Retained for the
      # cta_d_stat() contract; do not rename without updating accessors.
      overall_ess        = if (no_tree) NA_real_ else best_wess,
      has_weights        = any(w != 1),
      # Final selected tree training confusion (rows=actual, cols=predicted).
      # Captured at the exact moment the winning candidate is selected, using
      # the same scoring predictions.  NULL for no-tree fits.
      training_confusion = if (no_tree) NULL else best_confusion
    ),
    class = "cta_tree"
  )
}

# ---- Prediction -------------------------------------------------------------

#' Predict class labels from a fitted CTA tree
#'
#' Routes each observation down the tree to a leaf, returning that leaf's
#' majority class.  Observations with missing-coded attribute values at a
#' split node are routed to that node's majority class.
#'
#' @param object A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @param newdata Data frame or matrix matching training X column layout.
#' @param ... Unused.
#' @return Integer vector of predicted class labels, length \code{nrow(newdata)}.
#' @export
predict.cta_tree <- function(object, newdata,
                             missing_action = c("majority", "na"), ...) {
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  if (isTRUE(object$no_tree)) return(rep(NA_integer_, nrow(newdata)))
  missing_action <- match.arg(missing_action)
  n_new      <- nrow(newdata)
  miss_codes <- object$miss_codes

  predict_one <- function(i) {
    nid <- object$root_id
    repeat {
      nd <- object$nodes[[nid]]
      if (is.null(nd) || isTRUE(nd$leaf) || length(nd$child_ids) == 0L)
        return(nd$majority_class %||% NA_integer_)

      j     <- nd$attr_col
      x_val <- newdata[[j]][i]

      miss_here <- is.na(x_val)
      if (!is.null(miss_codes) && !miss_here) miss_here <- x_val %in% miss_codes
      if (miss_here) {
        if (missing_action == "na") return(NA_integer_)
        return(nd$majority_class)
      }

      rule  <- nd$rule
      if (!is.null(rule$type) &&
          rule$type %in% c("multiclass_ordered","multiclass_nominal")) {
        # Multiclass: y_hat is actual class label; look up in split_labels.
        y_hat <- as.integer(oda_rule_predict_multiclass(
          x_val, rule, boundary = rule$boundary %||% "megaoda_halfopen"))
        sl  <- nd$split_labels
        ic  <- which(sl == y_hat)
        if (length(ic) == 0L) return(nd$majority_class)
      } else {
        # Binary: y_hat is rule-side 0/1; route by position (0->child 1, 1->child 2).
        y_hat <- as.integer(oda_rule_predict(x_val, rule))
        if (!(y_hat %in% 0:1)) return(nd$majority_class)
        ic <- y_hat + 1L
      }
      nid <- nd$child_ids[ic[1L]]
      if (is.na(nid) || nid < 1L) return(nd$majority_class)
    }
  }

  vapply(seq_len(n_new), predict_one, integer(1L))
}

# ---- Print ------------------------------------------------------------------

#' Print a CTA tree
#'
#' @param x A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.cta_tree <- function(x, ...) {
  if (isTRUE(x$no_tree)) {
    cat(sprintf(
      "\nCTA Tree  alpha_split=%.3f  mindenom=%d  prune=%.3f  max_depth=%d  loo=%s\n\n",
      x$alpha_split, x$mindenom, x$prune_alpha, x$max_depth, x$loo))
    cat("No tree found (leaf-only): no valid split passed significance, LOO STABLE, and MINDENOM constraints.\n")
    return(invisible(x))
  }
  cat(sprintf(
    "\nCTA Tree  alpha_split=%.3f  mindenom=%d  prune=%.3f  max_depth=%d  loo=%s\n\n",
    x$alpha_split, x$mindenom, x$prune_alpha, x$max_depth, x$loo))

  hdr <- sprintf("%-14s %4s %4s %6s %7s %8s %8s %8s  MODEL\n",
                 "ATTRIBUTE", "NODE", "LEV", "OBS", "p", "ESS", "WESS", "LOO")
  cat(hdr)
  cat(strrep("-", 82), "\n")

  for (nid in seq_len(x$n_nodes)) {
    nd <- x$nodes[[nid]]
    if (is.null(nd) || isTRUE(nd$leaf)) next

    rule_str <- .cta_rule_string(nd$rule)
    loo_str  <- nd$loo_status %||% "-"
    p_str    <- if (!is.na(nd$p_mc))
      if (nd$p_mc < 0.001) ".000" else sprintf("%.3f", nd$p_mc)
    else "-"

    cat(sprintf("%-14s %4d %4d %6d %7s %7.2f%% %7.2f%% %8s  %s\n",
                nd$attribute %||% "",
                nd$node_id, nd$depth, nd$n_obs, p_str,
                nd$ess %||% 0, nd$ess_weighted %||% 0,
                loo_str, rule_str))

    if (!is.null(nd$confusion)) {
      cat("  Node-local split confusion (this rule only, observations at this node)\n")
      .print_node_confusion(nd$confusion)
    }
    cat("\n")
  }

  n_split <- sum(vapply(x$nodes, function(nd) !isTRUE(nd$leaf), logical(1)))
  n_leaf  <- x$n_nodes - n_split
  cat(sprintf("Nodes: %d total  (%d split  %d leaf)\n",
              x$n_nodes, n_split, n_leaf))

  # Terminal endpoints section — uses cta_endpoint_table() as the source.
  ep <- tryCatch(cta_endpoint_table(x), error = function(e) NULL)
  if (!is.null(ep) && nrow(ep) > 0L) {
    cat("\nTerminal endpoints (*):\n")
    for (i in seq_len(nrow(ep))) {
      ccr     <- ep$class_counts_raw[[i]]
      cnt_str <- if (!is.null(ccr))
        paste(sprintf("%s:%d", names(ccr), as.integer(ccr)), collapse = " ")
      else ""
      tp    <- ep$target_prop[i]
      tp_str <- if (!is.na(tp)) sprintf("  target_prop=%.1f%%", tp * 100) else ""
      cat(sprintf("* endpoint %d  node %d:  path=%s  n=%d%s  predicted=%d%s\n",
                  ep$endpoint_id[i], ep$leaf_node_id[i],
                  ep$path[i], ep$n[i],
                  if (nzchar(cnt_str)) sprintf("  counts=[%s]", cnt_str) else "",
                  ep$predicted_class[i], tp_str))
    }
  }

  # Compact footer: ESS/WESS, D, strata, min endpoint denominator.
  ess_val <- x$overall_ess
  ess_str <- if (!is.null(ess_val) && is.finite(ess_val))
    sprintf("%.2f%%", ess_val) else "NA"
  d_val   <- cta_d_stat(x)
  d_str   <- if (!is.na(d_val)) sprintf("%.4f", d_val) else "NA"
  s_val   <- cta_strata(x)
  s_str   <- if (!is.na(s_val)) as.character(s_val) else "NA"
  mtd_val <- cta_min_terminal_denom(x)
  mtd_str <- if (!is.na(mtd_val)) as.character(mtd_val) else "NA"
  ess_label <- if (isTRUE(x$has_weights)) "WESS" else "ESS"
  cat(sprintf("%s: %s  D: %s  strata: %s  min_denom: %s\n",
              ess_label, ess_str, d_str, s_str, mtd_str))

  invisible(x)
}

# ---- Node table -------------------------------------------------------------

#' Extract the CTA node summary as a data frame
#'
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @return Data frame with one row per node.
#' @export
cta_node_table <- function(tree) {
  rows <- lapply(tree$nodes, function(nd) {
    if (is.null(nd)) return(NULL)
    data.frame(
      node_id      = nd$node_id,
      parent_id    = nd$parent_id,
      level        = nd$depth,
      depth        = nd$depth,
      leaf         = isTRUE(nd$leaf),
      attribute    = nd$attribute    %||% NA_character_,
      attr_type    = nd$attr_type    %||% NA_character_,
      n_obs        = nd$n_obs,
      n_weighted   = round(nd$n_weighted, 4),
      p_mc         = nd$p_mc         %||% NA_real_,
      ess          = nd$ess          %||% NA_real_,
      ess_weighted = nd$ess_weighted %||% NA_real_,
      loo_status   = nd$loo_status   %||% NA_character_,
      loo_ess      = nd$loo_ess      %||% NA_real_,
      loo_p        = nd$loo_p        %||% NA_real_,
      model        = if (!isTRUE(nd$leaf))
                       tryCatch(.cta_model_string(tree, nd),
                                error = function(e) NA_character_)
                     else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  do.call(rbind, rows)
}

# ---- Internal helpers -------------------------------------------------------

# .cta_model_string(tree, nd)
#
# Build the CTA.exe-style MODEL field for a split node.  Format per branch:
#   branch_condition-->class,correct_n/total_n,pct[*]
# where `*` marks a terminal (leaf) branch.
#
# correct_n and total_n are derived from the split node's own confusion matrix
# (rows = actual class, cols = predicted class).  The column for each predicted
# class gives the branch totals; the diagonal cell gives the correct count.
# This is available for all split nodes regardless of tree depth.
.cta_model_string <- function(tree, nd) {
  if (is.null(nd) || isTRUE(nd$leaf)) return(NA_character_)
  if (length(nd$child_ids) == 0L)      return(.cta_rule_string(nd$rule))

  conf <- nd$confusion
  if (is.table(conf))  conf <- as.matrix(conf)
  if (!is.matrix(conf) || prod(dim(conf)) == 0L) conf <- NULL

  parts <- character(length(nd$child_ids))
  for (k in seq_along(nd$child_ids)) {
    child_id  <- nd$child_ids[k]
    pred_cls  <- as.integer(nd$split_labels[k])
    pc_str    <- as.character(pred_cls)

    # Branch label ("<=cut" / ">cut" etc.)
    branch_lbl <- tryCatch(.cta_branch_string(nd, k),
                           error = function(e) pc_str)

    # Per-branch counts from the split node's confusion column
    if (!is.null(conf) && pc_str %in% colnames(conf)) {
      col_v    <- conf[, pc_str]
      total_n  <- as.integer(round(sum(col_v)))
      corr_n   <- if (pc_str %in% rownames(conf))
        as.integer(round(conf[pc_str, pc_str])) else NA_integer_
      pct_str  <- if (!is.na(corr_n) && total_n > 0L)
        sprintf("%.2f%%", 100 * corr_n / total_n) else "?"
      cnt_str  <- if (!is.na(corr_n))
        sprintf("%d/%d", corr_n, total_n) else sprintf("?/%d", total_n)
    } else {
      cnt_str <- "?/?"
      pct_str <- "?"
    }

    # Terminal marker: `*` if child is a leaf
    child_nd <- tree$nodes[[child_id]]
    star     <- if (!is.null(child_nd) && isTRUE(child_nd$leaf)) "*" else ""
    parts[k] <- sprintf("%s-->%d,%s,%s%s",
                        branch_lbl, pred_cls, cnt_str, pct_str, star)
  }
  paste(parts, collapse = "; ")
}

.cta_rule_string <- function(rule) {
  if (is.null(rule)) return("")
  tryCatch({
    t <- rule$type
    if (identical(t, "ordered_cut")) {
      cv <- rule$cut_value
      if (rule$direction == "0->1")
        sprintf("<=%g-->0; >%g-->1", cv, cv)
      else
        sprintf("<=%g-->1; >%g-->0", cv, cv)
    } else if (identical(t, "multiclass_ordered")) {
      cuts <- rule$cut_values; segs <- rule$seg_classes; K <- length(segs)
      parts <- character(K)
      for (k in seq_len(K)) {
        if      (k == 1L) parts[k] <- sprintf("<=%g-->%d", cuts[1L], segs[k])
        else if (k == K)  parts[k] <- sprintf(">%g-->%d",  cuts[K-1L], segs[k])
        else               parts[k] <- sprintf("%g<%s<=%g-->%d",
                                               cuts[k-1L], "x", cuts[k], segs[k])
      }
      paste(parts, collapse = "; ")
    } else if (t %in% c("binary_map", "nominal_cut")) {
      sprintf("{%s}-->0; {%s}-->1",
              paste(rule$left_levels,  collapse = ","),
              paste(rule$right_levels, collapse = ","))
    } else ""
  }, error = function(e) "")
}

.print_node_confusion <- function(conf) {
  if (is.table(conf)) conf <- as.matrix(conf)
  if (!is.matrix(conf) || nrow(conf) < 2L) return(invisible(NULL))

  classes <- rownames(conf) %||% as.character(seq_len(nrow(conf)))
  cw      <- max(7L, max(nchar(classes)) + 2L)

  cat(sprintf("%*s", 12L, ""), paste(sprintf("%*s", cw, classes), collapse=""), "\n")
  cat(sprintf("%*s", 12L, ""), strrep("-", cw * ncol(conf)), "\n")

  for (i in seq_len(nrow(conf))) {
    rt  <- sum(conf[i, ])
    pac <- if (rt > 0) 100 * conf[i, i] / rt else 0
    row_str <- paste(sprintf("%*d", cw, as.integer(round(conf[i, ]))), collapse="")
    cat(sprintf("  %6s  | %s | %6.2f%%\n", classes[i], row_str, pac))
  }

  cat(sprintf("%*s", 12L, ""), strrep("-", cw * ncol(conf)), "\n")
  col_tots <- colSums(conf)
  cat(sprintf("  %6s  | %s\n", "NP",
              paste(sprintf("%*d", cw, as.integer(round(col_tots))), collapse="")))
  invisible(NULL)
}
