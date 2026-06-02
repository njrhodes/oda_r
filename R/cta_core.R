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

# .cta_mc_precomp_ord()
#
# Precompute the fixed per-node/attribute structure used by .cta_fast_scan_perm().
# Called once per (node, attribute) pair before the MC permutation loop.
#
# Returns a list with:
#   keep       — logical mask: which raw observations survive miss_code removal
#   w_clean    — case weights for kept observations (numeric)
#   block_idx  — integer 1..m: unique-value block membership for each kept obs
#   uvals      — sorted unique x values (length m)
#   m          — number of unique values
#   nv         — per-block observation counts (integer, length m)
#   Nv         — cumsum(nv): left-side observation count at each cut (integer, length m)
#   Nm         — total kept observations (integer)
#   mindenom   — stored for use inside .cta_fast_scan_perm()
.cta_mc_precomp_ord <- function(x_raw, w_raw, miss_codes, mindenom) {
  miss_mask <- is.na(x_raw)
  if (!is.null(miss_codes)) miss_mask <- miss_mask | (x_raw %in% miss_codes)
  keep      <- !miss_mask
  x_c       <- x_raw[keep]
  w_c       <- as.numeric(w_raw[keep])
  uvals     <- sort(unique(x_c))
  m         <- length(uvals)
  block_idx <- match(x_c, uvals)          # integer 1..m per kept obs
  nv        <- tabulate(block_idx, nbins = m)
  Nv        <- cumsum(nv)
  list(keep      = keep,
       w_clean   = w_c,
       block_idx = block_idx,
       uvals     = uvals,
       m         = m,
       nv        = as.integer(nv),
       Nv        = as.integer(Nv),
       Nm        = as.integer(sum(nv)),
       mindenom  = as.integer(mindenom))
}

# .cta_fast_scan_perm()
#
# Fast per-permutation CTA ordered scan.  Equivalent to .cta_ordered_scan() called
# with a permuted y, but uses precomputed block structure to eliminate the O(m)
# R-level for-loop over unique x values.
#
# Replaces:
#   scan_b <- .cta_ordered_scan(x, y_star, w, priors_on, miss_codes, mindenom)
# with:
#   scan_b <- .cta_fast_scan_perm(y_star, precomp, priors_on)
# where precomp = .cta_mc_precomp_ord(x, w, miss_codes, mindenom) built once.
#
# Canonical rule preserved:
#   "0->1": rightmost eligible j where sens_wa = 1 - C1a[j] > 0.5
#   "1->0": leftmost  eligible j where sens_wa = C1a[j]     > 0.5
#   direction ties -> "0->1" wins
#   mindenom applied to both nL and nR
#
# Returns the same structure as .cta_ordered_scan(): list(rule, ess, sens_wa, spec_wa)
# or NULL when no eligible cut exists.
.cta_fast_scan_perm <- function(y_star, precomp, priors_on) {
  if (precomp$m < 2L) return(NULL)              # no candidate cuts possible
  y_c <- as.integer(y_star)[precomp$keep]
  n0  <- sum(y_c == 0L); n1 <- sum(y_c == 1L)
  if (n0 == 0L || n1 == 0L) return(NULL)

  w_c <- precomp$w_clean
  if (isTRUE(priors_on)) {
    W0 <- sum(w_c[y_c == 0L]); W1 <- sum(w_c[y_c == 1L])
    if (W0 <= 0 || W1 <= 0) return(NULL)
    w0_adj <- w_c * (y_c == 0L) / W0
    w1_adj <- w_c * (y_c == 1L) / W1
  } else {
    w0_adj <- w_c * (y_c == 0L)
    w1_adj <- w_c * (y_c == 1L)
  }

  # Block aggregation via rowsum() (C-level) — replaces the O(m) R for-loop
  z0a <- as.numeric(rowsum(w0_adj, precomp$block_idx, reorder = TRUE))
  z1a <- as.numeric(rowsum(w1_adj, precomp$block_idx, reorder = TRUE))
  C0a <- cumsum(z0a); C1a <- cumsum(z1a)

  Nv       <- precomp$Nv; Nm <- precomp$Nm; m <- precomp$m
  uvals    <- precomp$uvals; mindenom <- precomp$mindenom
  j_seq    <- seq_len(m - 1L)
  md_ok    <- Nv[j_seq] >= mindenom & (Nm - Nv[j_seq]) >= mindenom

  # "0->1": rightmost eligible j with sens = 1 - C1a[j] > 0.5
  elig01 <- j_seq[md_ok & C1a[j_seq] < 0.5]
  r01 <- if (length(elig01) > 0L) {
    j     <- elig01[length(elig01)]          # rightmost
    sens  <- 1 - C1a[j]; spec <- C0a[j]
    list(rule     = list(type = "ordered_cut", direction = "0->1",
                         cut_value = (uvals[j] + uvals[j + 1L]) / 2),
         ess      = (spec + sens - 1) * 100,
         sens_wa  = sens, spec_wa  = spec)
  } else NULL

  # "1->0": leftmost eligible j with sens = C1a[j] > 0.5
  elig10 <- j_seq[md_ok & C1a[j_seq] > 0.5]
  r10 <- if (length(elig10) > 0L) {
    j     <- elig10[1L]                      # leftmost
    sens  <- C1a[j]; spec <- 1 - C0a[j]
    list(rule     = list(type = "ordered_cut", direction = "1->0",
                         cut_value = (uvals[j] + uvals[j + 1L]) / 2),
         ess      = (spec + sens - 1) * 100,
         sens_wa  = sens, spec_wa  = spec)
  } else NULL

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

  # Precompute fixed block structure once per (node, attribute) call
  precomp <- .cta_mc_precomp_ord(x, w, miss_codes, mindenom)

  for (b in seq_len(mc_iter)) {
    iter_used <- b
    y_star    <- sample(y_coded, size = n, replace = FALSE)
    scan_b    <- .cta_fast_scan_perm(y_star, precomp, priors_on)
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
#'   (MegaODA LOO STABLE; accept when |WESSL - WESS| <= 0.01 pp; reports
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
    mc_stopup   = NULL,
    mc_seed     = NULL,
    loo         = "off",
    attr_names  = NULL,
    K_segments  = NULL,
    verbose     = FALSE,
    diag_env    = NULL
) {
  # ---- Validate and prep ----------------------------------------------------
  # CTA STOP policy: the single STOP value governs both tails (canon: CTA.exe
  # has no separate STOPUP command).  NULL (omitted) -> mirror mc_stop.
  # Explicit NA -> bypass nonsignificance stopping.  Numeric -> use as-is.
  mc_stopup_eff <- if (is.null(mc_stopup)) mc_stop else mc_stopup

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
  #   "OFF"     -  LOO mode is "off" or LOO did not run
  #   "STABLE"  -  candidate passed the LOO STABLE gate (|WESSL - WESS| <= 0.01 pp)
  #   "PVALUE"  -  candidate passed the LOO p-value gate
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
      mc_stopup  = mc_stopup_eff,
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
    # Canon: MC significance gate first; LOO only for Signif T candidates.
    # EXE family logs confirm LOO STABLE is reported only for selected/significant
    # nodes, never for Signif F candidates.  Running LOO unconditionally inside
    # oda_fit() for all candidates is incorrect gate ordering and wastes n-fold
    # LOO work on candidates that will be rejected by MC anyway.
    .vmsg("[CTA node ", pos_id, "] MC: ", attr_names[j])
    fit <- tryCatch(
      oda_fit(x = x_j, y = y_n, w = w_n, priors_on = priors_on,
              miss_codes = miss_codes, K_segments = K_segments,
              mcarlo = TRUE, mc_iter = as.integer(mc_iter),
              mc_target = mc_target, mc_stop = mc_stop,
              mc_stopup = mc_stopup_eff, mc_seed = seed_j,
              loo = "off",
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
    # LOO gate: only for Signif T candidates that survived MC and ESS screens.
    if (!identical(loo_arg, "off")) {
      seed_loo <- if (is.null(mc_seed)) NULL else seed_j + 1L
      t0_loo   <- proc.time()[["elapsed"]]
      # Recode y_n to {0,1} to match oda_fit()'s internal coding used for fit$rule.
      # oda_fit() recodes y to {0,1} internally; fit$rule is in that {0,1} space.
      bin_labels_n <- sort(unique(as.integer(y_n[!is.na(y_n)])))
      y_coded_n    <- ifelse(as.integer(y_n) == bin_labels_n[1L], 0L, 1L)
      # Clean miss codes before LOO: oda_loo_for_rule must receive the same
      # cleaned data that oda_univariate_core used for training.  Without this
      # cleaning, miss-coded obs (e.g. -9) that are not NA get classified by
      # oda_rule_predict() in the binary_map fast path, contaminating ess_loo
      # and causing the STABLE check to fail spuriously.
      loo_clean    <- oda_clean_xy(x_j, y_coded_n, w_n, miss_codes = miss_codes)
      loo_out  <- tryCatch(
        oda_loo_for_rule(
          x          = loo_clean$x,
          y          = loo_clean$y,
          rule       = fit$rule,
          w          = loo_clean$w,
          chance_model = "class",
          k_attr     = fit$k_attr,
          attr_type  = fit$attr_type %||% "ordered",
          priors_on  = priors_on,
          miss_codes = miss_codes,
          mc_iter    = as.integer(mc_iter),
          mc_target  = mc_target,
          mc_stop    = mc_stop,
          mc_stopup  = mc_stopup_eff,
          mc_seed    = seed_loo
        ),
        error = function(e) NULL
      )
      loo_elapsed <- round(proc.time()[["elapsed"]] - t0_loo, 2)
      elapsed     <- elapsed + loo_elapsed
      if (!is.null(diag_env)) {
        wessl_d <- if (!is.null(loo_out)) loo_out$ess_loo %||% NA_real_ else NA_real_
        loo_stable_flag <- if (identical(loo_arg, "stable"))
          !is.na(wessl_d) && isTRUE(all.equal(ess, wessl_d, tolerance = 1e-12))
        else NA
        diag_env$loo_log[[length(diag_env$loo_log) + 1L]] <- list(
          path        = "generic_oda",
          attr_name   = attr_names[j],
          pos_id      = pos_id,
          row_hash    = rh_gen,
          n_obs       = length(idx),
          loo_mode    = if (is.null(loo_out)) "failed" else "n_fold_or_algebraic",
          ess_loo     = wessl_d,
          stable      = loo_stable_flag,
          elapsed_sec = loo_elapsed
        )
      }
      if (identical(loo_arg, "stable")) {
        if (is.null(loo_out) || !isTRUE(loo_out$allowed)) {
          .vmsg("  -> rejected: LOO not applicable (", elapsed, "s)")
          .log_gen(p_mc, ess, FALSE, "loo_unstable")
          return(NULL)
        }
        # STABLE: WESSL must equal WESS (tolerance matching unioda_core canon).
        if (!isTRUE(all.equal(ess, loo_out$ess_loo %||% NA_real_,
                              tolerance = 1e-12))) {
          .vmsg("  -> rejected: LOO not stable (", elapsed, "s)")
          .log_gen(p_mc, ess, FALSE, "loo_unstable")
          return(NULL)
        }
      } else if (is.numeric(loo)) {
        loo_pv <- if (!is.null(loo_out)) loo_out$p_value %||% NA_real_ else NA_real_
        if (is.na(loo_pv) || loo_pv >= loo) {
          .vmsg("  -> rejected: LOO p=", round(loo_pv %||% NA_real_, 4),
                " (", elapsed, "s)")
          .log_gen(p_mc, ess, FALSE, "loo_p_fail")
          return(NULL)
        }
      }
      fit$loo <- loo_out  # attach for .loo_info() in .split_nd()
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
  # Fills env$nodes (list indexed by canonical integer node_id).
  # node_id must follow Appendix C binary-tree geometry: left child of X = 2X,
  # right child of X = 2X+1.  Children are computed from node_id, not a counter.
  # pred_class: class assigned by the parent ODA rule for the branch leading here.
  .ho_grow_canon <- function(idx, depth, node_id, parent_id, env,
                             pred_class = NULL) {
    nid <- node_id
    nd  <- .leaf_nd(nid, parent_id, depth, idx, pred_class = pred_class)
    y_node <- y[idx]
    n_obs  <- length(idx)

    if (length(unique(y_node)) < 2L || n_obs < mindenom || depth >= max_depth) {
      env$nodes[[nid]] <- nd
      return(invisible(NULL))
    }

    fast <- .fast_screen(idx, pos_label = nid)
    if (length(fast) == 0L) {
      env$nodes[[nid]] <- nd
      return(invisible(NULL))
    }

    # HO-CTA greedy: accept the first (highest fast-ESS) candidate passing full MC+LOO
    best <- NULL
    for (fp in fast) {
      cand <- .full_fit_one(fp$j, idx, pos_id = nid)
      if (!is.null(cand)) { best <- cand; break }
    }

    if (is.null(best)) {
      env$nodes[[nid]] <- nd
      return(invisible(NULL))
    }

    appl <- .apply_cand(best, idx)
    nd   <- .split_nd(nid, parent_id, depth, idx, best, appl)

    # Canonical child IDs: left = 2*nid, right = 2*nid+1  (Appendix C, Figure C.1).
    # For ordered_cut rules (including normalized binary_map 0/1):
    #   left (2*nid) = x <= cut_value  (matches EXE INCLUDE x<=cut -> left-child).
    #   right (2*nid+1) = x > cut_value.
    # For other rule types: left = first class by sorted label, right = second.
    left_id    <- 2L * nid
    right_id   <- 2L * nid + 1L
    valid_mask <- !is.na(appl$y_pred)
    child_ids  <- integer(0L)

    if (identical(appl$rule$type, "ordered_cut")) {
      x_j_vals <- X[[best$j]][idx]
      sides    <- oda_rule_side(x_j_vals, appl$rule)   # 0 = x<=cut, 1 = x>cut
      for (sv in c(0L, 1L)) {
        child_idx <- idx[valid_mask & sides == sv]
        if (length(child_idx) == 0L) next
        cid_k    <- if (sv == 0L) left_id else right_id
        pred_cls <- appl$y_pred[valid_mask & sides == sv][1L]
        .ho_grow_canon(child_idx, depth + 1L, cid_k, nid, env,
                       pred_class = as.integer(pred_cls))
        child_ids <- c(child_ids, cid_k)
      }
    } else {
      sl <- sort(unique(appl$y_pred[valid_mask]))
      for (k in seq_along(sl)) {
        child_idx <- idx[valid_mask & appl$y_pred == sl[k]]
        if (length(child_idx) == 0L) next
        cid_k <- if (k == 1L) left_id else right_id
        .ho_grow_canon(child_idx, depth + 1L, cid_k, nid, env,
                       pred_class = as.integer(sl[k]))
        child_ids <- c(child_ids, cid_k)
      }
    }
    nd$child_ids <- child_ids
    env$nodes[[nid]] <- nd
    return(invisible(NULL))
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
          # Binary routing: for ordered_cut use rule_side (0=<=cut -> child 1,
          # 1=>cut -> child 2), matching EXE's INCLUDE x<=cut -> left-child.
          # For other rule types use oda_rule_predict (class-label order).
          y_hat <- if (!is.null(rule$type) && rule$type == "ordered_cut")
            as.integer(oda_rule_side(x_val, rule))
          else
            as.integer(oda_rule_predict(x_val, rule))
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
    if (prune_alpha >= 1.0) return(list(nodes = nodes_list, removed_ids = integer(0)))

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

    current       <- nodes_list
    removed_ids   <- integer(0)
    removed_attrs <- character(0)

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

      # Capture attribute BEFORE mutating the node to a leaf.
      attr_before <- current[[best_nid]]$attribute %||% NA_character_

      # BFS from best_nid to collect all descendant split IDs that will become
      # unreachable after this collapse.  Includes best_nid itself as position [1].
      all_under    <- .bfs_ids(current, best_nid)
      desc_splits  <- all_under[vapply(all_under, function(did) {
        did != best_nid && {
          nd <- current[[did]]
          !is.null(nd) && !isTRUE(nd$leaf) && length(nd$child_ids) > 0L
        }
      }, logical(1L))]

      current[[best_nid]]$leaf      <- TRUE
      current[[best_nid]]$child_ids <- integer(0)

      # Track: direct collapse root + any now-unreachable descendant splits.
      removed_ids   <- c(removed_ids, best_nid, desc_splits)
      removed_attrs <- c(removed_attrs, attr_before,
                         rep(NA_character_, length(desc_splits)))

      .vmsg("[PRUNE] node ", best_nid, " (", attr_before, "): WESS ",
            round(current_wess, 2), "% -> ", round(best_wess_p, 2), "%")
    }

    list(nodes = current, removed_ids = removed_ids, removed_attrs = removed_attrs)
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

  # ---- ENUMERATE: EO-CTA top-3-node enumeration --------------------------------
  # EO-CTA canonical (CTA.exe trace / MPE Ch.11):
  #   For each root A, build exactly ONE expanded candidate:
  #     B = best-ESS valid split on the left branch  (trace-selected node 2), or leaf
  #     C = best-ESS valid split on the right branch (trace-selected node 3), or leaf
  #   Grow HO-CTA below B and C.  Score the full pruned tree.
  #
  # Full A×B×C outer-product enumeration is NOT canonical.  The EXE trace for
  # CTA_DEMO MODEL8 Tree 4 (A=V6) proves this: V4 (ESS=31.86%) is evaluated at
  # node 2 but NOT selected; only V2 (ESS=48.65%) is used for the expanded
  # candidate.  Evaluating B=V4 produces a 69.99% tree EXE never reports.
  #
  # (B=leaf, C=leaf) is excluded from this expanded loop; handled by the
  # root-only stump phase below under path-local scoring.

  best_wess                  <- -Inf
  best_nodes                 <- NULL
  best_n_nodes               <- 0L
  best_confusion             <- NULL
  best_prune_removed_ids     <- integer(0)
  best_prune_removed_attrs   <- character(0)
  best_unpruned_n_splits     <- 0L
  best_unpruned_ess          <- -Inf

  .vmsg("[ENUMERATE] ", length(root_cands), " root candidates (A\u00d7B\u00d7C)")

  for (i in seq_along(root_cands)) {
    A_cand <- root_cands[[i]]
    .vmsg("[ENUMERATE ", i, "/", length(root_cands), "] root=", A_cand$name,
          " ESS=", round(A_cand$ess, 2), "%")

    appl_A  <- .apply_cand(A_cand, seq_len(n))
    valid_A <- !is.na(appl_A$y_pred)
    sl_A    <- sort(unique(appl_A$y_pred[valid_A]))
    if (length(sl_A) < 2L) next

    # EXE canonical: for ordered_cut left = x<=cut, right = x>cut.
    if (identical(appl_A$rule$type, "ordered_cut")) {
      sides_A   <- oda_rule_side(X[[A_cand$j]], appl_A$rule)
      left_idx  <- which(valid_A & sides_A == 0L)
      right_idx <- which(valid_A & sides_A == 1L)
    } else {
      left_idx  <- which(valid_A & appl_A$y_pred == sl_A[1L])
      right_idx <- which(valid_A & appl_A$y_pred == sl_A[2L])
    }

    # All valid split candidates for B (left) and C (right) positions.
    B_cands <- .all_cands(left_idx,  pos_id = 2L)
    C_cands <- .all_cands(right_idx, pos_id = 3L)

    # EO-CTA: one expanded candidate per root A.
    # B = best-ESS valid split on left branch, or leaf if none.
    # C = best-ESS valid split on right branch, or leaf if none.
    B_options <- if (length(B_cands) > 0L) list(B_cands[[1L]]) else list(NULL)
    C_options <- if (length(C_cands) > 0L) list(C_cands[[1L]]) else list(NULL)

    .vmsg("  B: ", if (length(B_cands) > 0L) B_cands[[1L]]$name else "leaf")
    .vmsg("  C: ", if (length(C_cands) > 0L) C_cands[[1L]]$name else "leaf")

    for (bi in seq_along(B_options)) {
      B_opt     <- B_options[[bi]]
      B_is_leaf <- is.null(B_opt)

      # ---- Build B sub-tree (cached for all C at this (A, B)) -----------------
      if (B_is_leaf) {
        nd2 <- .leaf_nd(2L, 1L, 2L, left_idx, pred_class = as.integer(sl_A[1L]))
      } else {
        B_cand <- B_opt
        appl_B <- .apply_cand(B_cand, left_idx)
        sl_B   <- sort(unique(appl_B$y_pred[!is.na(appl_B$y_pred)]))
        if (length(sl_B) < 2L) next   # degenerate B split - skip

        nd2 <- .split_nd(2L, 1L, 2L, left_idx, B_cand, appl_B)

        # Grow HO-CTA below B's children (depth 3+); node 2 is already placed.
        # Canonical: B is at node 2; left child = 4 (2×2), right child = 5 (2×2+1).
        B_sub_env       <- new.env(parent = emptyenv())
        B_sub_env$nodes <- list()

        if (identical(appl_B$rule$type, "ordered_cut")) {
          sides_B     <- oda_rule_side(X[[B_cand$j]][left_idx], appl_B$rule)
          B_left_idx  <- left_idx[!is.na(appl_B$y_pred) & sides_B == 0L]
          B_right_idx <- left_idx[!is.na(appl_B$y_pred) & sides_B == 1L]
        } else {
          B_left_idx  <- left_idx[!is.na(appl_B$y_pred) & appl_B$y_pred == sl_B[1L]]
          B_right_idx <- left_idx[!is.na(appl_B$y_pred) & appl_B$y_pred == sl_B[2L]]
        }

        .ho_grow_canon(B_left_idx,  3L, 4L, 2L, B_sub_env,
                       pred_class = as.integer(sl_B[1L]))
        .ho_grow_canon(B_right_idx, 3L, 5L, 2L, B_sub_env,
                       pred_class = as.integer(sl_B[2L]))
        nd2$child_ids <- c(4L, 5L)
      }

      # ---- Inner C loop -------------------------------------------------------
      for (ci in seq_along(C_options)) {
        C_opt     <- C_options[[ci]]
        C_is_leaf <- is.null(C_opt)

        # (B=leaf, C=leaf) handled by root-only stump phase - skip here
        if (B_is_leaf && C_is_leaf) next

        # Canonical: C is always the right child of the root (node 3 = 2×1+1).
        c_root_id <- 3L

        if (C_is_leaf) {
          nd_C <- .leaf_nd(3L, 1L, 2L, right_idx,
                           pred_class = as.integer(sl_A[2L]))
        } else {
          C_cand <- C_opt
          appl_C <- .apply_cand(C_cand, right_idx)
          sl_C   <- sort(unique(appl_C$y_pred[!is.na(appl_C$y_pred)]))
          if (length(sl_C) < 2L) next   # degenerate C split - skip

          nd_C <- .split_nd(3L, 1L, 2L, right_idx, C_cand, appl_C)

          # Grow HO-CTA below C's children; C is at node 3.
          # Canonical: left child = 6 (2×3), right child = 7 (2×3+1).
          C_sub_env       <- new.env(parent = emptyenv())
          C_sub_env$nodes <- list()

          if (identical(appl_C$rule$type, "ordered_cut")) {
            sides_C     <- oda_rule_side(X[[C_cand$j]][right_idx], appl_C$rule)
            C_left_idx  <- right_idx[!is.na(appl_C$y_pred) & sides_C == 0L]
            C_right_idx <- right_idx[!is.na(appl_C$y_pred) & sides_C == 1L]
          } else {
            C_left_idx  <- right_idx[!is.na(appl_C$y_pred) & appl_C$y_pred == sl_C[1L]]
            C_right_idx <- right_idx[!is.na(appl_C$y_pred) & appl_C$y_pred == sl_C[2L]]
          }

          .ho_grow_canon(C_left_idx,  3L, 6L, 3L, C_sub_env,
                         pred_class = as.integer(sl_C[1L]))
          .ho_grow_canon(C_right_idx, 3L, 7L, 3L, C_sub_env,
                         pred_class = as.integer(sl_C[2L]))
          nd_C$child_ids <- c(6L, 7L)
        }

        B_label <- if (B_is_leaf) "leaf" else B_opt$name
        C_label <- if (C_is_leaf) "leaf" else C_opt$name
        .vmsg("  [B=", B_label, " C=", C_label, "]")

        # ---- Assemble full candidate node list --------------------------------
        # Sparse list indexed by canonical node ID (Appendix C geometry).
        # root=1, B=2, C=3, B-subtree={4,5,8,9,...}, C-subtree={6,7,12,13,...}
        cand_nodes <- list()

        nd1 <- .split_nd(1L, 0L, 1L, seq_len(n), A_cand, appl_A)
        nd1$split_labels <- as.integer(sl_A)
        nd1$child_ids    <- c(2L, 3L)   # canonical: left=2, right=3
        cand_nodes[[1L]] <- nd1
        cand_nodes[[2L]] <- nd2
        cand_nodes[[3L]] <- nd_C

        if (!B_is_leaf) {
          for (bid in seq_along(B_sub_env$nodes)) {
            nb <- B_sub_env$nodes[[bid]]
            if (!is.null(nb)) cand_nodes[[bid]] <- nb
          }
        }

        if (!C_is_leaf) {
          for (cv in seq_along(C_sub_env$nodes)) {
            nc <- C_sub_env$nodes[[cv]]
            if (!is.null(nc)) cand_nodes[[cv]] <- nc
          }
        }

        # ---- Prune, score, compare -------------------------------------------
        # Capture actual split count from cand_nodes BEFORE pruning.
        unpruned_splits_cand <- sum(vapply(cand_nodes, function(nd)
          !is.null(nd) && !isTRUE(nd$leaf) && length(nd$child_ids) > 0L,
          logical(1L)))

        prune_result        <- .prune_tree(cand_nodes)
        pruned_nodes        <- prune_result$nodes
        removed_cand        <- prune_result$removed_ids
        removed_attrs_cand  <- prune_result$removed_attrs

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
          # Capture PRUNE provenance for the winning candidate.
          best_prune_removed_ids   <- removed_cand
          best_prune_removed_attrs <- removed_attrs_cand
          best_unpruned_n_splits   <- unpruned_splits_cand
          best_unpruned_ess <- tryCatch({
            p_before <- .predict_all(cand_nodes, root_id = 1L)
            .wess_classes(y, p_before, w)
          }, error = function(e) wess_cand)
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
    # EXE canonical: for ordered_cut left = x<=cut, right = x>cut.
    if (identical(appl_A$rule$type, "ordered_cut")) {
      sides_A   <- oda_rule_side(X[[A_cand$j]], appl_A$rule)
      left_idx  <- which(valid_A & sides_A == 0L)
      right_idx <- which(valid_A & sides_A == 1L)
    } else {
      left_idx  <- which(valid_A & appl_A$y_pred == sl_A[1L])
      right_idx <- which(valid_A & appl_A$y_pred == sl_A[2L])
    }
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
      # pred_class for each leaf: actual majority class on that side.
      # For ordered_cut, left = x<=cut and right = x>cut (canonical); the
      # predicted class for each side comes from the ODA rule applied to that idx.
      left_pred  <- appl_A$y_pred[left_idx[1L]]
      right_pred <- appl_A$y_pred[right_idx[1L]]
      stump_nodes[[2L]] <- .leaf_nd(2L, 1L, 2L, left_idx,
                                    pred_class = as.integer(left_pred))
      stump_nodes[[3L]] <- .leaf_nd(3L, 1L, 2L, right_idx,
                                    pred_class = as.integer(right_pred))

      best_wess    <- wess_stump
      best_nodes   <- stump_nodes
      best_n_nodes <- 3L
      # Capture final-tree confusion at the moment of selection.
      # stump_preds[ok] is path-local: missing-root obs already excluded.
      best_confusion             <- .make_training_conf(y[ok], stump_preds[ok])
      # Stumps are not pruned: unpruned == pruned; no removed nodes.
      best_prune_removed_ids     <- integer(0)
      best_prune_removed_attrs   <- character(0)
      best_unpruned_n_splits     <- 1L   # stump has exactly one split: the root
      best_unpruned_ess          <- wess_stump
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
      training_confusion = if (no_tree) NULL else best_confusion,
      # PRUNE provenance: what pruning did (or did not do) to the winning
      # candidate.  NULL for no-tree fits.  Fields:
      #   unpruned_n_splits  -- split-node count before pruning
      #   pruned_n_splits    -- split-node count after pruning (final tree)
      #   removed_node_ids   -- node IDs collapsed to leaves by PRUNE
      #   removed_attrs      -- attribute name at each removed node
      #   unpruned_ess       -- ESS/WESS of winning candidate before pruning
      #   pruned_ess         -- ESS/WESS of winning candidate after pruning
      #                         (equals overall_ess)
      prune_info = if (no_tree) NULL else {
        pruned_n_splits <- sum(vapply(best_nodes, function(nd)
          !is.null(nd) && !isTRUE(nd$leaf) && length(nd$child_ids) > 0L,
          logical(1L)))
        list(
          # unpruned_n_splits: actual split count captured from cand_nodes BEFORE
          # .prune_tree() ran.  For stumps this equals 1L.
          unpruned_n_splits = as.integer(best_unpruned_n_splits),
          pruned_n_splits   = as.integer(pruned_n_splits),
          # removed_node_ids: directly collapsed node IDs plus any descendant split
          # IDs that became unreachable (whole-subtree collapse).
          removed_node_ids  = best_prune_removed_ids,
          # removed_attrs: attribute names captured BEFORE mutation; NA for
          # descendant splits whose root was collapsed.
          removed_attrs     = best_prune_removed_attrs,
          unpruned_ess      = best_unpruned_ess,
          pruned_ess        = best_wess
        )
      }
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
        # Binary routing: use rule_side (0=x<=cut -> child 1, 1=x>cut -> child 2),
        # matching the tree-building convention in .ho_grow_canon which assigns
        # left_id to sv=0 (x<=cut) and right_id to sv=1 (x>cut), regardless of
        # rule direction.  oda_rule_predict would flip the index for "1->0" rules,
        # routing obs to the wrong branch.
        y_hat <- if (!is.null(rule$type) && rule$type == "ordered_cut")
          as.integer(oda_rule_side(x_val, rule))
        else
          as.integer(oda_rule_predict(x_val, rule))
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
    cat("No tree found (leaf-only): no valid split passed significance, LOO gate, and MINDENOM constraints.\n")
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

  # Terminal endpoints section  -  uses cta_endpoint_table() as the source.
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
