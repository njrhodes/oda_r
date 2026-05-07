###############################################################################
# R/cta_core.R — Classification Tree Analysis (CTA)
#
# Implements EO-CTA (Enumerated Optimal CTA), corresponding to MegaODA's
# ENUMERATE command.
#
# Algorithm (Yarnold & Soltysik 2016, Chapters 10-11):
#
#   ENUMERATE: Evaluate every valid (root A × left-child B × right-child C)
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
#' @param loo LOO mode: "off", "stable", or "pvalue".
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
    K_segments  = NULL
) {
  # ---- Validate and prep ----------------------------------------------------
  if (!is.data.frame(X)) X <- as.data.frame(X)
  n <- nrow(X)
  stopifnot(length(y) == n)
  y <- as.integer(y)
  if (is.null(w)) w <- rep(1.0, n) else {
    w <- as.numeric(w)
    stopifnot(length(w) == n, all(w >= 0))
  }
  if (is.null(attr_names))
    attr_names <- colnames(X) %||% paste0("V", seq_len(ncol(X)))
  n_attrs <- ncol(X)
  C       <- length(unique(y))
  loo_arg <- if (is.numeric(loo)) "pvalue" else as.character(loo)

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
    list(
      node_id = nid, parent_id = parent_id, depth = depth,
      n_obs = length(idx), n_weighted = sum(w[idx]),
      leaf = TRUE, majority_class = maj,
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

  # 2×2 binary confusion matrix
  .conf_binary <- function(y_node, y_pred_raw) {
    conf <- oda_confusion_binary(y_node, y_pred_raw)
    matrix(c(conf$TN, conf$FP, conf$FN, conf$TP), nrow = 2L, ncol = 2L,
           dimnames = list(actual = c("0","1"), predicted = c("0","1")))
  }

  # LOO metadata from a fit result
  .loo_info <- function(fit) {
    lr <- fit$loo
    list(
      status  = if (!is.null(lr) && isTRUE(lr$allowed)) "STABLE" else "OFF",
      ess_loo = if (!is.null(lr)) lr$ess_loo %||% NA_real_ else NA_real_,
      p_value = if (!is.null(lr)) lr$p_value %||% NA_real_ else NA_real_
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
      p_mc = cand$p_mc, loo_status = li$status,
      loo_ess = li$ess_loo, loo_p = li$p_value,
      confusion = conf,
      split_labels = as.integer(sl),
      child_ids = integer(0)
    )
  }

  # Fast screen: mcarlo=FALSE for all attributes at idx.
  # Returns list of list(j, ess), sorted by ess descending.
  .fast_screen <- function(idx) {
    if (length(idx) < mindenom) return(list())
    y_n <- y[idx]; w_n <- w[idx]
    if (length(unique(y_n)) < 2L) return(list())
    result <- list()
    for (j in seq_len(n_attrs)) {
      x_j <- X[[j]][idx]
      x_v <- x_j[!is.na(x_j)]
      if (!is.null(miss_codes)) x_v <- x_v[!(x_v %in% miss_codes)]
      if (length(x_v) < 2L || length(unique(x_v)) < 2L) next
      fit_f <- tryCatch(
        oda_fit(x = x_j, y = y_n, w = w_n, priors_on = priors_on,
                miss_codes = miss_codes, K_segments = K_segments,
                mcarlo = FALSE, loo = "off", mindenom = mindenom),
        error = function(e) list(ok = FALSE))
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
    seed_j <- if (is.null(mc_seed)) NULL else
      (as.integer(mc_seed) + as.integer(pos_id) * 100L + j) %% .Machine$integer.max
    fit <- tryCatch(
      oda_fit(x = x_j, y = y_n, w = w_n, priors_on = priors_on,
              miss_codes = miss_codes, K_segments = K_segments,
              mcarlo = TRUE, mc_iter = as.integer(mc_iter),
              mc_target = mc_target, mc_stop = mc_stop,
              mc_stopup = mc_stopup, mc_seed = seed_j,
              loo = loo_arg, eval_order = "loo_then_mc",
              mindenom = mindenom),
      error = function(e) list(ok = FALSE))
    if (!isTRUE(fit$ok))                                            return(NULL)
    p_mc <- fit$p_mc
    if (is.na(p_mc) || p_mc >= alpha_split || p_mc >= prune_alpha) return(NULL)
    ess <- .get_ess(fit)
    if (is.na(ess) || ess < ess_min)                               return(NULL)
    if (identical(loo_arg, "stable")) {
      if (is.null(fit$loo) || !isTRUE(fit$loo$allowed))           return(NULL)
    } else if (is.numeric(loo)) {
      loo_pv <- fit$loo$p_value %||% NA_real_
      if (is.na(loo_pv) || loo_pv >= loo)                         return(NULL)
    }
    list(j = j, name = attr_names[j], fit = fit, ess = ess, p_mc = p_mc)
  }

  # Get ALL valid candidates at a node position (fast screen then ALL full MC+LOO).
  # pos_id: unique position for seed derivation (1=root,2=left,3=right).
  .all_cands <- function(idx, pos_id) {
    fast <- .fast_screen(idx)
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

    fast <- .fast_screen(idx)
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
        if (miss_here) return(nd$majority_class)
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
          # Binary: y_hat is rule-side 0/1; route by position (0→child 1, 1→child 2).
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

  # ---- ENUMERATE: top-3-node enumeration ------------------------------------

  root_cands <- .all_cands(seq_len(n), pos_id = 1L)

  # No valid root → single-leaf tree
  if (length(root_cands) == 0L) {
    nodes_out <- list(.leaf_nd(1L, 0L, 1L, seq_len(n)))
    return(structure(
      list(nodes = nodes_out, root_id = 1L, n_nodes = 1L, n = n, C = C,
           attr_names = attr_names, n_attrs = n_attrs, miss_codes = miss_codes,
           priors_on = priors_on, alpha_split = alpha_split,
           mindenom = mindenom, prune_alpha = prune_alpha,
           max_depth = max_depth, ess_min = ess_min, loo = loo),
      class = "cta_tree"))
  }

  # ---- ENUMERATE: for each valid root candidate, grow full HO-CTA, pick best ----
  # Canonical algorithm: try each root attribute (sorted by root-level ESS desc),
  # grow a complete HO-CTA tree below each, compute full-tree WESS, return best.
  # Ties broken by enumeration order (highest root ESS first).

  best_wess    <- -Inf
  best_nodes   <- NULL
  best_n_nodes <- 0L

  for (A_cand in root_cands) {
    appl_A  <- .apply_cand(A_cand, seq_len(n))
    valid_A <- !is.na(appl_A$y_pred)
    sl_A    <- sort(unique(appl_A$y_pred[valid_A]))
    if (length(sl_A) < 2L) next

    left_idx  <- which(valid_A & appl_A$y_pred == sl_A[1L])
    right_idx <- which(valid_A & appl_A$y_pred == sl_A[2L])

    # Grow full HO-CTA on left branch; first id allocated = 2
    B_env         <- new.env(parent = emptyenv())
    B_env$nodes   <- vector("list", n + 4L)
    B_env$counter <- 1L
    cid2  <- .ho_grow(left_idx, 2L, 1L, B_env,
                      pred_class = as.integer(sl_A[1L]))
    B_max <- B_env$counter

    # Grow full HO-CTA on right branch; first id allocated = B_max + 1
    C_env         <- new.env(parent = emptyenv())
    C_env$nodes   <- vector("list", n + B_max + 4L)
    C_env$counter <- B_max
    cid3  <- .ho_grow(right_idx, 2L, 1L, C_env,
                      pred_class = as.integer(sl_A[2L]))
    C_max <- C_env$counter

    # Build candidate node list
    total_size <- C_max
    cand_nodes <- vector("list", total_size)

    nd1 <- .split_nd(1L, 0L, 1L, seq_len(n), A_cand, appl_A)
    nd1$split_labels <- as.integer(sl_A)
    nd1$child_ids    <- c(cid2, cid3)
    cand_nodes[[1L]] <- nd1

    for (bid in seq_len(B_max)[-1L]) {
      nd_b <- B_env$nodes[[bid]]
      if (!is.null(nd_b)) cand_nodes[[bid]] <- nd_b
    }
    for (cid_v in seq(cid3, C_max)) {
      nd_c <- C_env$nodes[[cid_v]]
      if (!is.null(nd_c)) cand_nodes[[cid_v]] <- nd_c
    }

    # Evaluate full-tree WESS
    wess_cand <- tryCatch({
      preds <- .predict_all(cand_nodes, root_id = 1L)
      .wess_classes(y, preds, w)
    }, error = function(e) -Inf)

    if (wess_cand > best_wess) {
      best_wess    <- wess_cand
      best_nodes   <- cand_nodes
      best_n_nodes <- sum(!vapply(cand_nodes, is.null, logical(1L)))
    }
  }  # end A_cand loop

  # Fallback: enumeration produced no valid tree
  if (is.null(best_nodes)) {
    best_nodes   <- list(.leaf_nd(1L, 0L, 1L, seq_len(n)))
    best_n_nodes <- 1L
  }

  structure(
    list(
      nodes       = best_nodes,
      root_id     = 1L,
      n_nodes     = best_n_nodes,
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
      loo         = loo
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
predict.cta_tree <- function(object, newdata, ...) {
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
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
      if (miss_here) return(nd$majority_class)

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
        # Binary: y_hat is rule-side 0/1; route by position (0→child 1, 1→child 2).
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

    if (!is.null(nd$confusion)) .print_node_confusion(nd$confusion)
    cat("\n")
  }

  n_split <- sum(vapply(x$nodes, function(nd) !isTRUE(nd$leaf), logical(1)))
  n_leaf  <- x$n_nodes - n_split
  cat(sprintf("Nodes: %d total  (%d split  %d leaf)\n",
              x$n_nodes, n_split, n_leaf))
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
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  do.call(rbind, rows)
}

# ---- Internal helpers -------------------------------------------------------

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
