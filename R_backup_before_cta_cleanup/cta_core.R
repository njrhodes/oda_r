###############################################################################
# R/cta_core.R — Classification Tree Analysis (CTA)
#
# CTA is recursive ODA. Each node:
#   1. Applies oda_fit() to every attribute on the node's population.
#   2. Selects the attribute with highest ESS that passes alpha_split and
#      has n_obs >= mindenom.  (ENUMERATE: all attributes evaluated every node)
#   3. Splits the population by predicted class; recurses into each child.
#
# MegaODA command → parameter mapping:
#   MINDENOM n    → mindenom      (raw obs count, not weighted)
#   PRUNE 0.05    → prune_alpha   (don't grow branches with p >= prune_alpha;
#                                  set to 1.0 for unpruned tree)
#   ENUMERATE     → always TRUE here (we always evaluate all attributes)
#   LOO STABLE    → loo = "stable"
#   WEIGHT V5     → w argument
#   MISSING (-9)  → miss_codes
#   MC ITER 5000 CUTOFF 0.05 STOP 99.9 → mc_iter, alpha_split, mc_stop
#
# Three "tree types" in MegaODA output (Unpruned / Pruned / Enumerated)
# all grow the same way here — differentiate only by prune_alpha:
#   Unpruned   = prune_alpha 1.0  (keep every significant split)
#   Pruned     = prune_alpha 0.05 (drop branches with p >= prune_alpha)
#   Enumerated = same as unpruned but every attribute evaluated (always true here)
#
# Key design decisions:
#   - oda_fit() handles LOO internally.  When loo="stable" and the LOO ESS
#     differs from training ESS, oda_fit() returns ok=FALSE with reason
#     "loo_not_stable". CTA treats that as a non-split at that attribute.
#   - mindenom is a RAW OBSERVATION COUNT (not weighted sum), matching
#     MegaODA's OBS column.
#   - Class labels in y are preserved as-is; oda_fit() recodes internally
#     per attribute fit.  CTA does NOT recode y globally.
#   - ESS for best-attribute selection: binary engine stores in $ess;
#     multiclass in $ess_pac.  Both normalised to [0,100] with 0 = chance.
#   - split_labels stored on each split node for prediction routing.
###############################################################################

#' Fit a Classification Tree Analysis (CTA) model
#'
#' Builds a classification tree by recursively applying ODA at each node.
#' Every attribute is evaluated at every node (ENUMERATE). The attribute
#' with the highest ESS that passes significance is selected for the split.
#'
#' @param X Data frame or matrix of attribute columns.
#' @param y Class variable vector.
#' @param w Optional numeric case weights (MegaODA WEIGHT). Default: unit weights.
#' @param priors_on Logical; use prior-odds weighting at each node. Default TRUE.
#' @param miss_codes Numeric vector of missing-value codes (MegaODA MISSING ALL).
#' @param alpha_split p threshold to split a node (MegaODA MC CUTOFF). Default 0.05.
#' @param mindenom Minimum raw observation count to attempt a split
#'   (MegaODA MINDENOM). Default 5.
#' @param prune_alpha Don't grow branches with p >= prune_alpha (MegaODA PRUNE).
#'   Use 1.0 (default) for unpruned; use alpha_split for pruned.
#' @param max_depth Maximum depth. Root is depth 1. Default 10.
#' @param ess_min Minimum ESS to split. Default 0.
#' @param mc_iter Maximum MC iterations per node fit. Default 25000.
#' @param mc_target MC significance threshold for early stopping. Default 0.05.
#' @param mc_stop Confidence level (percent) for STOP (lower tail). Default 99.9.
#' @param mc_stopup Confidence level (percent) for STOPUP. Default 20.
#' @param mc_seed Integer base seed. Each node*attribute pair gets a unique seed
#'   derived from this. Default NULL (unseeded).
#' @param loo LOO mode passed to oda_fit at every node: "off" (default),
#'   "stable" (MegaODA LOO STABLE — only split if LOO ESS = training ESS),
#'   or "pvalue".
#' @param attr_names Optional character vector of attribute names.
#' @param K_segments Segments for multiclass ordered splits. Default = C.
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

  # Preserve class labels exactly; oda_fit() recodes internally per fit
  y <- as.integer(y)
  if (is.null(w)) w <- rep(1.0, n) else {
    w <- as.numeric(w)
    stopifnot(length(w) == n, all(w >= 0))
  }

  if (is.null(attr_names)) {
    attr_names <- colnames(X) %||% paste0("V", seq_len(ncol(X)))
  }
  n_attrs <- ncol(X)
  C       <- length(unique(y))

  # ---- Node store -----------------------------------------------------------
  nodes        <- vector("list", 2L * n)  # pre-allocate (will trim)
  node_counter <- 0L

  # ---- Extract ESS from oda_fit result (handles binary and multiclass) ------
  .get_ess <- function(fit) {
    # Binary engine: $ess is the objective-space ESS (priors_on or not)
    # Multiclass engine: $ess_pac is ESS from mean PAC, $ess_pv from mean PV
    # Both are scaled 0-100 with 0 = chance baseline.
    if (!is.null(fit$ess)     && is.finite(fit$ess))     return(fit$ess)
    if (!is.null(fit$ess_pac) && is.finite(fit$ess_pac)) return(fit$ess_pac)
    NA_real_
  }

  # ---- Internal recursive node builder -------------------------------------
  fit_node <- function(idx, parent_id, depth) {
    node_counter <<- node_counter + 1L
    nid <- node_counter

    y_node <- y[idx]
    w_node <- w[idx]
    n_obs  <- length(idx)       # raw count — used for mindenom
    n_wt   <- sum(w_node)       # weighted count — reported as n_weighted

    # Weighted majority class
    ctab <- tapply(w_node, y_node, sum)
    majority <- as.integer(names(which.max(ctab)))

    # Base node record
    nd <- list(
      node_id       = nid,
      parent_id     = parent_id,
      depth         = depth,
      n_obs         = n_obs,
      n_weighted    = n_wt,
      leaf          = TRUE,
      majority_class = majority,
      attribute     = NA_character_,
      attr_col      = NA_integer_,
      attr_type     = NA_character_,
      rule          = NULL,
      ess           = NA_real_,
      ess_weighted  = NA_real_,
      p_mc          = NA_real_,
      loo_status    = NA_character_,
      loo_ess       = NA_real_,
      loo_p         = NA_real_,
      confusion     = NULL,
      split_labels  = integer(0),
      child_ids     = integer(0)
    )

    # ---- Stopping rules -------------------------------------------------------
    # Pure node, too small, or at max depth → leaf immediately
    if (length(unique(y_node)) == 1L ||
        n_obs < mindenom ||
        depth >= max_depth) {
      nodes[[nid]] <<- nd
      return(nid)
    }

    # ---- Attribute selection (ENUMERATE = always try all) -------------------
    best_ess  <- -Inf
    best      <- NULL

    for (j in seq_len(n_attrs)) {
      x_j <- X[[j]][idx]

      # Skip if this attribute is constant or all-missing in this node
      x_valid <- x_j[!is.na(x_j)]
      if (length(x_valid) < 2L || length(unique(x_valid)) < 2L) next

      # Unique seed per (node, attribute) for reproducibility
      seed_j <- if (is.null(mc_seed)) NULL else
        (as.integer(mc_seed) + nid * 100L + j) %% .Machine$integer.max

      # Run with loo=loo so LOO STABLE is enforced per-attribute.
      # If loo="stable" and the solution is not stable, oda_fit() returns
      # ok=FALSE (reason="loo_not_stable") and the attribute is skipped.
      fit_j <- tryCatch(
        oda_fit(
          x          = x_j,
          y          = y_node,
          w          = w_node,
          priors_on  = priors_on,
          miss_codes = miss_codes,
          K_segments = K_segments,
          mcarlo     = TRUE,
          mc_iter    = as.integer(mc_iter),
          mc_target  = mc_target,
          mc_stop    = mc_stop,
          mc_stopup  = mc_stopup,
          mc_seed    = seed_j,
          loo        = loo
        ),
        error = function(e) list(ok = FALSE, reason = conditionMessage(e))
      )

      if (!isTRUE(fit_j$ok)) next   # no valid rule, or LOO STABLE rejected

      ess_j <- .get_ess(fit_j)
      if (is.na(ess_j) || ess_j < ess_min) next

      if (ess_j > best_ess + 1e-10) {
        best_ess <- ess_j
        best <- list(
          j    = j,
          name = attr_names[j],
          fit  = fit_j,
          ess  = ess_j,
          p_mc = fit_j$p_mc
        )
      }
    }

    # No LOO-stable candidate found — leaf
    if (is.null(best)) {
      nodes[[nid]] <<- nd
      return(nid)
    }

    # Alpha gate on the max-ESS winner (applied after LOO STABLE + max ESS)
    if (is.na(best$p_mc) || best$p_mc >= alpha_split ||
        best$p_mc >= prune_alpha) {
      nodes[[nid]] <<- nd
      return(nid)
    }

    # ---- Apply the winning split --------------------------------------------
    rule    <- best$fit$rule
    engine  <- best$fit$engine %||% "binary"
    x_best  <- X[[best$j]][idx]

    y_pred_node <- if (engine == "multiclass") {
      as.integer(oda_rule_predict_multiclass(
        x_best, rule,
        boundary = rule$boundary %||% "megaoda_halfopen"))
    } else {
      as.integer(oda_rule_predict(x_best, rule))
    }

    # ---- LOO status string --------------------------------------------------
    # When loo="stable", oda_fit() only returns ok=TRUE if the solution is
    # LOO stable (ESS_LOO == ESS_train). Extract the LOO result from the winner.
    loo_res     <- best$fit$loo
    loo_status  <- if (!is.null(loo_res) && isTRUE(loo_res$allowed)) "STABLE" else "OFF"
    loo_ess_val <- if (!is.null(loo_res)) loo_res$ess_loo  %||% NA_real_ else NA_real_
    loo_p_val   <- if (!is.null(loo_res)) loo_res$p_value  %||% NA_real_ else NA_real_

    # ---- Build confusion for this split node --------------------------------
    # Use raw counts (unit weights) so the matrix matches MegaODA OBS counts
    conf_node <- oda_confusion_binary(y_node, y_pred_node)     # binary
    if (engine == "multiclass") {
      conf_node <- best$fit$confusion    # already raw counts from oda_fit fix
    } else {
      # Convert binary TP/TN/FP/FN list to 2x2 matrix
      tp <- conf_node$TP; tn <- conf_node$TN
      fp <- conf_node$FP; fn <- conf_node$FN
      conf_node <- matrix(c(tn, fp, fn, tp),
                          nrow = 2, ncol = 2,
                          dimnames = list(actual    = c("0","1"),
                                         predicted = c("0","1")))
    }

    # ---- Update node record --------------------------------------------------
    nd$leaf         <- FALSE
    nd$attribute    <- best$name
    nd$attr_col     <- best$j
    nd$attr_type    <- best$fit$attr_type %||% "ordered"
    nd$rule         <- rule
    nd$ess          <- best$ess
    nd$ess_weighted <- best$ess   # diverges from ess when w != 1; see note below
    nd$p_mc         <- best$p_mc
    nd$loo_status   <- loo_status
    nd$loo_ess      <- loo_ess_val
    nd$loo_p        <- loo_p_val
    nd$confusion    <- conf_node

    # Note on ess_weighted: MegaODA's WESS is the ESS in the case-weighted
    # objective space (WEIGHT command). With unit weights ess == ess_weighted.
    # With non-unit case weights the binary priors-adjustment incorporates w,
    # so best$ess already reflects case weights. The distinction between ESS
    # and WESS in MegaODA output arises when LOO is used on weighted data;
    # full replication requires per-fold weighted ESS computation (future work).

    nodes[[nid]] <<- nd   # store before recursion so IDs are allocated

    # ---- Recurse into children -----------------------------------------------
    split_labels <- sort(unique(y_pred_node))
    child_ids    <- integer(length(split_labels))

    for (k in seq_along(split_labels)) {
      child_idx <- idx[y_pred_node == split_labels[k]]
      if (length(child_idx) == 0L) next
      cid          <- fit_node(child_idx, parent_id = nid, depth = depth + 1L)
      child_ids[k] <- cid
    }

    nd$split_labels <- as.integer(split_labels)
    nd$child_ids    <- child_ids[child_ids > 0L]
    nodes[[nid]] <<- nd   # write back with child_ids and split_labels

    return(nid)
  }

  # ---- Grow tree from root --------------------------------------------------
  root_id <- fit_node(seq_len(n), parent_id = 0L, depth = 1L)

  # Trim pre-allocated nulls
  nodes <- nodes[seq_len(node_counter)]

  structure(
    list(
      nodes       = nodes,
      root_id     = root_id,
      n_nodes     = node_counter,
      n           = n,
      C           = C,
      attr_names  = attr_names,
      n_attrs     = n_attrs,
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
#' Routes each observation down the tree until it reaches a leaf, then
#' returns that leaf's majority class.
#'
#' @param object A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @param newdata Data frame or matrix with the same attribute columns as
#'   training X (same column order and names).
#' @param ... Unused.
#' @return Integer vector of predicted class labels, length \code{nrow(newdata)}.
#' @export
predict.cta_tree <- function(object, newdata, ...) {
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  n_new <- nrow(newdata)

  predict_one <- function(i) {
    nid <- object$root_id
    repeat {
      nd <- object$nodes[[nid]]
      if (is.null(nd) || isTRUE(nd$leaf) || length(nd$child_ids) == 0L)
        return(nd$majority_class %||% NA_integer_)

      j     <- nd$attr_col
      x_val <- newdata[[j]][i]

      # Missing attribute at this node → route to majority class of current node
      if (is.na(x_val)) return(nd$majority_class)

      # Apply the rule
      rule <- nd$rule
      y_hat <- if (!is.null(rule$type) &&
                   rule$type %in% c("multiclass_ordered","multiclass_nominal")) {
        as.integer(oda_rule_predict_multiclass(
          x_val, rule,
          boundary = rule$boundary %||% "megaoda_halfopen"))
      } else {
        as.integer(oda_rule_predict(x_val, rule))
      }

      # Find matching child
      split_labels <- nd$split_labels
      idx <- which(split_labels == y_hat)
      if (length(idx) == 0L) return(nd$majority_class)
      nid <- nd$child_ids[idx[1L]]
      if (is.na(nid) || nid < 1L) return(nd$majority_class)
    }
  }

  vapply(seq_len(n_new), predict_one, integer(1L))
}

# ---- Print (MegaODA node table format) -------------------------------------

#' Print a CTA tree
#'
#' Displays each split node in MegaODA node-table format: attribute, node id,
#' depth, raw obs count, p-value, ESS, LOO status, and rule string, followed
#' by the node confusion matrix.
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
                nd$node_id,
                nd$depth,
                nd$n_obs,
                p_str,
                nd$ess        %||% 0,
                nd$ess_weighted %||% 0,
                loo_str,
                rule_str))

    if (!is.null(nd$confusion)) .print_node_confusion(nd$confusion)
    cat("\n")
  }

  n_split <- sum(vapply(x$nodes, function(nd) !isTRUE(nd$leaf), logical(1)))
  n_leaf  <- x$n_nodes - n_split
  cat(sprintf("Nodes: %d total  (%d split  %d leaf)\n",
              x$n_nodes, n_split, n_leaf))
  invisible(x)
}

# ---- Node table as data frame -----------------------------------------------

#' Extract the CTA node summary as a data frame
#'
#' Returns a tidy data frame with one row per node, mirroring the MegaODA
#' node table (ATTRIBUTE, NODE, LEV, OBS, p, ESS, LOO, TYP columns).
#'
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @return Data frame with columns: node_id, parent_id, depth, leaf,
#'   attribute, attr_type, n_obs, n_weighted, p_mc, ess, ess_weighted,
#'   loo_status, loo_ess.
#' @export
cta_node_table <- function(tree) {
  rows <- lapply(tree$nodes, function(nd) {
    if (is.null(nd)) return(NULL)
    data.frame(
      node_id      = nd$node_id,
      parent_id    = nd$parent_id,
      depth        = nd$depth,
      leaf         = isTRUE(nd$leaf),
      attribute    = nd$attribute     %||% NA_character_,
      attr_type    = nd$attr_type     %||% NA_character_,
      n_obs        = nd$n_obs,
      n_weighted   = round(nd$n_weighted, 4),
      p_mc         = nd$p_mc          %||% NA_real_,
      ess          = nd$ess           %||% NA_real_,
      ess_weighted = nd$ess_weighted  %||% NA_real_,
      loo_status   = nd$loo_status    %||% NA_character_,
      loo_ess      = nd$loo_ess       %||% NA_real_,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  do.call(rbind, rows)
}

# ---- Internal helpers -------------------------------------------------------

# Format a rule as a string matching MegaODA MODEL column.
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

# Print a confusion matrix in MegaODA layout.
.print_node_confusion <- function(conf) {
  if (is.table(conf)) conf <- as.matrix(conf)
  if (!is.matrix(conf) || nrow(conf) < 2L) return(invisible(NULL))

  classes <- rownames(conf) %||% as.character(seq_len(nrow(conf)))
  cw <- max(7L, max(nchar(classes)) + 2L)

  # Header row
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
