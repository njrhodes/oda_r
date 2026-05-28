###############################################################################
# R/cta_s3.R — S3 methods and read-only report tables for cta_tree
#
# Classes produced here:
#   cta_tree_summary — structured output of summary.cta_tree()
#
# Public API:
#   summary.cta_tree()
#   print.cta_tree_summary()
#   cta_endpoint_table()
#
# These functions are read-only: no refitting, no prediction calls,
# no calls to predict.cta_tree(), no statistical behavior changes.
###############################################################################

# ---- Internal helpers ------------------------------------------------------- #

# Derive status string from tree fields (read-only).
.cta_tree_status <- function(tree) {
  if (isTRUE(tree$no_tree)) return("no_tree")
  s <- cta_strata(tree)
  if (is.na(s)) return("no_tree")
  if (s == 2L) return("stump")
  "valid_tree"
}

# ---- summary.cta_tree ------------------------------------------------------- #

#' Summarize a fitted CTA tree
#'
#' Returns a structured list with class \code{"cta_tree_summary"} capturing
#' tree-level metadata.  All fields are read directly from stored objects;
#' no refitting or prediction is performed.
#'
#' @param object A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @param ... Unused.
#' @return A list of class \code{"cta_tree_summary"} with fields:
#' \describe{
#'   \item{\code{status}}{Character: \code{"valid_tree"}, \code{"stump"},
#'     or \code{"no_tree"}.}
#'   \item{\code{no_tree}}{Logical; \code{TRUE} for leaf-only fits.}
#'   \item{\code{root_attribute}}{Character attribute name at the root split;
#'     \code{NA_character_} for no-tree fits.}
#'   \item{\code{n_nodes}}{Total number of nodes including leaves.}
#'   \item{\code{n_splits}}{Number of non-leaf (split) nodes.}
#'   \item{\code{n_leaves}}{Number of terminal leaf endpoints (= \code{strata}).}
#'   \item{\code{strata}}{Alias for \code{n_leaves}; \code{NA_integer_} for
#'     no-tree fits.}
#'   \item{\code{overall_ess}}{WESS when weights are active, ESS otherwise;
#'     \code{NA_real_} when absent.}
#'   \item{\code{d}}{D statistic (\code{NA_real_} for no-tree or ESS <= 0).}
#'   \item{\code{min_terminal_denom}}{Smallest leaf \code{n_obs};
#'     \code{NA_integer_} for no-tree fits.}
#'   \item{\code{endpoint_denominators}}{Named integer vector of leaf
#'     \code{n_obs}; \code{integer(0)} for no-tree fits.}
#'   \item{\code{has_weights}}{Logical; \code{TRUE} when case weights are
#'     active.}
#'   \item{\code{mindenom}}{MINDENOM used when fitting.}
#'   \item{\code{alpha_split}}{Significance threshold used when fitting.}
#'   \item{\code{prune_alpha}}{Pruning threshold used when fitting.}
#'   \item{\code{loo}}{LOO mode string used when fitting.}
#' }
#' @seealso \code{\link{oda_cta_fit}}, \code{\link{cta_node_table}},
#'   \code{\link{cta_strata}}, \code{\link{cta_d_stat}}
#' @export
summary.cta_tree <- function(object, ...) {
  tree <- object
  stopifnot(inherits(tree, "cta_tree"))

  no_tree <- isTRUE(tree$no_tree)
  status  <- .cta_tree_status(tree)

  # Root attribute — read from stored root node; NA for no-tree.
  root_attr <- NA_character_
  if (!no_tree && !is.null(tree$root_id)) {
    root_nd   <- tree$nodes[[tree$root_id]]
    root_attr <- if (!is.null(root_nd)) root_nd$attribute %||% NA_character_
                 else NA_character_
  }

  # Node counts: n_splits = non-leaf nodes; n_nodes from tree field.
  n_nodes  <- as.integer(tree$n_nodes %||% length(tree$nodes))
  n_splits <- if (!no_tree) {
    as.integer(sum(vapply(tree$nodes,
                          function(nd) !is.null(nd) && !isTRUE(nd$leaf),
                          logical(1L))))
  } else 0L

  strata          <- cta_strata(tree)
  ep_denoms       <- cta_endpoint_denominators(tree)
  min_term_denom  <- cta_min_terminal_denom(tree)
  d               <- cta_d_stat(tree)
  overall_ess     <- if (!is.null(tree$overall_ess)) tree$overall_ess
                     else NA_real_

  s <- list(
    status                = status,
    no_tree               = no_tree,
    root_attribute        = root_attr,
    n_nodes               = n_nodes,
    n_splits              = n_splits,
    n_leaves              = strata,
    strata                = strata,
    overall_ess           = overall_ess,
    d                     = d,
    min_terminal_denom    = min_term_denom,
    endpoint_denominators = ep_denoms,
    has_weights           = isTRUE(tree$has_weights),
    mindenom              = tree$mindenom    %||% NA_integer_,
    alpha_split           = tree$alpha_split %||% NA_real_,
    prune_alpha           = tree$prune_alpha %||% NA_real_,
    loo                   = tree$loo         %||% NA_character_
  )
  class(s) <- "cta_tree_summary"
  s
}

# ---- print.cta_tree_summary ------------------------------------------------- #

#' Print a CTA tree summary
#'
#' @param x A \code{cta_tree_summary} from \code{\link{summary.cta_tree}}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.cta_tree_summary <- function(x, ...) {
  wt_str <- if (isTRUE(x$has_weights)) "  weighted=TRUE" else ""
  cat(sprintf("\nCTA Tree Summary  status=%s%s\n", x$status, wt_str))

  md  <- if (!is.null(x$mindenom)    && !is.na(x$mindenom))    as.character(x$mindenom)    else "?"
  al  <- if (!is.null(x$alpha_split) && !is.na(x$alpha_split)) sprintf("%.3f", x$alpha_split) else "?"
  pr  <- if (!is.null(x$prune_alpha) && !is.na(x$prune_alpha)) sprintf("%.3f", x$prune_alpha) else "?"
  lo  <- if (!is.null(x$loo)         && !is.na(x$loo))         x$loo                          else "?"
  cat(sprintf("  mindenom=%s  alpha_split=%s  prune=%s  loo=%s\n", md, al, pr, lo))

  if (!isTRUE(x$no_tree)) {
    cat(sprintf("  root: %s\n", x$root_attribute %||% "?"))
    nl  <- if (!is.null(x$n_leaves)  && !is.na(x$n_leaves))  x$n_leaves  else NA_integer_
    nsp <- if (!is.null(x$n_splits)  && !is.na(x$n_splits))  x$n_splits  else NA_integer_
    cat(sprintf("  nodes: %d total  (%d split  %d leaf)\n",
                x$n_nodes %||% NA_integer_, nsp, nl))

    ess_str <- if (!is.null(x$overall_ess) && !is.na(x$overall_ess))
      sprintf("%.2f%%", x$overall_ess) else "NA"
    d_str   <- if (!is.null(x$d) && !is.na(x$d))
      sprintf("%.4f", x$d) else "NA"
    mtd_str <- if (!is.null(x$min_terminal_denom) && !is.na(x$min_terminal_denom))
      as.character(x$min_terminal_denom) else "NA"
    cat(sprintf("  overall_ess=%s  D=%s  min_denom=%s\n", ess_str, d_str, mtd_str))
  } else {
    cat("  No tree found (leaf-only)\n")
  }

  cat("\n")
  invisible(x)
}

# =============================================================================
# cta_endpoint_table() — terminal-leaf report table
# =============================================================================

# Internal: build a human-readable branch label for one branch of a split node.
#
# split_nd   — a split node list (has $attribute, $rule, $split_labels)
# branch_idx — 1-based index into split_nd$child_ids / split_nd$split_labels
#
# Returns a character string like "V14<=0.5" or "V14>0.5".
.cta_branch_string <- function(split_nd, branch_idx) {
  attr_nm <- split_nd$attribute %||% "?"
  rule    <- split_nd$rule
  cls     <- split_nd$split_labels[branch_idx]

  if (is.null(rule)) return(sprintf("%s[->%d]", attr_nm, as.integer(cls)))

  tryCatch({
    t <- rule$type
    if (identical(t, "ordered_cut")) {
      cv        <- rule$cut_value
      # direction "0->1": <=cut assigned class 0, >cut assigned class 1.
      # direction "1->0": <=cut assigned class 1, >cut assigned class 0.
      left_cls  <- if (identical(rule$direction, "0->1")) 0L else 1L
      if (as.integer(cls) == left_cls)
        sprintf("%s<=%g", attr_nm, cv)
      else
        sprintf("%s>%g", attr_nm, cv)
    } else if (identical(t, "multiclass_ordered")) {
      cuts <- rule$cut_values
      segs <- as.integer(rule$seg_classes)
      k    <- which(segs == as.integer(cls))[1L]
      K    <- length(segs)
      if (is.na(k))
        sprintf("%s[->%d]", attr_nm, as.integer(cls))
      else if (k == 1L)
        sprintf("%s<=%g", attr_nm, cuts[1L])
      else if (k == K)
        sprintf("%s>%g", attr_nm, cuts[K - 1L])
      else
        sprintf("%g<%s<=%g", cuts[k - 1L], attr_nm, cuts[k])
    } else if (t %in% c("binary_map", "nominal_cut")) {
      lvls <- if (as.integer(cls) == 0L) rule$left_levels else rule$right_levels
      sprintf("%s in {%s}", attr_nm, paste(lvls, collapse = ","))
    } else {
      sprintf("%s[->%d]", attr_nm, as.integer(cls))
    }
  }, error = function(e) sprintf("%s[->%d]", attr_nm, as.integer(cls)))
}

# Internal: reconstruct the root-to-leaf path as a character vector of branch
# labels, one per ancestor split.  Returns character(0) for the root leaf
# (no-tree case).
.cta_path_segments <- function(tree, leaf_nid) {
  segs   <- character(0)
  cur_id <- leaf_nid

  repeat {
    nd        <- tree$nodes[[cur_id]]
    if (is.null(nd)) break
    parent_id <- nd$parent_id %||% 0L
    if (parent_id == 0L) break          # reached root; no further ancestor

    parent_nd <- tree$nodes[[parent_id]]
    if (is.null(parent_nd)) break

    branch_k <- which(parent_nd$child_ids == cur_id)[1L]
    seg <- if (is.na(branch_k))
      sprintf("[node%d]", cur_id)
    else
      .cta_branch_string(parent_nd, branch_k)

    segs   <- c(seg, segs)              # prepend → builds root-to-leaf order
    cur_id <- parent_id
  }

  segs
}

#' Canonical terminal endpoint map for a fitted CTA tree
#'
#' Returns one row per terminal leaf (endpoint) of a \code{cta_tree}.  All
#' values are read directly from stored node fields; no refitting or prediction
#' is performed.  This is the canonical endpoint map for reporting, translation,
#' ORT, and staged workflows.
#'
#' Leaf class counts are stored on every terminal node at fit time
#' (\code{class_counts_raw}, \code{class_counts_weighted}).  \code{target_n}
#' and \code{target_prop} are derived from the stored counts.
#'
#' ESS, WESS, p, LOO status, LOO ESS/WESSL, and LOOp are canonical split-node
#' report metrics (see \code{\link{cta_node_table}}).  Terminal endpoints are
#' connected to those metrics through their parent split-node lineage.  The
#' \code{parent_split_*} columns expose the immediate parent split's canonical
#' metrics for auditability.  They are not recomputed ESS at the leaf.
#'
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @param target_class Integer class label to use as the target (positive)
#'   class for \code{target_n} and \code{target_prop}.  When \code{NULL}
#'   (default), the function auto-detects: for binary trees with classes
#'   \code{0} and \code{1}, class \code{1} is used; otherwise \code{target_n}
#'   and \code{target_prop} are \code{NA}.
#' @return A \code{data.frame} with one row per terminal leaf and columns:
#' \describe{
#'   \item{\code{endpoint_id}}{Integer sequential endpoint index 1..n.}
#'   \item{\code{leaf_node_id}}{Integer tree node identifier for this leaf.}
#'   \item{\code{terminal_marker}}{Character \code{"*"} on every row.}
#'   \item{\code{terminal}}{Logical \code{TRUE} on every row.}
#'   \item{\code{depth}}{Integer depth from root (root = 1).}
#'   \item{\code{parent_split_node_id}}{Integer parent split node identifier.}
#'   \item{\code{path}}{Character; AND-joined branch labels from root to this
#'     leaf (e.g. \code{"V14<=0.5 AND V15>0.5"}).  \code{"/"} for the root
#'     leaf.}
#'   \item{\code{n}}{Integer raw observation count at this endpoint.}
#'   \item{\code{class_counts_raw}}{List column; each element is a named
#'     integer vector of raw per-class counts, or \code{NULL}.}
#'   \item{\code{class_counts_weighted}}{List column; each element is a named
#'     numeric vector of weighted per-class counts, or \code{NULL}.}
#'   \item{\code{predicted_class}}{Integer class label assigned to this
#'     endpoint (stored leaf majority class).}
#'   \item{\code{target_n}}{Integer count of \code{target_class} observations
#'     at this endpoint (\code{NA} when not resolvable).}
#'   \item{\code{target_prop}}{Numeric proportion \code{target_n / n}
#'     (\code{NA} when not resolvable).}
#'   \item{\code{parent_split_attribute}}{Character attribute name of the
#'     immediate parent split node (\code{NA} if not available).}
#'   \item{\code{parent_split_ess}}{Numeric ESS of the parent split node
#'     (\code{NA} if not available).}
#'   \item{\code{parent_split_wess}}{Numeric WESS (weighted ESS) of the parent
#'     split node (\code{NA} if not available).}
#'   \item{\code{parent_split_loo_status}}{Character LOO status of the parent
#'     split node (\code{NA} if not available).}
#'   \item{\code{parent_split_loo_ess}}{Numeric LOO ESS/WESSL of the parent
#'     split node (\code{NA} if not available).}
#'   \item{\code{parent_split_p_mc}}{Numeric MC p-value of the parent split
#'     node (\code{NA} if not available).}
#' }
#' For a no-tree fit the returned data frame has zero rows but the correct
#' column structure and types.
#' @seealso \code{\link{oda_cta_fit}}, \code{\link{cta_node_table}},
#'   \code{\link{summary.cta_tree}}, \code{\link{cta_strata}},
#'   \code{\link{cta_endpoint_denominators}}, \code{\link{cta_endpoint_summary}},
#'   \code{\link{cta_endpoint_counts}}
#' @examples
#' data(mtcars)
#' X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
#' y    <- as.integer(mtcars$am)
#' tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
#' cta_endpoint_table(tree)
cta_endpoint_table <- function(tree, target_class = NULL) {
  stopifnot(inherits(tree, "cta_tree"))

  # ------------------------------------------------------------------
  # Empty data frame template — returned for no-tree or zero-leaf cases
  # ------------------------------------------------------------------
  empty_df <- function() {
    df <- data.frame(
      endpoint_id              = integer(0),
      leaf_node_id             = integer(0),
      terminal_marker          = character(0),
      terminal                 = logical(0),
      depth                    = integer(0),
      parent_split_node_id     = integer(0),
      path                     = character(0),
      n                        = integer(0),
      predicted_class          = integer(0),
      target_n                 = integer(0),
      target_prop              = numeric(0),
      parent_split_attribute   = character(0),
      parent_split_ess         = numeric(0),
      parent_split_wess        = numeric(0),
      parent_split_loo_status  = character(0),
      parent_split_loo_ess     = numeric(0),
      parent_split_p_mc        = numeric(0),
      stringsAsFactors         = FALSE
    )
    df$class_counts_raw      <- list()
    df$class_counts_weighted <- list()
    df
  }

  if (isTRUE(tree$no_tree)) return(empty_df())

  leaves <- Filter(function(nd) !is.null(nd) && isTRUE(nd$leaf), tree$nodes)
  if (length(leaves) == 0L) return(empty_df())

  # Sort leaves by node_id for stable ordering
  node_ids <- vapply(leaves, function(nd) nd$node_id, integer(1L))
  leaves   <- leaves[order(node_ids)]

  # ------------------------------------------------------------------
  # Build per-leaf rows
  # ------------------------------------------------------------------
  nl <- length(leaves)
  leaf_node_id_vec           <- integer(nl)
  parent_split_vec           <- integer(nl)
  depth_vec                  <- integer(nl)
  path_vec                   <- character(nl)
  n_vec                      <- integer(nl)
  predicted_class_vec        <- integer(nl)
  class_counts_raw_list      <- vector("list", nl)
  class_counts_weighted_list <- vector("list", nl)
  ps_attribute_vec           <- character(nl)
  ps_ess_vec                 <- numeric(nl)
  ps_wess_vec                <- numeric(nl)
  ps_loo_status_vec          <- character(nl)
  ps_loo_ess_vec             <- numeric(nl)
  ps_p_mc_vec                <- numeric(nl)

  for (i in seq_len(nl)) {
    nd   <- leaves[[i]]
    segs <- .cta_path_segments(tree, nd$node_id)

    leaf_node_id_vec[i]    <- nd$node_id
    parent_split_vec[i]    <- nd$parent_id       %||% NA_integer_
    depth_vec[i]           <- nd$depth           %||% NA_integer_
    path_vec[i]            <- if (length(segs) == 0L) "/"
                              else paste(segs, collapse = " AND ")
    n_vec[i]               <- nd$n_obs           %||% NA_integer_
    predicted_class_vec[i] <- nd$majority_class  %||% NA_integer_
    class_counts_raw_list[[i]]      <- nd$class_counts_raw      %||% NULL
    class_counts_weighted_list[[i]] <- nd$class_counts_weighted %||% NULL

    # Parent split lineage
    pid <- nd$parent_id %||% NA_integer_
    pnd <- if (!is.na(pid)) tree$nodes[[pid]] else NULL
    ps_attribute_vec[i]  <- if (!is.null(pnd)) pnd$attribute    %||% NA_character_ else NA_character_
    ps_ess_vec[i]        <- if (!is.null(pnd)) pnd$ess          %||% NA_real_      else NA_real_
    ps_wess_vec[i]       <- if (!is.null(pnd)) pnd$ess_weighted %||% NA_real_      else NA_real_
    ps_loo_status_vec[i] <- if (!is.null(pnd)) pnd$loo_status   %||% NA_character_ else NA_character_
    ps_loo_ess_vec[i]    <- if (!is.null(pnd)) pnd$loo_ess      %||% NA_real_      else NA_real_
    ps_p_mc_vec[i]       <- if (!is.null(pnd)) pnd$p_mc         %||% NA_real_      else NA_real_
  }

  # ------------------------------------------------------------------
  # Determine target_class (auto-detect binary 0/1 if not supplied)
  # ------------------------------------------------------------------
  if (is.null(target_class)) {
    first_counts <- Filter(Negate(is.null), class_counts_raw_list)
    if (length(first_counts) > 0L) {
      cls_nms <- names(first_counts[[1L]])
      if (!is.null(cls_nms) && length(cls_nms) == 2L &&
          all(sort(cls_nms) == c("0", "1"))) {
        target_class <- 1L
      }
    }
  }

  # ------------------------------------------------------------------
  # Compute target_n and target_prop
  # ------------------------------------------------------------------
  target_n_vec    <- rep(NA_integer_, nl)
  target_prop_vec <- rep(NA_real_,    nl)
  if (!is.null(target_class)) {
    tc_str <- as.character(target_class)
    for (i in seq_len(nl)) {
      ccr <- class_counts_raw_list[[i]]
      ni  <- n_vec[i]
      if (!is.null(ccr) && tc_str %in% names(ccr) &&
          !is.na(ni) && ni > 0L) {
        tn                <- as.integer(ccr[[tc_str]])
        target_n_vec[i]   <- tn
        target_prop_vec[i] <- tn / ni
      }
    }
  }

  df <- data.frame(
    endpoint_id              = seq_len(nl),
    leaf_node_id             = leaf_node_id_vec,
    terminal_marker          = rep("*", nl),
    terminal                 = rep(TRUE, nl),
    depth                    = depth_vec,
    parent_split_node_id     = parent_split_vec,
    path                     = path_vec,
    n                        = n_vec,
    predicted_class          = predicted_class_vec,
    target_n                 = target_n_vec,
    target_prop              = target_prop_vec,
    parent_split_attribute   = ps_attribute_vec,
    parent_split_ess         = ps_ess_vec,
    parent_split_wess        = ps_wess_vec,
    parent_split_loo_status  = ps_loo_status_vec,
    parent_split_loo_ess     = ps_loo_ess_vec,
    parent_split_p_mc        = ps_p_mc_vec,
    stringsAsFactors         = FALSE
  )
  df$class_counts_raw      <- class_counts_raw_list
  df$class_counts_weighted <- class_counts_weighted_list
  df
}

# =============================================================================
# cta_confusion_table() — final selected tree training confusion, tidy format
# =============================================================================

#' Final selected tree training confusion table
#'
#' Returns the stored full-tree training confusion matrix for the final
#' selected CTA model in tidy long format (one row per actual × predicted
#' class pair).
#'
#' The confusion matrix is captured at fit time at the exact moment the
#' winning candidate is selected, using the same scoring predictions.  For
#' the expanded ENUMERATE phase, predictions use majority-fallback for
#' missing attributes.  For the root-only stump phase, predictions are
#' path-local (observations whose root attribute is missing are excluded).
#'
#' This function does \strong{not} report split-node local confusion.
#' Split-node confusion reflects all observations at a node classified by
#' that node's rule alone; it is not the same as full-tree confusion for
#' trees with more than one split.  The two coincide incidentally for stumps
#' but the semantics here are always final-tree.
#'
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{\code{actual}}{Integer actual class label.}
#'   \item{\code{predicted}}{Integer predicted class label.}
#'   \item{\code{n}}{Integer raw count of observations with this actual ×
#'     predicted combination in the final selected tree.}
#' }
#' Rows are sorted by \code{actual} then \code{predicted}.
#' For a no-tree fit (or if \code{training_confusion} is absent), the
#' returned data frame has zero rows but the correct column structure.
#' @seealso \code{\link{oda_cta_fit}}, \code{\link{summary.cta_tree}},
#'   \code{\link{cta_endpoint_table}}, \code{\link{cta_node_table}}
#' @examples
#' data(mtcars)
#' X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
#' y    <- as.integer(mtcars$am)
#' tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
#' cta_confusion_table(tree)
cta_confusion_table <- function(tree) {
  stopifnot(inherits(tree, "cta_tree"))

  empty_df <- function() {
    data.frame(actual    = integer(0),
               predicted = integer(0),
               n         = integer(0),
               stringsAsFactors = FALSE)
  }

  conf <- tree$training_confusion
  if (is.null(conf) || !is.matrix(conf) || prod(dim(conf)) == 0L)
    return(empty_df())

  rn <- rownames(conf)
  cn <- colnames(conf)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(conf)) - 1L)
  if (is.null(cn)) cn <- as.character(seq_len(ncol(conf)) - 1L)

  nr <- nrow(conf); nc_m <- ncol(conf)
  actual    <- integer(nr * nc_m)
  predicted <- integer(nr * nc_m)
  n         <- integer(nr * nc_m)

  k <- 0L
  for (i in seq_len(nr)) {
    for (j in seq_len(nc_m)) {
      k <- k + 1L
      actual[k]    <- as.integer(rn[i])
      predicted[k] <- as.integer(cn[j])
      n[k]         <- as.integer(conf[i, j])
    }
  }

  df <- data.frame(actual = actual, predicted = predicted, n = n,
                   stringsAsFactors = FALSE)
  df[order(df$actual, df$predicted), , drop = FALSE]
}

# =============================================================================
# cta_endpoint_summary() — conservative endpoint reporting accessor
# =============================================================================

#' Endpoint reporting summary for a fitted CTA tree
#'
#' Returns one row per terminal leaf (endpoint) with stable endpoint
#' identifiers and stored node fields suitable for downstream reporting.
#' All values are read directly from stored node fields; no refitting,
#' no prediction, and no recomputation of tree metrics is performed.
#'
#' \strong{Scope:} This function reports structural endpoint fields only.
#' It does \emph{not} include endpoint class counts, target-class
#' proportions, event rates, odds, or staging order.  Per-endpoint class
#' counts are available via \code{\link{cta_endpoint_counts}}.
#' Staging-table and event-rate summaries are available via
#' \code{\link{cta_staging_table}}.
#'
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @return A \code{data.frame} with one row per terminal leaf and columns:
#' \describe{
#'   \item{\code{endpoint_id}}{Integer sequential index 1..n in node order.}
#'   \item{\code{endpoint_node_id}}{Integer tree node identifier for this
#'     leaf, corresponding to \code{node_id} in \code{\link{cta_endpoint_table}}.}
#'   \item{\code{path}}{Character; AND-joined branch labels from root to this
#'     leaf (e.g. \code{"V14<=0.5 AND V15>0.5"}).}
#'   \item{\code{depth}}{Integer depth from root (root = 1).}
#'   \item{\code{terminal_prediction}}{Integer class label assigned to this
#'     endpoint (stored leaf \code{majority_class}).}
#'   \item{\code{n_obs}}{Integer raw observation count at this endpoint.}
#'   \item{\code{n_weighted}}{Numeric weighted observation count.  Equals
#'     \code{n_obs} when case weights are not active (not \code{NA}).}
#'   \item{\code{denominator}}{Integer endpoint denominator (equal to
#'     \code{n_obs}); included to align with MPE/MDSA terminology.}
#' }
#' For a no-tree fit the returned data frame has zero rows but the correct
#' column structure and types.
#' @seealso \code{\link{oda_cta_fit}}, \code{\link{cta_endpoint_table}},
#'   \code{\link{cta_strata}}, \code{\link{cta_endpoint_denominators}},
#'   \code{\link{cta_endpoint_counts}}, \code{\link{cta_staging_table}}
#' @examples
#' data(mtcars)
#' X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
#' y    <- as.integer(mtcars$am)
#' tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
#' cta_endpoint_summary(tree)
cta_endpoint_summary <- function(tree) {
  stopifnot(inherits(tree, "cta_tree"))

  empty_df <- function() {
    data.frame(
      endpoint_id         = integer(0),
      endpoint_node_id    = integer(0),
      path                = character(0),
      depth               = integer(0),
      terminal_prediction = integer(0),
      n_obs               = integer(0),
      n_weighted          = numeric(0),
      denominator         = integer(0),
      stringsAsFactors    = FALSE
    )
  }

  if (isTRUE(tree$no_tree)) return(empty_df())

  leaves <- Filter(function(nd) !is.null(nd) && isTRUE(nd$leaf), tree$nodes)
  if (length(leaves) == 0L) return(empty_df())

  # Sort leaves by node_id for stable ordering (matches cta_endpoint_table).
  node_ids <- vapply(leaves, function(nd) nd$node_id, integer(1L))
  leaves   <- leaves[order(node_ids)]

  n <- length(leaves)
  endpoint_node_id    <- integer(n)
  path                <- character(n)
  depth               <- integer(n)
  terminal_prediction <- integer(n)
  n_obs               <- integer(n)
  n_weighted          <- numeric(n)

  for (i in seq_len(n)) {
    nd   <- leaves[[i]]
    segs <- .cta_path_segments(tree, nd$node_id)

    endpoint_node_id[i]    <- nd$node_id
    path[i]                <- if (length(segs) == 0L) "/"
                              else paste(segs, collapse = " AND ")
    depth[i]               <- nd$depth          %||% NA_integer_
    terminal_prediction[i] <- nd$majority_class  %||% NA_integer_
    n_obs[i]               <- nd$n_obs           %||% NA_integer_
    n_weighted[i]          <- nd$n_weighted      %||% NA_real_
  }

  data.frame(
    endpoint_id         = seq_len(n),
    endpoint_node_id    = endpoint_node_id,
    path                = path,
    depth               = depth,
    terminal_prediction = terminal_prediction,
    n_obs               = n_obs,
    n_weighted          = n_weighted,
    denominator         = n_obs,
    stringsAsFactors    = FALSE
  )
}

# =============================================================================
# cta_endpoint_counts() — per-endpoint × class count table
# =============================================================================

#' Per-endpoint class count table for a fitted CTA tree
#'
#' Returns one row per terminal endpoint (leaf) per actual class, read
#' directly from stored leaf node fields.  No refitting, no prediction,
#' and no recomputation from training data is performed.
#'
#' Class counts are stored at fit time by \code{\link{oda_cta_fit}} on
#' every terminal leaf.  Row order within each endpoint follows the order
#' of \code{names(leaf$class_counts_raw)}, which is ascending by class
#' label.  Endpoints are ordered by \code{node_id}, matching
#' \code{\link{cta_endpoint_summary}}.
#'
#' \strong{Scope:} This function exposes stored raw and weighted class
#' counts only.  It does \emph{not} include target-class proportions,
#' event rates, odds, or staging order.  Staging-table and event-rate
#' summaries are available via \code{\link{cta_staging_table}}.
#'
#' If any terminal leaf is missing the stored class counts (i.e., the
#' \code{cta_tree} was fitted by an earlier version of odacore that did
#' not store endpoint counts), the function stops with a clear error.
#'
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @return A \code{data.frame} with one row per terminal endpoint per
#'   actual class and columns:
#' \describe{
#'   \item{\code{endpoint_id}}{Integer sequential endpoint index 1..n in
#'     node order, matching \code{\link{cta_endpoint_summary}}.}
#'   \item{\code{endpoint_node_id}}{Integer tree node identifier for this
#'     endpoint leaf.}
#'   \item{\code{path}}{Character; AND-joined branch labels from root to
#'     this leaf (e.g. \code{"V14<=0.5 AND V15>0.5"}).}
#'   \item{\code{terminal_prediction}}{Integer class label assigned to
#'     this endpoint (stored leaf \code{majority_class}).}
#'   \item{\code{class}}{Character; actual class label for this row
#'     (e.g. \code{"0"}, \code{"1"}).}
#'   \item{\code{n_raw}}{Integer raw count of observations of this actual
#'     class reaching this endpoint.}
#'   \item{\code{n_weighted}}{Numeric weighted total for this actual class
#'     reaching this endpoint.  Equals \code{n_raw} when case weights are
#'     not active.}
#' }
#' For a no-tree fit the returned data frame has zero rows but the correct
#' column structure and types.
#' @seealso \code{\link{oda_cta_fit}}, \code{\link{cta_endpoint_summary}},
#'   \code{\link{cta_confusion_table}}, \code{\link{cta_endpoint_table}},
#'   \code{\link{cta_staging_table}}
#' @examples
#' data(mtcars)
#' X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
#' y    <- as.integer(mtcars$am)
#' tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
#' cta_endpoint_counts(tree)
cta_endpoint_counts <- function(tree) {
  stopifnot(inherits(tree, "cta_tree"))

  empty_df <- function() {
    data.frame(
      endpoint_id         = integer(0),
      endpoint_node_id    = integer(0),
      path                = character(0),
      terminal_prediction = integer(0),
      class               = character(0),
      n_raw               = integer(0),
      n_weighted          = numeric(0),
      stringsAsFactors    = FALSE
    )
  }

  if (isTRUE(tree$no_tree)) return(empty_df())

  leaves <- Filter(function(nd) !is.null(nd) && isTRUE(nd$leaf), tree$nodes)
  if (length(leaves) == 0L) return(empty_df())

  # Sort by node_id (matches cta_endpoint_summary ordering).
  node_ids <- vapply(leaves, function(nd) nd$node_id, integer(1L))
  leaves   <- leaves[order(node_ids)]

  # Guard: require both stored class-count fields on all leaves.
  missing_counts <- vapply(leaves,
                           function(nd) is.null(nd$class_counts_raw) ||
                                        is.null(nd$class_counts_weighted),
                           logical(1L))
  if (any(missing_counts)) {
    stop("endpoint class counts are unavailable for this cta_tree; ",
         "refit with a version of odacore that stores endpoint counts")
  }

  # Build output columns — one entry per leaf × class.
  out_endpoint_id         <- integer(0)
  out_endpoint_node_id    <- integer(0)
  out_path                <- character(0)
  out_terminal_prediction <- integer(0)
  out_class               <- character(0)
  out_n_raw               <- integer(0)
  out_n_weighted          <- numeric(0)

  for (i in seq_along(leaves)) {
    nd     <- leaves[[i]]
    segs   <- .cta_path_segments(tree, nd$node_id)
    path_s <- if (length(segs) == 0L) "/" else paste(segs, collapse = " AND ")
    pred   <- nd$majority_class %||% NA_integer_

    cls_names <- names(nd$class_counts_raw)
    m         <- length(cls_names)

    out_endpoint_id         <- c(out_endpoint_id,         rep(i,          m))
    out_endpoint_node_id    <- c(out_endpoint_node_id,    rep(nd$node_id, m))
    out_path                <- c(out_path,                rep(path_s,     m))
    out_terminal_prediction <- c(out_terminal_prediction, rep(pred,       m))
    out_class               <- c(out_class,               cls_names)
    out_n_raw               <- c(out_n_raw,               unname(nd$class_counts_raw))
    out_n_weighted          <- c(out_n_weighted,          unname(nd$class_counts_weighted))
  }

  data.frame(
    endpoint_id         = out_endpoint_id,
    endpoint_node_id    = out_endpoint_node_id,
    path                = out_path,
    terminal_prediction = out_terminal_prediction,
    class               = out_class,
    n_raw               = out_n_raw,
    n_weighted          = out_n_weighted,
    stringsAsFactors    = FALSE
  )
}

# =============================================================================
# cta_staging_table() — per-endpoint staging table ordered by target propensity
# =============================================================================

#' Staging table for a fitted CTA tree
#'
#' Returns one row per terminal endpoint ordered by ascending target-class
#' propensity (lowest to highest risk stratum).  Empirical counts,
#' proportions, and odds are computed from the stored leaf class counts.
#' When an endpoint is perfectly predicted (100 percent one class), the
#' empirical odds and proportion are undefined; the \code{adjust_perfect}
#' option adds one hypothetical misclassified observation to the undefined
#' profile so all endpoints can be ranked and compared — a canon remedy
#' anchored in Yarnold and Linden (2017).
#'
#' \strong{Scope:} The two-class case is handled automatically when
#' \code{target_class = NULL} (defaults to the numerically larger class
#' label, typically 1).  For trees with three or more classes
#' \code{target_class} must be supplied explicitly.
#'
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @param target_class Integer (or coercible); the class label treated as the
#'   target (positive / high-risk) class.  \code{NULL} (default) uses the
#'   numerically largest class label for binary trees, and stops for trees
#'   with three or more classes.
#' @param weighted Logical.  \code{FALSE} (default) uses raw observation
#'   counts; \code{TRUE} uses case-weighted counts.
#' @param adjust_perfect Logical.  \code{TRUE} (default) applies the
#'   one-hypothetical-misclassification adjustment to perfectly predicted
#'   endpoints so that all endpoints can be ordered by propensity.
#' @return A \code{data.frame} with one row per terminal endpoint, ordered
#'   by ascending target-class propensity (lowest to highest risk stratum),
#'   with columns:
#' \describe{
#'   \item{\code{stage}}{Integer rank 1..n, ascending by target proportion.}
#'   \item{\code{endpoint_id}}{Integer sequential endpoint index, matching
#'     \code{\link{cta_endpoint_summary}}.}
#'   \item{\code{endpoint_node_id}}{Integer tree node identifier.}
#'   \item{\code{path}}{Character; AND-joined branch labels from root.}
#'   \item{\code{terminal_prediction}}{Integer majority-class prediction.}
#'   \item{\code{target_class}}{Integer; the target class used for this
#'     table.}
#'   \item{\code{target_n}}{Numeric; raw (or weighted) count of target-class
#'     observations at this endpoint.}
#'   \item{\code{denominator}}{Numeric; total raw (or weighted) observations
#'     at this endpoint.}
#'   \item{\code{target_proportion}}{Numeric; empirical target-class
#'     proportion (\code{target_n / denominator}).}
#'   \item{\code{non_target_n}}{Numeric; denominator minus target_n.}
#'   \item{\code{odds}}{Numeric; empirical odds
#'     (\code{target_n / non_target_n}); \code{NA} when
#'     \code{perfectly_predicted} is \code{TRUE}.}
#'   \item{\code{perfectly_predicted}}{Logical; \code{TRUE} when the
#'     endpoint is 100 percent one class (\code{target_n == 0} or
#'     \code{non_target_n == 0}).}
#'   \item{\code{adjusted}}{Logical; \code{TRUE} when the
#'     one-hypothetical-misclassification adjustment has been applied.
#'     Always \code{FALSE} when \code{adjust_perfect = FALSE}.}
#'   \item{\code{adjusted_target_n}}{Numeric; target_n after adjustment.
#'     Equal to \code{target_n} when \code{adjusted} is \code{FALSE}.}
#'   \item{\code{adjusted_denominator}}{Numeric; denominator after
#'     adjustment.}
#'   \item{\code{adjusted_target_proportion}}{Numeric; adjusted proportion.}
#'   \item{\code{adjusted_non_target_n}}{Numeric; adjusted non-target
#'     count.}
#'   \item{\code{adjusted_odds}}{Numeric; adjusted odds.}
#'   \item{\code{weighted}}{Logical; the value of the \code{weighted}
#'     argument.}
#'   \item{\code{n_obs}}{Integer; raw observation count at this endpoint
#'     (from \code{\link{cta_endpoint_summary}}).}
#'   \item{\code{n_weighted}}{Numeric; weighted observation count.}
#' }
#' For a no-tree fit the returned data frame has zero rows but the correct
#' column structure and types.
#' @references
#' Yarnold PR, Linden A (2017). Computing propensity score weights for CTA
#' models involving perfectly predicted endpoints.
#' \emph{Optimal Data Analysis}, \strong{6}, 43-46.
#' @seealso \code{\link{oda_cta_fit}}, \code{\link{cta_endpoint_summary}},
#'   \code{\link{cta_endpoint_counts}}
#' @examples
#' data(mtcars)
#' X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
#' y    <- as.integer(mtcars$am)
#' tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
#' cta_staging_table(tree)
cta_staging_table <- function(tree, target_class = NULL, weighted = FALSE,
                               adjust_perfect = TRUE) {
  stopifnot(inherits(tree, "cta_tree"))

  empty_df <- function() {
    data.frame(
      stage                      = integer(0),
      endpoint_id                = integer(0),
      endpoint_node_id           = integer(0),
      path                       = character(0),
      terminal_prediction        = integer(0),
      target_class               = integer(0),
      target_n                   = numeric(0),
      denominator                = numeric(0),
      target_proportion          = numeric(0),
      non_target_n               = numeric(0),
      odds                       = numeric(0),
      perfectly_predicted        = logical(0),
      adjusted                   = logical(0),
      adjusted_target_n          = numeric(0),
      adjusted_denominator       = numeric(0),
      adjusted_target_proportion = numeric(0),
      adjusted_non_target_n      = numeric(0),
      adjusted_odds              = numeric(0),
      weighted                   = logical(0),
      n_obs                      = integer(0),
      n_weighted                 = numeric(0),
      stringsAsFactors           = FALSE
    )
  }

  if (isTRUE(tree$no_tree)) return(empty_df())

  # Consume endpoint tables (no refitting, no prediction).
  es <- cta_endpoint_summary(tree)
  ec <- cta_endpoint_counts(tree)

  if (nrow(es) == 0L) return(empty_df())

  # Resolve target_class.
  all_classes <- sort(unique(ec$class))          # character, ascending
  if (is.null(target_class)) {
    if (length(all_classes) != 2L)
      stop("target_class must be specified for trees with ",
           length(all_classes), " classes (found: ",
           paste(all_classes, collapse = ", "), ")")
    target_class_int <- as.integer(max(as.integer(all_classes)))
  } else {
    target_class_int <- as.integer(target_class)
  }
  target_class_chr <- as.character(target_class_int)
  if (!target_class_chr %in% all_classes)
    stop("target_class ", target_class_int,
         " not found in tree classes: ", paste(all_classes, collapse = ", "))

  # Select count column.
  count_col <- if (isTRUE(weighted)) "n_weighted" else "n_raw"

  # Per-endpoint denominator (sum over all classes).
  denom_by_ep <- tapply(ec[[count_col]], ec$endpoint_id, sum)

  # Per-endpoint target count.
  ec_target    <- ec[ec$class == target_class_chr, , drop = FALSE]
  target_by_ep <- setNames(ec_target[[count_col]],
                            as.character(ec_target$endpoint_id))

  # Align to es row order.
  ep_chr      <- as.character(es$endpoint_id)
  denominator <- as.numeric(denom_by_ep[ep_chr])
  target_n    <- as.numeric(target_by_ep[ep_chr])
  target_n[is.na(target_n)] <- 0.0      # guard: zero-count class absent from ec

  non_target_n      <- denominator - target_n
  target_proportion <- ifelse(denominator > 0, target_n / denominator, NA_real_)

  perfectly_predicted <- (denominator > 0) &
                         (target_n == 0 | non_target_n == 0)

  odds <- ifelse(perfectly_predicted | denominator == 0,
                 NA_real_,
                 target_n / non_target_n)

  # Canon adjustment: one hypothetical misclassified obs in undefined profile.
  adjusted         <- isTRUE(adjust_perfect) & perfectly_predicted
  adj_target_n     <- target_n
  adj_non_target_n <- non_target_n
  adj_denominator  <- denominator

  boost_target    <- adjusted & (target_n == 0)
  boost_nontarget <- adjusted & (non_target_n == 0)

  adj_target_n[boost_target]        <- adj_target_n[boost_target] + 1.0
  adj_denominator[boost_target]     <- adj_denominator[boost_target] + 1.0
  adj_non_target_n[boost_nontarget] <- adj_non_target_n[boost_nontarget] + 1.0
  adj_denominator[boost_nontarget]  <- adj_denominator[boost_nontarget] + 1.0

  adj_target_proportion <- ifelse(adj_denominator > 0,
                                  adj_target_n / adj_denominator,
                                  NA_real_)
  adj_odds <- ifelse(adj_non_target_n == 0 | adj_denominator == 0,
                     NA_real_,
                     adj_target_n / adj_non_target_n)

  # Sort: adjusted proportion when adjust_perfect, else empirical; tie by node_id.
  sort_prop <- if (isTRUE(adjust_perfect)) adj_target_proportion else target_proportion
  ord       <- order(sort_prop, es$endpoint_node_id, na.last = TRUE)
  n         <- nrow(es)

  data.frame(
    stage                      = seq_len(n),
    endpoint_id                = es$endpoint_id[ord],
    endpoint_node_id           = es$endpoint_node_id[ord],
    path                       = es$path[ord],
    terminal_prediction        = es$terminal_prediction[ord],
    target_class               = rep(target_class_int, n),
    target_n                   = target_n[ord],
    denominator                = denominator[ord],
    target_proportion          = target_proportion[ord],
    non_target_n               = non_target_n[ord],
    odds                       = odds[ord],
    perfectly_predicted        = perfectly_predicted[ord],
    adjusted                   = adjusted[ord],
    adjusted_target_n          = adj_target_n[ord],
    adjusted_denominator       = adj_denominator[ord],
    adjusted_target_proportion = adj_target_proportion[ord],
    adjusted_non_target_n      = adj_non_target_n[ord],
    adjusted_odds              = adj_odds[ord],
    weighted                   = rep(isTRUE(weighted), n),
    n_obs                      = es$n_obs[ord],
    n_weighted                 = es$n_weighted[ord],
    stringsAsFactors           = FALSE
  )
}

# =============================================================================
# cta_propensity_weights() — endpoint × class stabilized propensity weights
# =============================================================================

#' Endpoint-level propensity-score weights for a fitted CTA tree
#'
#' Returns one row per terminal endpoint per actual class, containing the
#' CTA-derived stabilized propensity-style weights described in Yarnold and
#' Linden (2017).  All values are computed on demand from the stored leaf
#' class counts; no refitting, no prediction, and no training-data
#' recomputation is performed.
#'
#' \strong{Formula:} For endpoint \eqn{s} and actual class \eqn{z},
#' \deqn{w_{s,z} = \frac{n_s \cdot \Pr(Z=z)}{n_{s,z}}}
#' where \eqn{n_s} is the endpoint denominator, \eqn{n_{s,z}} is the raw
#' count of class \eqn{z} observations at endpoint \eqn{s}, and
#' \eqn{\Pr(Z=z)} is the marginal class probability across the full
#' classified analytic sample.  This weighting makes each endpoint's
#' class distribution proportional to the global marginal, enabling
#' translational comparisons across strata.
#'
#' \strong{Perfect endpoints:} When \eqn{n_{s,z} = 0} for some class, the
#' empirical weight is undefined (\code{Inf}).  When \code{adjusted = TRUE}
#' (default), one hypothetical misclassified observation is added to the
#' absent class profile — and to the global marginal totals — so that all
#' endpoint × class cells yield finite adjusted weights.  This is the canon
#' remedy from Yarnold and Linden (2017).
#'
#' \strong{Scope:} Raw observation counts (\code{n_raw}) are used
#' exclusively.  The function does not return observation-level weights;
#' those would require endpoint membership per training observation, which
#' is not stored on the fitted tree.
#'
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @param target_class Integer (or coercible); annotation column only —
#'   does not filter output rows.  \code{NULL} (default) uses the
#'   numerically largest class label for binary trees, and stops for
#'   trees with three or more classes.
#' @param adjusted Logical.  \code{TRUE} (default) applies the
#'   one-hypothetical-misclassification adjustment so that all cells
#'   yield finite adjusted weights.  \code{FALSE} leaves undefined
#'   weights as \code{Inf} and adjusted columns equal to empirical.
#' @return A \code{data.frame} with one row per terminal endpoint per
#'   actual class, with columns:
#' \describe{
#'   \item{\code{endpoint_id}}{Integer sequential endpoint index.}
#'   \item{\code{endpoint_node_id}}{Integer tree node identifier.}
#'   \item{\code{path}}{Character; AND-joined branch labels from root.}
#'   \item{\code{terminal_prediction}}{Integer majority-class prediction.}
#'   \item{\code{class}}{Character; actual class label for this row.}
#'   \item{\code{target_class}}{Integer; design-annotation class label.}
#'   \item{\code{class_n}}{Integer; raw count of this class at this
#'     endpoint (empirical \eqn{n_{s,z}}).}
#'   \item{\code{endpoint_n}}{Integer; total raw observations at this
#'     endpoint (empirical \eqn{n_s}).}
#'   \item{\code{marginal_class_n}}{Integer; total raw observations of
#'     this class across all endpoints (empirical \eqn{N_z}).}
#'   \item{\code{marginal_total_n}}{Integer; total classified observations
#'     across all endpoints (empirical \eqn{N}).}
#'   \item{\code{marginal_class_probability}}{Numeric; empirical marginal
#'     class probability \eqn{\Pr(Z=z) = N_z / N}.}
#'   \item{\code{propensity_weight}}{Numeric; empirical stabilized weight
#'     \eqn{n_s \cdot \Pr(Z=z) / n_{s,z}}.  \code{Inf} when
#'     \code{class_n == 0}.}
#'   \item{\code{undefined_empirical}}{Logical; \code{TRUE} when
#'     \code{class_n == 0} (empirical weight is undefined).}
#'   \item{\code{perfectly_predicted_endpoint}}{Logical; \code{TRUE} when
#'     any class has \code{class_n == 0} at this endpoint.}
#'   \item{\code{adjusted}}{Logical; \code{TRUE} when the
#'     one-hypothetical-observation adjustment was applied to this row.}
#'   \item{\code{adjusted_class_n}}{Numeric; \code{class_n + 1} where
#'     \code{adjusted}, otherwise \code{class_n}.}
#'   \item{\code{adjusted_endpoint_n}}{Numeric; endpoint denominator
#'     after adjustment.}
#'   \item{\code{adjusted_marginal_class_n}}{Numeric; global class count
#'     after all hypothetical additions.}
#'   \item{\code{adjusted_marginal_total_n}}{Numeric; global total after
#'     all hypothetical additions.}
#'   \item{\code{adjusted_marginal_class_probability}}{Numeric; adjusted
#'     marginal class probability.}
#'   \item{\code{adjusted_propensity_weight}}{Numeric; adjusted weight
#'     \eqn{n_s^* \cdot \Pr^*(Z=z) / n_{s,z}^*}.  Finite whenever
#'     \code{adjusted_class_n > 0}.}
#' }
#' For a no-tree fit the returned data frame has zero rows but the correct
#' column structure and types.
#' @references
#' Yarnold PR, Linden A (2017). Computing propensity score weights for CTA
#' models involving perfectly predicted endpoints.
#' \emph{Optimal Data Analysis}, \strong{6}, 43-46.
#' @seealso \code{\link{oda_cta_fit}}, \code{\link{cta_endpoint_counts}},
#'   \code{\link{cta_staging_table}}
#' @examples
#' data(mtcars)
#' X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
#' y    <- as.integer(mtcars$am)
#' tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
#' cta_propensity_weights(tree)
cta_propensity_weights <- function(tree, target_class = NULL, adjusted = TRUE) {
  stopifnot(inherits(tree, "cta_tree"))

  empty_df <- function() {
    data.frame(
      endpoint_id                       = integer(0),
      endpoint_node_id                  = integer(0),
      path                              = character(0),
      terminal_prediction               = integer(0),
      class                             = character(0),
      target_class                      = integer(0),
      class_n                           = integer(0),
      endpoint_n                        = integer(0),
      marginal_class_n                  = integer(0),
      marginal_total_n                  = integer(0),
      marginal_class_probability        = numeric(0),
      propensity_weight                 = numeric(0),
      undefined_empirical               = logical(0),
      perfectly_predicted_endpoint      = logical(0),
      adjusted                          = logical(0),
      adjusted_class_n                  = numeric(0),
      adjusted_endpoint_n               = numeric(0),
      adjusted_marginal_class_n         = numeric(0),
      adjusted_marginal_total_n         = numeric(0),
      adjusted_marginal_class_probability = numeric(0),
      adjusted_propensity_weight        = numeric(0),
      stringsAsFactors                  = FALSE
    )
  }

  if (isTRUE(tree$no_tree)) return(empty_df())

  # Consume stored counts — no refitting, no prediction, no mutation of tree.
  ec <- cta_endpoint_counts(tree)
  if (nrow(ec) == 0L) return(empty_df())

  # ---- Resolve target_class (annotation only; does not filter rows) ---------
  all_classes <- sort(unique(ec$class))          # character, ascending
  if (is.null(target_class)) {
    if (length(all_classes) != 2L)
      stop("target_class must be specified for trees with ",
           length(all_classes), " classes (found: ",
           paste(all_classes, collapse = ", "), ")")
    target_class_int <- as.integer(max(as.integer(all_classes)))
  } else {
    target_class_int <- as.integer(target_class)
  }
  target_class_chr <- as.character(target_class_int)
  if (!target_class_chr %in% all_classes)
    stop("target_class ", target_class_int,
         " not found in tree classes: ", paste(all_classes, collapse = ", "))

  # ---- Empirical counts and marginals ---------------------------------------

  # Per-row class count (raw integer).
  class_n <- ec$n_raw                            # integer vector, length = nrow(ec)

  # Per-endpoint denominator aligned to ec rows.
  ep_n_by_id <- tapply(ec$n_raw, ec$endpoint_id, sum)
  endpoint_n <- as.integer(ep_n_by_id[as.character(ec$endpoint_id)])

  # Marginal class counts aligned to ec rows.
  marg_by_cls <- tapply(ec$n_raw, ec$class, sum)
  marginal_class_n <- as.integer(marg_by_cls[ec$class])

  marginal_total_n     <- as.integer(sum(ec$n_raw))
  marginal_class_prob  <- marginal_class_n / marginal_total_n   # numeric per row

  # ---- Empirical propensity weight -----------------------------------------
  undefined_empirical <- class_n == 0L
  propensity_weight   <- ifelse(
    undefined_empirical,
    Inf,
    as.numeric(endpoint_n) * marginal_class_prob / as.numeric(class_n)
  )

  # ---- Perfectly predicted endpoint flag -----------------------------------
  # An endpoint is perfectly predicted if any class has class_n == 0 there.
  pp_by_ep <- tapply(ec$n_raw == 0L, ec$endpoint_id, any)
  perfectly_predicted_endpoint <- as.logical(pp_by_ep[as.character(ec$endpoint_id)])

  # ---- Adjustment: one hypothetical obs added to absent class profiles -----
  needs_adjust <- isTRUE(adjusted) & undefined_empirical

  adj_class_n <- as.numeric(class_n)
  adj_class_n[needs_adjust] <- adj_class_n[needs_adjust] + 1.0

  # adj_endpoint_n is uniform within each endpoint: original n_s plus the
  # number of hypothetical additions made to that endpoint (one per absent class).
  adj_per_ep     <- tapply(as.integer(needs_adjust), ec$endpoint_id, sum)
  adj_endpoint_n <- as.numeric(endpoint_n) +
                    as.numeric(adj_per_ep[as.character(ec$endpoint_id)])

  # Adjusted global marginals: add one hypothetical per absent-class cell.
  # Keep names on marg_by_cls so adj_marg_cls_n can be indexed by class label.
  adj_additions  <- tapply(as.integer(needs_adjust), ec$class, sum)
  adj_marg_cls_n <- marg_by_cls + adj_additions[names(marg_by_cls)]
  adj_marg_total <- sum(adj_marg_cls_n)

  # Per-row adjusted marginals.
  adj_marg_cls_n_row  <- adj_marg_cls_n[ec$class]
  adj_marg_cls_prob   <- adj_marg_cls_n_row / adj_marg_total

  adj_propensity_weight <- adj_endpoint_n * adj_marg_cls_prob / adj_class_n
  row_adjusted          <- needs_adjust

  # ---- Assemble output ------------------------------------------------------
  data.frame(
    endpoint_id                         = ec$endpoint_id,
    endpoint_node_id                    = ec$endpoint_node_id,
    path                                = ec$path,
    terminal_prediction                 = ec$terminal_prediction,
    class                               = ec$class,
    target_class                        = rep(target_class_int, nrow(ec)),
    class_n                             = class_n,
    endpoint_n                          = endpoint_n,
    marginal_class_n                    = marginal_class_n,
    marginal_total_n                    = rep(marginal_total_n, nrow(ec)),
    marginal_class_probability          = marginal_class_prob,
    propensity_weight                   = propensity_weight,
    undefined_empirical                 = undefined_empirical,
    perfectly_predicted_endpoint        = perfectly_predicted_endpoint,
    adjusted                            = row_adjusted,
    adjusted_class_n                    = adj_class_n,
    adjusted_endpoint_n                 = adj_endpoint_n,
    adjusted_marginal_class_n           = adj_marg_cls_n_row,
    adjusted_marginal_total_n           = rep(adj_marg_total, nrow(ec)),
    adjusted_marginal_class_probability = adj_marg_cls_prob,
    adjusted_propensity_weight          = adj_propensity_weight,
    stringsAsFactors                    = FALSE
  )
}

# ---- cta_assign_endpoints --------------------------------------------------- #

#' Assign training (or new) observations to CTA terminal endpoints
#'
#' Traverses the fitted \code{cta_tree} for each row of \code{newdata} and
#' returns the terminal leaf reached, expressed as both its stored
#' \code{node_id} (\code{endpoint_node_id}) and its sequential endpoint index
#' (\code{endpoint_id}) matching \code{\link{cta_endpoint_summary}}.
#'
#' No endpoint membership is stored at fit time.  This function performs
#' the traversal on demand so the \code{cta_tree} object remains lean.  Use
#' the returned \code{endpoint_id} to join with
#' \code{\link{cta_propensity_weights}} and assign endpoint-level weights to
#' individual observations.
#'
#' \strong{Column order requirement:} \code{newdata} must have the same
#' attribute column order as the \code{X} matrix passed to
#' \code{\link{oda_cta_fit}}.  Traversal uses stored integer column positions
#' (\code{attr_col}), not column names.  If both \code{names(newdata)} and the
#' stored \code{tree$attr_names} are non-NULL, a warning is issued when they
#' disagree at the split attribute positions.
#'
#' \strong{Missingness:}
#' \itemize{
#'   \item \code{"na"} (default) — canonical path-local behaviour: when a
#'     split attribute is \code{NA} or a stored miss-code on the observation's
#'     actual traversal path, the row returns \code{NA} for both output
#'     columns.
#'   \item \code{"majority"} — routes the missing observation to the child
#'     subtree with the larger \code{n_obs}, then continues traversal until a
#'     terminal leaf is reached.  Ties are broken by picking the first child.
#' }
#'
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @param newdata A \code{data.frame} (or coercible object) with the same
#'   column order as the training \code{X}.
#' @param missing_action Character; one of \code{"na"} (default) or
#'   \code{"majority"}.  See Details.
#' @return A \code{data.frame} with one row per row of \code{newdata} and
#'   columns:
#' \describe{
#'   \item{\code{row_id}}{Integer; positional row index in \code{newdata}
#'     (1 to \code{nrow(newdata)}).}
#'   \item{\code{endpoint_node_id}}{Integer; \code{node_id} of the terminal
#'     leaf reached by traversal.  \code{NA_integer_} when the observation
#'     cannot be classified (missing split attribute with
#'     \code{missing_action = "na"}, or no-tree fit).}
#'   \item{\code{endpoint_id}}{Integer; sequential endpoint index matching
#'     \code{\link{cta_endpoint_summary}}.  \code{NA_integer_} under the same
#'     conditions as \code{endpoint_node_id}.}
#' }
#' For no-tree fits all rows have \code{endpoint_node_id = NA_integer_} and
#' \code{endpoint_id = NA_integer_}.
#' @seealso \code{\link{oda_cta_fit}}, \code{\link{cta_endpoint_summary}},
#'   \code{\link{cta_propensity_weights}}
#' @examples
#' data(mtcars)
#' X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
#' y    <- as.integer(mtcars$am)
#' tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
#' ep   <- cta_assign_endpoints(tree, X)
#' head(ep)
#' @export
cta_assign_endpoints <- function(tree, newdata,
                                  missing_action = c("na", "majority")) {
  stopifnot(inherits(tree, "cta_tree"))
  missing_action <- match.arg(missing_action)
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)

  n_new <- nrow(newdata)

  # --- schema for no-tree or zero-row result ---------------------------------
  empty_df <- function() {
    data.frame(
      row_id           = integer(0),
      endpoint_node_id = integer(0),
      endpoint_id      = integer(0)
    )
  }

  # No-tree: all NA.
  if (isTRUE(tree$no_tree)) {
    return(data.frame(
      row_id           = seq_len(n_new),
      endpoint_node_id = rep(NA_integer_, n_new),
      endpoint_id      = rep(NA_integer_, n_new)
    ))
  }

  # --- split nodes (used for validation and column mapping) -------------------
  split_nodes <- Filter(function(nd) !isTRUE(nd$leaf), tree$nodes)

  # --- column routing: positional (same-width) or name-based (wide/reordered) --
  #
  # Positions are canonical under the hood; names are a safe cross-reference
  # layer used only when newdata has a different column count than training X.
  #
  #  * same-width  (ncol(newdata) == n_attrs): positional routing.  Identity map.
  #    Works for named or unnamed newdata; existing callers are unaffected.
  #
  #  * wide/reordered (ncol(newdata) != n_attrs): name-based routing required.
  #    Both names(newdata) and tree$attr_names must be present.
  #    Required split variable names must all be found in names(newdata).
  #    Never silently fallback to position for a missing named variable.
  tree_names <- tree$attr_names
  nd_names   <- names(newdata)
  n_tr       <- tree$n_attrs %||% 0L
  n_nd       <- ncol(newdata)

  if (n_nd == n_tr) {
    # Same-width: positional (identity map).
    attr_col_in_newdata <- seq_len(n_nd)
    wide_routing        <- FALSE
  } else {
    # Wide / reordered: name-based mapping required — no silent positional fallback.
    if (is.null(nd_names) || length(nd_names) == 0L) {
      stop(sprintf(
        "cta_assign_endpoints: newdata has %d column(s) but training X had %d. ",
        n_nd, n_tr
      ), "Name-based routing is required for wide newdata, but names(newdata) is absent.")
    }
    if (is.null(tree_names) || length(tree_names) == 0L) {
      stop(sprintf(
        "cta_assign_endpoints: newdata has %d column(s) but training X had %d. ",
        n_nd, n_tr
      ), "Name-based routing is required, but the fitted tree has no attr_names.")
    }
    nd_name_pos <- setNames(seq_along(nd_names), nd_names)
    attr_col_in_newdata <- vapply(seq_along(tree_names), function(j) {
      nm <- tree_names[j]
      if (!is.na(nm) && nm %in% names(nd_name_pos)) {
        nd_name_pos[[nm]]
      } else {
        NA_integer_   # absent: caught by validation below
      }
    }, integer(1L))
    wide_routing <- TRUE
  }

  # --- column validation (runs even for zero-row newdata) ---------------------
  split_attr_cols <- unique(vapply(split_nodes,
                                   function(nd) nd$attr_col %||% NA_integer_,
                                   integer(1L)))
  split_attr_cols <- split_attr_cols[!is.na(split_attr_cols)]

  if (wide_routing) {
    # Name-based: every required split attribute must map to a column in newdata.
    valid_sac <- split_attr_cols[split_attr_cols >= 1L & split_attr_cols <= length(attr_col_in_newdata)]
    missing_idx <- valid_sac[is.na(attr_col_in_newdata[valid_sac])]
    if (length(missing_idx) > 0L) {
      missing_names <- tree_names[missing_idx]
      missing_names <- missing_names[!is.na(missing_names)]
      stop(sprintf(
        "cta_assign_endpoints: newdata is missing required split variable(s): %s.",
        paste(missing_names, collapse = ", ")
      ))
    }
  } else {
    # Positional: check newdata has enough columns.
    max_attr_col <- if (length(split_attr_cols) == 0L) 0L else max(split_attr_cols)
    if (ncol(newdata) < max_attr_col) {
      stop(sprintf(
        "cta_assign_endpoints: newdata has %d column(s) but tree requires at least %d.",
        ncol(newdata), max_attr_col
      ))
    }
  }

  if (n_new == 0L) return(empty_df())

  # --- build endpoint_node_id -> endpoint_id lookup --------------------------
  es <- cta_endpoint_summary(tree)
  ep_lookup <- integer(0)
  if (nrow(es) > 0L) {
    ep_lookup <- setNames(es$endpoint_id, as.character(es$endpoint_node_id))
  }

  miss_codes <- tree$miss_codes

  # --- traversal helper: returns terminal node_id for one observation --------
  traverse_one <- function(i) {
    nid <- tree$root_id
    repeat {
      nd <- tree$nodes[[nid]]
      if (is.null(nd) || isTRUE(nd$leaf) || length(nd$child_ids) == 0L)
        return(nd$node_id %||% NA_integer_)

      j        <- nd$attr_col
      j_actual <- if (!is.na(j) && j >= 1L && j <= length(attr_col_in_newdata))
                    attr_col_in_newdata[j] else j
      x_val    <- newdata[[j_actual]][i]

      miss_here <- is.na(x_val)
      if (!is.null(miss_codes) && !miss_here) miss_here <- x_val %in% miss_codes
      if (miss_here) {
        if (missing_action == "na") return(NA_integer_)
        # majority: route to child with largest n_obs; ties -> first child
        child_ns <- vapply(nd$child_ids, function(cid) {
          cn <- tree$nodes[[cid]]
          if (is.null(cn)) NA_integer_ else cn$n_obs %||% NA_integer_
        }, integer(1L))
        nid <- nd$child_ids[which.max(child_ns)]
        next
      }

      rule <- nd$rule
      if (!is.null(rule$type) &&
          rule$type %in% c("multiclass_ordered", "multiclass_nominal")) {
        y_hat <- as.integer(oda_rule_predict_multiclass(
          x_val, rule, boundary = rule$boundary %||% "megaoda_halfopen"))
        sl <- nd$split_labels
        ic <- which(sl == y_hat)
        if (length(ic) == 0L) {
          # unrecognised label: route to largest child
          child_ns <- vapply(nd$child_ids, function(cid) {
            cn <- tree$nodes[[cid]]
            if (is.null(cn)) NA_integer_ else cn$n_obs %||% NA_integer_
          }, integer(1L))
          nid <- nd$child_ids[which.max(child_ns)]
          next
        }
      } else {
        y_hat <- as.integer(oda_rule_predict(x_val, rule))
        if (!(y_hat %in% 0:1)) {
          child_ns <- vapply(nd$child_ids, function(cid) {
            cn <- tree$nodes[[cid]]
            if (is.null(cn)) NA_integer_ else cn$n_obs %||% NA_integer_
          }, integer(1L))
          nid <- nd$child_ids[which.max(child_ns)]
          next
        }
        ic <- y_hat + 1L
      }
      next_nid <- nd$child_ids[ic[1L]]
      if (is.na(next_nid) || next_nid < 1L) return(nd$node_id %||% NA_integer_)
      nid <- next_nid
    }
  }

  # --- traverse all rows -----------------------------------------------------
  endpoint_node_ids <- vapply(seq_len(n_new), traverse_one, integer(1L))
  endpoint_ids <- ifelse(
    is.na(endpoint_node_ids),
    NA_integer_,
    as.integer(ep_lookup[as.character(endpoint_node_ids)])
  )

  data.frame(
    row_id           = seq_len(n_new),
    endpoint_node_id = endpoint_node_ids,
    endpoint_id      = endpoint_ids,
    stringsAsFactors = FALSE
  )
}

###############################################################################
# cta_observation_weights()
# Convenience wrapper: assigns each observation to its endpoint via
# cta_assign_endpoints(), retrieves endpoint-level propensity weights via
# cta_propensity_weights(), and joins the two on (endpoint_id, class).
# Returns one row per row of newdata — never fewer, never more.
###############################################################################

#' Assign per-observation CTA propensity weights
#'
#' Convenience wrapper that calls \code{\link{cta_assign_endpoints}} and
#' \code{\link{cta_propensity_weights}} and returns a joined
#' observation-level data frame.  The \code{cta_tree} object is not mutated;
#' all computation is on-demand.
#'
#' \strong{Column order requirement:} Same as \code{\link{cta_assign_endpoints}}
#' — \code{newdata} must have the same attribute column order as the \code{X}
#' matrix passed to \code{\link{oda_cta_fit}}.
#'
#' Observations with \code{NA} endpoint (missing root split attribute under
#' \code{missing_action = "na"}) or \code{NA} class label receive
#' \code{assigned = FALSE} and \code{NA} for all weight columns.
#'
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @param newdata A \code{data.frame} (or coercible object) with the same
#'   column order as the training \code{X}.
#' @param y Class labels for each row of \code{newdata}.  Any type coercible
#'   to character; length must equal \code{nrow(newdata)}.
#' @param target_class Passed to \code{\link{cta_propensity_weights}} as an
#'   annotation parameter.  It identifies which class is treated as the design
#'   target (high-risk class) for the \code{target_class} output column; it does
#'   \emph{not} filter the endpoint × class rows used for the join.  Each
#'   observation is matched to its own \code{actual_class} regardless of this
#'   value.  \code{NULL} (default) lets \code{cta_propensity_weights} resolve
#'   the target class automatically (numerically largest class for binary trees;
#'   required explicitly for trees with three or more classes).
#' @param adjusted Passed to \code{\link{cta_propensity_weights}}.
#'   Default \code{TRUE}.
#' @param missing_action Passed to \code{\link{cta_assign_endpoints}}.
#'   One of \code{"na"} (default) or \code{"majority"}.
#' @return A \code{data.frame} with \code{nrow(newdata)} rows and columns:
#' \describe{
#'   \item{\code{row_id}}{Integer; positional row index (1 to
#'     \code{nrow(newdata)}).}
#'   \item{\code{actual_class}}{Character; class label from \code{y}.}
#'   \item{\code{endpoint_node_id}}{Integer; node ID of terminal leaf, or
#'     \code{NA_integer_} when unroutable.}
#'   \item{\code{endpoint_id}}{Integer; sequential endpoint index matching
#'     \code{\link{cta_endpoint_summary}}, or \code{NA_integer_}.}
#'   \item{\code{target_class}}{Integer; resolved design target class from
#'     \code{\link{cta_propensity_weights}}, or \code{NA_integer_} when
#'     unassigned.}
#'   \item{\code{propensity_weight}}{Numeric; unadjusted propensity weight for
#'     the observation's endpoint–class cell, or \code{NA}.}
#'   \item{\code{adjusted_propensity_weight}}{Numeric; adjusted propensity
#'     weight, or \code{NA}.}
#'   \item{\code{undefined_empirical}}{Logical; \code{TRUE} when the
#'     endpoint–class cell has zero observed frequency, or \code{NA}.}
#'   \item{\code{perfectly_predicted_endpoint}}{Logical; \code{TRUE} when the
#'     endpoint is perfectly predicted, or \code{NA}.}
#'   \item{\code{adjusted}}{Logical; \code{TRUE} when the adjusted weight was
#'     applied, or \code{NA}.}
#'   \item{\code{assigned}}{Logical; \code{TRUE} when a propensity weight was
#'     successfully matched.}
#' }
#' @seealso \code{\link{cta_assign_endpoints}},
#'   \code{\link{cta_propensity_weights}}, \code{\link{oda_cta_fit}}
#' @export
cta_observation_weights <- function(
    tree,
    newdata,
    y,
    target_class   = NULL,
    adjusted       = TRUE,
    missing_action = c("na", "majority")
) {
  stopifnot(inherits(tree, "cta_tree"))
  missing_action <- match.arg(missing_action)
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  n <- nrow(newdata)
  if (length(y) != n)
    stop("cta_observation_weights: length(y) must equal nrow(newdata).",
         call. = FALSE)

  actual_class <- as.character(y)

  ep <- cta_assign_endpoints(tree, newdata, missing_action = missing_action)
  pw <- cta_propensity_weights(tree, target_class = target_class,
                               adjusted = adjusted)

  # Build composite key (endpoint_id + "\001" + class) for order-preserving join
  pw_key <- if (nrow(pw) > 0L)
    paste(pw$endpoint_id, pw$class, sep = "\001")
  else
    character(0L)

  ep_key  <- paste(ep$endpoint_id, actual_class, sep = "\001")
  pw_idx  <- match(ep_key, pw_key)

  # NA endpoint or NA class → no match possible
  pw_idx[is.na(ep$endpoint_id)] <- NA_integer_
  pw_idx[is.na(y)]              <- NA_integer_

  # Warn for classified obs whose class is not represented in pw
  unmatched <- !is.na(ep$endpoint_id) & !is.na(y) & is.na(pw_idx)
  if (any(unmatched))
    warning(sprintf(
      paste0("cta_observation_weights: %d classified observation(s) could not",
             " be matched to endpoint-class propensity weights"),
      sum(unmatched)
    ), call. = FALSE)

  assigned <- !is.na(pw_idx)

  # Helper: index into pw column, preserving NA for unmatched rows
  get_col <- function(col) { v <- pw[[col]]; v[pw_idx] }

  data.frame(
    row_id                    = ep$row_id,
    actual_class              = actual_class,
    endpoint_node_id          = ep$endpoint_node_id,
    endpoint_id               = ep$endpoint_id,
    target_class              = get_col("target_class"),
    propensity_weight         = get_col("propensity_weight"),
    adjusted_propensity_weight = get_col("adjusted_propensity_weight"),
    undefined_empirical       = get_col("undefined_empirical"),
    perfectly_predicted_endpoint = get_col("perfectly_predicted_endpoint"),
    adjusted                  = get_col("adjusted"),
    assigned                  = assigned,
    stringsAsFactors          = FALSE
  )
}

# =============================================================================
# Internal helpers for cta_plot_data() / plot.cta_tree()
# =============================================================================

# Map a class integer to a display label string.
# class_labels: NULL → "class=C"; named vector → matched by name;
# positional vector → index cls+1.
.class_label <- function(cls, class_labels) {
  if (is.null(class_labels)) return(paste0("class=", cls))
  nm <- names(class_labels)
  if (!is.null(nm)) {
    idx <- match(as.character(cls), nm)
    if (!is.na(idx)) return(class_labels[[idx]])
  } else {
    idx <- as.integer(cls) + 1L
    if (idx >= 1L && idx <= length(class_labels)) return(class_labels[[idx]])
  }
  paste0("class=", cls)
}

# Validate and normalise endpoint_palette to a function(n) -> character(n).
# NULL   → default gradient (warm pink → pale yellow → soft green).
# fn     → used as-is.
# vector → interpolated via colorRampPalette.
.build_palette_fn <- function(endpoint_palette) {
  if (is.null(endpoint_palette))
    return(grDevices::colorRampPalette(c("#ffebee", "#fffde7", "#e8f5e9")))
  if (is.function(endpoint_palette))
    return(endpoint_palette)
  if (is.character(endpoint_palette))
    return(grDevices::colorRampPalette(endpoint_palette))
  stop("endpoint_palette must be NULL, a function, or a character vector",
       call. = FALSE)
}

# =============================================================================
# cta_plot_data() — pure layout-data contract for tree visualisation
# =============================================================================

#' Extract layout data for plotting a CTA tree
#'
#' Returns a pure-data list describing tree topology and layout coordinates.
#' No graphics are produced.  Use this as input to \code{\link{plot.cta_tree}}
#' or to custom rendering code.
#'
#' \strong{Layout:} leaves receive sequential integer x-positions in
#' depth-first (left-to-right) traversal order; internal nodes are centred
#' over their children.  \code{y = -depth} so the root sits at the top.
#'
#' \strong{Target-class enrichment:} when \code{target_class} is supplied,
#' each terminal leaf is joined to \code{\link{cta_staging_table}} and
#' annotated with target-class counts, proportions, and a continuous display
#' color derived from the endpoint's rank among all endpoints by ascending
#' target-class proportion.  Colors encode relative position within this
#' tree's endpoint distribution and do \emph{not} imply clinical thresholds
#' or categories.
#'
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @param target_class Integer (or \code{NULL}); the class label treated as
#'   the target when enriching endpoint annotations.  \code{NULL} (default)
#'   returns the structural layout only.  For binary trees with
#'   \code{target_class = NULL} the enrichment is skipped entirely — no
#'   endpoint columns are added.
#' @param class_labels Optional character vector of display names for class
#'   labels.  Supply as a \emph{named} vector, e.g.
#'   \code{c("0" = "Alive", "1" = "Deceased")}, or as a positional vector
#'   (index \eqn{k+1} maps to class \eqn{k}).  \code{NULL} (default) uses
#'   \code{"class=C"} formatting.
#' @param digits Integer number of decimal places for percentage formatting
#'   in \code{endpoint_label}.  Default \code{1}.
#' @param endpoint_palette Palette for endpoint fill colors, used only when
#'   \code{target_class} is supplied.  Accepts \code{NULL} (default gradient),
#'   a palette \code{function(n)} returning \code{n} color strings, or a
#'   character vector of colors interpolated via
#'   \code{\link[grDevices]{colorRampPalette}}.
#' @return When \code{target_class = NULL}: a list with elements
#'   \code{nodes}, \code{edges}, \code{no_tree}, \code{has_weights}.
#'   When \code{target_class} is supplied: the same list plus an
#'   \code{endpoints} data.frame and \code{target_class_used} integer.
#' \describe{
#'   \item{\code{nodes}}{A \code{data.frame} with one row per node.
#'     Always-present columns: \code{node_id} (integer), \code{parent_id}
#'     (integer), \code{depth} (integer), \code{x} (numeric), \code{y}
#'     (numeric), \code{leaf} (logical), \code{attribute} (character; NA
#'     for leaves), \code{n_obs} (integer), \code{majority_class} (integer),
#'     \code{ess} (numeric; NA for leaves), \code{label} (character display
#'     text).  Additional columns when \code{target_class} supplied (NA on
#'     split nodes): \code{endpoint_id} (integer), \code{stage} (integer),
#'     \code{target_class} (integer), \code{target_n} (numeric),
#'     \code{denominator} (numeric), \code{target_proportion} (numeric),
#'     \code{target_rank} (integer; ascending rank of proportion, ties
#'     broken by \code{ties.method = "first"}), \code{endpoint_fill_color}
#'     (character hex color), \code{predicted_label} (character),
#'     \code{target_label} (character), \code{endpoint_label} (character
#'     multi-line display text for this endpoint).}
#'   \item{\code{edges}}{A \code{data.frame} with one row per parent-to-child
#'     edge: \code{from_node_id}, \code{to_node_id} (integer),
#'     \code{x0}, \code{y0}, \code{x1}, \code{y1} (numeric), \code{label}
#'     (character branch condition, e.g. \code{"V14<=0.5"}).}
#'   \item{\code{endpoints}}{(target_class only) Staging table joined with
#'     leaf layout coordinates.  One row per endpoint, ordered by stage.}
#'   \item{\code{target_class_used}}{(target_class only) The integer
#'     target_class argument used.}
#'   \item{\code{no_tree}}{Logical; \code{TRUE} for leaf-only fits.}
#'   \item{\code{has_weights}}{Logical; \code{TRUE} when case weights active.}
#' }
#' @seealso \code{\link{plot.cta_tree}}, \code{\link{cta_staging_table}},
#'   \code{\link{oda_cta_fit}}
#' @export
cta_plot_data <- function(tree, target_class = NULL, class_labels = NULL,
                          digits = 1, endpoint_palette = NULL) {
  stopifnot(inherits(tree, "cta_tree"))

  no_tree     <- isTRUE(tree$no_tree)
  has_weights <- isTRUE(tree$has_weights)

  empty_nodes <- data.frame(
    node_id        = integer(0),
    parent_id      = integer(0),
    depth          = integer(0),
    x              = numeric(0),
    y              = numeric(0),
    leaf           = logical(0),
    attribute      = character(0),
    n_obs          = integer(0),
    majority_class = integer(0),
    ess            = numeric(0),
    label          = character(0),
    stringsAsFactors = FALSE
  )
  empty_edges <- data.frame(
    from_node_id = integer(0),
    to_node_id   = integer(0),
    x0           = numeric(0),
    y0           = numeric(0),
    x1           = numeric(0),
    y1           = numeric(0),
    label        = character(0),
    stringsAsFactors = FALSE
  )

  all_nodes <- Filter(Negate(is.null), tree$nodes)
  if (no_tree || length(all_nodes) == 0L)
    return(list(nodes = empty_nodes, edges = empty_edges,
                no_tree = no_tree, has_weights = has_weights))

  # ------------------------------------------------------------------
  # Assign layout x-positions via depth-first (left-to-right) traversal.
  # Leaves receive sequential integers 1, 2, 3, ...
  # Internal nodes receive the mean of their children's x-positions.
  # ------------------------------------------------------------------
  x_pos        <- numeric(0)   # named by node_id character
  leaf_counter <- 0L

  assign_x <- function(nid) {
    nd <- tree$nodes[[nid]]
    if (is.null(nd)) return(NA_real_)
    is_leaf <- isTRUE(nd$leaf) || length(nd$child_ids) == 0L
    if (is_leaf) {
      leaf_counter <<- leaf_counter + 1L
      x_pos[as.character(nid)] <<- as.numeric(leaf_counter)
      return(as.numeric(leaf_counter))
    }
    child_xs <- vapply(nd$child_ids, assign_x, numeric(1L))
    cx <- mean(child_xs, na.rm = TRUE)
    x_pos[as.character(nid)] <<- cx
    cx
  }

  assign_x(tree$root_id)

  # ------------------------------------------------------------------
  # Build nodes data.frame
  # ------------------------------------------------------------------
  n          <- length(all_nodes)
  nid_vec    <- integer(n)
  pid_vec    <- integer(n)
  depth_vec  <- integer(n)
  x_vec      <- numeric(n)
  y_vec      <- numeric(n)
  leaf_vec   <- logical(n)
  attr_vec   <- character(n)
  nobs_vec   <- integer(n)
  maj_vec    <- integer(n)
  ess_vec    <- numeric(n)
  lbl_vec    <- character(n)

  ess_label <- if (has_weights) "WESS" else "ESS"

  for (i in seq_along(all_nodes)) {
    nd  <- all_nodes[[i]]
    nid <- nd$node_id
    is_leaf  <- isTRUE(nd$leaf) || length(nd$child_ids) == 0L
    ess_val  <- if (is_leaf) NA_real_ else (nd$ess %||% NA_real_)
    x_i      <- x_pos[as.character(nid)] %||% NA_real_
    dep_i    <- nd$depth %||% NA_integer_

    nid_vec[i]   <- nid
    pid_vec[i]   <- nd$parent_id %||% NA_integer_
    depth_vec[i] <- dep_i
    x_vec[i]     <- x_i
    y_vec[i]     <- -dep_i
    leaf_vec[i]  <- is_leaf
    attr_vec[i]  <- if (!is_leaf) (nd$attribute %||% NA_character_) else NA_character_
    nobs_vec[i]  <- nd$n_obs %||% NA_integer_
    maj_vec[i]   <- nd$majority_class %||% NA_integer_
    ess_vec[i]   <- ess_val
    lbl_vec[i]   <- if (is_leaf) {
      sprintf("class=%d\nn=%d", nd$majority_class %||% NA_integer_,
              nd$n_obs %||% NA_integer_)
    } else if (!is.na(ess_val)) {
      sprintf("%s\n%s=%.1f%%\nn=%d",
              nd$attribute %||% "?", ess_label, ess_val,
              nd$n_obs %||% NA_integer_)
    } else {
      sprintf("%s\nn=%d", nd$attribute %||% "?", nd$n_obs %||% NA_integer_)
    }
  }

  nodes_df <- data.frame(
    node_id        = nid_vec,
    parent_id      = pid_vec,
    depth          = depth_vec,
    x              = x_vec,
    y              = y_vec,
    leaf           = leaf_vec,
    attribute      = attr_vec,
    n_obs          = nobs_vec,
    majority_class = maj_vec,
    ess            = ess_vec,
    label          = lbl_vec,
    stringsAsFactors = FALSE
  )

  # ------------------------------------------------------------------
  # Build edges data.frame
  # ------------------------------------------------------------------
  split_nds <- Filter(
    function(nd) !is.null(nd) && !isTRUE(nd$leaf) && length(nd$child_ids) > 0L,
    tree$nodes
  )

  n_edges    <- sum(vapply(split_nds, function(nd) length(nd$child_ids), integer(1L)))
  from_vec   <- integer(n_edges)
  to_vec     <- integer(n_edges)
  x0_vec     <- numeric(n_edges)
  y0_vec     <- numeric(n_edges)
  x1_vec     <- numeric(n_edges)
  y1_vec     <- numeric(n_edges)
  elbl_vec   <- character(n_edges)

  ei <- 0L
  for (nd in split_nds) {
    px <- x_pos[as.character(nd$node_id)] %||% NA_real_
    py <- -(nd$depth %||% 1L)
    for (k in seq_along(nd$child_ids)) {
      cid <- nd$child_ids[k]
      cnd <- tree$nodes[[cid]]
      if (is.null(cnd)) next
      ei <- ei + 1L
      cx <- x_pos[as.character(cid)] %||% NA_real_
      cy <- -(cnd$depth %||% 2L)
      from_vec[ei]  <- nd$node_id
      to_vec[ei]    <- cid
      x0_vec[ei]    <- px
      y0_vec[ei]    <- py
      x1_vec[ei]    <- cx
      y1_vec[ei]    <- cy
      elbl_vec[ei]  <- .cta_branch_string(nd, k)
    }
  }

  if (ei < n_edges) {
    keep     <- seq_len(ei)
    from_vec <- from_vec[keep]; to_vec   <- to_vec[keep]
    x0_vec   <- x0_vec[keep];   y0_vec   <- y0_vec[keep]
    x1_vec   <- x1_vec[keep];   y1_vec   <- y1_vec[keep]
    elbl_vec <- elbl_vec[keep]
  }

  edges_df <- data.frame(
    from_node_id = from_vec,
    to_node_id   = to_vec,
    x0           = x0_vec,
    y0           = y0_vec,
    x1           = x1_vec,
    y1           = y1_vec,
    label        = elbl_vec,
    stringsAsFactors = FALSE
  )

  # ------------------------------------------------------------------
  # Target-class enrichment (only when target_class is supplied)
  # ------------------------------------------------------------------
  if (!is.null(target_class)) {
    target_class_int <- as.integer(target_class)
    st   <- cta_staging_table(tree, target_class = target_class_int)
    tgt_label <- .class_label(target_class_int, class_labels)
    n_ep <- nrow(st)

    if (n_ep > 0L) {
      # Palette: build from user spec, then assign n_ep colors ordered by rank
      pal_fn <- .build_palette_fn(endpoint_palette)
      colors <- pal_fn(n_ep)
      if (length(colors) < n_ep) colors <- rep_len(colors, n_ep)

      # Rank endpoints by ascending target_proportion; ties.method="first"
      # (deterministic: resolved by enumeration order of staging table rows,
      # which are ordered by node_id through cta_endpoint_summary).
      st$target_rank         <- as.integer(rank(st$target_proportion,
                                                  ties.method = "first"))
      st$endpoint_fill_color <- colors[st$target_rank]
      st$predicted_label     <- vapply(st$terminal_prediction,
                                       function(mc) .class_label(mc, class_labels),
                                       character(1L))
      st$target_label        <- tgt_label

      # x/y layout coords for the endpoints data.frame
      xy        <- nodes_df[, c("node_id", "x", "y")]
      st$x      <- xy$x[match(st$endpoint_node_id, xy$node_id)]
      st$y      <- xy$y[match(st$endpoint_node_id, xy$node_id)]
      st$node_id <- st$endpoint_node_id

      ep_cols <- c("endpoint_id", "endpoint_node_id", "path", "stage",
                   "node_id", "target_class", "target_n", "denominator",
                   "target_proportion", "target_rank", "terminal_prediction",
                   "predicted_label", "target_label", "endpoint_fill_color",
                   "x", "y")
      endpoints_df <- st[order(st$stage), ep_cols]
      rownames(endpoints_df) <- NULL

      # Join target fields to nodes by node_id == endpoint_node_id
      li <- match(nodes_df$node_id, st$endpoint_node_id)

      nodes_df$endpoint_id         <- st$endpoint_id[li]
      nodes_df$stage               <- st$stage[li]

      tc_vec <- rep(NA_integer_, n)
      tc_vec[!is.na(li)] <- target_class_int
      nodes_df$target_class <- tc_vec

      nodes_df$target_n            <- st$target_n[li]
      nodes_df$denominator         <- st$denominator[li]
      nodes_df$target_proportion   <- st$target_proportion[li]
      nodes_df$target_rank         <- as.integer(st$target_rank[li])
      nodes_df$endpoint_fill_color <- st$endpoint_fill_color[li]

      nodes_df$predicted_label <- vapply(nodes_df$majority_class,
                                         function(mc) .class_label(mc, class_labels),
                                         character(1L))
      tl_vec <- rep(NA_character_, n)
      tl_vec[!is.na(li)] <- tgt_label
      nodes_df$target_label <- tl_vec

      # Assemble endpoint_label for leaf endpoints; NA for split nodes
      has_ep  <- !is.na(li) & nodes_df$leaf
      ep_lbl  <- rep(NA_character_, n)
      if (any(has_ep)) {
        ep_lbl[has_ep] <- sprintf(
          "%s\n%s%%, n=%d/%d\npred: %s\nStage %d",
          tgt_label,
          formatC(nodes_df$target_proportion[has_ep] * 100,
                  digits = digits, format = "f"),
          as.integer(round(nodes_df$target_n[has_ep])),
          as.integer(round(nodes_df$denominator[has_ep])),
          nodes_df$predicted_label[has_ep],
          as.integer(nodes_df$stage[has_ep])
        )
      }
      nodes_df$endpoint_label  <- ep_lbl
      nodes_df$label[has_ep]   <- ep_lbl[has_ep]

    } else {
      # Guard: valid tree should always have endpoints; add NA columns
      na_int <- rep(NA_integer_, n)
      na_dbl <- rep(NA_real_,    n)
      na_chr <- rep(NA_character_, n)
      nodes_df$endpoint_id         <- na_int
      nodes_df$stage               <- na_int
      nodes_df$target_class        <- na_int
      nodes_df$target_n            <- na_dbl
      nodes_df$denominator         <- na_dbl
      nodes_df$target_proportion   <- na_dbl
      nodes_df$target_rank         <- na_int
      nodes_df$endpoint_fill_color <- na_chr
      nodes_df$predicted_label     <- na_chr
      nodes_df$target_label        <- na_chr
      nodes_df$endpoint_label      <- na_chr
      endpoints_df <- data.frame(
        endpoint_id = integer(0), endpoint_node_id = integer(0),
        path = character(0), stage = integer(0), node_id = integer(0),
        target_class = integer(0), target_n = numeric(0),
        denominator = numeric(0), target_proportion = numeric(0),
        target_rank = integer(0), terminal_prediction = integer(0),
        predicted_label = character(0), target_label = character(0),
        endpoint_fill_color = character(0),
        x = numeric(0), y = numeric(0),
        stringsAsFactors = FALSE
      )
    }

    return(list(
      nodes             = nodes_df,
      edges             = edges_df,
      endpoints         = endpoints_df,
      no_tree           = no_tree,
      has_weights       = has_weights,
      target_class_used = target_class_int
    ))
  }

  list(nodes = nodes_df, edges = edges_df,
       no_tree = no_tree, has_weights = has_weights)
}

# =============================================================================
# plot.cta_tree() — base-R tree diagram
# =============================================================================

#' Plot a fitted CTA tree
#'
#' Produces a base-R tree diagram.  Calls \code{\link{cta_plot_data}} for
#' layout; uses only base graphics — no external package dependencies.
#'
#' Split (internal) nodes show the split attribute, node-level ESS or WESS,
#' and observation count.  Without \code{target_class}, leaf nodes show the
#' majority-class prediction and observation count.  With \code{target_class},
#' leaf nodes show the target-class count, percentage, predicted class, and
#' stage from \code{\link{cta_staging_table}}.  Edge labels show the branch
#' condition (e.g. \code{"V14<=0.5"}).
#'
#' \strong{Color note:} when \code{target_class} is supplied, endpoint fill
#' colors are assigned by ascending rank of each endpoint's target-class
#' proportion within this tree.  Colors encode relative position in the
#' endpoint distribution and do \emph{not} imply clinical thresholds or
#' categories.  Supply a custom palette via \code{endpoint_palette} to change
#' the color encoding.
#'
#' @param x A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @param target_class Integer target class for endpoint annotation; passed
#'   to \code{\link{cta_plot_data}}.  \code{NULL} (default) produces a
#'   structural plot without endpoint enrichment.
#' @param class_labels Optional display names; passed to
#'   \code{\link{cta_plot_data}}.
#' @param digits Decimal places for percentage labels; passed to
#'   \code{\link{cta_plot_data}}.  Default \code{1}.
#' @param main Character plot title.  Default \code{"CTA Tree"}.
#' @param show_counts Logical; include \code{n=a/b} counts in endpoint labels
#'   when \code{target_class} is supplied.  Default \code{TRUE}.
#' @param show_stage Logical; include \code{Stage s} line in endpoint labels
#'   when \code{target_class} is supplied.  Default \code{TRUE}.
#' @param endpoint_palette Palette for endpoint fill colors when
#'   \code{target_class} is supplied; passed to \code{\link{cta_plot_data}}.
#'   \code{NULL} uses the default gradient. Accepts a palette function or a
#'   character vector of colors.
#' @param endpoint_fill Default fill colour for leaf nodes when
#'   \code{target_class} is \code{NULL}.  Default \code{"#D9F7E6"}.
#' @param split_fill Fill colour for split (internal) nodes.
#'   Default \code{"#D9EAF7"}.
#' @param node_col_split Legacy alias for \code{split_fill}; overrides it
#'   when non-\code{NULL}.
#' @param node_col_leaf Legacy alias for \code{endpoint_fill}; overrides it
#'   when non-\code{NULL}.
#' @param edge_col Colour for edge lines and arrowheads.  Default \code{"grey40"}.
#' @param border_col Border colour for all nodes.  Default \code{"grey30"}.
#' @param text_col Text colour for node labels.  Default \code{"black"}.
#' @param arrow_col Arrowhead and line colour for directed edges.  \code{NULL}
#'   (default) uses \code{edge_col}.
#' @param show_caption Logical; if \code{TRUE} and \code{target_class} is
#'   supplied, adds a bottom caption: \emph{"Endpoint fill: relative
#'   target-class proportion within this tree. Not a clinical threshold."}
#'   Default \code{FALSE}.
#' @param cex Text expansion factor for node labels.  Default \code{0.75}.
#' @param ... Unused; included for S3 compatibility.
#' @return \code{invisible(pd)}, where \code{pd} is the \code{\link{cta_plot_data}}
#'   list used to render the plot.  The caller can inspect layout coordinates
#'   and endpoint annotations from the returned object.
#' @seealso \code{\link{cta_plot_data}}, \code{\link{cta_staging_table}},
#'   \code{\link{oda_cta_fit}}
#' @export
plot.cta_tree <- function(x,
                          target_class     = NULL,
                          class_labels     = NULL,
                          digits           = 1,
                          main             = "CTA Tree",
                          show_counts      = TRUE,
                          show_stage       = TRUE,
                          endpoint_palette = NULL,
                          endpoint_fill    = "#D9F7E6",
                          split_fill       = "#D9EAF7",
                          node_col_split   = NULL,
                          node_col_leaf    = NULL,
                          edge_col         = "grey40",
                          border_col       = "grey30",
                          text_col         = "black",
                          arrow_col        = NULL,
                          show_caption     = FALSE,
                          cex              = 0.75, ...) {
  # Legacy color-arg compat: non-NULL node_col_* overrides the new names
  split_fill_used   <- node_col_split %||% split_fill
  leaf_fill_default <- node_col_leaf  %||% endpoint_fill
  arrow_col_used    <- arrow_col %||% edge_col

  pd <- cta_plot_data(x, target_class = target_class,
                      class_labels = class_labels, digits = digits,
                      endpoint_palette = endpoint_palette)

  has_target <- !is.null(pd$target_class_used)

  if (pd$no_tree || nrow(pd$nodes) == 0L) {
    graphics::plot.new()
    graphics::title(main = main)
    graphics::text(0.5, 0.5, "no tree (leaf only)", cex = 1.2, col = "grey50")
    return(invisible(pd))
  }

  nd <- pd$nodes
  ed <- pd$edges

  xr  <- range(nd$x, na.rm = TRUE)
  yr  <- range(nd$y, na.rm = TRUE)
  if (diff(xr) < 1) xr <- xr + c(-0.5, 0.5)
  if (diff(yr) < 1) yr <- yr + c(-0.5, 0.5)
  xpad <- diff(xr) * 0.12 + 0.6
  ypad <- diff(yr) * 0.14 + 0.5

  hw <- 0.40
  hh <- if (has_target) 0.40 else 0.30

  # Ellipse polygon vertices for split (internal) nodes.
  .ellipse_xy <- function(cx, cy, a = hw, b = hh, n = 72L) {
    theta <- seq(0, 2 * pi, length.out = n)
    list(x = cx + a * cos(theta), y = cy + b * sin(theta))
  }

  op <- graphics::par(mar = c(2, 1, 3, 1))
  on.exit(graphics::par(op), add = TRUE)

  graphics::plot.new()
  graphics::plot.window(
    xlim = xr + c(-xpad, xpad),
    ylim = yr + c(-ypad, ypad)
  )
  graphics::title(main = main, cex.main = 1)

  # Directed arrows — drawn before nodes so shapes sit on top
  if (nrow(ed) > 0L) {
    for (i in seq_len(nrow(ed))) {
      graphics::arrows(
        ed$x0[i], ed$y0[i] - hh,
        ed$x1[i], ed$y1[i] + hh,
        length = 0.08, angle = 20,
        col = arrow_col_used, lwd = 1.2
      )
      mx <- (ed$x0[i] + ed$x1[i]) / 2
      my <- (ed$y0[i] - hh + ed$y1[i] + hh) / 2
      graphics::text(mx, my, labels = ed$label[i],
                     cex = cex * 0.82, col = "grey30", font = 1)
    }
  }

  # Nodes: ellipses for split nodes, rectangles for terminal leaves
  for (i in seq_len(nrow(nd))) {
    cx <- nd$x[i]
    cy <- nd$y[i]
    if (is.na(cx) || is.na(cy)) next

    if (nd$leaf[i]) {
      # Terminal endpoint: rectangle
      fc <- if (has_target && !is.na(nd$endpoint_fill_color[i]))
              nd$endpoint_fill_color[i]
            else
              leaf_fill_default

      display_lbl <- if (has_target && !is.na(nd$target_proportion[i])) {
        pct_str <- formatC(nd$target_proportion[i] * 100,
                           digits = digits, format = "f")
        parts <- c(
          nd$target_label[i],
          if (isTRUE(show_counts))
            sprintf("%s%%, n=%d/%d", pct_str,
                    as.integer(round(nd$target_n[i])),
                    as.integer(round(nd$denominator[i])))
          else
            paste0(pct_str, "%"),
          paste0("pred: ", nd$predicted_label[i])
        )
        if (isTRUE(show_stage) && !is.na(nd$stage[i]))
          parts <- c(parts, paste0("Stage ", as.integer(nd$stage[i])))
        paste(parts, collapse = "\n")
      } else {
        nd$label[i]
      }

      graphics::rect(cx - hw, cy - hh, cx + hw, cy + hh,
                     col = fc, border = border_col, lwd = 0.8)

    } else {
      # Split (internal) node: ellipse
      fc          <- split_fill_used
      display_lbl <- nd$label[i]

      ep <- .ellipse_xy(cx, cy)
      graphics::polygon(ep$x, ep$y,
                        col = fc, border = border_col, lwd = 0.8)
    }

    graphics::text(cx, cy, labels = display_lbl,
                   cex = cex, adj = c(0.5, 0.5), col = text_col)
  }

  # Optional caption explaining endpoint color semantics
  if (isTRUE(show_caption) && has_target) {
    graphics::mtext(
      paste0("Endpoint fill: relative target-class proportion within this tree.",
             " Not a clinical threshold."),
      side = 1, line = 0.8, cex = cex * 0.72, col = "grey40"
    )
  }

  invisible(pd)
}
