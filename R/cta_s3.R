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

#' Terminal-leaf report table for a fitted CTA tree
#'
#' Returns one row per terminal leaf (endpoint) of a \code{cta_tree}.  All
#' values are read directly from stored node fields; no refitting or prediction
#' is performed.
#'
#' Per-endpoint class-count and target-class summaries are deferred to the
#' future CTA confusion/endpoint-count table because leaf class counts are not
#' currently stored at fit time.
#'
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @return A \code{data.frame} with one row per terminal leaf and columns:
#' \describe{
#'   \item{\code{node_id}}{Integer node identifier.}
#'   \item{\code{parent_id}}{Integer parent node identifier
#'     (\code{0} for the root leaf in a no-tree fit).}
#'   \item{\code{depth}}{Integer depth from root (root = 1).}
#'   \item{\code{path}}{Character; AND-joined branch labels from root to this
#'     leaf (e.g. \code{"V14<=0.5 AND V15>0.5"}).  \code{"/"} for the root
#'     leaf.}
#'   \item{\code{majority_class}}{Integer class label assigned to this leaf.}
#'   \item{\code{n_obs}}{Integer raw observation count at this leaf.}
#'   \item{\code{n_weighted}}{Numeric weighted observation count
#'     (\code{NA} when weights not active).}
#'   \item{\code{ess}}{Numeric ESS stored on this leaf (\code{NA_real_} —
#'     ESS is a split-node metric and is always \code{NA} for leaf nodes).}
#'   \item{\code{ess_weighted}}{Numeric WESS stored on this leaf
#'     (\code{NA_real_}).}
#'   \item{\code{loo_status}}{Character LOO status (\code{NA} for leaf nodes).}
#'   \item{\code{loo_ess}}{Numeric LOO ESS (\code{NA} for leaf nodes).}
#' }
#' For a no-tree fit the returned data frame has zero rows but the correct
#' column structure.
#' @seealso \code{\link{oda_cta_fit}}, \code{\link{summary.cta_tree}},
#'   \code{\link{cta_strata}}, \code{\link{cta_endpoint_denominators}}
#' @examples
#' data(mtcars)
#' X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
#' y    <- as.integer(mtcars$am)
#' tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
#' cta_endpoint_table(tree)
cta_endpoint_table <- function(tree) {
  stopifnot(inherits(tree, "cta_tree"))

  # ------------------------------------------------------------------
  # Empty data frame template — returned for no-tree or zero-leaf cases
  # ------------------------------------------------------------------
  empty_df <- function() {
    data.frame(
      node_id        = integer(0),
      parent_id      = integer(0),
      depth          = integer(0),
      path           = character(0),
      majority_class = integer(0),
      n_obs          = integer(0),
      n_weighted     = numeric(0),
      ess            = numeric(0),
      ess_weighted   = numeric(0),
      loo_status     = character(0),
      loo_ess        = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  if (isTRUE(tree$no_tree)) return(empty_df())

  leaves <- Filter(function(nd) !is.null(nd) && isTRUE(nd$leaf), tree$nodes)
  if (length(leaves) == 0L) return(empty_df())

  # ------------------------------------------------------------------
  # Build per-leaf rows
  # ------------------------------------------------------------------
  n <- length(leaves)
  node_id        <- integer(n)
  parent_id      <- integer(n)
  depth          <- integer(n)
  path           <- character(n)
  majority_class <- integer(n)
  n_obs          <- integer(n)
  n_weighted     <- numeric(n)
  ess            <- numeric(n)
  ess_weighted   <- numeric(n)
  loo_status     <- character(n)
  loo_ess        <- numeric(n)

  for (i in seq_len(n)) {
    nd <- leaves[[i]]
    segs <- .cta_path_segments(tree, nd$node_id)

    node_id[i]        <- nd$node_id
    parent_id[i]      <- nd$parent_id %||% NA_integer_
    depth[i]          <- nd$depth     %||% NA_integer_
    path[i]           <- if (length(segs) == 0L) "/"
                         else paste(segs, collapse = " AND ")
    majority_class[i] <- nd$majority_class %||% NA_integer_
    n_obs[i]          <- nd$n_obs          %||% NA_integer_
    n_weighted[i]     <- nd$n_weighted     %||% NA_real_
    ess[i]            <- nd$ess            %||% NA_real_
    ess_weighted[i]   <- nd$ess_weighted   %||% NA_real_
    loo_status[i]     <- nd$loo_status     %||% NA_character_
    loo_ess[i]        <- nd$loo_ess        %||% NA_real_
  }

  df <- data.frame(
    node_id        = node_id,
    parent_id      = parent_id,
    depth          = depth,
    path           = path,
    majority_class = majority_class,
    n_obs          = n_obs,
    n_weighted     = n_weighted,
    ess            = ess,
    ess_weighted   = ess_weighted,
    loo_status     = loo_status,
    loo_ess        = loo_ess,
    stringsAsFactors = FALSE
  )

  df[order(df$node_id), , drop = FALSE]
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
#' proportions, event rates, odds, or staging order.  Staging-table and
#' event-rate reporting require per-endpoint class counts, which are not
#' currently stored at fit time.  That capability is deferred until
#' fit-time endpoint class-count storage is designed.
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
#'   \code{\link{cta_strata}}, \code{\link{cta_endpoint_denominators}}
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
