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
