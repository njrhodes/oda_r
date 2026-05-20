###############################################################################
# R/cta_s3.R — S3 methods for cta_tree summary objects
#
# Classes produced here:
#   cta_tree_summary — structured output of summary.cta_tree()
#
# Public API:
#   summary.cta_tree()
#   print.cta_tree_summary()
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
