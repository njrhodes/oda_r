# R/cta_family.R
# Read-only CTA endpoint accessors and MDSA family skeleton.
#
# These are reporting/comparison functions only.
# They read fields already stored on the cta_tree object (tree$nodes,
# tree$overall_ess, tree$has_weights).
# They must not call predict.cta_tree(), recompute confusion/ESS/WESS/D,
# or add new fields to cta_tree.

# ---- Leaf node helpers (internal) -------------------------------------------

# Return all terminal leaf nodes from a cta_tree$nodes list.
.cta_leaves <- function(tree) {
  Filter(function(nd) !is.null(nd) && isTRUE(nd$leaf), tree$nodes)
}

# ---- Exported read-only accessors -------------------------------------------

#' Number of terminal leaf endpoints in a CTA tree
#'
#' Counts the terminal (leaf) nodes in a fitted \code{cta_tree}.  Returns
#' \code{NA_integer_} for no-tree fits (where \code{tree$no_tree} is
#' \code{TRUE}).
#'
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @return Integer count of terminal leaf nodes, or \code{NA_integer_} for
#'   no-tree fits.
#' @seealso \code{\link{cta_endpoint_denominators}},
#'   \code{\link{cta_min_terminal_denom}}
#' @export
cta_strata <- function(tree) {
  stopifnot(inherits(tree, "cta_tree"))
  if (isTRUE(tree$no_tree)) return(NA_integer_)
  length(.cta_leaves(tree))
}

#' Terminal endpoint denominators of a CTA tree
#'
#' Returns the observation counts (\code{n_obs}) for each terminal leaf node,
#' named by node ID.  These are the raw row counts stored at fit time — they
#' are not recomputed from training data or predictions.
#'
#' Returns \code{integer(0)} for no-tree fits.
#'
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @return Named integer vector of leaf \code{n_obs} values, named by node ID
#'   (as character); \code{integer(0)} for no-tree fits.
#' @seealso \code{\link{cta_strata}}, \code{\link{cta_min_terminal_denom}}
#' @export
cta_endpoint_denominators <- function(tree) {
  stopifnot(inherits(tree, "cta_tree"))
  if (isTRUE(tree$no_tree)) return(integer(0))
  leaves <- .cta_leaves(tree)
  if (length(leaves) == 0L) return(integer(0))
  nobs <- vapply(leaves, function(nd) nd$n_obs,   integer(1L))
  ids  <- vapply(leaves, function(nd) nd$node_id, integer(1L))
  setNames(nobs, as.character(ids))
}

#' Minimum terminal endpoint denominator of a CTA tree
#'
#' Returns the smallest leaf \code{n_obs} across all terminal endpoints.
#' This value drives the next MINDENOM step in the MDSA descendant family:
#' \code{next_mindenom = cta_min_terminal_denom(tree) + 1L}.
#'
#' Returns \code{NA_integer_} for no-tree fits.
#'
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @return Minimum leaf \code{n_obs} as an integer, or \code{NA_integer_} for
#'   no-tree fits.
#' @seealso \code{\link{cta_strata}}, \code{\link{cta_endpoint_denominators}}
#' @export
cta_min_terminal_denom <- function(tree) {
  stopifnot(inherits(tree, "cta_tree"))
  if (isTRUE(tree$no_tree)) return(NA_integer_)
  denoms <- cta_endpoint_denominators(tree)
  if (length(denoms) == 0L) return(NA_integer_)
  min(denoms)
}

#' D statistic for a fitted CTA tree
#'
#' Computes the parsimony-normalized classification criterion:
#'
#' \deqn{D = \frac{100}{\text{ESS} / \text{strata}} - \text{strata}}
#'
#' where \code{strata} is the number of terminal leaf endpoints and
#' \code{ESS} is \code{tree$overall_ess} (WESS when case weights are active,
#' ESS otherwise).
#'
#' Returns \code{NA_real_} when:
#' \itemize{
#'   \item \code{tree$no_tree} is \code{TRUE};
#'   \item \code{tree$overall_ess} is missing, non-finite, or \eqn{\le 0};
#'   \item \code{strata < 2}.
#' }
#'
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @return Numeric scalar D, or \code{NA_real_}.
#' @seealso \code{\link{cta_strata}}, \code{\link{cta_min_terminal_denom}}
#' @export
cta_d_stat <- function(tree) {
  stopifnot(inherits(tree, "cta_tree"))
  if (isTRUE(tree$no_tree)) return(NA_real_)
  ess <- tree$overall_ess
  if (is.null(ess) || length(ess) != 1L || !is.finite(ess) || ess <= 0)
    return(NA_real_)
  s <- cta_strata(tree)
  if (is.na(s) || s < 2L) return(NA_real_)
  100 / (ess / s) - s
}

# ---- Internal MDSA family skeleton ------------------------------------------
# Internal constructors for the future cta_descendant_family() loop (Phase 2E).

#' Construct a single MDSA descendant family member (internal)
#'
#' Collects tree-level metadata for one fitted CTA tree in an MDSA chain.
#' \code{d} and \code{overall_ess} are read directly from the supplied
#' \code{tree} object; both will be \code{NA} for no-tree fits.
#'
#' @param mindenom Integer MINDENOM used to fit \code{tree}.
#' @param tree A \code{cta_tree} from \code{\link{oda_cta_fit}}.
#' @return A named list with fields: \code{mindenom}, \code{no_tree},
#'   \code{tree}, \code{strata}, \code{endpoint_denominators},
#'   \code{min_terminal_denom}, \code{has_weights}, \code{overall_ess},
#'   \code{d}, \code{next_mindenom}.
#' @keywords internal
new_cta_family_member <- function(mindenom, tree) {
  stopifnot(inherits(tree, "cta_tree"))
  mindenom  <- as.integer(mindenom)
  no_tree   <- isTRUE(tree$no_tree)
  strata    <- cta_strata(tree)
  ep_denoms <- cta_endpoint_denominators(tree)
  min_td    <- cta_min_terminal_denom(tree)
  next_md   <- if (is.na(min_td)) NA_integer_ else min_td + 1L
  list(
    mindenom              = mindenom,
    no_tree               = no_tree,
    tree                  = tree,
    strata                = strata,
    endpoint_denominators = ep_denoms,
    min_terminal_denom    = min_td,
    has_weights           = isTRUE(tree$has_weights),
    overall_ess           = tree$overall_ess %||% NA_real_,
    d                     = cta_d_stat(tree),
    next_mindenom         = next_md
  )
}

#' Construct a CTA descendant family object (internal)
#'
#' Skeleton container for MDSA descendant family members built by
#' \code{\link{new_cta_family_member}}.  The full search loop
#' (\code{cta_descendant_family}) is deferred to Phase 2E.
#'
#' @param members A list of objects from \code{\link{new_cta_family_member}}.
#' @return A list of class \code{cta_family}.
#' @keywords internal
new_cta_family <- function(members = list()) {
  stopifnot(is.list(members))
  structure(list(members = members), class = "cta_family")
}
