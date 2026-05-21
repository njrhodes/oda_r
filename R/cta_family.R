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
#' Container for MDSA descendant family members built by
#' \code{\link{new_cta_family_member}}.  Used internally by
#' \code{\link{cta_descendant_family}}.
#'
#' @param members A list of objects from \code{\link{new_cta_family_member}}.
#' @param mindenoms Integer vector of MINDENOM values, same length as
#'   \code{members}.
#' @param summary Data frame summary with one row per member.
#' @param min_d_idx Integer index of the feasible member with minimum D;
#'   \code{NA_integer_} if none.
#' @param terminated Logical; \code{TRUE} when the loop has exited.
#' @param termination_reason Character: one of \code{"no_tree"},
#'   \code{"max_steps"}, \code{"no_next_mindenom"}.
#' @return A list of class \code{cta_family}.
#' @keywords internal
new_cta_family <- function(members       = list(),
                           mindenoms     = integer(0L),
                           summary       = data.frame(),
                           min_d_idx     = NA_integer_,
                           terminated    = FALSE,
                           termination_reason = NA_character_) {
  stopifnot(is.list(members))
  structure(
    list(
      members            = members,
      mindenoms          = mindenoms,
      summary            = summary,
      min_d_idx          = as.integer(min_d_idx),
      terminated         = terminated,
      termination_reason = termination_reason
    ),
    class = "cta_family"
  )
}

# ---- Public MDSA descendant family function ---------------------------------

#' MDSA descendant family for CTA
#'
#' Traces the MDSA descendant family by fitting CTA models starting at
#' \code{start_mindenom} and stepping according to the novometric MDSA rule:
#' next MINDENOM = minimum terminal endpoint denominator + 1.  The family
#' terminates when a no-tree fit is produced or \code{max_steps} is reached.
#'
#' @param X Data frame of predictor attributes; passed to
#'   \code{\link{oda_cta_fit}}.
#' @param y Integer class vector; passed to \code{\link{oda_cta_fit}}.
#' @param w Optional numeric case-weight vector; passed to
#'   \code{\link{oda_cta_fit}}.
#' @param ... Additional arguments forwarded to \code{\link{oda_cta_fit}}
#'   (e.g. \code{alpha_split}, \code{prune_alpha}, \code{mc_iter},
#'   \code{mc_seed}, \code{loo}, \code{miss_codes}, \code{verbose}).
#' @param start_mindenom Integer MINDENOM for the first family member.
#'   Defaults to \code{1L}.
#' @param max_steps Integer safety cap on the number of CTA fits; prevents
#'   unbounded loops.  Defaults to \code{20L}.
#'
#' @return A list of class \code{cta_family} with fields:
#' \describe{
#'   \item{members}{List of \code{new_cta_family_member} objects in order,
#'     including the terminal no-tree member.}
#'   \item{mindenoms}{Integer vector of MINDENOM values tried.}
#'   \item{summary}{Data frame with one row per member: \code{mindenom},
#'     \code{status} (\code{"valid_tree"}, \code{"stump"}, or
#'     \code{"no_tree"}), \code{strata}, \code{min_terminal_denom},
#'     \code{overall_ess}, \code{d}, \code{no_tree}.}
#'   \item{min_d_idx}{Integer index of the feasible (non-no-tree) member with
#'     minimum D; \code{NA_integer_} if no feasible member exists.}
#'   \item{terminated}{Logical; always \code{TRUE}.}
#'   \item{termination_reason}{Character: one of \code{"no_tree"},
#'     \code{"max_steps"}, \code{"no_next_mindenom"}.}
#' }
#' @seealso \code{\link{oda_cta_fit}}, \code{\link{cta_d_stat}},
#'   \code{\link{cta_min_terminal_denom}}, \code{\link{cta_strata}}
#' @export
cta_descendant_family <- function(
    X,
    y,
    w              = NULL,
    ...,
    start_mindenom = 1L,
    max_steps      = 20L
) {
  start_mindenom <- as.integer(start_mindenom)
  max_steps      <- as.integer(max_steps)
  stopifnot(
    is.data.frame(X),
    length(y) == nrow(X),
    !is.na(start_mindenom) && start_mindenom >= 1L,
    !is.na(max_steps)      && max_steps      >= 1L
  )

  members     <- vector("list", 0L)
  current_md  <- start_mindenom
  term_reason <- NA_character_

  for (step in seq_len(max_steps)) {
    tree   <- oda_cta_fit(X, y, w = w, mindenom = current_md, ...)
    member <- new_cta_family_member(current_md, tree)
    members <- c(members, list(member))

    if (isTRUE(tree$no_tree)) {
      term_reason <- "no_tree"
      break
    }

    next_md <- member$next_mindenom    # cta_min_terminal_denom(tree) + 1L
    if (is.na(next_md) || next_md <= current_md) {
      term_reason <- "no_next_mindenom"
      break
    }

    current_md <- next_md
  }

  if (is.na(term_reason)) term_reason <- "max_steps"

  # ---- Assemble return fields ------------------------------------------------
  mindenoms <- vapply(members, `[[`, integer(1L),  "mindenom")
  d_vals    <- vapply(members, `[[`, double(1L),   "d")
  no_trees  <- vapply(members, `[[`, logical(1L),  "no_tree")
  strata_v  <- vapply(members, function(m) {
    s <- m[["strata"]]; if (is.null(s)) NA_integer_ else as.integer(s)
  }, integer(1L))
  mintd_v   <- vapply(members, function(m) {
    v <- m[["min_terminal_denom"]]; if (is.null(v)) NA_integer_ else as.integer(v)
  }, integer(1L))
  ess_v     <- vapply(members, function(m) {
    v <- m[["overall_ess"]]; if (is.null(v)) NA_real_ else as.double(v)
  }, double(1L))

  status_v <- ifelse(
    no_trees,
    "no_tree",
    ifelse(!is.na(strata_v) & strata_v == 2L, "stump", "valid_tree")
  )

  summary_df <- data.frame(
    mindenom           = mindenoms,
    status             = status_v,
    strata             = strata_v,
    min_terminal_denom = mintd_v,
    overall_ess        = ess_v,
    d                  = d_vals,
    no_tree            = no_trees,
    stringsAsFactors   = FALSE
  )

  feasible  <- which(!no_trees)
  min_d_idx <- if (length(feasible) == 0L) {
    NA_integer_
  } else {
    as.integer(feasible[which.min(d_vals[feasible])])
  }

  new_cta_family(
    members            = members,
    mindenoms          = mindenoms,
    summary            = summary_df,
    min_d_idx          = min_d_idx,
    terminated         = TRUE,
    termination_reason = term_reason
  )
}

# ---- Public CTA family reporting surface ------------------------------------

#' Tidy table of a CTA descendant family
#'
#' Returns a \code{data.frame} with one row per family member.  All values are
#' read from the stored \code{cta_family} object — no refitting or recomputation
#' is performed.
#'
#' @param family A \code{cta_family} from \code{\link{cta_descendant_family}}.
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{\code{index}}{Integer position of the member in the chain.}
#'   \item{\code{mindenom}}{Integer MINDENOM used for this fit.}
#'   \item{\code{status}}{Character: \code{"valid_tree"}, \code{"stump"}, or
#'     \code{"no_tree"}.}
#'   \item{\code{no_tree}}{Logical; \code{TRUE} for the terminal no-tree member.}
#'   \item{\code{strata}}{Integer number of terminal leaf endpoints;
#'     \code{NA} for no-tree members.}
#'   \item{\code{min_terminal_denom}}{Integer minimum leaf \code{n_obs};
#'     \code{NA} for no-tree members.}
#'   \item{\code{next_mindenom}}{Integer MINDENOM for the next chain step
#'     (\code{min_terminal_denom + 1}); \code{NA} for no-tree members.}
#'   \item{\code{overall_ess}}{Numeric overall ESS or WESS; \code{NA} for
#'     no-tree members.}
#'   \item{\code{has_weights}}{Logical; \code{TRUE} when case weights were
#'     active for this fit.}
#'   \item{\code{d}}{Numeric D statistic; \code{NA} for no-tree members.}
#'   \item{\code{selected_min_d}}{Logical; \code{TRUE} for the feasible member
#'     with minimum D (i.e., at index \code{family$min_d_idx}).  All
#'     \code{FALSE} when no feasible member exists.}
#' }
#' @seealso \code{\link{cta_descendant_family}}, \code{\link{summary.cta_family}}
#' @export
cta_family_table <- function(family) {
  stopifnot(inherits(family, "cta_family"))
  n <- length(family$members)
  if (n == 0L) {
    return(data.frame(
      index              = integer(0),
      mindenom           = integer(0),
      status             = character(0),
      no_tree            = logical(0),
      strata             = integer(0),
      min_terminal_denom = integer(0),
      next_mindenom      = integer(0),
      overall_ess        = double(0),
      has_weights        = logical(0),
      d                  = double(0),
      selected_min_d     = logical(0),
      stringsAsFactors   = FALSE
    ))
  }
  s <- family$summary
  next_md_v <- vapply(family$members, function(m) {
    v <- m[["next_mindenom"]]
    if (is.null(v)) NA_integer_ else as.integer(v)
  }, integer(1L))
  has_wt_v <- vapply(family$members, function(m) {
    isTRUE(m[["has_weights"]])
  }, logical(1L))
  sel <- rep(FALSE, n)
  mid <- family$min_d_idx
  if (!is.na(mid) && mid >= 1L && mid <= n) sel[mid] <- TRUE
  data.frame(
    index              = seq_len(n),
    mindenom           = s$mindenom,
    status             = s$status,
    no_tree            = s$no_tree,
    strata             = s$strata,
    min_terminal_denom = s$min_terminal_denom,
    next_mindenom      = next_md_v,
    overall_ess        = s$overall_ess,
    has_weights        = has_wt_v,
    d                  = s$d,
    selected_min_d     = sel,
    stringsAsFactors   = FALSE
  )
}

#' Summary of a CTA descendant family
#'
#' Returns a structured S3 summary object for a \code{cta_family}.  All values
#' are read from stored fields — no refitting or recomputation is performed.
#'
#' @param object A \code{cta_family} from \code{\link{cta_descendant_family}}.
#' @param ... Unused; included for S3 compatibility.
#' @return A list of class \code{c("cta_family_summary", "list")} with fields:
#' \describe{
#'   \item{\code{n_members}}{Integer number of family members.}
#'   \item{\code{min_d_idx}}{Integer index of the feasible member with minimum D;
#'     \code{NA_integer_} if none.}
#'   \item{\code{terminated}}{Logical; always \code{TRUE} for a completed chain.}
#'   \item{\code{termination_reason}}{Character: one of \code{"no_tree"},
#'     \code{"max_steps"}, \code{"no_next_mindenom"}.}
#'   \item{\code{has_weights}}{Logical; \code{TRUE} when any family member used
#'     case weights.}
#'   \item{\code{table}}{A \code{data.frame} from \code{\link{cta_family_table}}.}
#' }
#' @seealso \code{\link{cta_descendant_family}}, \code{\link{cta_family_table}}
#' @export
summary.cta_family <- function(object, ...) {
  stopifnot(inherits(object, "cta_family"))
  has_wt <- any(vapply(object$members, function(m) isTRUE(m[["has_weights"]]),
                       logical(1L)))
  structure(
    list(
      n_members          = length(object$members),
      min_d_idx          = object$min_d_idx,
      terminated         = object$terminated,
      termination_reason = object$termination_reason,
      has_weights        = has_wt,
      table              = cta_family_table(object)
    ),
    class = c("cta_family_summary", "list")
  )
}

#' Print a CTA family summary
#'
#' Compact display of the \code{cta_family_summary} object returned by
#' \code{\link{summary.cta_family}}.
#'
#' @param x A \code{cta_family_summary} object.
#' @param ... Unused; included for S3 compatibility.
#' @return \code{invisible(x)}.
#' @seealso \code{\link{summary.cta_family}}, \code{\link{cta_family_table}}
#' @export
print.cta_family_summary <- function(x, ...) {
  ess_label <- if (isTRUE(x$has_weights)) "WESS" else "ESS"
  cat("CTA descendant family\n")
  cat(sprintf("  members: %d  |  termination: %s  |  min-D index: %s\n",
              x$n_members,
              x$termination_reason,
              if (is.na(x$min_d_idx)) "none" else as.character(x$min_d_idx)))
  tbl <- x$table
  names(tbl)[names(tbl) == "mindenom"]    <- "MINDENOM"
  names(tbl)[names(tbl) == "overall_ess"] <- ess_label
  names(tbl)[names(tbl) == "d"]           <- "D"
  print(tbl, row.names = FALSE)
  invisible(x)
}

#' Print a CTA descendant family
#'
#' Calls \code{\link{summary.cta_family}} and prints the result.
#'
#' @param x A \code{cta_family} from \code{\link{cta_descendant_family}}.
#' @param ... Passed to \code{\link{print.cta_family_summary}}.
#' @return \code{invisible(x)}.
#' @seealso \code{\link{summary.cta_family}}, \code{\link{cta_family_table}}
#' @export
print.cta_family <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}
