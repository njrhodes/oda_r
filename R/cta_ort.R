###############################################################################
# R/cta_ort.R -- Optimal Recursive Tree (ORT) engine and S3 methods
#
# Entry point: cta_fit(X, y, ..., recursive = TRUE) in R/oda_fit.R
# Internal engine: .cta_ort_fit_internal()
#
# Classes:
#   cta_ort -- composite recursive tree; also inherits cta_tree for backward compat
#
# Lean-fit invariant extends to ORT:
#   - No training X/y stored at any ORT node
#   - Each ort_node$model is a lean cta_tree
#   - $strata table is computed at fit time and stored
#   - No row indices, no per-stratum X/y, no cached membership vectors
###############################################################################

# ---- Internal helpers ------------------------------------------------------- #

# Weighted majority class from integer y and optional numeric w.
.ort_majority_class <- function(y, w) {
  cls <- sort(unique(y[!is.na(y)]))
  if (length(cls) == 0L) return(NA_integer_)
  if (is.null(w)) {
    ct <- vapply(cls, function(c) sum(y == c, na.rm = TRUE), integer(1L))
  } else {
    ct <- vapply(cls, function(c) sum(w[!is.na(y) & y == c], na.rm = TRUE), double(1L))
  }
  as.integer(cls[which.max(ct)])
}

# Named integer vector of raw class counts (by class label).
.ort_class_counts <- function(y) {
  cls <- sort(unique(y[!is.na(y)]))
  if (length(cls) == 0L) return(integer(0))
  ct  <- vapply(cls, function(c) sum(y == c, na.rm = TRUE), integer(1L))
  setNames(ct, as.character(cls))
}

# ---- .cta_ort_fit_internal -------------------------------------------------- #

# Recursive ORT engine.
#
# counter_env fields (mutable):
#   $next_node_id  -- next available ORT node_id (starts at 1L)
#   $total_nodes   -- total nodes allocated so far (for max_nodes guard)
#   $all_nodes     -- list of completed ort_node lists, keyed by character node_id
#
# Returns: node_id (integer) of the node just created.
.cta_ort_fit_internal <- function(
  X,                 # data.frame, training rows at this ORT node
  y,                 # integer class vector (same length as nrow(X))
  w,                 # numeric weight vector or NULL
  path_conditions,   # character vector of conditions from root to here
  path_ess,          # numeric vector of ESS values from root (empty at root)
  path_mindenom,     # integer vector of MINDENOM values from root (empty at root)
  depth,             # integer recursion depth (root = 0L)
  counter_env,       # mutable environment for node_id tracking
  mc_iter, alpha_split, prune_alpha, loo,
  min_n, max_depth, max_nodes, verbose
  # mc_seed intentionally absent: seed is set once in .cta_ort_fit(); child
  # fits consume the RNG stream in deterministic R->L traversal order.
) {
  # Allocate this node's ID (DFS pre-order)
  my_node_id <- counter_env$next_node_id
  counter_env$next_node_id <- counter_env$next_node_id + 1L
  counter_env$total_nodes  <- counter_env$total_nodes  + 1L

  n     <- nrow(X)
  y_int <- as.integer(y)

  path_str <- if (length(path_conditions) == 0L) "/"
              else paste(path_conditions, collapse = " AND ")

  # ---- Guard checks ----------------------------------------------------------
  stop_reason <- NULL

  if (n < min_n) {
    stop_reason <- "min_n"
  } else if (depth >= max_depth) {
    stop_reason <- "max_depth"
  } else if (counter_env$total_nodes > max_nodes) {
    stop_reason <- "max_nodes"
  } else if (length(unique(y_int)) < 2L) {
    # Pure node: only one class present -- no discrimination possible
    stop_reason <- "no_tree"
  }

  if (!is.null(stop_reason)) {
    nd <- list(
      node_id         = my_node_id,
      depth           = depth,
      n               = n,
      class_counts    = .ort_class_counts(y_int),
      path_conditions = path_conditions,
      path            = path_str,
      path_ess        = path_ess,
      path_mindenom   = path_mindenom,
      is_terminal     = TRUE,
      stop_reason     = stop_reason,
      model           = NULL,
      level_mindenom  = NA_integer_,
      level_ess       = NA_real_,
      level_d         = NA_real_,
      terminal_class  = .ort_majority_class(y_int, w),
      child_ids       = integer(0),
      right_child_id  = NA_integer_,
      left_child_id   = NA_integer_
    )
    counter_env$all_nodes[[as.character(my_node_id)]] <- nd
    return(my_node_id)
  }

  # ---- MDSA family scan ------------------------------------------------------
  if (verbose) {
    message(sprintf(
      "[ORT] depth=%d node=%d n=%d -- running MDSA family scan",
      depth, my_node_id, n
    ))
  }

  fam <- cta_descendant_family(
    X, y_int, w = w,
    mc_iter     = mc_iter,
    alpha_split = alpha_split,
    prune_alpha = prune_alpha,
    loo         = loo
    # mc_seed not passed: RNG stream flows from single top-level seed
  )

  min_d_idx <- fam$min_d_idx
  # Guard: length 0 or NA both mean no valid min-D member
  if (length(min_d_idx) != 1L) min_d_idx <- NA_integer_

  if (is.na(min_d_idx)) {
    # All family members are no-tree -> terminal
    last_tree <- fam$members[[length(fam$members)]]$tree
    nd <- list(
      node_id         = my_node_id,
      depth           = depth,
      n               = n,
      class_counts    = .ort_class_counts(y_int),
      path_conditions = path_conditions,
      path            = path_str,
      path_ess        = path_ess,
      path_mindenom   = path_mindenom,
      is_terminal     = TRUE,
      stop_reason     = "no_tree",
      model           = last_tree,
      level_mindenom  = NA_integer_,
      level_ess       = NA_real_,
      level_d         = NA_real_,
      terminal_class  = .ort_majority_class(y_int, w),
      child_ids       = integer(0),
      right_child_id  = NA_integer_,
      left_child_id   = NA_integer_
    )
    counter_env$all_nodes[[as.character(my_node_id)]] <- nd
    return(my_node_id)
  }

  # ---- Extract min-D member --------------------------------------------------
  best_member  <- fam$members[[min_d_idx]]
  best_tree    <- best_member$tree
  best_md      <- best_member$mindenom
  best_ess     <- best_member$overall_ess %||% NA_real_
  best_d       <- best_member$d           %||% NA_real_

  new_path_ess      <- c(path_ess,      best_ess)
  new_path_mindenom <- c(path_mindenom, best_md)

  if (verbose) {
    message(sprintf(
      "[ORT] depth=%d node=%d -- MINDENOM=%d, ESS=%.2f%%, D=%.4f",
      depth, my_node_id, best_md,
      if (is.na(best_ess)) 0 else best_ess,
      if (is.na(best_d))   0 else best_d
    ))
  }

  # ---- Get endpoint assignments for X subset ---------------------------------
  ep_summary <- cta_endpoint_summary(best_tree)
  n_ep       <- nrow(ep_summary)
  ep_df      <- cta_assign_endpoints(best_tree, X, missing_action = "na")

  # ---- Recurse into each endpoint -- right first (reverse endpoint_id order) --
  # ep_summary is ordered by ascending node_id (ep_id 1 = lowest node = left)
  # Right-first: process ep_id = n_ep, n_ep-1, ..., 1
  child_ids <- integer(n_ep)
  ep_order  <- rev(seq_len(n_ep))

  for (ep_idx in ep_order) {
    ep_node_id  <- ep_summary$endpoint_node_id[ep_idx]

    # Condition string(s) for this endpoint's path within the sub-tree
    ep_segs <- .cta_path_segments(best_tree, ep_node_id)
    ep_cond <- if (length(ep_segs) == 0L) sprintf("[ep%d]", ep_idx)
               else paste(ep_segs, collapse = " AND ")

    # Row subset belonging to this endpoint
    in_ep  <- which(!is.na(ep_df$endpoint_id) & ep_df$endpoint_id == ep_idx)
    X_sub  <- X[in_ep, , drop = FALSE]
    y_sub  <- y_int[in_ep]
    w_sub  <- if (!is.null(w)) w[in_ep] else NULL

    child_id <- .cta_ort_fit_internal(
      X               = X_sub,
      y               = y_sub,
      w               = w_sub,
      path_conditions = c(path_conditions, ep_cond),
      path_ess        = new_path_ess,
      path_mindenom   = new_path_mindenom,
      depth           = depth + 1L,
      counter_env     = counter_env,
      mc_iter         = mc_iter,
      alpha_split     = alpha_split,
      prune_alpha     = prune_alpha,
      loo             = loo,
      min_n           = min_n,
      max_depth       = max_depth,
      max_nodes       = max_nodes,
      verbose         = verbose
    )

    child_ids[ep_idx] <- child_id
  }

  # ---- Assemble and store this node ------------------------------------------
  nd <- list(
    node_id         = my_node_id,
    depth           = depth,
    n               = n,
    class_counts    = .ort_class_counts(y_int),
    path_conditions = path_conditions,
    path            = path_str,
    path_ess        = path_ess,
    path_mindenom   = path_mindenom,
    is_terminal     = FALSE,
    stop_reason     = NA_character_,
    model           = best_tree,
    level_mindenom  = best_md,
    level_ess       = best_ess,
    level_d         = best_d,
    terminal_class  = NA_integer_,
    child_ids       = child_ids,
    right_child_id  = child_ids[n_ep],    # last ep (highest node_id = "right")
    left_child_id   = child_ids[1L]       # first ep (lowest node_id = "left")
  )
  counter_env$all_nodes[[as.character(my_node_id)]] <- nd
  return(my_node_id)
}

# ---- .cta_ort_build_strata -------------------------------------------------- #

.ort_empty_strata <- function() {
  df <- data.frame(
    stratum_id     = integer(0),
    node_id        = integer(0),
    path           = character(0),
    depth          = integer(0),
    n              = integer(0),
    prop_class1    = double(0),
    odds_class1    = double(0),
    terminal_class = integer(0),
    stop_reason    = character(0),
    path_ess       = character(0),
    path_mindenom  = character(0),
    stringsAsFactors = FALSE
  )
  df$class_counts <- list()
  df
}

# Build the flat strata table from completed ort_nodes.
# Terminal nodes sorted ascending by prop_class1 (Stage 1 = lowest risk).
.cta_ort_build_strata <- function(ort_nodes) {
  term_keys <- Filter(function(k) isTRUE(ort_nodes[[k]]$is_terminal),
                      names(ort_nodes))
  if (length(term_keys) == 0L) return(.ort_empty_strata())

  # Build row data for each terminal node
  node_id_v       <- integer(length(term_keys))
  path_v          <- character(length(term_keys))
  depth_v         <- integer(length(term_keys))
  n_v             <- integer(length(term_keys))
  prop_v          <- double(length(term_keys))
  odds_v          <- double(length(term_keys))
  term_cls_v      <- integer(length(term_keys))
  stop_v          <- character(length(term_keys))
  path_ess_v      <- character(length(term_keys))
  path_md_v       <- character(length(term_keys))
  class_counts_l  <- vector("list", length(term_keys))

  for (i in seq_along(term_keys)) {
    nd  <- ort_nodes[[term_keys[i]]]
    cc  <- nd$class_counts        # named integer vector, e.g. c("0"=18, "1"=2)
    n   <- nd$n

    cls_names <- names(cc)
    idx1 <- which(cls_names == "1")
    if (length(idx1) == 0L) idx1 <- length(cls_names)  # fallback: last class
    n1   <- if (length(idx1) > 0L) as.integer(cc[idx1[1L]]) else 0L

    prop <- if (n > 0L && !is.na(n1)) n1 / n else NA_real_
    odds <- if (!is.na(prop) && is.finite(prop) && prop > 0 && prop < 1)
              prop / (1 - prop)
            else if (!is.na(prop) && prop == 1) Inf
            else if (!is.na(prop) && prop == 0) 0
            else NA_real_

    ess_vec <- nd$path_ess
    md_vec  <- nd$path_mindenom
    path_ess_str <- if (length(ess_vec) == 0L) NA_character_
                    else paste(round(ess_vec, 2L), collapse = " -> ")
    path_md_str  <- if (length(md_vec) == 0L) NA_character_
                    else paste(md_vec, collapse = " -> ")

    node_id_v[i]      <- nd$node_id
    path_v[i]         <- nd$path
    depth_v[i]        <- nd$depth
    n_v[i]            <- n
    prop_v[i]         <- prop
    odds_v[i]         <- odds
    term_cls_v[i]     <- nd$terminal_class %||% NA_integer_
    stop_v[i]         <- nd$stop_reason    %||% NA_character_
    path_ess_v[i]     <- path_ess_str
    path_md_v[i]      <- path_md_str
    class_counts_l[[i]] <- cc
  }

  df <- data.frame(
    node_id        = node_id_v,
    path           = path_v,
    depth          = depth_v,
    n              = n_v,
    prop_class1    = prop_v,
    odds_class1    = odds_v,
    terminal_class = term_cls_v,
    stop_reason    = stop_v,
    path_ess       = path_ess_v,
    path_mindenom  = path_md_v,
    stringsAsFactors = FALSE
  )
  df$class_counts <- class_counts_l

  # Sort ascending by prop_class1, then node_id for ties
  ord <- order(prop_v, node_id_v, na.last = TRUE)
  df  <- df[ord, ]
  df$stratum_id <- seq_len(nrow(df))
  df  <- df[, c("stratum_id", setdiff(names(df), "stratum_id"))]
  rownames(df) <- NULL
  df
}

# ---- .cta_ort_fit (entry point called from cta_fit) ------------------------- #

# Called from cta_fit() when recursive = TRUE.
# Returns a cta_ort / cta_tree dual-tagged object.
.cta_ort_fit <- function(
  X,
  y,
  w           = NULL,
  mc_seed     = 42L,
  mc_iter     = 5000L,
  alpha_split = 0.05,
  prune_alpha = 0.05,
  loo         = "stable",
  min_n       = 30L,
  max_depth   = 8L,
  max_nodes   = 31L,
  verbose     = FALSE
) {
  X         <- as.data.frame(X)
  y         <- as.integer(y)
  n         <- nrow(X)
  min_n     <- as.integer(min_n)
  max_depth <- as.integer(max_depth)
  max_nodes <- as.integer(max_nodes)

  # Set seed once for the entire recursive run.
  # Child fits consume the RNG stream in deterministic right-then-left
  # traversal order; seeds are not reset per node.
  if (!is.null(mc_seed)) set.seed(as.integer(mc_seed))

  # Mutable state for the recursion
  counter_env <- new.env(parent = emptyenv())
  counter_env$next_node_id <- 1L
  counter_env$total_nodes  <- 0L
  counter_env$all_nodes    <- list()

  # Run the recursive engine (mc_seed not passed -- RNG stream flows through)
  .cta_ort_fit_internal(
    X               = X,
    y               = y,
    w               = w,
    path_conditions = character(0),
    path_ess        = numeric(0),
    path_mindenom   = integer(0),
    depth           = 0L,
    counter_env     = counter_env,
    mc_iter         = mc_iter,
    alpha_split     = alpha_split,
    prune_alpha     = prune_alpha,
    loo             = loo,
    min_n           = min_n,
    max_depth       = max_depth,
    max_nodes       = max_nodes,
    verbose         = verbose
  )

  ort_nodes  <- counter_env$all_nodes
  root_node  <- ort_nodes[["1"]]
  strata     <- .cta_ort_build_strata(ort_nodes)

  # Root-level cta_tree model for backward compatibility.
  # Non-terminal root: model is the MDSA min-D tree.
  # Terminal root (guard fired): synthesize a no-tree cta_tree.
  root_model <- root_node$model
  if (is.null(root_model)) {
    # Forced terminal at root (min_n / max_depth / max_nodes guard before any fit)
    root_model <- oda_cta_fit(
      X, y, w = w,
      mindenom = n + 1L,
      mc_iter  = 1L
      # mc_seed not passed: no MC runs at this mindenom, seed irrelevant
    )
  }

  # Build the dual-tagged ORT object (extends the root cta_tree)
  obj <- structure(
    modifyList(root_model, list(
      recursive    = TRUE,
      ort_nodes    = ort_nodes,
      ort_root_id  = 1L,
      strata       = strata,
      n_strata     = nrow(strata),
      ort_settings = list(
        mc_seed     = mc_seed,
        mc_iter     = mc_iter,
        alpha_split = alpha_split,
        prune_alpha = prune_alpha,
        loo         = loo,
        min_n       = min_n,
        max_depth   = max_depth,
        max_nodes   = max_nodes
      )
    )),
    class = c("cta_ort", "cta_tree")
  )

  # Strata consistency check: route training X through the composite tree
  # and verify per-stratum row counts match stored n values.
  if (!is.null(strata) && nrow(strata) > 0L) {
    tryCatch({
      train_pred <- predict(obj, X, type = "all", missing_action = "na")
      predicted_n <- vapply(strata$stratum_id, function(sid) {
        sum(!is.na(train_pred$stratum_id) & train_pred$stratum_id == sid,
            na.rm = TRUE)
      }, integer(1L))
      obj$strata_check_passed <- isTRUE(all(predicted_n == strata$n))
    }, error = function(e) {
      obj$strata_check_passed <<- FALSE
    })
  } else {
    obj$strata_check_passed <- (nrow(strata) == 0L)
  }

  obj
}

# ---- predict.cta_ort -------------------------------------------------------- #

#' Predict method for Optimal Recursive Tree
#'
#' Routes each row of \code{newdata} down the composite ORT by recursively
#' applying each node's \code{cta_tree} model via
#' \code{\link{cta_assign_endpoints}}.
#'
#' @param object A \code{cta_ort} from \code{cta_fit(..., recursive = TRUE)}.
#' @param newdata Data frame or matrix matching the training X column layout.
#' @param type Character; one of \code{"class"} (default), \code{"stratum"},
#'   \code{"path"}, or \code{"all"}.
#' @param missing_action Passed to each node-level
#'   \code{\link{cta_assign_endpoints}} call.  \code{"na"} (default):
#'   observations with a missing split attribute return \code{NA}.
#'   \code{"majority"}: route to the majority-class child.
#' @param ... Unused.
#' @return For \code{type = "class"}: integer vector of predicted class labels
#'   (length \code{nrow(newdata)}).  For \code{type = "stratum"}: integer
#'   stratum_id vector.  For \code{type = "path"}: character path vector.
#'   For \code{type = "all"}: data.frame with columns
#'   \code{predicted_class}, \code{stratum_id}, \code{path},
#'   \code{prop_class1}, \code{stop_reason}.
#' @export
predict.cta_ort <- function(object, newdata,
                             type           = c("class", "stratum", "path", "all"),
                             missing_action = c("na", "majority"),
                             ...) {
  type           <- match.arg(type)
  missing_action <- match.arg(missing_action)
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)

  n_new     <- nrow(newdata)
  ort_nodes <- object$ort_nodes %||% list()
  strata    <- object$strata

  # Build lookup maps from terminal ORT node_id to strata fields
  if (!is.null(strata) && nrow(strata) > 0L) {
    nid_chr    <- as.character(strata$node_id)
    strat_sid  <- setNames(strata$stratum_id,    nid_chr)
    strat_prop <- setNames(strata$prop_class1,   nid_chr)
    strat_stop <- setNames(strata$stop_reason,   nid_chr)
    strat_path <- setNames(strata$path,          nid_chr)
    strat_cls  <- setNames(strata$terminal_class, nid_chr)
  } else {
    strat_sid  <- character(0)
    strat_prop <- numeric(0)
    strat_stop <- character(0)
    strat_path <- character(0)
    strat_cls  <- integer(0)
  }

  # Route one observation (row i) down the ORT.
  route_one <- function(i) {
    nd_id <- 1L
    repeat {
      nd_key <- as.character(nd_id)
      nd     <- ort_nodes[[nd_key]]
      if (is.null(nd) || isTRUE(nd$is_terminal)) break

      model <- nd$model
      if (is.null(model) || isTRUE(model$no_tree)) break

      row_df <- newdata[i, , drop = FALSE]
      ep_df  <- cta_assign_endpoints(model, row_df,
                                     missing_action = missing_action)
      ep_id  <- ep_df$endpoint_id[1L]

      if (is.na(ep_id)) {
        return(list(cls = NA_integer_, sid = NA_integer_,
                    path = NA_character_, prop = NA_real_,
                    stop = NA_character_))
      }

      cids <- nd$child_ids
      if (length(cids) == 0L || ep_id < 1L || ep_id > length(cids)) break
      nd_id <- cids[ep_id]
      if (is.na(nd_id) || nd_id < 1L) break
    }

    nd_key_t <- as.character(nd_id)
    list(
      cls  = as.integer(strat_cls[nd_key_t]  %||% NA_integer_),
      sid  = as.integer(strat_sid[nd_key_t]  %||% NA_integer_),
      path = as.character(strat_path[nd_key_t] %||% NA_character_),
      prop = as.double(strat_prop[nd_key_t]  %||% NA_real_),
      stop = as.character(strat_stop[nd_key_t] %||% NA_character_)
    )
  }

  results <- lapply(seq_len(n_new), route_one)

  if (type == "class") {
    return(vapply(results, function(r) r$cls, integer(1L)))
  }
  if (type == "stratum") {
    return(vapply(results, function(r) r$sid, integer(1L)))
  }
  if (type == "path") {
    return(vapply(results, function(r) r$path, character(1L)))
  }
  # type == "all"
  data.frame(
    predicted_class  = vapply(results, function(r) r$cls,  integer(1L)),
    stratum_id       = vapply(results, function(r) r$sid,  integer(1L)),
    path             = vapply(results, function(r) r$path, character(1L)),
    prop_class1      = vapply(results, function(r) r$prop, double(1L)),
    stop_reason      = vapply(results, function(r) r$stop, character(1L)),
    stringsAsFactors = FALSE
  )
}

# ---- print.cta_ort ---------------------------------------------------------- #

#' Print method for Optimal Recursive Tree
#'
#' @param x A \code{cta_ort} object.
#' @param ... Unused.
#' @return \code{invisible(x)}.
#' @export
print.cta_ort <- function(x, ...) {
  strata <- x$strata
  ns     <- x$n_strata   %||% 0L
  s      <- x$ort_settings %||% list()
  cat("Optimal Recursive Tree (cta_ort)\n")
  cat(sprintf("  Strata: %d terminal strata\n", ns))
  if (!is.null(s$min_n))
    cat(sprintf("  Guards: min_n=%d, max_depth=%d, max_nodes=%d\n",
                s$min_n, s$max_depth, s$max_nodes))
  if (!is.null(s$mc_seed))
    cat(sprintf("  mc_seed=%d, mc_iter=%d\n", s$mc_seed, s$mc_iter))
  chk <- x$strata_check_passed
  if (!is.null(chk))
    cat(sprintf("  Strata consistency check: %s\n",
                if (isTRUE(chk)) "PASSED" else "FAILED"))
  if (!is.null(strata) && nrow(strata) > 0L) {
    cat("\nTerminal strata (ascending class-1 proportion):\n")
    for (i in seq_len(nrow(strata))) {
      r <- strata[i, ]
      cat(sprintf("  Stage %-2d  n=%-6d  prop=%.4f  stop=%-10s  %s\n",
                  r$stratum_id, r$n,
                  r$prop_class1  %||% NA_real_,
                  r$stop_reason  %||% "",
                  r$path         %||% ""))
    }
  }
  invisible(x)
}

# ---- summary.cta_ort -------------------------------------------------------- #

#' Summary method for Optimal Recursive Tree
#'
#' Returns a structured list of class \code{"cta_ort_summary"} capturing
#' tree-level metadata for the composite ORT.
#'
#' @param object A \code{cta_ort} object.
#' @param ... Unused.
#' @return A list of class \code{"cta_ort_summary"}.
#' @export
summary.cta_ort <- function(object, ...) {
  ort_nodes <- object$ort_nodes %||% list()
  n_total   <- length(ort_nodes)
  is_term_v <- vapply(ort_nodes, function(nd) isTRUE(nd$is_terminal), logical(1L))
  n_terminal <- sum(is_term_v)
  n_split    <- n_total - n_terminal
  depths     <- vapply(ort_nodes, function(nd) nd$depth %||% NA_integer_, integer(1L))
  max_d_obs  <- if (length(depths) > 0L) max(depths, na.rm = TRUE) else NA_integer_

  out <- list(
    n_strata            = object$n_strata %||% 0L,
    n_nodes_total       = n_total,
    n_split_nodes       = n_split,
    n_terminal_nodes    = n_terminal,
    max_depth_observed  = max_d_obs,
    strata              = object$strata,
    ort_settings        = object$ort_settings,
    strata_check_passed = object$strata_check_passed
  )
  class(out) <- "cta_ort_summary"
  out
}

#' Print method for cta_ort_summary
#'
#' @param x A \code{cta_ort_summary} object.
#' @param ... Unused.
#' @return \code{invisible(x)}.
#' @export
print.cta_ort_summary <- function(x, ...) {
  cat("ORT Summary\n")
  cat(sprintf("  %d terminal strata, %d total nodes (%d split + %d terminal)\n",
              x$n_strata, x$n_nodes_total, x$n_split_nodes, x$n_terminal_nodes))
  cat(sprintf("  Max depth observed: %d\n", x$max_depth_observed %||% NA_integer_))
  if (!is.null(x$strata_check_passed))
    cat(sprintf("  Strata consistency: %s\n",
                if (isTRUE(x$strata_check_passed)) "PASSED" else "FAILED"))
  if (!is.null(x$ort_settings)) {
    s <- x$ort_settings
    cat(sprintf("  Settings: min_n=%d, max_depth=%d, max_nodes=%d, mc_iter=%d\n",
                s$min_n %||% NA_integer_, s$max_depth %||% NA_integer_,
                s$max_nodes %||% NA_integer_, s$mc_iter %||% NA_integer_))
  }
  if (!is.null(x$strata) && nrow(x$strata) > 0L) {
    cat("\nStrata table:\n")
    st <- x$strata
    cols <- intersect(c("stratum_id","n","prop_class1","stop_reason","path"),
                      names(st))
    print(st[, cols, drop = FALSE], row.names = FALSE)
  }
  invisible(x)
}

# ---- ort_plot_data ---------------------------------------------------------- #

#' Renderer-independent layout data for an ORT composite tree
#'
#' Computes node positions and edge metadata for \code{\link{plot.cta_ort}}.
#' Terminal nodes receive integer x-slot positions (left-to-right in DFS
#' right-first order); internal nodes are centered over their children.
#'
#' @param object A \code{cta_ort} from \code{cta_fit(..., recursive = TRUE)}.
#' @param target_class Integer target class for terminal node annotation, or
#'   \code{NULL} for a structural plot.
#' @param class_labels Optional named character vector of class display names.
#' @param digits Integer decimal places for proportion labels.  Default 1.
#' @return A list with elements:
#' \describe{
#'   \item{\code{nodes}}{data.frame: \code{node_id}, \code{depth}, \code{x},
#'     \code{y}, \code{is_terminal}, \code{label}, \code{n},
#'     \code{stop_reason}.}
#'   \item{\code{edges}}{data.frame: \code{from_id}, \code{to_id},
#'     \code{x0}, \code{y0}, \code{x1}, \code{y1}, \code{label}.}
#'   \item{\code{strata}}{The strata table from the ORT object.}
#' }
#' @export
ort_plot_data <- function(object, target_class = NULL,
                          class_labels = NULL, digits = 1L) {
  stopifnot(inherits(object, "cta_ort"))
  ort_nodes <- object$ort_nodes %||% list()
  strata    <- object$strata

  empty_nodes <- function() {
    data.frame(node_id = integer(0), depth = integer(0),
               x = double(0), y = double(0), is_terminal = logical(0),
               label = character(0), n = integer(0),
               stop_reason = character(0), stringsAsFactors = FALSE)
  }
  empty_edges <- function() {
    data.frame(from_id = integer(0), to_id = integer(0),
               x0 = double(0), y0 = double(0), x1 = double(0), y1 = double(0),
               label = character(0), stringsAsFactors = FALSE)
  }

  if (length(ort_nodes) == 0L) {
    return(list(nodes = empty_nodes(), edges = empty_edges(), strata = strata))
  }

  # Assign x positions via DFS right-first leaf enumeration.
  # Terminal nodes get integer slots 1..n_terminal.
  # Internal nodes center over their children.
  x_pos <- setNames(rep(NA_real_, length(ort_nodes)), names(ort_nodes))
  leaf_counter <- 0L

  assign_x <- function(nid_key) {
    nd <- ort_nodes[[nid_key]]
    if (is.null(nd)) return(NA_real_)
    if (isTRUE(nd$is_terminal)) {
      leaf_counter <<- leaf_counter + 1L
      x_pos[nid_key] <<- as.double(leaf_counter)
      return(as.double(leaf_counter))
    }
    # Recurse right-first (reverse order of child_ids)
    cids     <- nd$child_ids
    child_xs <- vapply(rev(as.character(cids)), assign_x, double(1L))
    cx <- mean(child_xs, na.rm = TRUE)
    x_pos[nid_key] <<- cx
    return(cx)
  }
  assign_x("1")

  # Build strata lookup for terminal node label enrichment
  if (!is.null(strata) && nrow(strata) > 0L) {
    strat_sid  <- setNames(strata$stratum_id,  as.character(strata$node_id))
    strat_prop <- setNames(strata$prop_class1, as.character(strata$node_id))
  } else {
    strat_sid  <- character(0)
    strat_prop <- numeric(0)
  }

  # Build node data frame
  n_nodes     <- length(ort_nodes)
  node_ids_v  <- vapply(ort_nodes, `[[`, integer(1L),  "node_id")
  depths_v    <- vapply(ort_nodes, `[[`, integer(1L),  "depth")
  is_term_v   <- vapply(ort_nodes, function(nd) isTRUE(nd$is_terminal), logical(1L))
  n_v         <- vapply(ort_nodes, `[[`, integer(1L),  "n")
  stop_v      <- vapply(ort_nodes, function(nd) nd$stop_reason %||% NA_character_,
                        character(1L))

  node_df <- data.frame(
    node_id     = node_ids_v,
    depth       = depths_v,
    x           = x_pos[names(ort_nodes)],
    y           = -depths_v,
    is_terminal = is_term_v,
    n           = n_v,
    stop_reason = stop_v,
    stringsAsFactors = FALSE
  )

  # Build node labels
  node_df$label <- vapply(seq_len(nrow(node_df)), function(i) {
    nd_key <- as.character(node_df$node_id[i])
    nd     <- ort_nodes[[nd_key]]

    if (isTRUE(nd$is_terminal)) {
      n_here <- nd$n
      maj    <- nd$terminal_class %||% NA_integer_
      lbl_cls <- if (!is.null(class_labels)) {
        class_labels[as.character(maj)] %||% as.character(maj)
      } else as.character(maj %||% "?")

      if (!is.null(target_class)) {
        prop_val <- strat_prop[nd_key] %||% NA_real_
        sid_val  <- strat_sid[nd_key]  %||% NA_integer_
        if (!is.na(prop_val) && !is.na(sid_val)) {
          sprintf("Stage %d\nn=%d\n%.*f%%", sid_val, n_here,
                  as.integer(digits), 100 * prop_val)
        } else {
          sprintf("n=%d\n%s", n_here, lbl_cls)
        }
      } else {
        sprintf("n=%d\n%s", n_here, lbl_cls)
      }
    } else {
      # Split node: root attribute of sub-tree, ESS, MINDENOM, n
      model     <- nd$model
      root_attr <- if (!is.null(model) && !isTRUE(model$no_tree)) {
        root_nd <- model$nodes[[model$root_id]]
        root_nd$attribute %||% "?"
      } else "?"
      ess_lbl <- if (!is.null(nd$level_ess) && !is.na(nd$level_ess))
                   sprintf("ESS=%.1f%%", nd$level_ess) else ""
      md_lbl  <- if (!is.null(nd$level_mindenom) && !is.na(nd$level_mindenom))
                   sprintf("MD=%d", nd$level_mindenom) else ""
      parts <- c(root_attr, ess_lbl, md_lbl, sprintf("n=%d", nd$n))
      parts <- parts[nchar(parts) > 0L]
      paste(parts, collapse = "\n")
    }
  }, character(1L))

  # Build edges
  edge_rows <- list()
  for (nd_key in names(ort_nodes)) {
    nd <- ort_nodes[[nd_key]]
    if (isTRUE(nd$is_terminal) || length(nd$child_ids) == 0L) next

    i_par <- which(node_df$node_id == nd$node_id)[1L]
    if (is.na(i_par)) next
    px <- node_df$x[i_par]
    py <- node_df$y[i_par]

    for (ep_idx in seq_along(nd$child_ids)) {
      cid    <- nd$child_ids[ep_idx]
      if (is.na(cid) || cid < 1L) next
      c_nd   <- ort_nodes[[as.character(cid)]]
      if (is.null(c_nd)) next

      i_ch   <- which(node_df$node_id == cid)[1L]
      if (is.na(i_ch)) next
      cx <- node_df$x[i_ch]
      cy <- node_df$y[i_ch]

      # Edge label = last element of child's path_conditions
      child_conds <- c_nd$path_conditions
      edge_lbl    <- if (length(child_conds) > 0L)
                       child_conds[length(child_conds)]
                     else ""

      edge_rows[[length(edge_rows) + 1L]] <- list(
        from_id = nd$node_id,
        to_id   = cid,
        x0 = px, y0 = py,
        x1 = cx, y1 = cy,
        label = edge_lbl
      )
    }
  }

  if (length(edge_rows) > 0L) {
    edge_df <- data.frame(
      from_id = vapply(edge_rows, `[[`, integer(1L),   "from_id"),
      to_id   = vapply(edge_rows, `[[`, integer(1L),   "to_id"),
      x0      = vapply(edge_rows, `[[`, double(1L),    "x0"),
      y0      = vapply(edge_rows, `[[`, double(1L),    "y0"),
      x1      = vapply(edge_rows, `[[`, double(1L),    "x1"),
      y1      = vapply(edge_rows, `[[`, double(1L),    "y1"),
      label   = vapply(edge_rows, `[[`, character(1L), "label"),
      stringsAsFactors = FALSE
    )
  } else {
    edge_df <- empty_edges()
  }

  list(nodes = node_df, edges = edge_df, strata = strata)
}

# ---- plot.cta_ort ----------------------------------------------------------- #

#' Plot method for Optimal Recursive Tree
#'
#' Renders the composite ORT using G1 base-R conventions: ellipses for split
#' nodes, rectangles for terminal nodes, directed arrows for edges.
#'
#' @param x A \code{cta_ort} object.
#' @param target_class Integer target class for terminal node annotation;
#'   \code{NULL} for a structural plot.
#' @param class_labels Optional named character vector of class display names.
#' @param digits Decimal places for proportion labels.  Default \code{1}.
#' @param main Plot title.  Default \code{"ORT"}.
#' @param split_fill Fill color for split (internal) ellipse nodes.
#' @param endpoint_fill Default fill for terminal rectangle nodes.
#' @param endpoint_palette Palette for terminal nodes when \code{target_class}
#'   is supplied; a \code{function(n)} or character color vector.
#' @param border_col Border color for all nodes.  Default \code{"grey30"}.
#' @param text_col Text color for node labels.  Default \code{"black"}.
#' @param edge_col Color for directed edge arrows.  Default \code{"grey40"}.
#' @param arrow_col Arrow color; \code{NULL} (default) uses \code{edge_col}.
#' @param show_caption Logical; add color-encoding caption when
#'   \code{target_class} is supplied.  Default \code{FALSE}.
#' @param cex Text expansion factor.  Default \code{0.75}.
#' @param ... Unused.
#' @return \code{invisible(pd)}, the layout list from \code{\link{ort_plot_data}}.
#' @seealso \code{\link{ort_plot_data}}, \code{\link{plot.cta_tree}}
#' @export
plot.cta_ort <- function(x,
  target_class     = NULL,
  class_labels     = NULL,
  digits           = 1L,
  main             = "ORT",
  split_fill       = "#D9EAF7",
  endpoint_fill    = "#D9F7E6",
  endpoint_palette = NULL,
  border_col       = "grey30",
  text_col         = "black",
  edge_col         = "grey40",
  arrow_col        = NULL,
  show_caption     = FALSE,
  cex              = 0.75,
  ...
) {
  arrow_col <- arrow_col %||% edge_col
  pd        <- ort_plot_data(x, target_class = target_class,
                             class_labels = class_labels, digits = digits)
  nodes_df  <- pd$nodes
  edges_df  <- pd$edges
  strata    <- pd$strata

  if (nrow(nodes_df) == 0L) {
    plot.new()
    title(main = main)
    return(invisible(pd))
  }

  # Determine terminal fill colors
  term_mask <- nodes_df$is_terminal
  n_term    <- sum(term_mask)

  if (!is.null(target_class) && n_term > 0L) {
    if (!is.null(endpoint_palette)) {
      ep_colors <- if (is.function(endpoint_palette)) endpoint_palette(n_term)
                   else colorRampPalette(endpoint_palette)(n_term)
    } else {
      ep_colors <- colorRampPalette(c("#FFFFFF", "#1A9641"))(n_term)
    }
    # Map colors by ascending prop_class1 rank
    term_rows  <- nodes_df[term_mask, ]
    strat_prop <- if (!is.null(strata) && nrow(strata) > 0L)
                    setNames(strata$prop_class1, as.character(strata$node_id))
                  else numeric(0)
    props  <- strat_prop[as.character(term_rows$node_id)]
    props[is.na(props)] <- 0
    rnk    <- rank(props, ties.method = "first")
    fill_map <- setNames(ep_colors[rnk], as.character(term_rows$node_id))
  } else {
    term_rows <- nodes_df[term_mask, ]
    fill_map  <- setNames(rep(endpoint_fill, nrow(term_rows)),
                          as.character(term_rows$node_id))
  }

  x_range <- range(nodes_df$x, na.rm = TRUE)
  y_range <- range(nodes_df$y, na.rm = TRUE)
  x_pad   <- max(0.6, diff(x_range) * 0.15 + 0.4)
  y_pad   <- 0.6

  bot_margin <- if (isTRUE(show_caption) && !is.null(target_class)) 3 else 1
  op <- par(mar = c(bot_margin, 1, 3, 1))
  on.exit(par(op), add = TRUE)

  plot.new()
  plot.window(
    xlim = c(x_range[1] - x_pad, x_range[2] + x_pad),
    ylim = c(y_range[1] - y_pad, y_range[2] + y_pad)
  )
  title(main = main)

  # Draw edges first (below nodes)
  if (nrow(edges_df) > 0L) {
    for (i in seq_len(nrow(edges_df))) {
      x0  <- edges_df$x0[i]; y0 <- edges_df$y0[i]
      x1  <- edges_df$x1[i]; y1 <- edges_df$y1[i]
      lbl <- edges_df$label[i]
      dx  <- x1 - x0; dy <- y1 - y0; len <- sqrt(dx^2 + dy^2)
      if (len > 0) {
        shrink <- 0.22
        x0s <- x0 + dx/len * shrink; y0s <- y0 + dy/len * shrink
        x1s <- x1 - dx/len * shrink; y1s <- y1 - dy/len * shrink
      } else {
        x0s <- x0; y0s <- y0; x1s <- x1; y1s <- y1
      }
      segments(x0s, y0s, x1s, y1s, col = edge_col, lwd = 1.2)
      arrows(x0s, y0s, x1s, y1s,
             length = 0.08, angle = 20, col = arrow_col, lwd = 1.2, code = 2)
      if (!is.na(lbl) && nchar(lbl) > 0L) {
        mx <- (x0 + x1) / 2; my <- (y0 + y1) / 2
        text(mx, my, labels = lbl, cex = cex * 0.85,
             col = text_col, adj = c(0.5, -0.2))
      }
    }
  }

  # Draw nodes
  ew <- 0.38; eh <- 0.22   # ellipse half-widths
  rw <- 0.36; rh <- 0.20   # rect half-widths

  for (i in seq_len(nrow(nodes_df))) {
    row    <- nodes_df[i, ]
    nid_c  <- as.character(row$node_id)
    cx     <- row$x; cy <- row$y
    lbl    <- row$label %||% ""

    if (isTRUE(row$is_terminal)) {
      fill_c <- fill_map[nid_c] %||% endpoint_fill
      rect(cx - rw, cy - rh, cx + rw, cy + rh,
           col = fill_c, border = border_col)
    } else {
      theta <- seq(0, 2*pi, length.out = 60)
      polygon(cx + ew * cos(theta), cy + eh * sin(theta),
              col = split_fill, border = border_col)
    }

    # Multi-line label
    lines_v <- strsplit(lbl, "\n", fixed = TRUE)[[1]]
    n_ln    <- length(lines_v)
    y_top   <- cy + (n_ln - 1) * 0.065
    for (li in seq_along(lines_v)) {
      text(cx, y_top - (li - 1L) * 0.13, lines_v[li],
           cex = cex, col = text_col, adj = c(0.5, 0.5))
    }
  }

  if (isTRUE(show_caption) && !is.null(target_class)) {
    mtext(
      paste0("Endpoint fill: relative target-class proportion within this tree.",
             " Not a clinical threshold."),
      side = 1, line = 1.5, cex = cex * 0.8, col = "grey40"
    )
  }

  invisible(pd)
}
