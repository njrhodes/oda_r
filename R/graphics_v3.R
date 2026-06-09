###############################################################################
# R/graphics_v3.R  -  Graphics v3: direct ggplot2 tree renderers
#
# Public API:
#   v3C1 - tree renderers:
#     plot_cta_tree()            -  ggplot CTA tree diagram
#     plot_lort_tree()           -  ggplot LORT (recursive CTA) tree diagram
#   v3C2 - balance renderers:
#     plot_oda_balance()         -  ODA covariate balance dot plot
#     plot_smd_balance()         -  SMD absolute-value covariate balance dot plot
#     plot_balance_love()        -  Love-plot wrapper (calls plot_smd_balance)
#     plot_cta_balance()         -  CTA multivariate balance: tree or message panel
#   v3C3 - evidence-interval balance renderers:
#     plot_oda_balance_effects() -  forest plot for oda_balance_effect_table
#     plot_cta_balance_effects() -  evidence card for cta_balance_effect_summary
#
# Rules:
#   - All functions return ggplot objects.
#   - No model fitting inside any plot function.
#   - ggplot2 is in Suggests; missing -> clear error via .require_ggplot2().
#   - base plot methods (plot.cta_tree, plot.cta_ort) are unchanged.
#   - Use ggplot2:: prefix; do not attach ggplot2.
###############################################################################

# Suppress R CMD CHECK "no visible binding" notes for ggplot2 aes() column names
utils::globalVariables(c(
  # v3C1  -  tree rendering
  "x", "y", "label",
  "x0", "y0_adj", "x1", "y1_adj", "xmid", "ymid",
  "xmin", "xmax", "ymin", "ymax",
  "fill_value", "pred_chr",
  # v3C1 depth bands (LORT)
  "fill_band",
  # v3C2  -  balance rendering
  "attr_fac", "ess_display_bar",
  "sig_color", "significance_label",
  "abs_smd", "bal_col",
  # v4 simple renderer
  "x0_r", "y0_r", "x1_r", "y1_r", "xmid_r", "ymid_r",
  # v3C3  -  evidence-interval balance renderers
  "attribute", "estimate", "boot_lo", "boot_hi",
  "chance_lo", "chance_hi", "analysis"
))

# ---- Internal helpers ------------------------------------------------------- #

.require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop(
      "Install ggplot2 to use Graphics v3 plotting functions:\n",
      "  install.packages(\"ggplot2\")",
      call. = FALSE
    )
  invisible(TRUE)
}

# ---- v4 label builders --------------------------------------------------- #

# Build a canonical node label for a CTA split or terminal node.
# node_row: single-row data.frame from the normalized nodes frame.
# label_detail: "simple" (default) | "full"
# show_node_ess: if TRUE, append ESS/WESS to split-node labels.
# ess_label: "ESS" or "WESS" (passed from the renderer).
.build_cta_node_label <- function(node_row,
                                   label_detail  = c("simple", "full"),
                                   show_node_ess = FALSE,
                                   ess_label     = "ESS",
                                   show_p        = TRUE,
                                   show_loo      = TRUE) {
  label_detail <- match.arg(label_detail)
  is_lf        <- isTRUE(node_row$is_leaf) || isTRUE(node_row$leaf)

  if (is_lf) {
    # Terminal node  --------------------------------------------------------
    # When target enrichment columns are present, use them.
    stg  <- node_row$stage      %||% NA_integer_
    denom <- node_row$denominator %||% NA_real_
    prop  <- node_row$target_proportion %||% NA_real_
    pred  <- node_row$predicted_label %||% NA_character_
    if (is.na(pred))
      pred <- paste0("class ", node_row$majority_class %||% NA_integer_)

    if (!is.na(stg) && !is.na(denom) && !is.na(prop)) {
      rate_pct <- round(prop * 100)
      return(sprintf("Stage %d\n%s\nn = %d\nrate = %d%%",
                     as.integer(stg), pred,
                     as.integer(round(denom)), rate_pct))
    } else {
      return(sprintf("%s\nn = %d", pred,
                     as.integer(node_row$n_obs %||% NA_integer_)))
    }
  }

  # Split node  ---------------------------------------------------------------
  # Use split_cond if available (first outgoing edge label = full condition).
  # Fall back to attribute name only.
  cond  <- node_row$split_cond %||% NA_character_
  attr_nm <- node_row$attribute %||% "?"
  n_val   <- as.integer(node_row$n_obs %||% NA_integer_)

  if (!is.na(cond) && nchar(cond) > 0L) {
    base <- sprintf("%s\nn = %d", cond, n_val)
  } else {
    base <- sprintf("%s\nn = %d", attr_nm, n_val)
  }

  if (isTRUE(show_node_ess)) {
    ess_v <- node_row$ess %||% NA_real_
    if (!is.na(ess_v) && is.finite(ess_v))
      base <- sprintf("%s\n%s=%.1f%%", base, ess_label, ess_v)
  }

  # -- p-value and LOO status lines at base of split-node label -------------
  # show_p    = TRUE: append "MC p = X.XXX" when p_mc is available.
  # show_loo  = TRUE: append "LOO: STABLE", "LOO p = X.XXX", or nothing
  #   depending on loo_status:
  #     "STABLE" -> "LOO: STABLE"
  #     "PVALUE" -> "LOO p = X.XXX" (loo_p value)
  #     "OFF" / NA -> nothing
  if (isTRUE(show_p)) {
    pmc_v <- node_row$p_mc %||% NA_real_
    if (!is.na(pmc_v) && is.finite(pmc_v))
      base <- sprintf("%s\nMC p = %.3f", base, pmc_v)
  }

  if (isTRUE(show_loo)) {
    lstatus <- node_row$loo_status %||% NA_character_
    if (!is.na(lstatus) && lstatus == "STABLE") {
      base <- paste0(base, "\nLOO: STABLE")
    } else if (!is.na(lstatus) && lstatus == "PVALUE") {
      loo_p_v <- node_row$loo_p %||% NA_real_
      if (!is.na(loo_p_v) && is.finite(loo_p_v))
        base <- sprintf("%s\nLOO p = %.3f", base, loo_p_v)
    }
  }

  base
}

# Build a canonical node label for a LORT stratum node.
# Uses the same structure as .build_cta_node_label but adapts for ORT fields.
.build_lort_stage_label <- function(node_row,
                                     label_detail  = c("simple", "full"),
                                     show_node_ess = FALSE,
                                     ess_label     = "ESS") {
  label_detail <- match.arg(label_detail)
  is_term      <- isTRUE(node_row$is_terminal) || isTRUE(node_row$is_leaf)

  n_val <- as.integer(node_row$n %||% node_row$n_obs %||% NA_integer_)
  nid   <- as.integer(node_row$node_id %||% NA_integer_)

  if (is_term) {
    prop  <- node_row$fill_value %||% NA_real_
    cls   <- node_row$predicted_class %||% NA_integer_
    if (!is.na(prop) && !is.na(cls)) {
      return(sprintf("Stage %d\nclass %d\nn = %d\nrate = %d%%",
                     nid, as.integer(cls), n_val, round(prop * 100)))
    } else if (!is.na(cls)) {
      return(sprintf("Stage %d\nclass %d\nn = %d", nid, as.integer(cls), n_val))
    } else {
      return(sprintf("Stage %d\nn = %d", nid, n_val))
    }
  }

  # Non-terminal: show stratum node id and n
  base <- sprintf("Node %d\nn = %d", nid, n_val)
  base
}

# ---- v4 simple tree renderer ----------------------------------------------- #

# New canonical tree renderer.  Owns layout completely (leaf-first algorithm).
# Accepts:
#   nodes : normalized nodes data.frame (node_id, parent_id, depth, leaf/is_leaf,
#           n_obs, attribute, ess, majority_class, + optional endpoint fields)
#   edges : data.frame (from_node_id, to_node_id, label)
#   pal   : .tree_pal() result
#   color_by : "none" | "target_rate" | "prediction"
#   opts  : list of rendering parameters
.render_simple_tree_gg <- function(nodes, edges, pal, color_by, opts) {
  .require_ggplot2()

  # ---- Topology -------------------------------------------------------
  node_ids   <- nodes$node_id
  parent_ids <- nodes$parent_id %||% rep(NA_integer_, nrow(nodes))
  is_leaf    <- nodes$is_leaf

  # Build children list keyed by character node_id
  ch <- setNames(vector("list", length(node_ids)), as.character(node_ids))
  for (i in seq_along(node_ids)) {
    pid <- parent_ids[i]
    if (!is.na(pid) && pid > 0L) {
      pk <- as.character(pid)
      ch[[pk]] <- c(ch[[pk]], node_ids[i])
    }
  }

  # Find root: parent_id == 0 or NA or minimum depth
  root_id <- node_ids[which(is.na(parent_ids) | parent_ids == 0L)]
  if (length(root_id) == 0L) {
    dep_col <- nodes$depth %||% abs(nodes$y)
    root_id <- node_ids[which.min(dep_col)]
  }
  root_id <- root_id[1L]

  # DFS to get leaf order (left-to-right inorder traversal)
  get_leaves_lr <- function(nid) {
    kids <- ch[[as.character(nid)]]
    if (length(kids) == 0L) return(nid)
    unlist(lapply(kids, get_leaves_lr))
  }
  leaves_lr <- get_leaves_lr(root_id)

  # ---- Box dimensions ------------------------------------------------
  max_line_len <- max(vapply(nodes$label %||% character(nrow(nodes)), function(l) {
    segs <- strsplit(as.character(l), "\n", fixed = TRUE)[[1L]]
    max(nchar(segs))
  }, integer(1L)), na.rm = TRUE)
  box_w <- max(1.8, min(3.6, max_line_len * 0.115))
  half_w <- box_w / 2

  max_lines <- max(vapply(nodes$label %||% character(nrow(nodes)), function(l) {
    length(strsplit(as.character(l), "\n", fixed = TRUE)[[1L]])
  }, integer(1L)), na.rm = TRUE)
  half_h <- max(0.30, 0.12 + 0.09 * max(1L, max_lines))

  level_gap <- max(2.0, half_h * 2 + 0.8)

  # ---- Leaf-first x layout -------------------------------------------
  leaf_gap   <- box_w + 0.8
  x_asgn     <- setNames(numeric(length(node_ids)), as.character(node_ids))

  do_assign <- function(nid) {
    kids <- ch[[as.character(nid)]]
    if (length(kids) == 0L) {
      leaf_idx <- which(leaves_lr == nid) - 1L
      x_asgn[as.character(nid)] <<- leaf_idx * leaf_gap
      return(invisible(NULL))
    }
    for (kid in kids) do_assign(kid)
    kxs <- x_asgn[as.character(kids)]
    x_asgn[as.character(nid)] <<- mean(kxs, na.rm = TRUE)
  }
  do_assign(root_id)

  # Check and fix sibling overlap (up to 10 iterations)
  for (.it in seq_len(10L)) {
    overlap <- FALSE
    for (nid in node_ids[!is_leaf]) {
      kids <- ch[[as.character(nid)]]
      if (length(kids) < 2L) next
      xs <- sort(x_asgn[as.character(kids)])
      if (any(diff(xs) < (box_w + 0.4))) { overlap <- TRUE; break }
    }
    if (!overlap) break
    leaf_gap <- leaf_gap + 0.4
    do_assign(root_id)
  }

  nodes$x <- x_asgn[as.character(nodes$node_id)]
  nodes$y <- -(nodes$depth * level_gap)

  # ---- Attach split condition to split nodes -------------------------
  # Use first outgoing edge label (left/lower branch) as the node condition.
  if (nrow(edges) > 0L && "from_node_id" %in% names(edges)) {
    first_e <- edges[!duplicated(edges$from_node_id), c("from_node_id", "label"),
                     drop = FALSE]
    nodes$split_cond <- first_e$label[match(nodes$node_id, first_e$from_node_id)]
  } else {
    nodes$split_cond <- NA_character_
  }

  # ---- Rebuild node labels -------------------------------------------
  ess_lbl_str <- opts$ess_label %||% "ESS"
  nodes$label <- vapply(seq_len(nrow(nodes)), function(i) {
    .build_cta_node_label(
      nodes[i, , drop = FALSE],
      label_detail  = opts$label_detail  %||% "simple",
      show_node_ess = isTRUE(opts$show_node_ess),
      ess_label     = ess_lbl_str,
      show_p        = isTRUE(opts$show_p),
      show_loo      = isTRUE(opts$show_loo)
    )
  }, character(1L))

  # ---- Edge coordinates from new node positions ----------------------
  if (nrow(edges) > 0L && "from_node_id" %in% names(edges)) {
    edges$x0_r  <- nodes$x[match(edges$from_node_id, nodes$node_id)]
    edges$y0_r  <- nodes$y[match(edges$from_node_id, nodes$node_id)] - half_h
    edges$x1_r  <- nodes$x[match(edges$to_node_id,   nodes$node_id)]
    edges$y1_r  <- nodes$y[match(edges$to_node_id,   nodes$node_id)] + half_h
    edges$xmid_r <- (edges$x0_r + edges$x1_r) / 2
    edges$ymid_r <- (edges$y0_r + edges$y1_r) / 2

    # Strip attribute prefix from edge labels (keep just the condition)
    if (isTRUE(opts$short_edge_labels)) {
      edges$label <- sub("^[A-Za-z_][A-Za-z0-9_.]*\\s*(<=?|>=?)",
                         "\\1", edges$label)
      edges$label <- sub("^[A-Za-z_][A-Za-z0-9_.]*\\s+in\\b",
                         "in", edges$label)
    }

    edges <- edges[!is.na(edges$x0_r) & !is.na(edges$x1_r), , drop = FALSE]
  }

  # ---- Box bounds ----------------------------------------------------
  nodes$xmin <- nodes$x - half_w
  nodes$xmax <- nodes$x + half_w
  nodes$ymin <- nodes$y - half_h
  nodes$ymax <- nodes$y + half_h

  nd_int <- nodes[!is_leaf, , drop = FALSE]
  nd_lf  <- nodes[ is_leaf,  , drop = FALSE]

  # ---- Build plot ----------------------------------------------------
  p <- ggplot2::ggplot()

  # Edges
  if (nrow(edges) > 0L && "x0_r" %in% names(edges)) {
    p <- p + ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = x0_r, y = y0_r, xend = x1_r, yend = y1_r),
      color     = pal$edge_color,
      linewidth = 0.55,
      arrow     = ggplot2::arrow(length = ggplot2::unit(0.09, "inches"),
                                 type = "closed", ends = "last")
    )
    if (isTRUE(opts$show_rule) && any(nchar(edges$label) > 0L)) {
      el <- edges[nchar(edges$label) > 0L, , drop = FALSE]
      p <- p + ggplot2::geom_label(
        data = el,
        ggplot2::aes(x = xmid_r, y = ymid_r, label = label),
        size          = opts$edge_text_size %||% 3.2,
        linewidth     = 0,
        label.padding = ggplot2::unit(0.08, "lines"),
        fill          = "white",
        color         = "#333333",
        lineheight    = 0.9
      )
    }
  }

  # Internal (split) nodes -- ellipse shape (MPE CTA canon: circles for splits)
  if (nrow(nd_int) > 0L) {
    ell_pts2 <- do.call(rbind, lapply(seq_len(nrow(nd_int)), function(i) {
      nd    <- nd_int[i, , drop = FALSE]
      theta <- seq(0, 2 * pi, length.out = 41L)[-41L]
      data.frame(
        x    = nd$x + half_w * cos(theta),
        y    = nd$y + half_h * sin(theta),
        .eid = nd$node_id,
        stringsAsFactors = FALSE
      )
    }))
    p <- p +
      ggplot2::geom_polygon(
        data = ell_pts2,
        ggplot2::aes(x = x, y = y, group = .eid),
        fill      = "white",
        color     = pal$internal_color,
        linewidth = 0.7
      ) +
      ggplot2::geom_text(
        data = nd_int,
        ggplot2::aes(x = x, y = y, label = label),
        size       = opts$node_text_size %||% 3.5,
        lineheight = 0.95,
        vjust      = 0.5
      )
  }

  # Terminal (leaf) nodes -- fill varies by color_by
  if (nrow(nd_lf) > 0L) {
    has_fill <- "fill_value" %in% names(nd_lf) && !all(is.na(nd_lf$fill_value))
    has_pred <- "predicted_class" %in% names(nd_lf) && !all(is.na(nd_lf$predicted_class))

    if (color_by == "target_rate" && has_fill) {
      p <- p +
        ggplot2::geom_rect(
          data = nd_lf,
          ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                       fill = fill_value),
          color = "#333333", linewidth = 1.0
        ) +
        ggplot2::scale_fill_gradient(
          low      = pal$gradient_low,
          high     = pal$gradient_high,
          na.value = "white",
          name     = "Target\nrate",
          labels   = function(v) paste0(round(v * 100), "%")
        )
    } else if (color_by == "prediction" && has_pred) {
      nd_lf$pred_chr <- as.character(nd_lf$predicted_class)
      p <- p +
        ggplot2::geom_rect(
          data = nd_lf,
          ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                       fill = pred_chr),
          color = "#333333", linewidth = 1.0
        ) +
        ggplot2::scale_fill_manual(
          values   = pal$pred_colors,
          na.value = "white",
          name     = "Prediction"
        )
    } else {
      p <- p + ggplot2::geom_rect(
        data = nd_lf,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill      = "white",
        color     = "#333333",
        linewidth = 1.0
      )
    }

    p <- p + ggplot2::geom_text(
      data = nd_lf,
      ggplot2::aes(x = x, y = y, label = label),
      size       = opts$node_text_size %||% 3.5,
      lineheight = 0.95,
      vjust      = 0.5
    )
  }

  # Theme
  p <- p +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = if (color_by == "none") "none" else "right",
      plot.title      = ggplot2::element_text(face = "bold", size = 12,
                                              hjust = 0.5),
      plot.subtitle   = ggplot2::element_text(size = 10, hjust = 0.5,
                                              color = "#555555")
    )

  if (!is.null(opts$main))
    p <- p + ggplot2::ggtitle(label = opts$main, subtitle = opts$subtitle)

  p
}

# Return a ggplot message panel (no_tree, empty ORT, error fallback).
# Used for generic errors and structural empties.
.gg_message_panel <- function(msg, title = NULL, subtitle = NULL) {
  .require_ggplot2()
  p <- ggplot2::ggplot(
    data.frame(x = 0.5, y = 0.5, label = msg),
    ggplot2::aes(x = x, y = y, label = label)
  ) +
    ggplot2::geom_text(size = 4.2, hjust = 0.5, vjust = 0.5,
                       color = "#444444", lineheight = 1.4) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::theme_void()
  if (!is.null(title))
    p <- p + ggplot2::ggtitle(title, subtitle = subtitle)
  p
}

# Styled result card for favorable balance / no-tree outcomes.
# Produces a bordered panel with a positive interpretation rather than a
# blank void canvas (which looks like an error screen).
.gg_result_card <- function(header, body, footer = NULL,
                             title = NULL, subtitle = NULL,
                             header_color = "#1a7a4a") {
  .require_ggplot2()
  body_y   <- if (!is.null(footer)) 0.48 else 0.50
  footer_y <- 0.20
  box_df   <- data.frame(xmin = 0.08, xmax = 0.92, ymin = 0.12, ymax = 0.90)

  p <- ggplot2::ggplot() +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 10,  hjust = 0.5,
                                            color = "#555555")
    ) +
    # Card border
    ggplot2::geom_rect(
      data = box_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "#f9fdf9", color = "#6abd8e", linewidth = 1.2
    ) +
    # Header band
    ggplot2::geom_rect(
      data = data.frame(xmin = 0.08, xmax = 0.92, ymin = 0.72, ymax = 0.90),
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "#e8f6ef", color = NA
    ) +
    # Header text
    ggplot2::annotate("text", x = 0.50, y = 0.81,
                      label = header, hjust = 0.5, vjust = 0.5,
                      size = 4.8, fontface = "bold", color = header_color,
                      lineheight = 1.1) +
    # Body text
    ggplot2::annotate("text", x = 0.50, y = body_y,
                      label = body, hjust = 0.5, vjust = 0.5,
                      size = 3.8, color = "#333333", lineheight = 1.4)

  if (!is.null(footer))
    p <- p + ggplot2::annotate("text", x = 0.50, y = footer_y,
                                label = footer, hjust = 0.5, vjust = 0.5,
                                size = 3.2, color = "#777777", lineheight = 1.2)

  if (!is.null(title))
    p <- p + ggplot2::ggtitle(title, subtitle = subtitle)
  p
}

# Default color palette for tree diagrams.
.tree_pal <- function(palette = NULL) {
  list(
    internal_fill  = if (!is.null(palette[["internal"]])) palette[["internal"]] else "#d9eaf7",
    internal_color = "#3a6ea5",
    leaf_color     = "#3a6ea5",
    gradient_low   = if (!is.null(palette[["low"]]))  palette[["low"]]  else "#f5f5f5",
    gradient_high  = if (!is.null(palette[["high"]])) palette[["high"]] else "#b82929",
    pred_colors    = c("0" = "#7ec8e3", "1" = "#e07b5a"),
    edge_color     = "#555555"
  )
}

# Coerce input to cta_plot_data output.
# Accepts: cta_tree OR list(nodes, edges, ...) from cta_plot_data().
.coerce_cta_plot_data <- function(x, target_class, digits) {
  if (inherits(x, "cta_tree"))
    return(cta_plot_data(x, target_class = as.integer(target_class),
                         digits = as.numeric(digits)))
  if (is.list(x) && all(c("nodes", "edges") %in% names(x)))
    return(x)
  stop(
    "plot_cta_tree: 'x' must be a cta_tree object or cta_plot_data() output.",
    call. = FALSE
  )
}

# Coerce input to ort_plot_data output.
# Accepts: cta_ort OR list(nodes, edges, ...) from ort_plot_data().
.coerce_lort_plot_data <- function(x, target_class, digits) {
  if (inherits(x, "cta_ort"))
    return(ort_plot_data(x, target_class = as.integer(target_class),
                         digits = as.integer(digits)))
  if (is.list(x) && all(c("nodes", "edges") %in% names(x)))
    return(x)
  stop(
    "plot_lort_tree: 'x' must be a cta_ort object or ort_plot_data() output.",
    call. = FALSE
  )
}

# Normalize CTA plot-data nodes to the shared rendering schema.
# Returns nodes data.frame with: x, y, is_leaf, label, n,
#   fill_value (numeric 0-1 for target_rate), predicted_class (integer).
.normalize_cta_nodes <- function(nodes) {
  nodes$is_leaf    <- isTRUE(nodes$leaf) | nodes$leaf
  nodes$n          <- nodes$n_obs
  # fill_value: target_proportion when available (target_class was supplied)
  if ("target_proportion" %in% names(nodes)) {
    nodes$fill_value <- nodes$target_proportion
  } else {
    nodes$fill_value <- NA_real_
  }
  nodes$predicted_class <- nodes$majority_class
  nodes
}

# Normalize ORT plot-data nodes to the shared rendering schema.
.normalize_lort_nodes <- function(nodes, strata) {
  nodes$is_leaf         <- nodes$is_terminal
  # n is already nodes$n for ORT
  nodes$fill_value      <- NA_real_
  nodes$predicted_class <- NA_integer_

  if (!is.null(strata) && nrow(strata) > 0L) {
    m <- match(nodes$node_id, strata$node_id)
    ok <- !is.na(m)
    nodes$fill_value[ok]      <- strata$prop_class1[m[ok]]
    nodes$predicted_class[ok] <- strata$terminal_class[m[ok]]
  }
  nodes
}

# Core ggplot tree renderer.
# nodes: normalized data.frame (is_leaf, x, y, label, fill_value, predicted_class)
# edges: data.frame (x0, y0, x1, y1, label)
# pal: .tree_pal() result
# color_by: "target_rate" | "prediction" | "none"
# opts: list of rendering parameters
.render_tree_gg <- function(nodes, edges, pal, color_by, opts) {
  .require_ggplot2()

  # Node box dimensions in data units (y-axis: depth * -1, so 1 unit = 1 level)
  # Adaptive height: scale with the tallest label across all nodes (uniform height
  # keeps siblings aligned; per-node height would misalign siblings).
  max_lines <- max(vapply(nodes$label, function(l) {
    length(strsplit(as.character(l), "\n", fixed = TRUE)[[1L]])
  }, integer(1L)), na.rm = TRUE)
  half_h <- max(0.25, 0.12 + 0.09 * max(1L, max_lines))

  # Auto-compute node half-width from longest single line in all labels
  nw <- opts$node_width
  if (is.null(nw) || is.na(nw)) {
    max_line <- max(vapply(nodes$label, function(l) {
      segs <- strsplit(as.character(l), "\n", fixed = TRUE)[[1L]]
      max(nchar(segs))
    }, integer(1L)), na.rm = TRUE)
    nw <- max(1.4, min(3.8, max_line * 0.105))
  }
  half_w <- nw / 2

  # Short edge labels: strip attribute-name prefix, keeping operator + value.
  # "V14<=0.5" -> "<= 0.5",  "x3>0.5" -> "> 0.5".  Non-matching labels pass through.
  if (isTRUE(opts$short_edge_labels) && nrow(edges) > 0L) {
    edges$label <- sub("^[A-Za-z_][A-Za-z0-9_.]*(<=?|>=?)", "\\1 ", edges$label)
  }

  # Add rect corners
  nodes$xmin <- nodes$x - half_w
  nodes$xmax <- nodes$x + half_w
  nodes$ymin <- nodes$y - half_h
  nodes$ymax <- nodes$y + half_h

  # Adjust edges to attach to box edges (not node centers)
  edges$y0_adj <- edges$y0 - half_h   # parent bottom
  edges$y1_adj <- edges$y1 + half_h   # child top
  edges$xmid   <- (edges$x0 + edges$x1) / 2
  edges$ymid   <- (edges$y0_adj + edges$y1_adj) / 2

  nd_int <- nodes[!nodes$is_leaf, , drop = FALSE]
  nd_lf  <- nodes[nodes$is_leaf,  , drop = FALSE]

  p <- ggplot2::ggplot()

  # Depth bands (LORT mode): shaded horizontal stripes at each depth level to
  # visually separate stratum layers.  Only drawn when opts$depth_bands == TRUE.
  if (isTRUE(opts$depth_bands) && nrow(nodes) > 0L) {
    depths  <- sort(unique(nodes$depth %||% abs(nodes$y)))
    x_range <- range(c(nodes$xmin, nodes$xmax), na.rm = TRUE)
    band_df <- do.call(rbind, lapply(seq_along(depths), function(i) {
      d <- depths[i]
      y_ctr   <- -d
      data.frame(
        xmin = x_range[1L] - 0.1,
        xmax = x_range[2L] + 0.1,
        ymin = y_ctr - 0.5,
        ymax = y_ctr + 0.5,
        fill_band = if (i %% 2L == 0L) "#f0f4f8" else "#fafafa",
        stringsAsFactors = FALSE
      )
    }))
    p <- p + ggplot2::geom_rect(
      data = band_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = band_df$fill_band, color = NA, alpha = 0.6
    )
  }

  # ---- Edges ------------------------------------------------------------------
  if (nrow(edges) > 0L) {
    p <- p + ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = x0, y = y0_adj, xend = x1, yend = y1_adj),
      color     = pal$edge_color,
      linewidth = 0.55,
      arrow     = ggplot2::arrow(length = ggplot2::unit(0.09, "inches"),
                                 type = "closed", ends = "last")
    )
    if (isTRUE(opts$show_rule) && any(nchar(edges$label) > 0L)) {
      edges_lbl <- edges[nchar(edges$label) > 0L, , drop = FALSE]
      p <- p + ggplot2::geom_label(
        data = edges_lbl,
        ggplot2::aes(x = xmid, y = ymid, label = label),
        size          = opts$edge_text_size,
        linewidth     = 0,
        label.padding = ggplot2::unit(0.08, "lines"),
        fill          = "white",
        color         = "#333333",
        lineheight    = 0.9
      )
    }
  }

  # ---- Internal (split) nodes  -  ellipse shape (MPE canon) ------------------
  # Internal/split nodes are drawn as ellipses (circles) per MPE CTA figure
  # convention. Terminal/leaf nodes remain rectangles.
  if (nrow(nd_int) > 0L) {
    ell_pts <- do.call(rbind, lapply(seq_len(nrow(nd_int)), function(i) {
      nd    <- nd_int[i, , drop = FALSE]
      theta <- seq(0, 2 * pi, length.out = 41L)[-41L]
      data.frame(
        x      = nd$x + half_w * cos(theta),
        y      = nd$y + half_h * sin(theta),
        .eid   = nd$node_id,
        stringsAsFactors = FALSE
      )
    }))
    p <- p +
      ggplot2::geom_polygon(
        data = ell_pts,
        ggplot2::aes(x = x, y = y, group = .eid),
        fill      = pal$internal_fill,
        color     = pal$internal_color,
        linewidth = 0.5
      ) +
      ggplot2::geom_text(
        data = nd_int,
        ggplot2::aes(x = x, y = y, label = label),
        size       = opts$node_text_size,
        lineheight = 0.9,
        vjust      = 0.5
      )
  }

  # ---- Leaf/terminal nodes  -  fill varies by color_by -------------------------
  if (nrow(nd_lf) > 0L) {
    has_fill  <- "fill_value" %in% names(nd_lf) && !all(is.na(nd_lf$fill_value))
    has_pred  <- "predicted_class" %in% names(nd_lf) &&
                 !all(is.na(nd_lf$predicted_class))

    if (color_by == "target_rate" && has_fill) {
      p <- p +
        ggplot2::geom_rect(
          data = nd_lf,
          ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                       fill = fill_value),
          color = pal$leaf_color, linewidth = 0.5
        ) +
        ggplot2::scale_fill_gradient(
          low      = pal$gradient_low,
          high     = pal$gradient_high,
          na.value = "white",
          name     = "Target\nrate",
          labels   = function(v) paste0(round(v * 100), "%")
        )
    } else if (color_by == "prediction" && has_pred) {
      nd_lf$pred_chr <- as.character(nd_lf$predicted_class)
      p <- p +
        ggplot2::geom_rect(
          data = nd_lf,
          ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                       fill = pred_chr),
          color = pal$leaf_color, linewidth = 0.5
        ) +
        ggplot2::scale_fill_manual(
          values   = pal$pred_colors,
          na.value = "white",
          name     = "Prediction"
        )
    } else {
      # "none" or fallback
      p <- p + ggplot2::geom_rect(
        data = nd_lf,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill      = "white",
        color     = pal$leaf_color,
        linewidth = 0.5
      )
    }

    # LORT terminal nodes: overlay a dashed border to visually distinguish
    # final strata from intra-stratum split nodes.
    if (isTRUE(opts$lort_terminal_dash) && nrow(nd_lf) > 0L) {
      p <- p + ggplot2::geom_rect(
        data = nd_lf,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill      = NA,
        color     = "#2c6e8a",
        linewidth = 0.85,
        linetype  = "dashed"
      )
    }

    p <- p + ggplot2::geom_text(
      data = nd_lf,
      ggplot2::aes(x = x, y = y, label = label),
      size       = opts$node_text_size,
      lineheight = 0.9,
      vjust      = 0.5
    )
  }

  # ---- Theme ------------------------------------------------------------------
  base_theme <- switch(opts$theme,
    "minimal" = ggplot2::theme_minimal(),
    ggplot2::theme_void()   # "clean" and default
  )

  p <- p +
    base_theme +
    ggplot2::theme(
      axis.text        = ggplot2::element_blank(),
      axis.title       = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank(),
      panel.grid       = ggplot2::element_blank(),
      legend.position  = if (color_by == "none") "none" else "right",
      plot.title       = ggplot2::element_text(face = "bold", size = 12,
                                               hjust = 0.5),
      plot.subtitle    = ggplot2::element_text(size = 10, hjust = 0.5,
                                               color = "#555555")
    )

  if (!is.null(opts$main))
    p <- p + ggplot2::ggtitle(label    = opts$main,
                               subtitle = opts$subtitle)

  p
}

# ---- plot_cta_tree ---------------------------------------------------------- #

#' Plot a CTA tree using ggplot2
#'
#' Renders a publication-quality tree diagram for a fitted CTA tree.  Requires
#' the \pkg{ggplot2} package (listed in \code{Suggests}); if unavailable, a
#' clear error is raised.
#'
#' @param x A \code{cta_tree} object from \code{\link{cta_fit}} /
#'   \code{\link{oda_cta_fit}}, or the list returned by
#'   \code{\link{cta_plot_data}}.
#' @param target_class Integer; target class for endpoint coloring and
#'   target-rate annotation (default \code{1L}).  Ignored when \code{x} is
#'   already \code{cta_plot_data} output.
#' @param color_by Character; controls leaf-node fill color.
#'   \code{"none"} (default): white fill (B/W publication default);
#'   \code{"target_rate"}: continuous gradient by target-class proportion;
#'   \code{"prediction"}: discrete fill by predicted class.
#' @param label_detail Character; node label verbosity.  \code{"simple"}
#'   (default): canonical MPE-style labels (attribute + condition + n for
#'   split nodes; Stage/class/n/rate for terminal nodes).  \code{"full"}:
#'   same content (reserved for future extension).
#' @param show_node_ess Logical; if \code{TRUE}, append the node-level
#'   ESS/WESS to split-node labels.  Default \code{FALSE}.
#' @param show_p Logical; if \code{TRUE} (default), append \code{MC p = X.XXX}
#'   to each split-node label when the MC permutation p-value is available.
#' @param show_loo Logical; if \code{TRUE} (default), append the LOO result
#'   to each split-node label: \code{LOO: STABLE} when \code{loo = "stable"},
#'   \code{LOO p = X.XXX} when \code{loo = "pvalue"}.  Nothing is shown when
#'   \code{loo = "off"}.
#' @param main Character; plot title.  Default: auto-generated from tree
#'   structure (n, endpoints, ESS/D).
#' @param subtitle Character; plot subtitle.
#' @param show_rule Logical; show branch condition labels on edges.
#'   Default \code{TRUE}.
#' @param show_metrics Logical; if \code{TRUE}, appends an ESS/WESS and D
#'   line to the plot subtitle.  Default \code{FALSE}.
#' @param short_edge_labels Logical; if \code{TRUE} (default), strip the
#'   attribute-name prefix from edge labels so that \code{"x1 <= 24.5"}
#'   renders as \code{"<= 24.5"}.
#' @param node_text_size Numeric; ggplot text size for node labels.
#'   Default \code{3.5}.
#' @param edge_text_size Numeric; ggplot text size for edge labels.
#'   Default \code{3.2}.
#' @param palette Named list for color overrides: \code{internal}, \code{low},
#'   \code{high}.  \code{NULL} uses defaults.
#' @return A \code{\link[ggplot2]{ggplot}} object.  Print it, modify it, or
#'   save with \code{ggplot2::ggsave()}.
#' @seealso \code{\link{cta_fit}}, \code{\link{cta_plot_data}},
#'   \code{\link{plot_lort_tree}}, \code{\link{plot_cta_family}},
#'   \code{\link[ggplot2]{ggsave}}
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   X <- data.frame(x1 = c(1,2,3,4,5,6,7,8),
#'                   x2 = c(0L,0L,1L,0L,1L,1L,0L,1L))
#'   y <- c(1L,1L,1L,1L,2L,2L,2L,2L)
#'   tree <- cta_fit(X, y, mindenom=1L, mc_iter=500L, mc_seed=42L, loo="off")
#'   p <- plot_cta_tree(tree)
#'   print(p)
#' }
#' }
#' @export
plot_cta_tree <- function(x,
                           target_class       = 1L,
                           color_by           = c("none", "target_rate", "prediction"),
                           label_detail       = c("simple", "full"),
                           show_node_ess      = FALSE,
                           show_p             = TRUE,
                           show_loo           = TRUE,
                           main               = NULL,
                           subtitle           = NULL,
                           show_rule          = TRUE,
                           show_metrics       = FALSE,
                           short_edge_labels  = TRUE,
                           node_text_size     = 3.5,
                           edge_text_size     = 3.2,
                           palette            = NULL) {
  .require_ggplot2()
  color_by     <- match.arg(color_by)
  label_detail <- match.arg(label_detail)

  # Coerce to plot-data
  pd <- .coerce_cta_plot_data(x, target_class = target_class, digits = 1)

  # no_tree: return a styled result card (not a black void)
  if (isTRUE(pd$no_tree)) {
    n_val  <- pd$training_n %||% NA_integer_
    footer <- if (!is.na(n_val)) sprintf("n = %d", n_val) else NULL
    return(.gg_result_card(
      header       = "No CTA Tree Found",
      body         = "No valid split found above the significance threshold.\nAll candidate attributes failed the minimum denominator or LOO gate.",
      footer       = footer,
      title        = main %||% "CTA Tree",
      subtitle     = subtitle,
      header_color = "#7a3a1a"
    ))
  }

  nodes <- pd$nodes
  edges <- pd$edges

  if (nrow(nodes) == 0L)
    return(.gg_message_panel("No nodes to display.", title = main, subtitle = subtitle))

  nodes <- .normalize_cta_nodes(nodes)

  pal     <- .tree_pal(palette)
  ess_lbl <- pd$ess_label %||% "ESS"

  # Auto-title
  if (is.null(main)) {
    ess_val <- pd$overall_ess %||% NA_real_
    d_val   <- pd$d           %||% NA_real_
    n_val   <- pd$training_n  %||% NA_integer_
    n_ep    <- sum(nodes$is_leaf, na.rm = TRUE)
    parts   <- "CTA"
    if (!is.na(n_val))   parts <- paste0(parts, " \u2014 n=", n_val)
    if (n_ep > 0L)       parts <- paste0(parts, ", ", n_ep, " endpoint",
                                          if (n_ep != 1L) "s" else "")
    if (!is.na(ess_val)) parts <- paste0(parts, ", ", ess_lbl, "=",
                                          sprintf("%.1f%%", ess_val))
    if (!is.na(d_val) && is.finite(d_val))
      parts <- paste0(parts, ", D=", sprintf("%.4f", d_val))
    main <- parts
  }

  effective_subtitle <- subtitle
  if (isTRUE(show_metrics)) {
    ess_val <- pd$overall_ess %||% NA_real_
    d_val   <- pd$d           %||% NA_real_
    if (!is.na(ess_val)) {
      d_str      <- if (is.finite(d_val)) sprintf("%.4f", d_val) else "NA"
      metrics_ln <- sprintf("%s: %.2f%%  |  D: %s", ess_lbl, ess_val, d_str)
      effective_subtitle <- if (!is.null(subtitle))
        paste0(subtitle, "\n", metrics_ln)
      else
        metrics_ln
    }
  }

  opts <- list(
    main              = main,
    subtitle          = effective_subtitle,
    label_detail      = label_detail,
    show_node_ess     = isTRUE(show_node_ess),
    show_p            = isTRUE(show_p),
    show_loo          = isTRUE(show_loo),
    ess_label         = ess_lbl,
    show_rule         = isTRUE(show_rule),
    short_edge_labels = isTRUE(short_edge_labels),
    node_text_size    = node_text_size,
    edge_text_size    = edge_text_size
  )

  .render_simple_tree_gg(nodes, edges, pal, color_by, opts)
}

# ---- plot_lort_tree --------------------------------------------------------- #

#' Plot a LORT (Locally Optimal Recursive Tree) using ggplot2
#'
#' Renders a publication-quality CTA tree diagram for a single sub-tree within
#' a LORT object (indexed inspection), or a named list of plots for all sub-trees
#' (\code{show_all = TRUE}).  Requires the \pkg{ggplot2} package.
#'
#' @param x A \code{cta_ort} object from \code{\link{lort_fit}}.
#' @param index Integer or character; which LORT node (sub-tree) to render.
#'   Default \code{1L} (root sub-tree).  Ignored when \code{show_all = TRUE}.
#' @param show_all Logical; if \code{TRUE}, return a named list of ggplot
#'   objects, one per LORT node, in node-id order.  Default \code{FALSE}.
#' @param target_class Integer; target class for endpoint coloring and
#'   target-rate annotation (default \code{1L}).
#' @param color_by Character; controls leaf-node fill color.
#'   \code{"none"} (default): white fill; \code{"target_rate"}: gradient;
#'   \code{"prediction"}: discrete fill by predicted class.
#' @param label_detail Character; \code{"simple"} (default) or \code{"full"}.
#' @param show_node_ess Logical; append node-level ESS to split labels.
#'   Default \code{FALSE}.
#' @param show_p Logical; append \code{MC p = X.XXX} to split-node labels.
#'   Default \code{TRUE}.
#' @param show_loo Logical; append \code{LOO: STABLE} or \code{LOO p = X.XXX}
#'   to split-node labels when LOO was active.  Default \code{TRUE}.
#' @param main Character; plot title.  Default: auto-generated.  When
#'   \code{show_all = TRUE} each plot's title includes the node id.
#' @param subtitle Character; plot subtitle.
#' @param show_rule Logical; show branch condition labels on edges.
#' @param show_metrics Logical; append ESS/D to subtitle.  Default \code{FALSE}.
#' @param short_edge_labels Logical; strip attribute-name prefix from edge labels.
#'   Default \code{TRUE}.
#' @param node_text_size Numeric; text size for node labels.  Default \code{3.5}.
#' @param edge_text_size Numeric; text size for edge labels.  Default \code{3.2}.
#' @param palette Named list for color overrides.  \code{NULL} uses defaults.
#' @param show_path Logical; if \code{TRUE}, delegates to
#'   \code{\link{plot_lort_path}} to render the full path from root to
#'   \code{index}.  Default \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link{plot_lort_path}}
#'   when \code{show_path = TRUE}; otherwise ignored.
#' @return A \code{\link[ggplot2]{ggplot}} object, or (when \code{show_all =
#'   TRUE}) a named list of ggplot objects.
#' @seealso \code{\link{lort_fit}}, \code{\link{plot_cta_tree}},
#'   \code{\link{plot_cta_family}}, \code{\link[ggplot2]{ggsave}}
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   X <- data.frame(
#'     A = c(rep(0L,20), rep(1L,20), rep(1L,20)),
#'     B = c(rep(0L,20), rep(0L,20), rep(1L,20))
#'   )
#'   y <- c(rep(0L,40), rep(1L,20))
#'   lort <- lort_fit(X, y, mc_iter=100L, mc_seed=42L, loo="off", min_n=5L)
#'   p <- plot_lort_tree(lort, index=1L)
#'   print(p)
#' }
#' }
#' @export
plot_lort_tree <- function(x,
                            index              = 1L,
                            show_all           = FALSE,
                            show_path          = FALSE,
                            target_class       = 1L,
                            color_by           = c("none", "target_rate", "prediction"),
                            label_detail       = c("simple", "full"),
                            show_node_ess      = FALSE,
                            show_p             = TRUE,
                            show_loo           = TRUE,
                            main               = NULL,
                            subtitle           = NULL,
                            show_rule          = TRUE,
                            show_metrics       = FALSE,
                            short_edge_labels  = TRUE,
                            node_text_size     = 3.5,
                            edge_text_size     = 3.2,
                            palette            = NULL,
                            ...) {
  .require_ggplot2()
  color_by     <- match.arg(color_by)
  label_detail <- match.arg(label_detail)

  # Require a cta_ort object for indexed inspection
  if (!inherits(x, "cta_ort"))
    stop("plot_lort_tree: 'x' must be a cta_ort object.", call. = FALSE)

  ort_nodes <- x$ort_nodes
  if (is.null(ort_nodes) || length(ort_nodes) == 0L)
    return(.gg_message_panel("No LORT nodes to display.",
                              title = main %||% "LORT", subtitle = subtitle))

  node_keys <- names(ort_nodes)

  # When show_path = TRUE, delegate to plot_lort_path()
  if (isTRUE(show_path))
    return(plot_lort_path(x, index = index,
                          target_class      = target_class,
                          color_by          = color_by,
                          label_detail      = label_detail,
                          show_node_ess     = show_node_ess,
                          show_rule         = show_rule,
                          show_metrics      = show_metrics,
                          short_edge_labels = short_edge_labels,
                          node_text_size    = node_text_size,
                          edge_text_size    = edge_text_size,
                          palette           = palette))

  # Helper: render full local CTA at key k
  .render_one <- function(k, title_override = NULL) {
    nd  <- ort_nodes[[k]]
    if (is.null(nd) || is.null(nd$model))
      return(.gg_message_panel(
        sprintf("LORT index %s: no model (stop_reason=%s).",
                k, nd$stop_reason %||% "?"),
        title    = title_override %||% sprintf("LORT index %s", k),
        subtitle = subtitle
      ))

    auto_title <- title_override %||% main %||%
      sprintf("LORT index %s -- full local CTA (n=%d)", k,
              nd$model$n %||% NA_integer_)

    plot_cta_tree(
      nd$model,
      target_class      = target_class,
      color_by          = color_by,
      label_detail      = label_detail,
      show_node_ess     = show_node_ess,
      show_p            = show_p,
      show_loo          = show_loo,
      main              = auto_title,
      subtitle          = subtitle,
      show_rule         = show_rule,
      show_metrics      = show_metrics,
      short_edge_labels = short_edge_labels,
      node_text_size    = node_text_size,
      edge_text_size    = edge_text_size,
      palette           = palette
    )
  }

  if (isTRUE(show_all)) {
    plots <- setNames(
      lapply(node_keys, .render_one),
      paste0("index_", node_keys)
    )
    return(plots)
  }

  # Single indexed sub-tree: full local CTA model
  k <- as.character(index)
  if (!k %in% node_keys)
    stop(sprintf("plot_lort_tree: index '%s' not found in LORT nodes (%s).",
                 k, paste(node_keys, collapse = ", ")),
         call. = FALSE)

  .render_one(k)
}

# ---- plot_lort_path ---------------------------------------------------------- #

#' Plot the full local CTA models along a LORT recursion path
#'
#' Returns a named list of ggplot objects, one per LORT node on the path from
#' the root to the requested \code{index}.  Each panel shows the \emph{full}
#' local CTA model embedded at that LORT node -- not a stump summary.
#'
#' The list is named \code{index_1}, \code{index_2}, etc. (one name per LORT
#' node on the path).  Terminal nodes with no model get a message panel.
#'
#' @param x A \code{cta_ort} object from \code{\link{lort_fit}}.
#' @param index Integer; target LORT node index (end of path).
#' @param layout Character; \code{"multipanel"} (default) returns a single
#'   arranged figure with one panel per path node -- requires the
#'   \pkg{patchwork} package.  \code{"list"} returns the current named list of
#'   ggplot objects.
#' @param ncol Integer; number of columns in the multipanel layout.  Default
#'   \code{1L} (vertical stack).
#' @param target_class Integer; target class for node coloring.  Default
#'   \code{1L}.
#' @param color_by Character; leaf fill mode.  Default \code{"none"}.
#' @param label_detail Character; \code{"simple"} (default) or \code{"full"}.
#' @param show_node_ess Logical.  Default \code{FALSE}.
#' @param show_p Logical; append \code{MC p = X.XXX} to split-node labels.
#'   Default \code{TRUE}.
#' @param show_loo Logical; append LOO status/p to split-node labels.
#'   Default \code{TRUE}.
#' @param show_rule Logical.  Default \code{TRUE}.
#' @param show_metrics Logical.  Default \code{FALSE}.
#' @param short_edge_labels Logical.  Default \code{TRUE}.
#' @param node_text_size Numeric.  Default \code{3.5}.
#' @param edge_text_size Numeric.  Default \code{3.2}.
#' @param palette Named list; color overrides.
#' @param ... Ignored; reserved.
#' @return With \code{layout = "multipanel"}: a single \code{patchwork}/ggplot
#'   object containing all path panels.  With \code{layout = "list"}: a named
#'   list of ggplot objects.
#' @seealso \code{\link{lort_index_path}}, \code{\link{lort_local_tree}},
#'   \code{\link{lort_path_table}}, \code{\link{plot_lort_tree}}
#' @export
plot_lort_path <- function(x,
                            index             = 1L,
                            layout            = c("multipanel", "list"),
                            ncol              = 1L,
                            target_class      = 1L,
                            color_by          = c("none", "target_rate", "prediction"),
                            label_detail      = c("simple", "full"),
                            show_node_ess     = FALSE,
                            show_p            = TRUE,
                            show_loo          = TRUE,
                            show_rule         = TRUE,
                            show_metrics      = FALSE,
                            short_edge_labels = TRUE,
                            node_text_size    = 3.5,
                            edge_text_size    = 3.2,
                            palette           = NULL,
                            ...) {
  .require_ggplot2()
  layout       <- match.arg(layout)
  color_by     <- match.arg(color_by)
  label_detail <- match.arg(label_detail)

  if (!inherits(x, "cta_ort"))
    stop("plot_lort_path: 'x' must be a cta_ort object.", call. = FALSE)

  path      <- lort_index_path(x, index)
  ort_nodes <- x$ort_nodes %||% list()

  plots <- setNames(vector("list", nrow(path)),
                    paste0("index_", path$lort_index))

  for (i in seq_len(nrow(path))) {
    row <- path[i, ]
    nid <- row$lort_index
    nd  <- ort_nodes[[as.character(nid)]]

    # Panel title: show recursion position + key stats
    if (nid == 1L) {
      path_desc <- "root"
    } else {
      path_desc <- sprintf("via ep%d: %s",
                           row$incoming_endpoint_id %||% 0L,
                           row$incoming_path_condition %||% "?")
    }

    ess_str <- if (!is.na(row$local_ess)) sprintf(", ESS=%.1f%%", row$local_ess) else ""
    d_str   <- if (!is.na(row$local_d))   sprintf(", D=%.3f",     row$local_d)   else ""

    auto_title <- sprintf("LORT index %d | n=%d | %s%s%s",
                          nid, row$n, path_desc, ess_str, d_str)

    is_no_tree <- is.null(nd) || is.null(nd$model) || isTRUE(nd$model$no_tree)
    if (is_no_tree) {
      # No local CTA tree: forced-terminal (null model) or model returned no_tree
      stop_rsn <- row$stop_reason %||% "no valid split found"
      plots[[i]] <- .gg_result_card(
        header       = "No Local CTA Tree",
        body         = sprintf("Stop reason: %s", stop_rsn),
        footer       = sprintf("LORT index %d | n=%d | %s", nid, row$n, path_desc),
        title        = auto_title,
        header_color = "#7a3a1a"
      )
    } else {
      plots[[i]] <- plot_cta_tree(
        nd$model,
        target_class      = target_class,
        color_by          = color_by,
        label_detail      = label_detail,
        show_node_ess     = show_node_ess,
        show_p            = show_p,
        show_loo          = show_loo,
        main              = auto_title,
        show_rule         = show_rule,
        show_metrics      = show_metrics,
        short_edge_labels = short_edge_labels,
        node_text_size    = node_text_size,
        edge_text_size    = edge_text_size,
        palette           = palette
      )
    }
  }

  if (layout == "list")
    return(plots)

  # layout == "multipanel": arrange into a single figure via patchwork
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("plot_lort_path(layout = 'multipanel') requires the 'patchwork' package.\n",
         "Install it with install.packages('patchwork'), ",
         "or use layout = 'list' for the individual panels.",
         call. = FALSE)

  patchwork::wrap_plots(plots, ncol = as.integer(ncol))
}

# ---- plot_cta_family -------------------------------------------------------- #

#' Plot a CTA descendant family member using ggplot2
#'
#' Renders a publication-quality CTA tree diagram for a single member of a
#' \code{cta_family} object (indexed inspection), or a named list of plots for
#' all members (\code{show_all = TRUE}).  Requires the \pkg{ggplot2} package.
#'
#' @param family A \code{cta_family} object from
#'   \code{\link{cta_descendant_family}}.
#' @param index Integer or \code{"min_d"}; which family member to render.
#'   Default \code{1L}.  Use \code{"min_d"} to render the member with minimum
#'   D-statistic.  Ignored when \code{show_all = TRUE}.
#' @param min_d Logical; convenience shorthand for \code{index = "min_d"}.
#'   When \code{TRUE}, renders the minimum-D family member regardless of
#'   \code{index}.  Default \code{FALSE}.
#' @param show_all Logical; if \code{TRUE}, render all family members.
#'   Output format is controlled by \code{layout}.  Default \code{FALSE}.
#' @param layout Character; \code{"multipanel"} (default, requires
#'   \pkg{patchwork}) returns a single combined figure; \code{"list"} returns
#'   a named list of ggplot objects.  Only relevant when \code{show_all = TRUE}.
#' @param ncol Integer; number of columns in the multipanel grid.  Default
#'   \code{1L}.  Only used when \code{show_all = TRUE} and \code{layout =
#'   "multipanel"}.
#' @param target_class Integer; target class for endpoint coloring (default
#'   \code{1L}).
#' @param color_by Character; leaf-node fill.  \code{"none"} (default),
#'   \code{"target_rate"}, or \code{"prediction"}.
#' @param label_detail Character; \code{"simple"} (default) or \code{"full"}.
#' @param show_node_ess Logical; append node ESS to split labels.
#'   Default \code{FALSE}.
#' @param show_p Logical; append \code{MC p = X.XXX} to split-node labels.
#'   Default \code{TRUE}.
#' @param show_loo Logical; append LOO status/p to split-node labels.
#'   Default \code{TRUE}.
#' @param main Character; plot title.  Default: auto-generated with MINDENOM
#'   and D.
#' @param subtitle Character; plot subtitle.
#' @param show_rule Logical; show edge condition labels.  Default \code{TRUE}.
#' @param show_metrics Logical; append ESS/D to subtitle.  Default \code{FALSE}.
#' @param short_edge_labels Logical; strip attribute prefix from edge labels.
#'   Default \code{TRUE}.
#' @param node_text_size Numeric; text size for node labels.  Default \code{3.5}.
#' @param edge_text_size Numeric; text size for edge labels.  Default \code{3.2}.
#' @param palette Named list for color overrides.
#' @return A \code{\link[ggplot2]{ggplot}} object (single member or multipanel),
#'   or (when \code{show_all = TRUE} and \code{layout = "list"}) a named list
#'   of ggplot objects.
#' @seealso \code{\link{cta_descendant_family}}, \code{\link{plot_cta_tree}},
#'   \code{\link{plot_lort_tree}}, \code{\link[ggplot2]{ggsave}}
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   X <- data.frame(x1 = c(rep(0L,20), rep(1L,20)),
#'                   x2 = c(rep(0L,10), rep(1L,10), rep(0L,10), rep(1L,10)))
#'   y <- c(rep(0L,30), rep(1L,10))
#'   fam <- cta_descendant_family(X, y, mc_iter=200L, mc_seed=42L, loo="off")
#'   p   <- plot_cta_family(fam, index=1L)
#'   print(p)
#' }
#' }
#' @export
plot_cta_family <- function(family,
                             index              = 1L,
                             min_d              = FALSE,
                             show_all           = FALSE,
                             layout             = c("multipanel", "list"),
                             ncol               = 1L,
                             target_class       = 1L,
                             color_by           = c("none", "target_rate", "prediction"),
                             label_detail       = c("simple", "full"),
                             show_node_ess      = FALSE,
                             show_p             = TRUE,
                             show_loo           = TRUE,
                             main               = NULL,
                             subtitle           = NULL,
                             show_rule          = TRUE,
                             show_metrics       = FALSE,
                             short_edge_labels  = TRUE,
                             node_text_size     = 3.5,
                             edge_text_size     = 3.2,
                             palette            = NULL) {
  .require_ggplot2()
  if (!inherits(family, "cta_family"))
    stop("plot_cta_family: 'family' must be a cta_family object.", call. = FALSE)
  layout       <- match.arg(layout)
  color_by     <- match.arg(color_by)
  label_detail <- match.arg(label_detail)
  # min_d = TRUE is a convenience shorthand for index = "min_d"
  if (isTRUE(min_d)) index <- "min_d"

  members <- family$members
  if (length(members) == 0L)
    return(.gg_message_panel("No family members to display.",
                              title = main %||% "CTA Family", subtitle = subtitle))

  # Helper: render one family member
  .render_member <- function(k, title_override = NULL) {
    mem <- members[[k]]
    if (is.null(mem) || is.null(mem$tree))
      return(.gg_message_panel(
        sprintf("Family member %d: no tree.", k),
        title = title_override %||% sprintf("CTA Family -- MINDENOM=%d", k),
        subtitle = subtitle
      ))

    md_val <- mem$mindenom %||% k
    d_val  <- mem$d        %||% NA_real_
    no_tr  <- isTRUE(mem$no_tree)

    auto_title <- title_override %||% main %||% {
      lbl <- sprintf("CTA Family -- MINDENOM=%d", as.integer(md_val))
      ess_val   <- mem$overall_ess %||% NA_real_
      ess_label <- if (isTRUE(mem$has_weights)) "WESS" else "ESS"
      if (!no_tr && !is.na(ess_val) && is.finite(ess_val))
        lbl <- sprintf("%s, %s=%.2f%%", lbl, ess_label, ess_val)
      if (!no_tr && !is.na(d_val) && is.finite(d_val))
        lbl <- sprintf("%s, D=%.4f", lbl, d_val)
      if (no_tr) lbl <- paste0(lbl, " [no tree]")
      lbl
    }

    plot_cta_tree(
      mem$tree,
      target_class      = target_class,
      color_by          = color_by,
      label_detail      = label_detail,
      show_node_ess     = show_node_ess,
      show_p            = show_p,
      show_loo          = show_loo,
      main              = auto_title,
      subtitle          = subtitle,
      show_rule         = show_rule,
      show_metrics      = show_metrics,
      short_edge_labels = short_edge_labels,
      node_text_size    = node_text_size,
      edge_text_size    = edge_text_size,
      palette           = palette
    )
  }

  if (isTRUE(show_all)) {
    plots <- setNames(
      lapply(seq_along(members), .render_member),
      paste0("member_", seq_along(members))
    )
    if (identical(layout, "list")) return(plots)
    # layout == "multipanel": combine via patchwork
    if (!requireNamespace("patchwork", quietly = TRUE))
      stop("plot_cta_family(layout = 'multipanel') requires the 'patchwork' package.\n",
           "  install.packages(\"patchwork\")", call. = FALSE)
    return(patchwork::wrap_plots(plots, ncol = as.integer(ncol)))
  }

  # Resolve index
  k <- if (identical(index, "min_d")) {
    mid <- family$min_d_idx %||% NA_integer_
    if (is.na(mid))
      stop("plot_cta_family: no feasible min-D member (all no-tree).", call. = FALSE)
    as.integer(mid)
  } else {
    as.integer(index)
  }

  if (k < 1L || k > length(members))
    stop(sprintf("plot_cta_family: index %d out of range (1-%d).",
                 k, length(members)), call. = FALSE)

  .render_member(k)
}

###############################################################################
# v3C2  -  Balance renderers
#
# Public API:
#   plot_oda_balance()     -  ODA ESS/WESS covariate balance dot plot
#   plot_smd_balance()     -  SMD absolute-value covariate balance dot plot
#   plot_balance_love()    -  Love-plot wrapper (calls plot_smd_balance)
#   plot_cta_balance()     -  CTA multivariate balance: tree or message panel
#
# Rules:
#   - Pure renderers: no fitting, no recomputation.  Read pre-computed
#     plot-data objects only.
#   - plot_oda_balance renders only what is in oda_balance_plot_data$rows.
#     If abs_smd is absent from the plot-data it is not plotted.
#   - Accepts the parent table type as a convenience and calls the
#     corresponding plot_data() transform; never calls the fit function.
###############################################################################

# ---- Internal balance helpers ------------------------------------------------ #

# Coerce to oda_balance_plot_data.
# Accepts oda_balance_plot_data directly, or oda_balance_table (calls
# oda_balance_plot_data -- the pure transform, NOT oda_balance_table).
.coerce_oda_balance_pd <- function(x, p_col, rank_by) {
  if (inherits(x, "oda_balance_plot_data")) return(x)
  if (inherits(x, "oda_balance_table"))
    return(oda_balance_plot_data(x, p_col = p_col, rank_by = rank_by))
  stop(
    "plot_oda_balance: 'x' must be an oda_balance_plot_data or oda_balance_table object.",
    call. = FALSE
  )
}

# Coerce to cta_balance_plot_data.
.coerce_cta_balance_pd <- function(x, target_class, digits) {
  if (inherits(x, "cta_balance_plot_data")) return(x)
  if (inherits(x, "cta_balance_table"))
    return(cta_balance_plot_data(x, target_class = as.integer(target_class),
                                  digits = as.integer(digits)))
  stop(
    "plot_cta_balance: 'x' must be a cta_balance_plot_data or cta_balance_table object.",
    call. = FALSE
  )
}

# Balance color palette.
.balance_pal <- function(palette = NULL) {
  list(
    imbalanced   = palette[["imbalanced"]]   %||% "#c0392b",
    balanced     = palette[["balanced"]]     %||% "#2980b9",
    unclassified = palette[["unclassified"]] %||% "#95a5a6"
  )
}

# Shared theme finishing for balance plots.
.balance_theme <- function(p, base_theme, main, subtitle) {
  p +
    base_theme +
    ggplot2::theme(
      legend.position    = "right",
      plot.title         = ggplot2::element_text(face = "bold", size = 12,
                                                  hjust = 0.5),
      plot.subtitle      = ggplot2::element_text(size = 10, hjust = 0.5,
                                                  color = "#555555"),
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = main, subtitle = subtitle)
}

# ---- plot_oda_balance -------------------------------------------------------- #

#' Plot ODA covariate balance
#'
#' Renders a horizontal dot-plot of ODA-based covariate balance diagnostics.
#' Each covariate is shown as a point; the x-axis is ESS or WESS (0-100 \%),
#' and point color reflects significance status.  The function is a pure
#' renderer: it does not fit any ODA models and does not accept \code{group}
#' or \code{X} arguments.  If \code{abs_smd} is absent from the plot-data it
#' is not plotted.
#'
#' @param x An \code{"oda_balance_plot_data"} object from
#'   \code{\link{oda_balance_plot_data}}, or an \code{"oda_balance_table"}
#'   object from \code{\link{oda_balance_table}} (coerced internally via
#'   \code{\link{oda_balance_plot_data}}; never calls the fitting function).
#' @param p_col Character; which p-value column drives significance colour when
#'   coercing from an \code{oda_balance_table}.  One of \code{"p_mc"}
#'   (default), \code{"p_sidak"}, \code{"p_bonferroni"}.  Ignored when
#'   \code{x} is already \code{oda_balance_plot_data}.
#' @param rank_by Character; sort order when coercing from
#'   \code{oda_balance_table}: \code{"abs_ess"} (default), \code{"p"},
#'   \code{"abs_smd"}.
#' @param main Character; plot title.  Default: auto-generated summary.
#' @param subtitle Character; plot subtitle.
#' @param show_significance Logical; annotate significantly imbalanced
#'   covariates with a \code{"*"} label.  Default \code{TRUE}.
#' @param palette Named list for color overrides: \code{imbalanced},
#'   \code{balanced}, \code{unclassified}.
#' @param theme Character; \code{"clean"} (default, \code{theme_bw} base) or
#'   \code{"minimal"}.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @seealso \code{\link{oda_balance_plot_data}}, \code{\link{oda_balance_table}}
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   group <- c(rep(0L, 20), rep(1L, 20))
#'   X     <- data.frame(A = c(rep(0L,20), rep(1L,20)),
#'                        B = rnorm(40))
#'   bt  <- oda_balance_table(group, X, mcarlo = FALSE, mc_iter = 100L)
#'   pd  <- oda_balance_plot_data(bt)
#'   p   <- plot_oda_balance(pd)
#'   print(p)
#' }
#' }
#' @export
plot_oda_balance <- function(x,
                              p_col             = "p_mc",
                              rank_by           = "abs_ess",
                              main              = NULL,
                              subtitle          = NULL,
                              show_significance = TRUE,
                              palette           = NULL,
                              theme             = c("clean", "minimal")) {
  .require_ggplot2()
  theme <- match.arg(theme)

  pd   <- .coerce_oda_balance_pd(x, p_col = p_col, rank_by = rank_by)
  rows <- pd$rows

  if (nrow(rows) == 0L)
    return(.gg_message_panel("No covariates to display.",
                              title    = main %||% "ODA Covariate Balance",
                              subtitle = subtitle))

  pal <- .balance_pal(palette)

  # Sort by rank (1 = most imbalanced); reverse levels so top of Y = rank 1
  rows <- rows[order(rows$rank), ]
  rows$attr_fac  <- factor(rows$attribute, levels = rev(rows$attribute))
  rows$sig_color <- ifelse(is.na(rows$significant), "unknown",
                    ifelse(rows$significant, "imbalanced", "balanced"))

  ess_xlab   <- paste0(pd$ess_label %||% "ESS", " (%)")
  n_sig      <- pd$n_significant %||% 0L
  n_cov      <- pd$n_covariates  %||% nrow(rows)
  main_title <- main %||% sprintf(
    "ODA Covariate Balance  (%d of %d covariates significant)",
    n_sig, n_cov
  )

  base_theme <- if (identical(theme, "minimal")) ggplot2::theme_minimal()
                else ggplot2::theme_bw()

  p <- ggplot2::ggplot(rows,
         ggplot2::aes(x = ess_display_bar, y = attr_fac)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = ess_display_bar,
                   y = attr_fac, yend = attr_fac),
      color = "#cccccc", linewidth = 0.6, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = sig_color),
      size = 3.5, na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(
      values   = c(imbalanced   = pal$imbalanced,
                   balanced     = pal$balanced,
                   unknown      = pal$unclassified),
      labels   = c(imbalanced   = "Significant",
                   balanced     = "Not significant",
                   unknown      = "p not computed"),
      name     = NULL,
      na.value = pal$unclassified,
      drop     = FALSE
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 100),
      expand = ggplot2::expansion(mult = c(0, 0.04))
    ) +
    ggplot2::labs(x = ess_xlab, y = NULL)

  # Significance labels ("*") next to imbalanced points  -  pure render of
  # whatever significance_label is in the plot-data; no recomputation.
  if (isTRUE(show_significance) && "significance_label" %in% names(rows)) {
    sig_pts <- rows[nchar(as.character(rows$significance_label)) > 0L, ,
                    drop = FALSE]
    if (nrow(sig_pts) > 0L) {
      p <- p + ggplot2::geom_text(
        data = sig_pts,
        ggplot2::aes(x = ess_display_bar, y = attr_fac,
                     label = significance_label),
        hjust = -0.7, vjust = 0.4, size = 4.5, color = pal$imbalanced
      )
    }
  }

  .balance_theme(p, base_theme, main_title, subtitle)
}

# ---- plot_smd_balance -------------------------------------------------------- #

#' Plot SMD covariate balance
#'
#' Renders a horizontal dot-plot of absolute standardized mean differences
#' (|SMD|) for each covariate.  Vertical reference lines at 0.10 (and
#' optionally 0.20) mark conventional balance thresholds.  Points are colored
#' by whether |SMD| < 0.10.
#'
#' @param x A \code{"smd_balance_table"} object from
#'   \code{\link{smd_balance_table}}.
#' @param ref_010 Logical; draw a dashed reference line at |SMD| = 0.10.
#'   Default \code{TRUE}.
#' @param ref_020 Logical; draw a dotted reference line at |SMD| = 0.20.
#'   Default \code{FALSE}.
#' @param main Character; plot title.  Default \code{"SMD Covariate Balance"}.
#' @param subtitle Character; plot subtitle.
#' @param palette Named list for color overrides: \code{imbalanced},
#'   \code{balanced}, \code{unclassified}.
#' @param theme Character; \code{"clean"} (default) or \code{"minimal"}.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @seealso \code{\link{smd_balance_table}}, \code{\link{plot_balance_love}}
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   group <- c(rep(0L, 20), rep(1L, 20))
#'   X     <- data.frame(A = c(rep(0L,20), rep(1L,20)),
#'                        B = rnorm(40))
#'   smd <- smd_balance_table(group, X)
#'   p   <- plot_smd_balance(smd)
#'   print(p)
#' }
#' }
#' @export
plot_smd_balance <- function(x,
                              ref_010  = TRUE,
                              ref_020  = FALSE,
                              main     = NULL,
                              subtitle = NULL,
                              palette  = NULL,
                              theme    = c("clean", "minimal")) {
  .require_ggplot2()
  if (!inherits(x, "smd_balance_table"))
    stop("plot_smd_balance: 'x' must be a smd_balance_table object.",
         call. = FALSE)
  theme <- match.arg(theme)
  pal   <- .balance_pal(palette)
  tbl   <- as.data.frame(x)

  ok <- !is.na(tbl$abs_smd)
  if (!any(ok))
    return(.gg_message_panel("No abs_smd values available to plot.",
                              title    = main %||% "SMD Covariate Balance",
                              subtitle = subtitle))

  # Sort by abs_smd descending; most imbalanced at top of Y axis
  tbl <- tbl[order(tbl$abs_smd, decreasing = TRUE, na.last = TRUE), ]
  tbl$attr_fac <- factor(tbl$attribute, levels = rev(tbl$attribute))
  tbl$bal_col  <- ifelse(is.na(tbl$balanced_010), "unknown",
                  ifelse(tbl$balanced_010, "balanced", "imbalanced"))

  x_max      <- max(c(tbl$abs_smd, 0.22), na.rm = TRUE)
  base_theme <- if (identical(theme, "minimal")) ggplot2::theme_minimal()
                else ggplot2::theme_bw()

  # Number of attributes for y-axis scale expansion
  n_attr  <- nrow(tbl[ok, ])
  y_top   <- if (n_attr > 0L) as.numeric(n_attr) + 0.6 else 1.5

  p <- ggplot2::ggplot(tbl,
         ggplot2::aes(x = abs_smd, y = attr_fac)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = abs_smd,
                   y = attr_fac, yend = attr_fac),
      color = "#cccccc", linewidth = 0.7, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = bal_col),
      size = 4.5, na.rm = TRUE      # larger points for journal print
    ) +
    ggplot2::scale_color_manual(
      values   = c(imbalanced   = pal$imbalanced,
                   balanced     = pal$balanced,
                   unknown      = pal$unclassified),
      labels   = c(imbalanced   = "|SMD| \u2265 0.10",
                   balanced     = "|SMD| < 0.10",
                   unknown      = "Unknown"),
      name     = NULL,
      na.value = pal$unclassified,
      drop     = FALSE
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, x_max * 1.10),   # extra space for inline threshold label
      expand = ggplot2::expansion(mult = c(0, 0.04))
    ) +
    ggplot2::labs(x = "|SMD|", y = NULL)

  # Reference lines with inline text annotations (not just legend)
  if (isTRUE(ref_010)) {
    p <- p +
      ggplot2::geom_vline(xintercept = 0.10, linetype = "dashed",
                           color = "#555555", linewidth = 0.6) +
      ggplot2::annotate("text", x = 0.10, y = y_top,
                         label = "|SMD| = 0.10", hjust = -0.08, vjust = 1,
                         size = 3.0, color = "#555555", fontface = "italic")
  }
  if (isTRUE(ref_020)) {
    p <- p +
      ggplot2::geom_vline(xintercept = 0.20, linetype = "dotted",
                           color = "#888888", linewidth = 0.5) +
      ggplot2::annotate("text", x = 0.20, y = y_top,
                         label = "|SMD| = 0.20", hjust = -0.08, vjust = 1,
                         size = 3.0, color = "#888888", fontface = "italic")
  }

  # Apply shared balance theme first, then override legend to "bottom" so it
  # does not compete with the x-axis annotation of the threshold lines.
  # Order matters: .balance_theme() sets legend.position = "right"; the +
  # theme() below wins because ggplot2 themes are last-write-wins.
  .balance_theme(p, base_theme, main %||% "SMD Covariate Balance", subtitle) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.key.size  = ggplot2::unit(0.5, "lines"),
      legend.text      = ggplot2::element_text(size = 9),
      axis.text.y      = ggplot2::element_text(size = 10, face = "plain"),
      axis.title.x     = ggplot2::element_text(size = 9, color = "#555555")
    )
}

# ---- plot_balance_love ------------------------------------------------------- #

#' Love plot for covariate balance (SMD)
#'
#' A direct alias for \code{\link{plot_smd_balance}}.  Produces a
#' Cleveland-style Love plot of absolute SMD with conventional threshold
#' reference lines.
#'
#' @param x A \code{"smd_balance_table"} object from
#'   \code{\link{smd_balance_table}}.
#' @param ... Arguments forwarded to \code{\link{plot_smd_balance}}.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @seealso \code{\link{plot_smd_balance}}, \code{\link{smd_balance_table}}
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   group <- c(rep(0L, 20), rep(1L, 20))
#'   X     <- data.frame(A = c(rep(0L,20), rep(1L,20)),
#'                        B = rnorm(40))
#'   smd <- smd_balance_table(group, X)
#'   p   <- plot_balance_love(smd)
#'   print(p)
#' }
#' }
#' @export
plot_balance_love <- function(x, ...) {
  plot_smd_balance(x, ...)
}

# ---- plot_cta_balance -------------------------------------------------------- #

#' Plot CTA multivariate covariate balance
#'
#' Renders the CTA covariate balance result.  When no discriminating tree was
#' found (\code{status = "no_tree"}), a message panel confirms favorable
#' evidence of multivariable balance under the declared constraints.  When a
#' valid tree or stump was found, the tree diagram is rendered via
#' \code{\link{plot_cta_tree}}.
#'
#' This function is a pure renderer.  It does not fit any CTA models and does
#' not accept \code{group} or \code{X} arguments.
#'
#' @param x A \code{"cta_balance_plot_data"} object from
#'   \code{\link{cta_balance_plot_data}}, or a \code{"cta_balance_table"}
#'   object from \code{\link{cta_balance_table}} (coerced internally via
#'   \code{\link{cta_balance_plot_data}}; never calls the fitting function).
#' @param target_class Integer; target class for leaf-node coloring.
#'   Default \code{1L}.
#' @param color_by Character; leaf-node fill: \code{"target_rate"} (default),
#'   \code{"prediction"}, \code{"none"}.
#' @param main Character; plot title.  Default: auto-generated from ESS/WESS.
#' @param subtitle Character; plot subtitle.
#' @param ... Additional arguments forwarded to \code{\link{plot_cta_tree}}
#'   when a tree is rendered.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @seealso \code{\link{cta_balance_plot_data}}, \code{\link{cta_balance_table}},
#'   \code{\link{plot_cta_tree}}
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   X <- data.frame(
#'     A = c(rep(0L,20), rep(1L,20), rep(1L,20)),
#'     B = c(rep(0L,20), rep(0L,20), rep(1L,20))
#'   )
#'   group <- c(rep(0L, 40), rep(1L, 20))
#'   ct  <- cta_balance_table(group, X, mindenom = 5L,
#'                             mc_iter = 200L, mc_seed = 42L)
#'   cpd <- cta_balance_plot_data(ct)
#'   p   <- plot_cta_balance(cpd)
#'   print(p)
#' }
#' }
#' @export
plot_cta_balance <- function(x,
                              target_class = 1L,
                              color_by     = c("target_rate", "prediction", "none"),
                              main         = NULL,
                              subtitle     = NULL,
                              ...) {
  .require_ggplot2()
  color_by <- match.arg(color_by)

  pd     <- .coerce_cta_balance_pd(x, target_class = target_class, digits = 1L)
  status <- pd$status

  # no_tree: favorable balance -- styled result card (positive outcome, not error)
  if (identical(status, "no_tree")) {
    raw_msg  <- pd$no_tree_message %||% ""
    # Extract MINDENOM / alpha from message if present, else generic footer
    footer   <- if (nchar(raw_msg) > 0L) raw_msg else NULL
    return(.gg_result_card(
      header   = "No Discriminating Tree Found",
      body     = "No combination of covariates predicted group membership\nabove the declared significance threshold.\n\nThis is favorable evidence of multivariable covariate balance.",
      footer   = footer,
      title    = main %||% "CTA Covariate Balance",
      subtitle = subtitle
    ))
  }

  # fit_error: error message
  if (identical(status, "fit_error")) {
    msg <- pd$no_tree_message %||% "CTA fitting error."
    return(.gg_message_panel(msg,
                              title    = main %||% "CTA Covariate Balance \u2014 Error",
                              subtitle = subtitle))
  }

  # Tree/stump: render via plot_cta_tree (pure pass-through to existing renderer)
  if (is.null(pd$cta_pd)) {
    return(.gg_message_panel("CTA plot data unavailable.",
                              title    = main %||% "CTA Covariate Balance",
                              subtitle = subtitle))
  }

  ess_lbl    <- pd$ess_label %||% "ESS"
  ess_val    <- pd$ess_display
  main_title <- main %||% if (!is.na(ess_val %||% NA_real_)) {
    sprintf("CTA Covariate Balance  \u2014  %s = %.1f%%", ess_lbl, ess_val)
  } else {
    "CTA Covariate Balance"
  }

  plot_cta_tree(pd$cta_pd,
                target_class = as.integer(target_class),
                color_by     = color_by,
                main         = main_title,
                subtitle     = subtitle,
                ...)
}

###############################################################################
# v3C3  -  Evidence-interval balance renderers
#
# plot_oda_balance_effects()  -  forest plot for oda_balance_effect_table
# plot_cta_balance_effects()  -  evidence card for cta_balance_effect_summary
#
# Both are pure renderers.  All bootstrap/chance computation lives in the
# builder layer (balance.R).  These functions accept the pre-computed objects
# only; they do NOT accept group/X/w arguments.
###############################################################################

# ---- plot_oda_balance_effects ----------------------------------------------- #

#' Forest plot of ODA covariate balance evidence intervals
#'
#' Renders a forest plot from an \code{\link{oda_balance_effect_table}} object.
#' Each covariate is displayed as one row.  A thin gray segment shows the
#' chance (null) confidence interval; a thick black segment shows the
#' bootstrap model CI; a point shows the observed ESS/WESS.  A vertical
#' dashed line marks the chance upper bound (chance_hi) as a visual reference.
#'
#' When the object contains multiple analysis scales (e.g.,
#' \code{compare_weights = TRUE}), the plot is faceted by \code{analysis}.
#'
#' \strong{This function does not fit any models.}  Pass a pre-computed
#' \code{oda_balance_effect_table} from \code{\link{oda_balance_effect_table}}.
#'
#' @param x An \code{"oda_balance_effect_table"} object.
#' @param main Optional character; plot title.  Defaults to
#'   \code{"ODA Covariate Balance -- Evidence Intervals"}.
#' @param subtitle Optional character; plot subtitle.
#' @param x_label Optional character; x-axis label.  Defaults to the metric
#'   label from the data (\code{"ESS (\%)"} or \code{"WESS (\%)"}).
#' @param xlim Optional numeric(2); x-axis limits.  Auto-computed when
#'   \code{NULL}.
#' @param ... Ignored; reserved for future use.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{oda_balance_effect_table}}
#' @examples
#' \donttest{
#' group <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
#' X     <- data.frame(v1 = c(1, 2, 3, 4, 5, 6, 7, 8),
#'                     v2 = c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L))
#' et <- oda_balance_effect_table(group, X,
#'                                 nboot = 50L, chance_iter = 50L,
#'                                 mc_iter = 200L, mc_seed = 1L)
#' plot_oda_balance_effects(et)
#' }
#' @export
plot_oda_balance_effects <- function(x,
                                      main     = NULL,
                                      subtitle = NULL,
                                      x_label  = NULL,
                                      xlim     = NULL,
                                      ...) {
  .require_ggplot2()

  if (!inherits(x, "oda_balance_effect_table"))
    stop("x must be an 'oda_balance_effect_table' object from oda_balance_effect_table().",
         call. = FALSE)

  rows <- x$rows
  meta <- x$meta

  if (nrow(rows) == 0L)
    return(.gg_message_panel("No covariate rows to display.",
                              title = main %||% "ODA Covariate Balance -- Evidence Intervals"))

  # Reverse factor so first covariate appears at the top.
  rows$attribute <- factor(rows$attribute, levels = rev(unique(rows$attribute)))

  # Determine x-axis metric label.
  metrics  <- unique(rows$metric)
  met_lbl  <- x_label %||% if (length(metrics) == 1L)
                              paste0(metrics[[1L]], " (%)")
                           else "ESS/WESS (%)"

  # Auto x-axis limits: 0 to ceiling of max observed/boot_hi/chance_hi + 5.
  if (is.null(xlim)) {
    x_max <- max(c(rows$estimate, rows$boot_hi, rows$chance_hi),
                 na.rm = TRUE)
    if (!is.finite(x_max)) x_max <- 50
    xlim <- c(0, ceiling(x_max / 5) * 5 + 5)
  }

  # Mean chance_hi across the scale for the reference line.
  ref_chance_hi <- mean(rows$chance_hi, na.rm = TRUE)

  multi_scale <- length(unique(rows$analysis)) > 1L

  # Build the ggplot.
  p <- ggplot2::ggplot(rows, ggplot2::aes(y = attribute)) +
    # Chance interval: thin light gray
    ggplot2::geom_segment(
      ggplot2::aes(x = chance_lo, xend = chance_hi, yend = attribute),
      colour = "gray70", linewidth = 0.6, na.rm = TRUE
    ) +
    # Bootstrap CI: thick black
    ggplot2::geom_segment(
      ggplot2::aes(x = boot_lo, xend = boot_hi, yend = attribute),
      colour = "black", linewidth = 1.8, na.rm = TRUE
    ) +
    # Observed point
    ggplot2::geom_point(
      ggplot2::aes(x = estimate),
      shape = 21L, fill = "white", colour = "black", size = 2.5, na.rm = TRUE
    ) +
    # Reference line at chance_hi (mean across rows)
    ggplot2::geom_vline(
      xintercept = ref_chance_hi,
      linetype = "dashed", colour = "gray50", linewidth = 0.5
    ) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::scale_x_continuous(
      name   = met_lbl,
      breaks = pretty(seq(xlim[1L], xlim[2L], length.out = 6L))
    ) +
    ggplot2::scale_y_discrete(name = NULL) +
    ggplot2::labs(
      title    = main %||% "ODA Covariate Balance -- Evidence Intervals",
      subtitle = subtitle %||%
        paste0("Thick black = bootstrap CI; thin gray = chance CI; ",
               "point = observed ", if (length(metrics) == 1L) metrics[[1L]] else "ESS/WESS")
    ) +
    ggplot2::theme_bw(base_size = 11L) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.text.y        = ggplot2::element_text(size = 9L),
      plot.title         = ggplot2::element_text(size = 11L, face = "bold"),
      plot.subtitle      = ggplot2::element_text(size = 8L, colour = "gray40"),
      strip.background   = ggplot2::element_rect(fill = "gray90"),
      strip.text         = ggplot2::element_text(size = 9L, face = "bold")
    )

  if (multi_scale)
    p <- p + ggplot2::facet_wrap(~ analysis, ncol = 1L)

  p
}

# ---- plot_cta_balance_effects ----------------------------------------------- #

#' Evidence card for CTA multivariate covariate balance
#'
#' Renders an evidence-interval card from a
#' \code{\link{cta_balance_effect_summary}} object.  Each row of the card
#' corresponds to one analysis scale.  The plot uses the same interval
#' encoding as \code{\link{plot_oda_balance_effects}}: thick black = bootstrap
#' CI, thin gray = chance CI, open circle = observed ESS/WESS.
#'
#' When \code{status = "no_tree"} for all rows, a favorable-balance message
#' panel is returned instead of an interval plot.
#'
#' \strong{This function does not fit any models.}
#'
#' @param x A \code{"cta_balance_effect_summary"} object from
#'   \code{\link{cta_balance_effect_summary}}.
#' @param main Optional character; plot title.
#' @param subtitle Optional character; plot subtitle.
#' @param xlim Optional numeric(2); x-axis limits.
#' @param ... Ignored; reserved for future use.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{cta_balance_effect_summary}}
#' @examples
#' \donttest{
#' group <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
#' X     <- data.frame(v1 = c(1, 2, 3, 4, 5, 6, 7, 8),
#'                     v2 = c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L))
#' ces <- cta_balance_effect_summary(group, X, mindenom = 5L,
#'                                    mc_iter = 200L, mc_seed = 42L,
#'                                    nboot = 20L, chance_iter = 20L)
#' plot_cta_balance_effects(ces)
#' }
#' @export
plot_cta_balance_effects <- function(x,
                                      main     = NULL,
                                      subtitle = NULL,
                                      xlim     = NULL,
                                      ...) {
  .require_ggplot2()

  if (!inherits(x, "cta_balance_effect_summary"))
    stop("x must be a 'cta_balance_effect_summary' object from cta_balance_effect_summary().",
         call. = FALSE)

  rows <- x$rows
  meta <- x$meta

  if (nrow(rows) == 0L)
    return(.gg_message_panel("No rows to display.",
                              title = main %||% "CTA Covariate Balance -- Evidence"))

  # All no_tree: favorable balance result card.
  all_no_tree <- all(rows$status == "no_tree", na.rm = TRUE)
  if (all_no_tree) {
    return(.gg_result_card(
      header   = "No Discriminating Tree Found",
      body     = paste0("CTA found no combination of baseline covariates ",
                        "that predicted group membership above the declared ",
                        "significance threshold.\n\n",
                        "This is favorable evidence of multivariable covariate balance."),
      title    = main %||% "CTA Covariate Balance -- Evidence",
      subtitle = subtitle
    ))
  }

  # Factor: reverse so first analysis (usually "unweighted") is at top.
  rows$analysis <- factor(rows$analysis, levels = rev(unique(rows$analysis)))

  metrics  <- unique(rows$metric)
  met_lbl  <- paste0(if (length(metrics) == 1L) metrics[[1L]] else "ESS/WESS", " (%)")

  if (is.null(xlim)) {
    x_max <- max(c(rows$estimate, rows$boot_hi, rows$chance_hi),
                 na.rm = TRUE)
    if (!is.finite(x_max)) x_max <- 50
    xlim <- c(0, ceiling(x_max / 5) * 5 + 5)
  }

  ref_chance_hi <- mean(rows$chance_hi, na.rm = TRUE)

  # Build status annotation string for each row.
  rows$status_lbl <- vapply(seq_len(nrow(rows)), function(i) {
    st <- rows$status[i] %||% ""
    root <- rows$root_attribute[i] %||% ""
    n_ep <- rows$n_endpoints[i]
    if (identical(st, "no_tree")) return("no tree")
    if (identical(st, "fit_error")) return("fit error")
    lbl <- st
    if (!is.na(root) && nchar(root) > 0L)
      lbl <- paste0(lbl, ", root=", root)
    if (!is.na(n_ep))
      lbl <- paste0(lbl, ", ep=", n_ep)
    lbl
  }, character(1L))

  p <- ggplot2::ggplot(rows, ggplot2::aes(y = analysis)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = chance_lo, xend = chance_hi, yend = analysis),
      colour = "gray70", linewidth = 0.6, na.rm = TRUE
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = boot_lo, xend = boot_hi, yend = analysis),
      colour = "black", linewidth = 2.5, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = estimate),
      shape = 21L, fill = "white", colour = "black", size = 3.5, na.rm = TRUE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = xlim[2L], label = status_lbl),
      hjust = 1, size = 2.8, colour = "gray40", na.rm = TRUE
    ) +
    ggplot2::geom_vline(
      xintercept = ref_chance_hi,
      linetype = "dashed", colour = "gray50", linewidth = 0.5
    ) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::scale_x_continuous(
      name   = met_lbl,
      breaks = pretty(seq(xlim[1L], xlim[2L], length.out = 6L))
    ) +
    ggplot2::scale_y_discrete(name = "Analysis") +
    ggplot2::labs(
      title    = main %||% "CTA Covariate Balance -- Evidence",
      subtitle = subtitle %||%
        paste0("Thick black = bootstrap CI; thin gray = chance CI; ",
               "point = observed ",
               if (length(metrics) == 1L) metrics[[1L]] else "ESS/WESS")
    ) +
    ggplot2::theme_bw(base_size = 11L) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.text.y        = ggplot2::element_text(size = 10L, face = "bold"),
      plot.title         = ggplot2::element_text(size = 11L, face = "bold"),
      plot.subtitle      = ggplot2::element_text(size = 8L, colour = "gray40")
    )

  p
}
