###############################################################################
# R/graphics_v3.R  -  Graphics v3: direct ggplot2 tree renderers
#
# Public API (v3C1):
#   plot_cta_tree()     -  ggplot CTA tree diagram
#   plot_lort_tree()    -  ggplot LORT (recursive CTA) tree diagram
#
# v3C2 balance renderers are deferred to the next slice.
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
  # v3C2  -  balance rendering
  "attr_fac", "ess_display_bar",
  "sig_color", "significance_label",
  "abs_smd", "bal_col"
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

# Return a ggplot message panel (no_tree, empty ORT, errors).
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
  half_h <- 0.28
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

  # ---- Internal (split) nodes  -  fixed fill -----------------------------------
  if (nrow(nd_int) > 0L) {
    p <- p +
      ggplot2::geom_rect(
        data = nd_int,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
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
#'   \code{"target_rate"} (default): continuous gradient by target-class
#'   proportion; \code{"prediction"}: discrete fill by predicted class;
#'   \code{"none"}: all white.
#' @param main Character; plot title.  Default: auto-generated from tree
#'   structure.
#' @param subtitle Character; plot subtitle.
#' @param show_n Logical; include observation count in node labels.  Currently
#'   uses the labels pre-computed by \code{\link{cta_plot_data}}.
#' @param show_percent Logical; include ESS/WESS in split-node labels.
#'   Currently uses pre-computed labels.
#' @param show_rule Logical; show branch condition labels on edges.  Default
#'   \code{TRUE}.
#' @param show_metrics Logical; if \code{TRUE}, appends an ESS/WESS and D
#'   line to the plot subtitle.  Default \code{FALSE}.
#' @param wrap_width Integer; maximum characters per line for long labels.
#'   Default \code{26L}.  Applied to node labels via \code{strwrap()}.
#' @param node_text_size Numeric; ggplot text size for node labels.
#'   Default \code{3.5}.
#' @param edge_text_size Numeric; ggplot text size for edge labels.
#'   Default \code{3.2}.
#' @param node_width Numeric; node box width in data units.  \code{NULL}
#'   (default) auto-sizes based on longest label.
#' @param palette Named list for color overrides: \code{internal}, \code{low},
#'   \code{high}.  \code{NULL} uses defaults.
#' @param theme Character; ggplot theme: \code{"clean"} (default, theme_void
#'   base) or \code{"minimal"} (theme_minimal base).
#' @return A \code{\link[ggplot2]{ggplot}} object.  Print it, modify it, or
#'   save with \code{ggplot2::ggsave()}.
#' @seealso \code{\link{cta_fit}}, \code{\link{cta_plot_data}},
#'   \code{\link{plot_lort_tree}}, \code{\link[ggplot2]{ggsave}}
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   X <- data.frame(x1 = c(1,2,3,4,5,6,7,8),
#'                   x2 = c(0L,0L,1L,0L,1L,1L,0L,1L))
#'   y <- c(1L,1L,1L,1L,2L,2L,2L,2L)
#'   tree <- cta_fit(X, y, mindenom=1L, mc_iter=500L, mc_seed=42L, loo="off")
#'   p <- plot_cta_tree(tree, color_by="prediction")
#'   print(p)
#' }
#' }
#' @export
plot_cta_tree <- function(x,
                           target_class   = 1L,
                           color_by       = c("target_rate", "prediction", "none"),
                           main           = NULL,
                           subtitle       = NULL,
                           show_n         = TRUE,
                           show_percent   = TRUE,
                           show_rule      = TRUE,
                           show_metrics   = FALSE,
                           wrap_width     = 26L,
                           node_text_size = 3.5,
                           edge_text_size = 3.2,
                           node_width     = NULL,
                           palette        = NULL,
                           theme          = c("clean", "minimal")) {
  .require_ggplot2()
  color_by <- match.arg(color_by)
  theme    <- match.arg(theme)

  # Coerce to plot-data
  pd <- .coerce_cta_plot_data(x, target_class = target_class, digits = 1)

  # no_tree: return message panel
  if (isTRUE(pd$no_tree)) {
    msg <- "No tree found (leaf-only result).\nCTA found no valid split above the significance threshold."
    return(.gg_message_panel(msg, title = main %||% "CTA Tree", subtitle = subtitle))
  }

  nodes <- pd$nodes
  edges <- pd$edges

  if (nrow(nodes) == 0L)
    return(.gg_message_panel("No nodes to display.", title = main, subtitle = subtitle))

  nodes <- .normalize_cta_nodes(nodes)
  pal   <- .tree_pal(palette)

  effective_subtitle <- subtitle
  if (isTRUE(show_metrics)) {
    ess_val <- pd$overall_ess %||% NA_real_
    d_val   <- pd$d           %||% NA_real_
    ess_lbl <- pd$ess_label   %||% "ESS"
    if (!is.na(ess_val)) {
      d_str      <- if (is.finite(d_val)) sprintf("%.4f", d_val) else "NA"
      metrics_ln <- sprintf("%s: %.2f%%  |  D: %s", ess_lbl, ess_val, d_str)
      effective_subtitle <- if (!is.null(subtitle))
        paste0(subtitle, "\n", metrics_ln)
      else
        metrics_ln
    }
  }

  opts  <- list(
    main           = main %||% "CTA Tree",
    subtitle       = effective_subtitle,
    show_rule      = isTRUE(show_rule),
    node_text_size = node_text_size,
    edge_text_size = edge_text_size,
    node_width     = node_width,
    theme          = theme
  )

  .render_tree_gg(nodes, edges, pal, color_by, opts)
}

# ---- plot_lort_tree --------------------------------------------------------- #

#' Plot a LORT (Locally Optimal Recursive Tree) using ggplot2
#'
#' Renders a publication-quality tree diagram for a fitted LORT object.
#' Requires the \pkg{ggplot2} package (listed in \code{Suggests}).
#'
#' @param x A \code{cta_ort} object from \code{\link{lort_fit}} (or
#'   \code{cta_fit(..., recursive=TRUE)}), or the list returned by
#'   \code{\link{ort_plot_data}}.
#' @param target_class Integer; target class for endpoint coloring and
#'   target-rate annotation (default \code{1L}).
#' @param color_by Character; controls terminal-node fill color.
#'   \code{"target_rate"} (default): gradient by target-class proportion;
#'   \code{"prediction"}: discrete fill by predicted class;
#'   \code{"none"}: all white.
#' @param main Character; plot title.  Default \code{"LORT"}.
#' @param subtitle Character; plot subtitle.
#' @param show_n Logical; include observation count in node labels.
#' @param show_percent Logical; include ESS in split-node labels.
#' @param show_rule Logical; show branch condition labels on edges.
#' @param show_metrics Logical; if \code{TRUE}, appends an ESS/WESS and D
#'   line to the plot subtitle.  Default \code{FALSE}.
#' @param wrap_width Integer; maximum characters per line.  Default \code{26L}.
#' @param node_text_size Numeric; text size for node labels.  Default \code{3.5}.
#' @param edge_text_size Numeric; text size for edge labels.  Default \code{3.2}.
#' @param node_width Numeric or \code{NULL}; node box width (auto if \code{NULL}).
#' @param palette Named list for color overrides.  \code{NULL} uses defaults.
#' @param theme Character; \code{"clean"} (default) or \code{"minimal"}.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @seealso \code{\link{lort_fit}}, \code{\link{ort_plot_data}},
#'   \code{\link{plot_cta_tree}}, \code{\link[ggplot2]{ggsave}}
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   X <- data.frame(
#'     A = c(rep(0L,20), rep(1L,20), rep(1L,20)),
#'     B = c(rep(0L,20), rep(0L,20), rep(1L,20))
#'   )
#'   y <- c(rep(0L,40), rep(1L,20))
#'   lort <- lort_fit(X, y, mc_iter=100L, mc_seed=42L, loo="off", min_n=5L)
#'   p <- plot_lort_tree(lort, color_by="prediction")
#'   print(p)
#' }
#' }
#' @export
plot_lort_tree <- function(x,
                            target_class   = 1L,
                            color_by       = c("target_rate", "prediction", "none"),
                            main           = NULL,
                            subtitle       = NULL,
                            show_n         = TRUE,
                            show_percent   = TRUE,
                            show_rule      = TRUE,
                            show_metrics   = FALSE,
                            wrap_width     = 26L,
                            node_text_size = 3.5,
                            edge_text_size = 3.2,
                            node_width     = NULL,
                            palette        = NULL,
                            theme          = c("clean", "minimal")) {
  .require_ggplot2()
  color_by <- match.arg(color_by)
  theme    <- match.arg(theme)

  # Coerce to plot-data
  pd <- .coerce_lort_plot_data(x, target_class = target_class, digits = 1L)

  nodes <- pd$nodes
  edges <- pd$edges

  if (nrow(nodes) == 0L)
    return(.gg_message_panel("No ORT nodes to display.",
                              title = main %||% "LORT", subtitle = subtitle))

  nodes <- .normalize_lort_nodes(nodes, pd$strata)

  # Standardize edge column names (ORT uses from_id/to_id; renderer expects x0,y0,x1,y1,label)
  if (!"x0" %in% names(edges) && "x0" %in% names(edges)) {
    # already correct  -  no-op
  }
  # edges from ort_plot_data already have x0,y0,x1,y1,label  -  pass through.

  effective_subtitle <- subtitle
  if (isTRUE(show_metrics)) {
    ess_val <- pd$overall_ess %||% NA_real_
    d_val   <- pd$d           %||% NA_real_
    ess_lbl <- pd$ess_label   %||% "ESS"
    if (!is.na(ess_val)) {
      d_str      <- if (is.finite(d_val)) sprintf("%.4f", d_val) else "NA"
      metrics_ln <- sprintf("%s: %.2f%%  |  D: %s", ess_lbl, ess_val, d_str)
      effective_subtitle <- if (!is.null(subtitle))
        paste0(subtitle, "\n", metrics_ln)
      else
        metrics_ln
    }
  }

  pal  <- .tree_pal(palette)
  opts <- list(
    main           = main %||% "LORT",
    subtitle       = effective_subtitle,
    show_rule      = isTRUE(show_rule),
    node_text_size = node_text_size,
    edge_text_size = edge_text_size,
    node_width     = node_width,
    theme          = theme
  )

  .render_tree_gg(nodes, edges, pal, color_by, opts)
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

  p <- ggplot2::ggplot(tbl,
         ggplot2::aes(x = abs_smd, y = attr_fac)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = abs_smd,
                   y = attr_fac, yend = attr_fac),
      color = "#cccccc", linewidth = 0.6, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = bal_col),
      size = 3.5, na.rm = TRUE
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
      limits = c(0, x_max * 1.05),
      expand = ggplot2::expansion(mult = c(0, 0.04))
    ) +
    ggplot2::labs(x = "|SMD|", y = NULL)

  if (isTRUE(ref_010))
    p <- p + ggplot2::geom_vline(xintercept = 0.10, linetype = "dashed",
                                  color = "#555555", linewidth = 0.5)
  if (isTRUE(ref_020))
    p <- p + ggplot2::geom_vline(xintercept = 0.20, linetype = "dotted",
                                  color = "#888888", linewidth = 0.5)

  .balance_theme(p, base_theme, main %||% "SMD Covariate Balance", subtitle)
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

  # no_tree: favorable balance message
  if (identical(status, "no_tree")) {
    msg <- pd$no_tree_message %||%
      "No discriminating tree found.\nThis is favorable evidence of covariate balance."
    return(.gg_message_panel(msg,
                              title    = main %||% "CTA Covariate Balance",
                              subtitle = subtitle))
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
