\name{plot.cta_tree}
\alias{plot.cta_tree}
\title{Plot a fitted CTA tree}
\description{
Native base-R CTA visualization.  Calls \code{\link{cta_plot_data}} for
layout; uses only base graphics  -  no external package dependencies.

Split (internal) nodes are drawn as \strong{ellipses}; terminal endpoint
nodes are drawn as \strong{rectangles}; edges are \strong{directed arrows}.
Split nodes show the split attribute, node-level ESS or WESS, and
observation count.  Without \code{target_class}, leaf nodes show the
majority-class prediction and observation count.  With \code{target_class},
leaf nodes show the target-class count, percentage, predicted class, and
stage from \code{\link{cta_staging_table}}.  Edge labels show the branch
condition (e.g. \code{"V14<=0.5"}).

\strong{Color note:} when \code{target_class} is supplied, endpoint fill
colors are assigned by ascending rank of each endpoint's target-class
proportion within this tree.  Colors encode relative position in the
endpoint distribution and do \emph{not} imply clinical thresholds or
categories.  Supply a custom palette via \code{endpoint_palette} to change
the color encoding.  Use \code{show_caption = TRUE} to render an explicit
note on the plot.

\code{\link{cta_plot_data}} is the renderer-independent data contract.
This function (\code{plot.cta_tree}) is the current native base-R renderer.
}
\usage{
\method{plot}{cta_tree}(x,
     target_class = NULL, class_labels = NULL, digits = 1,
     main = "CTA Tree", show_counts = TRUE, show_stage = TRUE,
     endpoint_palette = NULL, endpoint_fill = "#D9F7E6",
     split_fill = "#D9EAF7", node_col_split = NULL,
     node_col_leaf = NULL, edge_col = "grey40",
     border_col = "grey30", text_col = "black",
     arrow_col = NULL, show_caption = FALSE,
     cex = 0.75, \dots)
}
\arguments{
  \item{x}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
  \item{target_class}{Integer target class for endpoint annotation; passed
    to \code{\link{cta_plot_data}}.  \code{NULL} (default) produces a
    structural plot without endpoint enrichment.}
  \item{class_labels}{Optional display names for class labels; passed to
    \code{\link{cta_plot_data}}.  See \code{\link{cta_plot_data}} for
    accepted formats.}
  \item{digits}{Decimal places for percentage labels in enriched endpoint
    nodes; passed to \code{\link{cta_plot_data}}.  Default \code{1}.}
  \item{main}{Character plot title.  Default \code{"CTA Tree"}.}
  \item{show_counts}{Logical; include \code{n=a/b} raw counts in endpoint
    labels when \code{target_class} is supplied.  Default \code{TRUE}.}
  \item{show_stage}{Logical; include \code{Stage s} line in endpoint labels
    when \code{target_class} is supplied.  Default \code{TRUE}.}
  \item{endpoint_palette}{Palette for endpoint fill colors when
    \code{target_class} is supplied; passed to \code{\link{cta_plot_data}}.
    \code{NULL} uses the default gradient.  Accepts a palette
    \code{function(n)} or a character vector of colors.}
  \item{endpoint_fill}{Default fill colour for leaf (terminal) nodes when
    \code{target_class} is \code{NULL}.  Default \code{"#D9F7E6"}
    (light green).}
  \item{split_fill}{Fill colour for split (internal) ellipse nodes.
    Default \code{"#D9EAF7"} (light blue).}
  \item{node_col_split}{Legacy alias for \code{split_fill}; overrides it
    when non-\code{NULL}.}
  \item{node_col_leaf}{Legacy alias for \code{endpoint_fill}; overrides it
    when non-\code{NULL}.}
  \item{edge_col}{Colour for directed edge arrows.  Default \code{"grey40"}.}
  \item{border_col}{Border colour for all nodes.  Default \code{"grey30"}.}
  \item{text_col}{Text colour for node labels.  Default \code{"black"}.}
  \item{arrow_col}{Arrow colour for directed edges.  \code{NULL} (default)
    uses \code{edge_col}.}
  \item{show_caption}{Logical; if \code{TRUE} and \code{target_class} is
    supplied, adds a bottom caption: \emph{"Endpoint fill: relative
    target-class proportion within this tree. Not a clinical threshold."}
    Default \code{FALSE}.}
  \item{cex}{Text expansion factor for node labels.  Default \code{0.75}.}
  \item{\dots}{Unused; included for S3 compatibility.}
}
\value{
\code{invisible(pd)}, where \code{pd} is the \code{\link{cta_plot_data}}
list used to render the plot.  The caller can inspect layout coordinates,
enrichment columns, and endpoint annotations from the returned object.
}
\seealso{
\code{\link{cta_plot_data}}, \code{\link{cta_staging_table}},
\code{\link{oda_cta_fit}}
}
\examples{
data(mtcars)
X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
y    <- as.integer(mtcars$am)
tree <- suppressMessages(
  oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L,
              loo = "off")
)

# Structural plot
plot(tree)

# Target-class enriched plot with custom labels
plot(tree, target_class = 1L,
     class_labels = c("0" = "Manual", "1" = "Auto"))

# Custom palette (white to dark red)
plot(tree, target_class = 1L,
     endpoint_palette = c("#ffffff", "#c62828"))
}
