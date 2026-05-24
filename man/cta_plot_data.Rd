\name{cta_plot_data}
\alias{cta_plot_data}
\title{Extract layout data for plotting a CTA tree}
\description{
Returns a pure-data list describing tree topology and layout coordinates.
No graphics are produced.  Use this as input to \code{\link{plot.cta_tree}}
or to custom rendering code.

Layout algorithm: leaves receive sequential integer x-positions in
depth-first (left-to-right) order; internal nodes are centred over their
children.  \code{y = -depth} so the root sits at the top.

\strong{Target-class enrichment:} when \code{target_class} is supplied,
each terminal leaf is joined to \code{\link{cta_staging_table}} and
annotated with target-class counts, proportions, and a continuous display
color derived from the endpoint's rank among all endpoints by ascending
target-class proportion.  Colors encode relative position within this
tree's endpoint distribution and do \emph{not} imply clinical thresholds
or categories.
}
\usage{
cta_plot_data(tree, target_class = NULL, class_labels = NULL,
              digits = 1, endpoint_palette = NULL)
}
\arguments{
  \item{tree}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
  \item{target_class}{Integer (or \code{NULL}); the class label treated as
    the target when enriching endpoint annotations.  \code{NULL} (default)
    returns the structural layout only — no endpoint enrichment columns are
    added.}
  \item{class_labels}{Optional character vector of display names for class
    labels.  Supply as a \emph{named} vector, e.g.
    \code{c("0" = "Alive", "1" = "Deceased")}, or as a positional vector
    (index \eqn{k+1} maps to class \eqn{k}).  \code{NULL} (default) uses
    \code{"class=C"} formatting.}
  \item{digits}{Integer number of decimal places for percentage formatting
    in \code{endpoint_label}.  Default \code{1}.}
  \item{endpoint_palette}{Palette for endpoint fill colors, used only when
    \code{target_class} is supplied.  Accepts \code{NULL} (default gradient:
    warm pink to pale yellow to soft green), a palette \code{function(n)}
    returning \code{n} color strings, or a character vector of colors
    interpolated via \code{\link[grDevices]{colorRampPalette}}.}
}
\value{
When \code{target_class = NULL}: a list with elements
\code{nodes}, \code{edges}, \code{no_tree}, \code{has_weights}.

When \code{target_class} is supplied: the same list plus
\code{endpoints} (a staging data.frame with layout coordinates) and
\code{target_class_used} (the integer target class used).

\describe{
  \item{\code{nodes}}{A \code{data.frame} with one row per node.
    Always-present columns: \code{node_id} (integer), \code{parent_id}
    (integer), \code{depth} (integer), \code{x} (numeric), \code{y}
    (numeric), \code{leaf} (logical), \code{attribute} (character;
    \code{NA} for leaves), \code{n_obs} (integer), \code{majority_class}
    (integer), \code{ess} (numeric; \code{NA} for leaves), \code{label}
    (character multi-line display text).

    Additional columns present when \code{target_class} is supplied
    (values are \code{NA} on split nodes): \code{endpoint_id} (integer),
    \code{stage} (integer), \code{target_class} (integer),
    \code{target_n} (numeric), \code{denominator} (numeric),
    \code{target_proportion} (numeric; raw continuous proportion, not
    binned), \code{target_rank} (integer; ascending rank of proportion,
    ties broken by \code{ties.method = "first"}),
    \code{endpoint_fill_color} (character hex color assigned by rank
    within this tree — does not imply clinical thresholds or categories),
    \code{predicted_label} (character), \code{target_label} (character),
    \code{endpoint_label} (character multi-line display text).}
  \item{\code{edges}}{A \code{data.frame} with one row per parent-to-child
    edge and columns \code{from_node_id} (integer), \code{to_node_id}
    (integer), \code{x0}, \code{y0}, \code{x1}, \code{y1} (numeric
    centre-to-centre coordinates), \code{label} (character branch
    condition, e.g. \code{"V14<=0.5"}).}
  \item{\code{endpoints}}{(\code{target_class} only) A \code{data.frame}
    with one row per endpoint, ordered by ascending stage.  Columns include
    all staging fields from \code{\link{cta_staging_table}} plus layout
    coordinates \code{x}, \code{y}, display columns
    \code{predicted_label}, \code{target_label}, \code{endpoint_fill_color},
    and integer \code{target_rank}.}
  \item{\code{target_class_used}}{(\code{target_class} only) The integer
    \code{target_class} argument used for enrichment.}
  \item{\code{no_tree}}{Logical; \code{TRUE} for leaf-only fits.}
  \item{\code{has_weights}}{Logical; \code{TRUE} when case weights are
    active.}
}
}
\seealso{
\code{\link{plot.cta_tree}}, \code{\link{cta_staging_table}},
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

# Structural layout only
pd <- cta_plot_data(tree)
head(pd$nodes)
pd$edges

# Target-class enrichment
pd2 <- cta_plot_data(tree, target_class = 1L,
                     class_labels = c("0" = "Manual", "1" = "Auto"))
pd2$endpoints[, c("stage", "target_proportion", "endpoint_fill_color")]
}
