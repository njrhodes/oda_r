\name{cta_confusion_table}
\alias{cta_confusion_table}
\title{Final selected tree training confusion table}
\description{
Returns the stored full-tree training confusion matrix for the final selected
CTA model in tidy long format (one row per actual x predicted class pair).

The confusion matrix is captured at fit time at the exact moment the winning
candidate is selected, using the same scoring predictions.  For the expanded
ENUMERATE phase, predictions use majority-fallback for missing attributes.
For the root-only stump phase, predictions are path-local (observations whose
root attribute is missing are excluded).

This function does \strong{not} report split-node local confusion.
Split-node confusion reflects all observations at a node classified by that
node's rule alone; it is not the same as full-tree confusion for trees with
more than one split.  The two coincide incidentally for stumps but the
semantics here are always final-tree.
}
\usage{cta_confusion_table(tree)}
\arguments{
  \item{tree}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
}
\value{
A \code{data.frame} with columns:
\describe{
  \item{\code{actual}}{Integer actual class label.}
  \item{\code{predicted}}{Integer predicted class label.}
  \item{\code{n}}{Integer raw count of observations with this actual x
    predicted combination in the final selected tree.}
}
Rows are sorted by \code{actual} then \code{predicted}.
For a no-tree fit (or if \code{training_confusion} is absent), the returned
data frame has zero rows but the correct column structure.
}
\seealso{
  \code{\link{oda_cta_fit}}, \code{\link{summary.cta_tree}},
  \code{\link{cta_endpoint_table}}, \code{\link{cta_node_table}}
}
\examples{
data(mtcars)
X    <- mtcars[, c("cyl", "disp", "hp", "wt")]
y    <- as.integer(mtcars$am)
tree <- oda_cta_fit(X, y, mindenom = 5L, mc_iter = 500L, mc_seed = 42L)
cta_confusion_table(tree)
}
