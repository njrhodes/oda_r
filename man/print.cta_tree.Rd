\name{print.cta_tree}
\alias{print.cta_tree}
\title{Print a CTA tree in MegaODA node table format}
\description{Displays each split node with its attribute, depth, n, p-value,
ESS, LOO status, and rule string, followed by the node confusion matrix.}
\usage{\method{print}{cta_tree}(x, ...)}
\arguments{
  \item{x}{A \code{cta_tree} from \code{\link{oda_cta_fit}}.}
  \item{...}{Unused.}
}
\value{Invisibly returns \code{x}.}
