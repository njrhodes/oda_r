\name{cta_fit}
\alias{cta_fit}
\title{Fit a Classification Tree Analysis (CTA) model (public wrapper)}
\description{
Public entry point for CTA. Currently supports binary (two-class) outcome
variables only. Validates the class variable and delegates to
\code{\link{oda_cta_fit}}.
}
\usage{
cta_fit(X, y, verbose = FALSE, ...)
}
\arguments{
  \item{X}{Data frame or matrix of attribute columns.}
  \item{y}{Integer class variable vector. Must have exactly two distinct values.}
  \item{verbose}{Logical; if \code{TRUE}, emit \code{[CTA]} progress messages
    via \code{message()} at each major stage (fast screen, MC+LOO fit,
    ENUMERATE scoring, final selection). Default \code{FALSE}.}
  \item{...}{Additional arguments passed to \code{\link{oda_cta_fit}}.}
}
\value{
A \code{cta_tree} object. Use \code{predict()} to classify new observations
and \code{print()} for a node summary.
}
\seealso{
\code{\link{oda_cta_fit}}, \code{\link{oda_fit}}
}
