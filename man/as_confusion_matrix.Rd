\name{as_confusion_matrix}
\alias{as_confusion_matrix}
\title{Convert a tidy confusion data frame to a 2x2 integer matrix}
\description{
Converts the \code{data.frame} returned by \code{\link{cta_confusion_table}}
(columns \code{actual}, \code{predicted}, \code{n}) to a 2x2 integer matrix
suitable for \code{\link{novo_boot_ci}}.
}
\usage{
as_confusion_matrix(df)
}
\arguments{
  \item{df}{A \code{data.frame} with integer columns \code{actual},
    \code{predicted}, and \code{n}.  Must represent a binary (2-class)
    classification with class labels 0 and 1.}
}
\value{
A 2x2 integer matrix with rows = actual class (0/1) and columns =
predicted class (0/1), matching the \code{training_confusion} convention
used throughout \code{odacore}.  Row and column names are \code{"0"} and
\code{"1"}.
}
\examples{
# From raw data frame:
df <- data.frame(
  actual    = c(0L, 0L, 1L, 1L),
  predicted = c(0L, 1L, 0L, 1L),
  n         = c(146L, 40L, 36L, 33L)
)
m <- as_confusion_matrix(df)
novo_boot_ci(m, nboot = 200L, seed = 1L)

# From a fitted tree:
\dontrun{
ct <- cta_confusion_table(tree)
m  <- as_confusion_matrix(ct)
novo_boot_ci(m, nboot = 5000L, seed = 42L)
}
}
\seealso{
\code{\link{cta_confusion_table}}, \code{\link{novo_boot_ci}}
}
