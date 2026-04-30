\name{oda_confusion_binary}
\alias{oda_confusion_binary}
\title{Binary confusion table}
\description{
Compute a weighted binary confusion table from actual and predicted labels.
}
\usage{oda_confusion_binary(y, y_pred, w = NULL)}
\arguments{
  \item{y}{Actual class labels (0/1 integer).}
  \item{y_pred}{Predicted class labels (0/1 integer).}
  \item{w}{Optional numeric weights. Default: unit weights.}
}
\value{
Named list with integer count fields \code{TP}, \code{TN}, \code{FP},
\code{FN} (weighted sums), and rate fields \code{sensitivity},
\code{specificity} (proportions in [0, 1]), and \code{mean_pac}
(proportion in [0, 1]).

\strong{Note:} With unit weights these are raw integer counts. With
prior-odds weights (from \code{\link{oda_univariate_core}} with
\code{priors_on = TRUE}) they are weighted counts, not raw integers.
}
