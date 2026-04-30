\name{oda_confusion_multiclass}
\alias{oda_confusion_multiclass}
\title{Multiclass confusion matrix}
\description{
Compute a weighted multiclass confusion matrix with PAC and PV summaries.
}
\usage{oda_confusion_multiclass(y, y_pred, w = NULL)}
\arguments{
  \item{y}{Actual integer class labels.}
  \item{y_pred}{Predicted integer class labels.}
  \item{w}{Optional numeric weights. Default: unit weights.}
}
\value{
Named list:
\describe{
  \item{\code{confusion}}{C x C numeric matrix of weighted counts. With
    unit weights these are raw integer observation counts. Rows are actual
    classes; columns are predicted classes.}
  \item{\code{correct}}{Total weighted count of correct classifications.}
  \item{\code{overall_acc}}{Overall accuracy as a proportion [0, 1].}
  \item{\code{pac_by_class}}{Per-class sensitivity as proportions [0, 1].}
  \item{\code{mean_pac}}{Mean sensitivity across classes, proportion [0,1].}
  \item{\code{pv_by_class}}{Per-class predictive value, proportions [0, 1].}
  \item{\code{mean_pv}}{Mean predictive value, proportion [0, 1].}
}
}
