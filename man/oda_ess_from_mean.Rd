\name{oda_ess_from_mean}
\alias{oda_ess_from_mean}
\title{ESS from mean metric for a C-class problem}
\description{
Compute Effect Strength for Sensitivity from mean PAC or mean PV for a
problem with C classes.
}
\usage{oda_ess_from_mean(mean_metric, C)}
\arguments{
  \item{mean_metric}{Mean PAC or mean PV as a proportion [0, 1].}
  \item{C}{Number of classes.}
}
\value{ESS as a percentage [0, 100]. Chance baseline is 1/C.}
