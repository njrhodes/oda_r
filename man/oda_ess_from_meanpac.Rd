\name{oda_ess_from_meanpac}
\alias{oda_ess_from_meanpac}
\title{Effect Strength for Sensitivity from mean PAC}
\description{
Compute ESS (Effect Strength for Sensitivity) in percent, scaled against
a chance baseline.
}
\usage{oda_ess_from_meanpac(mean_pac, chance)}
\arguments{
  \item{mean_pac}{Mean PAC as a proportion [0, 1].}
  \item{chance}{Chance baseline as a proportion (e.g. 0.5 for 2-class).}
}
\value{ESS as a percentage [0, 100].}
