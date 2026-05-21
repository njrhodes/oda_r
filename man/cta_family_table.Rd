\name{cta_family_table}
\alias{cta_family_table}
\title{Tidy table of a CTA descendant family}
\description{
Returns a \code{data.frame} with one row per family member, reading all values
from the stored \code{cta_family} object.  No refitting or recomputation is
performed.
}
\usage{cta_family_table(family)}
\arguments{
  \item{family}{A \code{cta_family} from \code{\link{cta_descendant_family}}.}
}
\value{
A \code{data.frame} with columns:
\describe{
  \item{\code{index}}{Integer position of the member in the chain.}
  \item{\code{mindenom}}{Integer MINDENOM used for this fit.}
  \item{\code{status}}{Character: \code{"valid_tree"}, \code{"stump"}, or
    \code{"no_tree"}.}
  \item{\code{no_tree}}{Logical; \code{TRUE} for the terminal no-tree member.}
  \item{\code{strata}}{Integer number of terminal leaf endpoints;
    \code{NA} for no-tree members.}
  \item{\code{min_terminal_denom}}{Integer minimum leaf \code{n_obs};
    \code{NA} for no-tree members.}
  \item{\code{next_mindenom}}{Integer MINDENOM for the next chain step
    (\code{min_terminal_denom + 1}); \code{NA} for no-tree members.}
  \item{\code{overall_ess}}{Numeric overall ESS or WESS stored at fit time;
    \code{NA} for no-tree members.}
  \item{\code{has_weights}}{Logical; \code{TRUE} when case weights were active
    for this fit.}
  \item{\code{d}}{Numeric D statistic (\code{100 / (ESS / strata) - strata});
    \code{NA} for no-tree members.}
  \item{\code{selected_min_d}}{Logical; \code{TRUE} for the feasible member
    with minimum D (index \code{family$min_d_idx}).  All \code{FALSE} when no
    feasible member exists.}
}
}
\seealso{
  \code{\link{cta_descendant_family}}, \code{\link{summary.cta_family}}
}
\examples{
data(mtcars)
X <- mtcars[, c("cyl", "disp", "hp", "wt")]
y <- as.integer(mtcars$am)
fam <- suppressMessages(
  cta_descendant_family(X, y, start_mindenom = 1L, mc_iter = 200L,
                        mc_seed = 42L, loo = "off")
)
cta_family_table(fam)
}
