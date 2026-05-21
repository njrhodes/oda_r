\name{summary.cta_family}
\alias{summary.cta_family}
\alias{print.cta_family_summary}
\alias{print.cta_family}
\title{Summary and print methods for CTA descendant family}
\description{
\code{summary.cta_family} returns a structured S3 object summarising the CTA
descendant family.  \code{print.cta_family_summary} displays it compactly.
\code{print.cta_family} calls \code{summary} then prints.

All values are read from stored fields — no refitting or recomputation is
performed.
}
\usage{
\method{summary}{cta_family}(object, \dots)
\method{print}{cta_family_summary}(x, \dots)
\method{print}{cta_family}(x, \dots)
}
\arguments{
  \item{object}{A \code{cta_family} from \code{\link{cta_descendant_family}}.}
  \item{x}{A \code{cta_family_summary} (for \code{print.cta_family_summary})
    or \code{cta_family} (for \code{print.cta_family}).}
  \item{\dots}{Unused; included for S3 compatibility.}
}
\value{
\code{summary.cta_family} returns a list of class
\code{c("cta_family_summary", "list")} with fields:
\describe{
  \item{\code{n_members}}{Integer number of family members.}
  \item{\code{min_d_idx}}{Integer index of the feasible member with minimum D;
    \code{NA_integer_} if no feasible member exists.}
  \item{\code{terminated}}{Logical; always \code{TRUE} for a completed chain.}
  \item{\code{termination_reason}}{Character: one of \code{"no_tree"},
    \code{"max_steps"}, \code{"no_next_mindenom"}.}
  \item{\code{has_weights}}{Logical; \code{TRUE} when any family member used
    case weights.}
  \item{\code{table}}{A \code{data.frame} from \code{\link{cta_family_table}}.}
}
\code{print.cta_family_summary} and \code{print.cta_family} return their
first argument invisibly.
}
\seealso{
  \code{\link{cta_descendant_family}}, \code{\link{cta_family_table}}
}
\examples{
data(mtcars)
X <- mtcars[, c("cyl", "disp", "hp", "wt")]
y <- as.integer(mtcars$am)
fam <- suppressMessages(
  cta_descendant_family(X, y, start_mindenom = 1L, mc_iter = 200L,
                        mc_seed = 42L, loo = "off")
)
s <- summary(fam)
print(s)
print(fam)
}
