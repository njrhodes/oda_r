###############################################################################
# R/oda_power.R
# ODA power and sample size analysis.
#
# Canon references:
#   Rhodes, N. J. (2020). Statistical power analysis in ODA, CTA and
#     Novometrics. Optimal Data Analysis, 9.
#     https://odajournal.files.wordpress.com/2020/02/v9a5.pdf
#   Yarnold PR, Soltysik RC (2005). Optimal Data Analysis: A Guidebook with
#     Software for Windows. Washington, DC: APA Books.
#
# Exported:
#   oda_power()        — power of a 2-class ODA model via simulation.
#   oda_sample_size()  — minimum n achieving a target power (bisection).
#   print.oda_power
#   print.oda_sample_size
###############################################################################

# ---- internal helpers --------------------------------------------------------

# Two-sided Fisher exact p-value via hypergeometric distribution.
# Equivalent to stats::fisher.test(alternative = "two.sided")$p.value but
# avoids S3/C dispatch overhead in tight simulation loops.
# a  = cell[1,1], r1 = row-1 total, r2 = row-2 total, c1 = col-1 total.
.fisher_pval_2x2 <- function(a, r1, r2, c1) {
  lo       <- max(0L, c1 - r2)
  hi       <- min(r1, c1)
  if (lo > hi) return(1.0)
  prob_all <- stats::dhyper(lo:hi, r1, r2, c1)
  prob_obs <- stats::dhyper(a,     r1, r2, c1)
  # Tolerance matches fisher.test source (workspace = 200000)
  sum(prob_all[prob_all <= prob_obs * (1 + 1e-7)])
}

# Convert ESS (%) to (p1, p2) under the symmetric balanced convention.
# Convention: p2 = accuracy of class 1 = ESS/200 + 0.5; p1 = 1 - p2.
# ESS = 100 * (Mean PAC - 0.5) / 0.5; with p1 = 1 - p2, Mean PAC = p2,
# so ESS = 100 * (p2 - 0.5) / 0.5.
.ess_to_proportions <- function(ess) {
  if (!is.numeric(ess) || length(ess) != 1L || ess <= 0 || ess >= 100)
    stop("'ess' must be a single numeric in (0, 100).", call. = FALSE)
  p2 <- ess / 200 + 0.5
  list(p1 = 1 - p2, p2 = p2)
}

# Core simulation: returns a vector of nsim Fisher p-values.
# Seed management is the caller's responsibility.
.power_pvals <- function(nsim, r1, r2, p1, p2) {
  draws1 <- stats::rbinom(nsim, r1, p1)
  draws2 <- stats::rbinom(nsim, r2, p2)
  c1_vec <- draws1 + draws2          # col-1 marginal for each replicate
  mapply(.fisher_pval_2x2,
         a  = draws1,
         c1 = c1_vec,
         MoreArgs = list(r1 = r1L <- as.integer(r1),
                         r2 = r2L <- as.integer(r2)))
}

# ---- oda_power ---------------------------------------------------------------

#' ODA power analysis via simulation
#'
#' Estimates planning power for unit-weighted binary 2\ifelse{html}{\out{&times;}}{\eqn{\times}}2 ODA-equivalent
#' designs.  The implemented design assumes fixed group sizes \code{n1}/\code{n2}
#' and binomial outcome probabilities \code{p1}/\code{p2}, then evaluates whether
#' the resulting 2\ifelse{html}{\out{&times;}}{\eqn{\times}}2 table is significant by Fisher's exact test at the
#' (optionally Sidak-adjusted) alpha.
#'
#' This is the binary lowest-measurement planning case discussed by Rhodes (2020).
#' See also Yarnold and Soltysik (2005) for the underlying ODA/Fisher isomorphism.
#' \strong{Scope:} unit-weighted, binary class, binary (2-level) attribute only.
#' This is not a general CTA, LORT, SDA, weighted, or multiclass power method.
#'
#' \strong{Method:} For each Monte Carlo replicate, binomial draws are generated
#' under (\code{p1}, \code{p2}) with fixed group sizes (\code{n1}, \code{n2}).
#' The resulting 2\ifelse{html}{\out{&times;}}{\eqn{\times}}2 table is tested by Fisher's exact test; power is the
#' proportion of replicates in which the null is rejected.  The prospective
#' sampling treats group sizes as fixed and outcomes as binomial within each
#' group; the Fisher test is then applied to the generated table with its
#' realized marginals.  This is the standard simulation-based power approach
#' for 2\ifelse{html}{\out{&times;}}{\eqn{\times}}2 contingency analyses.
#'
#' \strong{Effect-size input:}
#' Specify the effect either as per-group proportions \code{p1} and \code{p2}
#' directly, or as \code{ess} (Effect Strength for Sensitivity, percent) under
#' the symmetric balanced convention:
#' \eqn{p_2 = \text{ESS}/200 + 0.5}, \eqn{p_1 = 1 - p_2}.
#'
#' \strong{Sidak correction:}
#' When \code{comp > 1}, the working \eqn{\alpha} is Sidak-adjusted:
#' \eqn{\alpha_{\text{adj}} = 1 - (1 - \alpha)^{1/\text{comp}}}.
#'
#' @param n1 Integer (or integer vector) giving the per-group sample size for
#'   class 0.  When a vector is supplied, power is estimated at each element.
#' @param n2 Integer (or integer vector) giving the per-group sample size for
#'   class 1.  Defaults to \code{n1} (balanced design).  If scalar and
#'   \code{n1} is a vector, \code{n2} is recycled.
#' @param p1 Probability of the event in class 0.  Ignored when \code{ess} is
#'   supplied.
#' @param p2 Probability of the event in class 1.  Ignored when \code{ess} is
#'   supplied.
#' @param ess Effect Strength for Sensitivity (percent, \eqn{0 < \text{ESS} < 100})
#'   under the symmetric balanced convention.  Mutually exclusive with
#'   \code{p1}/\code{p2}.
#' @param alpha Nominal significance level.  Default 0.05.  May be a vector
#'   to evaluate power at multiple alpha levels simultaneously.
#' @param comp Number of comparisons for Sidak multiple-comparison correction.
#'   Default 1 (no correction).  Must be a single positive integer.
#' @param nsim Number of Monte Carlo replications per cell.  Default 10000.
#' @param mc_seed Integer seed passed to \code{set.seed()} once before all
#'   simulations, or \code{NULL} to use the current RNG state.
#'
#' @return An object of class \code{"oda_power"}, a list with elements:
#' \describe{
#'   \item{\code{power}}{Numeric matrix (rows = n1, cols = alpha_adj) of
#'     estimated power.  Simplified to a named vector if one dimension is
#'     scalar, or to a scalar if both are.}
#'   \item{\code{n1}, \code{n2}}{Per-group sample sizes.}
#'   \item{\code{p1}, \code{p2}}{Per-group event rates used.}
#'   \item{\code{ess_input}}{ESS supplied, or \code{NA} if \code{p1}/\code{p2} used.}
#'   \item{\code{alpha}, \code{alpha_adj}}{Input and Sidak-adjusted alpha.}
#'   \item{\code{comp}, \code{nsim}, \code{mc_seed}}{Input parameters.}
#' }
#'
#' @references
#' Rhodes, N. J. (2020). Statistical power analysis in ODA, CTA and
#'   Novometrics. \emph{Optimal Data Analysis}, 9.
#'   \url{https://odajournal.files.wordpress.com/2020/02/v9a5.pdf}
#'
#' Yarnold PR, Soltysik RC (2005). \emph{Optimal Data Analysis: A Guidebook
#'   with Software for Windows.} Washington, DC: APA Books.
#'
#' @examples
#' # Power for ESS = 48%, n = 50 per group (CRAN-safe nsim; use 10000L for publication)
#' oda_power(n1 = 50, ess = 48, nsim = 500L, mc_seed = 42L)
#'
#' # Power curve across a range of n
#' oda_power(n1 = c(30, 50, 80), ess = 48, nsim = 500L, mc_seed = 42L)
#'
#' # Direct proportions (p1 = 0.26, p2 = 0.74)
#' oda_power(n1 = 50, p1 = 0.26, p2 = 0.74, nsim = 500L, mc_seed = 42L)
#'
#' # Sidak correction for 3 comparisons
#' oda_power(n1 = 80, ess = 48, comp = 3L, nsim = 500L, mc_seed = 42L)
#'
#' @export
oda_power <- function(n1,
                      n2      = n1,
                      p1      = NULL,
                      p2      = NULL,
                      ess     = NULL,
                      alpha   = 0.05,
                      comp    = 1L,
                      nsim    = 10000L,
                      mc_seed = NULL) {

  # ---- effect size -----------------------------------------------------------
  if (!is.null(ess) && (!is.null(p1) || !is.null(p2)))
    stop("Supply either 'ess' or ('p1','p2'), not both.", call. = FALSE)
  if (is.null(ess) && (is.null(p1) || is.null(p2)))
    stop("Supply either 'ess' or both 'p1' and 'p2'.", call. = FALSE)

  if (!is.null(ess)) {
    props     <- .ess_to_proportions(ess)
    p1        <- props$p1
    p2        <- props$p2
    ess_input <- ess
  } else {
    if (!is.numeric(p1) || length(p1) != 1L || p1 <= 0 || p1 >= 1)
      stop("'p1' must be a single numeric in (0, 1).", call. = FALSE)
    if (!is.numeric(p2) || length(p2) != 1L || p2 <= 0 || p2 >= 1)
      stop("'p2' must be a single numeric in (0, 1).", call. = FALSE)
    ess_input <- NA_real_
  }

  # ---- sample sizes ----------------------------------------------------------
  n1 <- as.integer(n1)
  n2 <- as.integer(n2)
  if (length(n2) == 1L && length(n1) > 1L) n2 <- rep(n2, length(n1))
  if (length(n1) != length(n2))
    stop("'n1' and 'n2' must have the same length.", call. = FALSE)
  if (any(n1 < 2L) || any(n2 < 2L))
    stop("'n1' and 'n2' must be >= 2.", call. = FALSE)

  # ---- alpha / comp ----------------------------------------------------------
  if (!is.numeric(alpha) || any(alpha <= 0) || any(alpha >= 1))
    stop("'alpha' must be in (0, 1).", call. = FALSE)
  comp <- as.integer(comp)
  if (length(comp) != 1L || comp < 1L)
    stop("'comp' must be a single positive integer.", call. = FALSE)

  alpha_adj <- if (comp > 1L) 1 - (1 - alpha)^(1 / comp) else alpha

  nsim <- as.integer(nsim)
  if (length(nsim) != 1L || nsim < 1L)
    stop("'nsim' must be a single positive integer.", call. = FALSE)

  # ---- simulation ------------------------------------------------------------
  if (!is.null(mc_seed)) set.seed(mc_seed)

  n_sizes <- length(n1)
  n_alpha <- length(alpha_adj)
  pwr_mat <- matrix(NA_real_, nrow = n_sizes, ncol = n_alpha,
                    dimnames = list(paste0("n=", n1),
                                   paste0("alpha=", round(alpha_adj, 6))))

  for (i in seq_len(n_sizes)) {
    r1 <- n1[i]; r2 <- n2[i]
    # Draw all replicates for this n-size in one pass (seed already set above)
    draws1 <- stats::rbinom(nsim, r1, p1)
    draws2 <- stats::rbinom(nsim, r2, p2)
    c1_vec <- draws1 + draws2   # col-1 marginal varies across replicates

    pvals <- mapply(.fisher_pval_2x2,
                    a  = draws1,
                    c1 = c1_vec,
                    MoreArgs = list(r1 = r1, r2 = r2))

    for (j in seq_len(n_alpha)) {
      pwr_mat[i, j] <- mean(pvals < alpha_adj[j])
    }
  }

  # Simplify output dimensions
  pwr_out <- if (n_sizes == 1L && n_alpha == 1L) {
    pwr_mat[1L, 1L]
  } else if (n_alpha == 1L) {
    setNames(pwr_mat[, 1L], paste0("n=", n1))
  } else {
    pwr_mat
  }

  structure(
    list(
      power     = pwr_out,
      n1        = n1,
      n2        = n2,
      p1        = p1,
      p2        = p2,
      ess_input = ess_input,
      alpha     = alpha,
      alpha_adj = alpha_adj,
      comp      = comp,
      nsim      = nsim,
      mc_seed   = mc_seed
    ),
    class = "oda_power"
  )
}

# ---- print.oda_power ---------------------------------------------------------

#' @export
print.oda_power <- function(x, digits = 3L, ...) {
  cat("ODA Power Analysis\n")
  cat(sprintf("  Effect:  p1 = %.4f, p2 = %.4f", x$p1, x$p2))
  if (!is.na(x$ess_input))
    cat(sprintf("  (ESS = %.2f%%)", x$ess_input))
  cat("\n")
  if (x$comp > 1L) {
    cat(sprintf("  alpha = %.4f  (Sidak-adjusted for %d comparisons: %.6f)\n",
                x$alpha, x$comp, x$alpha_adj))
  } else {
    cat(sprintf("  alpha = %.4f\n", x$alpha))
  }
  cat(sprintf("  nsim  = %d\n", x$nsim))
  cat("\n")

  pwr <- x$power
  if (is.matrix(pwr)) {
    cat("Power estimates:\n")
    print(round(pwr, digits), ...)
  } else if (length(pwr) > 1L) {
    cat("Power estimates:\n")
    print(round(pwr, digits), ...)
  } else {
    cat(sprintf("Power = %.3f\n", round(pwr, digits)))
  }
  invisible(x)
}

# ---- oda_sample_size ---------------------------------------------------------

#' ODA minimum sample size via bisection
#'
#' Finds the minimum per-group sample size \eqn{n} (balanced design) at which
#' power reaches or exceeds \code{power_target}.  Uses bisection over
#' \code{oda_power()} with a fixed RNG seed for stable search.
#'
#' \strong{Scope:} unit-weighted, binary class, binary (2-level) attribute only.
#' This is not a general CTA, LORT, SDA, weighted, or multiclass sample-size
#' method.  For unbalanced designs, call \code{oda_power()} directly across a
#' candidate grid.
#'
#' @param power_target Target power.  Default 0.80.
#' @param p1 Probability of the event in class 0.  Ignored when \code{ess} is
#'   supplied.
#' @param p2 Probability of the event in class 1.  Ignored when \code{ess} is
#'   supplied.
#' @param ess Effect Strength for Sensitivity (percent, \eqn{0 < \text{ESS} < 100}).
#'   Mutually exclusive with \code{p1}/\code{p2}.
#' @param alpha Nominal significance level.  Default 0.05.
#' @param comp Number of comparisons for Sidak correction.  Default 1.
#' @param nsim Number of Monte Carlo replications per candidate \eqn{n}.
#'   Default 10000.
#' @param mc_seed Integer seed used for every \code{oda_power()} call during
#'   the bisection.  Using a fixed seed yields stable, reproducible results.
#'   Default 42L.
#' @param n_min Minimum \eqn{n} to search.  Default 2.
#' @param n_max Maximum \eqn{n} to search.  If power at \code{n_max} is still
#'   below \code{power_target}, an error is raised.  Default 2000.
#'
#' @return An object of class \code{"oda_sample_size"}, a list with elements:
#' \describe{
#'   \item{\code{n}}{Minimum per-group sample size achieving \code{power_target}.}
#'   \item{\code{power_achieved}}{Estimated power at \code{n}.}
#'   \item{\code{power_target}}{Input target power.}
#'   \item{\code{p1}, \code{p2}, \code{ess_input}}{Effect-size inputs.}
#'   \item{\code{alpha}, \code{alpha_adj}, \code{comp}}{Alpha parameters.}
#'   \item{\code{nsim}, \code{mc_seed}}{Simulation parameters.}
#' }
#'
#' @references
#' Rhodes, N. J. (2020). Statistical power analysis in ODA, CTA and
#'   Novometrics. \emph{Optimal Data Analysis}, 9.
#'   \url{https://odajournal.files.wordpress.com/2020/02/v9a5.pdf}
#'
#' Yarnold PR, Soltysik RC (2005). \emph{Optimal Data Analysis: A Guidebook
#'   with Software for Windows.} Washington, DC: APA Books.
#'
#' @examples
#' # Minimum n for ESS = 48%, 80% power (CRAN-safe nsim; use 10000L for publication)
#' oda_sample_size(ess = 48, nsim = 500L, mc_seed = 42L)
#'
#' # 90% power target
#' oda_sample_size(ess = 48, power_target = 0.90, nsim = 500L, mc_seed = 42L)
#'
#' @export
oda_sample_size <- function(power_target = 0.80,
                            p1      = NULL,
                            p2      = NULL,
                            ess     = NULL,
                            alpha   = 0.05,
                            comp    = 1L,
                            nsim    = 10000L,
                            mc_seed = 42L,
                            n_min   = 2L,
                            n_max   = 2000L) {

  # ---- effect size (delegate validation to oda_power) ------------------------
  n_min <- as.integer(n_min)
  n_max <- as.integer(n_max)
  if (n_min < 2L) stop("'n_min' must be >= 2.", call. = FALSE)
  if (n_max <= n_min) stop("'n_max' must be > 'n_min'.", call. = FALSE)
  if (power_target <= 0 || power_target >= 1)
    stop("'power_target' must be in (0, 1).", call. = FALSE)

  # Helper: power at a scalar n (re-seeds on every call for stable bisection)
  .pwr_at_n <- function(n) {
    fit <- oda_power(n1 = n, p1 = p1, p2 = p2, ess = ess,
                     alpha = alpha, comp = comp,
                     nsim = nsim, mc_seed = mc_seed)
    as.numeric(fit$power)
  }

  # Check bounds
  pwr_min <- .pwr_at_n(n_min)
  if (pwr_min >= power_target) {
    n_found <- n_min
    pwr_found <- pwr_min
  } else {
    pwr_max <- .pwr_at_n(n_max)
    if (pwr_max < power_target)
      stop(sprintf(
        "Power at n_max = %d is %.3f < target %.3f. Increase 'n_max'.",
        n_max, pwr_max, power_target), call. = FALSE)

    # Bisection
    lo <- n_min; hi <- n_max
    while (hi - lo > 1L) {
      mid     <- (lo + hi) %/% 2L
      pwr_mid <- .pwr_at_n(mid)
      if (pwr_mid >= power_target) hi <- mid else lo <- mid
    }
    n_found   <- hi
    pwr_found <- .pwr_at_n(hi)
  }

  # Retrieve alpha_adj and effect-size fields from a final oda_power call
  fit_final <- oda_power(n1 = n_found, p1 = p1, p2 = p2, ess = ess,
                         alpha = alpha, comp = comp,
                         nsim = nsim, mc_seed = mc_seed)

  structure(
    list(
      n               = n_found,
      power_achieved  = pwr_found,
      power_target    = power_target,
      p1              = fit_final$p1,
      p2              = fit_final$p2,
      ess_input       = fit_final$ess_input,
      alpha           = alpha,
      alpha_adj       = fit_final$alpha_adj,
      comp            = as.integer(comp),
      nsim            = as.integer(nsim),
      mc_seed         = mc_seed
    ),
    class = "oda_sample_size"
  )
}

# ---- print.oda_sample_size ---------------------------------------------------

#' @export
print.oda_sample_size <- function(x, digits = 3L, ...) {
  cat("ODA Sample Size Analysis\n")
  cat(sprintf("  Effect:  p1 = %.4f, p2 = %.4f", x$p1, x$p2))
  if (!is.na(x$ess_input))
    cat(sprintf("  (ESS = %.2f%%)", x$ess_input))
  cat("\n")
  if (x$comp > 1L) {
    cat(sprintf("  alpha = %.4f  (Sidak-adjusted for %d comparisons: %.6f)\n",
                x$alpha, x$comp, x$alpha_adj))
  } else {
    cat(sprintf("  alpha = %.4f\n", x$alpha))
  }
  cat(sprintf("  Target power:    %.3f\n", x$power_target))
  cat(sprintf("  Minimum n/group: %d\n",   x$n))
  cat(sprintf("  Power achieved:  %.3f\n", round(x$power_achieved, digits)))
  cat(sprintf("  nsim = %d\n", x$nsim))
  invisible(x)
}
