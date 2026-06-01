# test-oda-mc-fast-scorer.R
#
# Acceptance tests for the fast MC permutation scorer introduced in
# oda_mc_p_value():
#   .oda_mc_precomp()     -- precompute x/w structure once
#   .oda_fast_ess_perm()  -- per-permutation max-ESS scalar
#
# These tests verify:
#   1. Scalar parity: fast scorer matches oda_univariate_core ESS to 1e-10
#      for binary and ordered attributes (unweighted + weighted, priors on/off).
#   2. MC ge_count parity: oda_mc_p_value() with fast path gives identical
#      ge_count/iter_used to a manual slow loop using oda_univariate_core.
#   3. Performance (smoke tier only): fast path is materially faster.
#
# EXE myeloma parity (MODEL1/MODEL30/MODEL56) is covered by test-cta.R at
# smoke tier and is not duplicated here.

# --------------------------------------------------------------------------- #
# Helper: slow MC loop (replicates oda_mc_p_value before the fast path).      #
# --------------------------------------------------------------------------- #
.mc_slow_ref <- function(x, y, w, attr_type, priors_on, primary = NULL,
                         secondary = NULL, miss_codes = NULL,
                         chance_model = "class", mc_iter, mc_target = 0.05,
                         mc_stop = 99.9, mc_stopup = NA, seed, ess_obs,
                         direction = "off") {
  y <- as.integer(y)
  n <- length(y)
  if (is.null(w)) w <- rep(1, n) else w <- as.numeric(w)
  set.seed(seed)
  conf_stop   <- if (!is.na(mc_stop)   && mc_stop   > 1) mc_stop   / 100 else mc_stop
  conf_stopup <- if (!is.na(mc_stopup) && mc_stopup > 1) mc_stopup / 100 else mc_stopup
  ge_count <- 0L; iter_used <- 0L
  min_check <- 50L; check_every <- 50L
  for (b in seq_len(mc_iter)) {
    iter_used <- b
    y_star <- sample(y, n, replace = FALSE)
    fit_b  <- oda_univariate_core(
      x = x, y = y_star, w = w,
      attr_type    = attr_type,
      priors_on    = priors_on,
      primary      = primary,
      secondary    = secondary,
      miss_codes   = miss_codes,
      loo          = "off",
      mcarlo       = FALSE,
      mc_iter      = 0L,
      mc_target    = mc_target,
      mc_stop      = mc_stop,
      mc_stopup    = mc_stopup,
      mc_adjust    = FALSE,
      mc_seed      = NULL,
      chance_model = chance_model,
      direction    = direction
    )
    ess_b <- if (!isTRUE(fit_b$ok)) 0 else fit_b$ess
    if (ess_b >= ess_obs - 1e-12) ge_count <- ge_count + 1L
    if (b >= min_check && (b %% check_every == 0L)) {
      if (ge_count == 0L) {
        lower <- 0
        upper <- if (!is.na(conf_stop))
          stats::qbeta(conf_stop, 1, b) else NA_real_
      } else if (ge_count == b) {
        upper <- 1
        lower <- if (!is.na(conf_stopup))
          stats::qbeta(1 - conf_stopup, b, 1) else NA_real_
      } else {
        upper <- if (!is.na(conf_stop))
          stats::qbeta(conf_stop, ge_count + 1, b - ge_count) else NA_real_
        lower <- if (!is.na(conf_stopup))
          stats::qbeta(1 - conf_stopup, ge_count, b - ge_count + 1) else NA_real_
      }
      if (!is.na(mc_target) && !is.na(upper) && upper < mc_target) break
      if (!is.na(mc_target) && !is.na(lower) && lower > mc_target) break
    }
  }
  list(ge_count  = ge_count,
       iter_used = iter_used,
       p_mc      = if (iter_used <= 0L) NA_real_
                   else if (ge_count == 0L) 0.0
                   else ge_count / iter_used)
}

# --------------------------------------------------------------------------- #
# Test 1 -- Scalar parity: binary attribute                                    #
# --------------------------------------------------------------------------- #
test_that("fast scorer matches oda_univariate_core ESS for binary x (unweighted)", {
  set.seed(1L)
  n  <- 40L
  x  <- sample(c(0L, 1L), n, replace = TRUE)
  y  <- sample(c(0L, 1L), n, replace = TRUE)
  w  <- rep(1, n)

  pc <- odacore:::.oda_mc_precomp(x, w, "binary", mindenom = 1L)
  expect_false(is.null(pc), label = "precomp not NULL for binary")

  set.seed(99L)
  for (i in seq_len(200L)) {
    y_perm <- sample(y, n, replace = FALSE)
    fast   <- odacore:::.oda_fast_ess_perm(y_perm, pc, "off")
    slow   <- {
      fit <- oda_univariate_core(x = x, y = y_perm, w = w,
               attr_type = "binary", priors_on = TRUE,
               primary = NULL, secondary = NULL,
               loo = "off", mcarlo = FALSE, mc_iter = 0L)
      if (!isTRUE(fit$ok)) 0 else fit$ess
    }
    expect_equal(fast, slow, tolerance = 1e-10,
      label = paste("binary unweighted perm", i))
  }
})

test_that("fast scorer matches oda_univariate_core ESS for binary x (weighted, priors on)", {
  set.seed(2L)
  n  <- 50L
  x  <- sample(c(0L, 1L), n, replace = TRUE)
  y  <- sample(c(0L, 1L), n, replace = TRUE)
  w  <- runif(n, 0.5, 2.5)

  pc <- odacore:::.oda_mc_precomp(x, w, "binary", mindenom = 1L)
  expect_false(is.null(pc))

  set.seed(77L)
  for (i in seq_len(200L)) {
    y_perm <- sample(y, n, replace = FALSE)
    fast   <- odacore:::.oda_fast_ess_perm(y_perm, pc, "off")
    slow   <- {
      fit <- oda_univariate_core(x = x, y = y_perm, w = w,
               attr_type = "binary", priors_on = TRUE,
               primary = NULL, secondary = NULL,
               loo = "off", mcarlo = FALSE, mc_iter = 0L)
      if (!isTRUE(fit$ok)) 0 else fit$ess
    }
    expect_equal(fast, slow, tolerance = 1e-10,
      label = paste("binary weighted perm", i))
  }
})

# --------------------------------------------------------------------------- #
# Test 2 -- Scalar parity: ordered attribute (multi-cut)                      #
# --------------------------------------------------------------------------- #
test_that("fast scorer matches oda_univariate_core ESS for ordered x (unweighted)", {
  set.seed(3L)
  n <- 60L
  x <- rnorm(n)
  y <- sample(c(0L, 1L), n, replace = TRUE)
  w <- rep(1, n)

  pc <- odacore:::.oda_mc_precomp(x, w, "ordered", mindenom = 1L)
  expect_false(is.null(pc))
  expect_gt(length(pc$adm_j), 1L)  # multiple cut positions

  set.seed(55L)
  for (i in seq_len(200L)) {
    y_perm <- sample(y, n, replace = FALSE)
    fast   <- odacore:::.oda_fast_ess_perm(y_perm, pc, "off")
    slow   <- {
      fit <- oda_univariate_core(x = x, y = y_perm, w = w,
               attr_type = "ordered", priors_on = TRUE,
               primary = NULL, secondary = NULL,
               loo = "off", mcarlo = FALSE, mc_iter = 0L)
      if (!isTRUE(fit$ok)) 0 else fit$ess
    }
    expect_equal(fast, slow, tolerance = 1e-10,
      label = paste("ordered unweighted perm", i))
  }
})

test_that("fast scorer matches oda_univariate_core ESS for ordered x (weighted, priors on/off)", {
  set.seed(4L)
  n <- 80L
  x <- round(rnorm(n), 1)   # ties present
  y <- sample(c(0L, 1L), n, replace = TRUE)
  w <- runif(n, 0.1, 3.0)

  for (priors in c(TRUE, FALSE)) {
    pc <- odacore:::.oda_mc_precomp(x, w, "ordered", mindenom = 1L)
    expect_false(is.null(pc))
    set.seed(33L)
    for (i in seq_len(200L)) {
      y_perm <- sample(y, n, replace = FALSE)
      fast   <- odacore:::.oda_fast_ess_perm(y_perm, pc, "off")
      slow   <- {
        fit <- oda_univariate_core(x = x, y = y_perm, w = w,
                 attr_type = "ordered", priors_on = priors,
                 primary = NULL, secondary = NULL,
                 loo = "off", mcarlo = FALSE, mc_iter = 0L)
        if (!isTRUE(fit$ok)) 0 else fit$ess
      }
      expect_equal(fast, slow, tolerance = 1e-10,
        label = paste("ordered weighted priors_on=", priors, "perm", i))
    }
  }
})

# --------------------------------------------------------------------------- #
# Test 3 -- MC ge_count parity: fast path == slow path for same seed          #
# --------------------------------------------------------------------------- #
test_that("oda_mc_p_value fast path gives identical ge_count/iter_used to slow reference", {
  # Binary weighted case (exercises the fast binary path)
  set.seed(7L)
  n  <- 60L
  x  <- sample(c(0L, 1L), n, replace = TRUE)
  y  <- sample(c(0L, 1L), n, replace = TRUE)
  w  <- runif(n, 0.5, 2.0)

  obs_fit <- oda_univariate_core(x = x, y = y, w = w, attr_type = "binary",
               priors_on = TRUE, loo = "off", mcarlo = FALSE, mc_iter = 0L)
  expect_true(obs_fit$ok)

  SEED <- 42L
  fast_res <- oda_mc_p_value(
    x = x, y = y, w = w,
    attr_type    = "binary",
    priors_on    = TRUE,
    primary      = NULL,
    secondary    = NULL,
    chance_model = "class",
    mc_iter      = 500L,
    mc_target    = 0.05,
    mc_stop      = 99.9,
    mc_stopup    = NA,
    seed         = SEED,
    ess_obs      = obs_fit$ess,
    direction    = "off"
  )
  slow_res <- .mc_slow_ref(
    x = x, y = y, w = w, attr_type = "binary", priors_on = TRUE,
    chance_model = "class", mc_iter = 500L, mc_target = 0.05,
    mc_stop = 99.9, mc_stopup = NA, seed = SEED,
    ess_obs = obs_fit$ess, direction = "off"
  )
  expect_identical(fast_res$ge_count,  slow_res$ge_count,
    label = "binary weighted: ge_count")
  expect_identical(fast_res$iter_used, slow_res$iter_used,
    label = "binary weighted: iter_used")

  # Ordered weighted case
  set.seed(8L)
  n2 <- 80L
  x2 <- rnorm(n2)
  y2 <- sample(c(0L, 1L), n2, replace = TRUE)
  w2 <- runif(n2, 0.5, 2.0)

  obs_fit2 <- oda_univariate_core(x = x2, y = y2, w = w2, attr_type = "ordered",
                priors_on = TRUE, loo = "off", mcarlo = FALSE, mc_iter = 0L)
  expect_true(obs_fit2$ok)

  fast_res2 <- oda_mc_p_value(
    x = x2, y = y2, w = w2,
    attr_type    = "ordered",
    priors_on    = TRUE,
    primary      = NULL,
    secondary    = NULL,
    chance_model = "class",
    mc_iter      = 500L,
    mc_target    = 0.05,
    mc_stop      = 99.9,
    mc_stopup    = NA,
    seed         = SEED,
    ess_obs      = obs_fit2$ess,
    direction    = "off"
  )
  slow_res2 <- .mc_slow_ref(
    x = x2, y = y2, w = w2, attr_type = "ordered", priors_on = TRUE,
    chance_model = "class", mc_iter = 500L, mc_target = 0.05,
    mc_stop = 99.9, mc_stopup = NA, seed = SEED,
    ess_obs = obs_fit2$ess, direction = "off"
  )
  expect_identical(fast_res2$ge_count,  slow_res2$ge_count,
    label = "ordered weighted: ge_count")
  expect_identical(fast_res2$iter_used, slow_res2$iter_used,
    label = "ordered weighted: iter_used")
})

# --------------------------------------------------------------------------- #
# Test 4 -- Performance: fast path materially faster (smoke tier)             #
# --------------------------------------------------------------------------- #
test_that("oda_mc_p_value fast path is materially faster than slow reference (smoke)", {
  skip_if_not_smoke("fast-scorer performance test")

  # Use myeloma data to reproduce the actual bottleneck scenario.
  # data.txt is space-delimited, no header; columns V1..V19.
  # V1 = class (0/1), V2 = case weight, V14 = binary predictor (0/1/-9).
  myeloma_dir <- test_path("fixtures", "myeloma")
  data_path   <- file.path(myeloma_dir, "data.txt")
  skip_if_not(file.exists(data_path), "myeloma fixture missing")

  dat <- read.table(data_path, header = FALSE)
  names(dat) <- paste0("V", seq_len(ncol(dat)))
  dat <- dat[dat$V2 != 0, ]   # EX V2=0

  # V14 is binary (0/1 after -9 removal) -- a primary bottleneck attribute
  miss9 <- -9
  x_v14 <- dat$V14
  y_v1  <- as.integer(dat$V1)
  w_v2  <- dat$V2
  # Clean miss codes as oda_univariate_core would
  clean <- odacore:::oda_clean_xy(x_v14, y_v1, w_v2, miss_codes = miss9)
  x_c   <- clean$x; y_c <- as.integer(clean$y); w_c <- clean$w

  obs_fit <- oda_univariate_core(x = x_c, y = y_c, w = w_c,
               attr_type = "binary", priors_on = TRUE,
               loo = "off", mcarlo = FALSE, mc_iter = 0L)
  skip_if_not(isTRUE(obs_fit$ok), "V14 not admissible in test data")

  SEED <- 1L; ITER <- 500L

  t_fast <- system.time(
    oda_mc_p_value(x = x_c, y = y_c, w = w_c,
      attr_type = "binary", priors_on = TRUE,
      primary = NULL, secondary = NULL,
      chance_model = "class",
      mc_iter = ITER, mc_target = 0.05,
      mc_stop = 99.9, mc_stopup = NA,
      seed = SEED, ess_obs = obs_fit$ess, direction = "off")
  )["elapsed"]

  t_slow <- system.time(
    .mc_slow_ref(x = x_c, y = y_c, w = w_c, attr_type = "binary",
      priors_on = TRUE, chance_model = "class",
      mc_iter = ITER, mc_target = 0.05, mc_stop = 99.9,
      mc_stopup = NA, seed = SEED, ess_obs = obs_fit$ess, direction = "off")
  )["elapsed"]

  message(sprintf(
    "[perf] binary V14 (n=%d, iter=%d): fast=%.3fs  slow=%.3fs  speedup=%.1fx",
    length(x_c), ITER, t_fast, t_slow, t_slow / max(t_fast, 0.001)))

  # Assert at least 3x speedup (typically 5-10x expected).
  expect_gt(t_slow / max(t_fast, 0.001), 3,
    label = "fast path should be at least 3x faster than slow reference")
})
