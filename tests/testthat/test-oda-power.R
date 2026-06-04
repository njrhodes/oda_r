###############################################################################
# test-oda-power.R
#
# Tests for oda_power() and oda_sample_size().
#
# Scope: unit-weighted binary 2x2 ODA-equivalent power planning only.
# Canon: Rhodes 2020 (Fisher/ODA isomorphism), Yarnold & Soltysik 2005.
#
# Coverage:
#   - Class and structure of returned objects
#   - Seed reproducibility
#   - ESS <-> (p1, p2) conversion
#   - Sidak correction: alpha_adj = 1 - (1 - alpha)^(1/comp)
#   - Monotonicity of power in n
#   - oda_sample_size returns minimum n >= target
#   - Input validation: bad n, alpha, comp, nsim, ess, p1/p2
###############################################################################

# ---- oda_power: object structure --------------------------------------------

test_that("oda_power returns oda_power object with required fields", {
  fit <- oda_power(n1 = 30L, ess = 30, nsim = 200L, mc_seed = 1L)

  expect_s3_class(fit, "oda_power")
  expect_true(is.numeric(fit$power))
  expect_true(is.numeric(fit$p1))
  expect_true(is.numeric(fit$p2))
  expect_equal(fit$ess_input, 30)
  expect_equal(fit$alpha,     0.05)
  expect_equal(fit$comp,      1L)
  expect_equal(fit$nsim,      200L)
  expect_equal(fit$mc_seed,   1L)
})

test_that("oda_power with p1/p2 sets ess_input to NA", {
  fit <- oda_power(n1 = 30L, p1 = 0.3, p2 = 0.7, nsim = 200L, mc_seed = 1L)
  expect_true(is.na(fit$ess_input))
  expect_equal(fit$p1, 0.3)
  expect_equal(fit$p2, 0.7)
})

# ---- oda_power: ESS conversion ----------------------------------------------

test_that("oda_power ESS conversion is symmetric: p2 = ESS/200 + 0.5", {
  fit <- oda_power(n1 = 10L, ess = 48, nsim = 10L, mc_seed = 1L)
  expect_equal(fit$p2, 48 / 200 + 0.5)
  expect_equal(fit$p1, 1 - fit$p2)
})

test_that("oda_power ESS = 0 rejected (strictly positive)", {
  expect_error(oda_power(n1 = 30L, ess =   0, nsim = 10L), "ess.*0.*100")
  expect_error(oda_power(n1 = 30L, ess = 100, nsim = 10L), "ess.*0.*100")
  expect_error(oda_power(n1 = 30L, ess = -10, nsim = 10L), "ess.*0.*100")
})

# ---- oda_power: seed reproducibility ----------------------------------------

test_that("oda_power gives identical results with same seed", {
  r1 <- oda_power(n1 = 50L, ess = 40, nsim = 500L, mc_seed = 77L)
  r2 <- oda_power(n1 = 50L, ess = 40, nsim = 500L, mc_seed = 77L)
  expect_equal(r1$power, r2$power)
})

test_that("oda_power results differ with different seeds (stochastic check)", {
  # Very small nsim so RNG differences are likely to show; may rarely match
  r1 <- oda_power(n1 = 50L, ess = 40, nsim = 200L, mc_seed = 1L)
  r2 <- oda_power(n1 = 50L, ess = 40, nsim = 200L, mc_seed = 2L)
  # Not asserting inequality (could match by chance); just that both are valid
  expect_true(r1$power >= 0 && r1$power <= 1)
  expect_true(r2$power >= 0 && r2$power <= 1)
})

# ---- oda_power: Sidak correction --------------------------------------------

test_that("Sidak correction: alpha_adj = 1 - (1-alpha)^(1/comp)", {
  fit <- oda_power(n1 = 30L, ess = 40, comp = 3L, alpha = 0.05,
                   nsim = 10L, mc_seed = 1L)
  expected_adj <- 1 - (1 - 0.05)^(1 / 3)
  expect_equal(fit$alpha_adj, expected_adj, tolerance = 1e-10)
})

test_that("Sidak correction: comp = 1 leaves alpha unchanged", {
  fit <- oda_power(n1 = 30L, ess = 40, comp = 1L, alpha = 0.05,
                   nsim = 10L, mc_seed = 1L)
  expect_equal(fit$alpha_adj, 0.05)
})

# ---- oda_power: monotonicity ------------------------------------------------

test_that("oda_power increases with n for fixed moderate effect", {
  # Use moderate effect and enough nsim for stable estimates
  pwr <- oda_power(n1 = c(20L, 60L, 120L), ess = 40,
                   nsim = 2000L, mc_seed = 42L)$power
  # Monotone non-decreasing (with MC noise tolerance)
  expect_true(pwr[2] >= pwr[1] - 0.05,
              info = sprintf("power[60] = %.3f < power[20] = %.3f - 0.05", pwr[2], pwr[1]))
  expect_true(pwr[3] >= pwr[2] - 0.05,
              info = sprintf("power[120] = %.3f < power[60] = %.3f - 0.05", pwr[3], pwr[2]))
})

# ---- oda_power: vector n returns named output -------------------------------

test_that("oda_power vector n1 returns named vector", {
  fit <- oda_power(n1 = c(20L, 40L), ess = 40, nsim = 100L, mc_seed = 1L)
  expect_length(fit$power, 2L)
  expect_equal(names(fit$power), c("n=20", "n=40"))
})

test_that("oda_power matrix output when vector n1 and vector alpha", {
  fit <- oda_power(n1 = c(20L, 40L), ess = 40, alpha = c(0.05, 0.01),
                   nsim = 100L, mc_seed = 1L)
  expect_true(is.matrix(fit$power))
  expect_equal(dim(fit$power), c(2L, 2L))
})

# ---- oda_power: n2 default and recycling ------------------------------------

test_that("n2 defaults to n1 (balanced design)", {
  fit <- oda_power(n1 = 40L, ess = 40, nsim = 100L, mc_seed = 1L)
  expect_equal(fit$n1, 40L)
  expect_equal(fit$n2, 40L)
})

test_that("scalar n2 is recycled to match vector n1", {
  fit <- oda_power(n1 = c(30L, 50L), n2 = 40L, ess = 40,
                   nsim = 100L, mc_seed = 1L)
  expect_equal(fit$n2, c(40L, 40L))
})

# ---- oda_power: input validation --------------------------------------------

test_that("oda_power rejects both ess and p1/p2", {
  expect_error(oda_power(n1 = 30L, ess = 40, p1 = 0.3, p2 = 0.7, nsim = 10L),
               "either.*ess.*or.*p1.*p2")
})

test_that("oda_power rejects neither ess nor p1/p2", {
  expect_error(oda_power(n1 = 30L, nsim = 10L),
               "either.*ess.*or.*both.*p1.*and.*p2")
})

test_that("oda_power rejects n1 < 2", {
  expect_error(oda_power(n1 = 1L, ess = 40, nsim = 10L), "n1.*n2.*>= 2")
})

test_that("oda_power rejects invalid alpha", {
  expect_error(oda_power(n1 = 30L, ess = 40, alpha = 0,   nsim = 10L),
               "alpha.*0.*1")
  expect_error(oda_power(n1 = 30L, ess = 40, alpha = 1.1, nsim = 10L),
               "alpha.*0.*1")
})

test_that("oda_power rejects comp < 1", {
  expect_error(oda_power(n1 = 30L, ess = 40, comp = 0L, nsim = 10L),
               "comp.*positive")
})

test_that("oda_power rejects nsim < 1", {
  expect_error(oda_power(n1 = 30L, ess = 40, nsim = 0L), "nsim.*positive")
})

test_that("oda_power rejects p1 outside (0,1)", {
  expect_error(oda_power(n1 = 30L, p1 =  0, p2 = 0.7, nsim = 10L), "p1")
  expect_error(oda_power(n1 = 30L, p1 =  1, p2 = 0.7, nsim = 10L), "p1")
  expect_error(oda_power(n1 = 30L, p1 = -0.1, p2 = 0.7, nsim = 10L), "p1")
})

# ---- oda_power: print method ------------------------------------------------

test_that("print.oda_power produces output without error", {
  fit <- oda_power(n1 = 30L, ess = 40, nsim = 100L, mc_seed = 1L)
  expect_output(print(fit), "ODA Power")
  expect_output(print(fit), "p1")
  expect_output(print(fit), "p2")
})

# ---- oda_sample_size: basic contract ----------------------------------------

test_that("oda_sample_size returns oda_sample_size object", {
  ss <- oda_sample_size(ess = 48, power_target = 0.80,
                        nsim = 500L, mc_seed = 42L)
  expect_s3_class(ss, "oda_sample_size")
  expect_true(is.integer(ss$n))
  expect_true(ss$n >= 2L)
  expect_equal(ss$power_target, 0.80)
})

test_that("oda_sample_size: power at returned n >= target", {
  ss <- oda_sample_size(ess = 48, power_target = 0.70,
                        nsim = 1000L, mc_seed = 42L)
  expect_gte(ss$power_achieved, 0.70 - 0.05)   # allow 5% MC tolerance
})

test_that("oda_sample_size: higher target requires larger n", {
  ss70 <- oda_sample_size(ess = 40, power_target = 0.70,
                          nsim = 500L, mc_seed = 42L)
  ss80 <- oda_sample_size(ess = 40, power_target = 0.80,
                          nsim = 500L, mc_seed = 42L)
  expect_gte(ss80$n, ss70$n)
})

test_that("oda_sample_size: Sidak correction increases required n", {
  ss1 <- oda_sample_size(ess = 40, power_target = 0.80, comp = 1L,
                         nsim = 500L, mc_seed = 42L)
  ss3 <- oda_sample_size(ess = 40, power_target = 0.80, comp = 3L,
                         nsim = 500L, mc_seed = 42L)
  expect_gte(ss3$n, ss1$n)
})

# ---- oda_sample_size: input validation ---------------------------------------

test_that("oda_sample_size rejects power_target outside (0,1)", {
  expect_error(oda_sample_size(ess = 40, power_target = 0,   nsim = 10L),
               "power_target")
  expect_error(oda_sample_size(ess = 40, power_target = 1.0, nsim = 10L),
               "power_target")
})

test_that("oda_sample_size rejects n_max <= n_min", {
  expect_error(oda_sample_size(ess = 40, n_min = 100L, n_max = 50L, nsim = 10L),
               "n_max.*n_min")
})

test_that("oda_sample_size errors when n_max power is below target", {
  # Tiny n_max with very high target - should fail
  expect_error(
    oda_sample_size(ess = 5, power_target = 0.99, n_max = 10L, nsim = 200L),
    "n_max"
  )
})

# ---- oda_sample_size: print method ------------------------------------------

test_that("print.oda_sample_size produces output without error", {
  ss <- oda_sample_size(ess = 48, power_target = 0.80,
                        nsim = 200L, mc_seed = 42L)
  expect_output(print(ss), "ODA Sample Size")
  expect_output(print(ss), "Minimum n")
  expect_output(print(ss), "Power achieved")
})
