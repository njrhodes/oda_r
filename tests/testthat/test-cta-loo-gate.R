###############################################################################
# test-cta-loo-gate.R
#
# Regression test: CTA LOO gate correctly rejects an ordered-cut attribute
# with higher training ESS but LOO instability in favour of a lower-ESS
# LOO-stable attribute.
#
# Synthetic dataset (n=43, binary class):
#   Trap   (V2): ordered, ESS~49.21% at cut<=7.5, LOO UNSTABLE
#                (tie boundary at Trap=7; dropping a class-2 obs there
#                 shifts the winning cut from 7.5 to 6.5, changing
#                 the held-out prediction)
#   Stable (V3): binary,  ESS~42.86% at cut<=0.5, LOO STABLE
#
# CTA.exe canonical output (ENUMERATE; LOO STABLE; MC ITER 5000):
#   V3  Node 1  Tree 3  ESS 42.86%  Best ESS 42.86%
#   INCLUDE V3>0.5 (i.e. cut_value=0.5, Stable=1 → class 1)
#   Selected None (no further splits)
#
# R result must agree exactly.
###############################################################################

.make_loo_gate_data <- function() {
  add <- function(cls, trap, stable, n)
    data.frame(Class = rep(cls, n), Trap = rep(trap, n), Stable = rep(stable, n))

  rbind(
    # Trap <= 6.5 region: c1=1, c2=22
    add(1L,  6, 0, 1), add(2L, 1, 0, 4), add(2L, 2, 0, 4),
    add(2L,  3, 0, 4), add(2L, 4, 0, 4), add(2L, 5, 0, 3),
    add(2L,  6, 0, 3),
    # Trap=7: tie-boundary bin that makes Trap LOO unstable (c1=1, c2=6)
    add(1L,  7, 0, 1), add(2L, 7, 0, 6),
    # Trap=8: c1=1, c2=0, Stable=1
    add(1L,  8, 1, 1),
    # Trap>8: c1=4, c2=8
    add(1L,  9, 1, 2), add(2L,  9, 0, 4),
    add(1L, 10, 0, 2), add(2L, 10, 0, 4)
  )
}

.loo_gate_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      d <- .make_loo_gate_data()
      fit <<- cta_fit(
        X           = d[, c("Trap", "Stable")],
        y           = d$Class,
        priors_on   = TRUE,
        alpha_split = 0.05,
        mindenom    = 1L,
        prune_alpha = 0.05,
        max_depth   = 3L,
        ess_min     = 0,
        mc_iter     = 5000L,
        mc_target   = 0.05,
        mc_stop     = 99.9,
        mc_stopup   = 99.9,
        mc_seed     = NULL,
        loo         = "stable",
        attr_names  = c("Trap", "Stable")
      )
    }
    fit
  }
})

# ---- root selection ---------------------------------------------------------

test_that("loo-gate: root is Stable, not LOO-unstable Trap", {
  skip_if_not_smoke("cta-loo-gate")
  tree <- .loo_gate_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_false(isTRUE(root$leaf))
  expect_equal(root$attribute, "Stable")
})

test_that("loo-gate: root cut value is 0.5 (matches CTA.exe V3>0.5)", {
  skip_if_not_smoke("cta-loo-gate")
  tree <- .loo_gate_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$rule$cut_value, 0.5, tolerance = 1e-9)
})

test_that("loo-gate: root ESS ~42.86% (matches CTA.exe)", {
  skip_if_not_smoke("cta-loo-gate")
  tree <- .loo_gate_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(round(root$ess, 2), 42.86)
})

test_that("loo-gate: root LOO status is STABLE", {
  skip_if_not_smoke("cta-loo-gate")
  tree <- .loo_gate_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$loo_status, "STABLE")
})

# ---- tree structure (CTA.exe: "Selected None" = no further splits) ----------

test_that("loo-gate: only 1 split node — no further splits below root", {
  skip_if_not_smoke("cta-loo-gate")
  tree <- .loo_gate_fit()
  tbl  <- cta_node_table(tree)
  expect_equal(sum(!tbl$leaf), 1L)
})

test_that("loo-gate: root n_obs = 43", {
  skip_if_not_smoke("cta-loo-gate")
  tree <- .loo_gate_fit()
  root <- tree$nodes[[tree$root_id]]
  expect_equal(root$n_obs, 43L)
})
