###############################################################################
# tests/testthat/helper-oda.R
# Loaded automatically by testthat before every test file.
#
# Strategy: always make harness functions available to tests by defining them
# here directly. The oda package is loaded via the normal test mechanism;
# these helpers fill in the unpublished test-support functions.
###############################################################################

# ---- Ensure package functions are available ---------------------------------
# devtools::test() loads the package; R CMD check installs then loads it.
# Either way, oda_* functions come from the package namespace.

# ---- Test support functions (not exported from package) ---------------------
# These duplicate the harness_utils.R implementations so tests work under
# R CMD check without requiring harness exports.

confusion_raw <- function(y, y_pred, C) {
  y <- as.integer(y); y_pred <- as.integer(y_pred)
  m <- matrix(0L, nrow = as.integer(C), ncol = as.integer(C))
  for (i in seq_along(y)) {
    yi <- y[i]; pi <- y_pred[i]
    if (!is.na(yi) && !is.na(pi) && yi >= 1L && yi <= C && pi >= 1L && pi <= C)
      m[yi, pi] <- m[yi, pi] + 1L
  }
  m
}

pac_overall <- function(conf) {
  100 * sum(diag(conf)) / sum(conf)
}

ess_pac_from_conf <- function(conf) {
  C  <- nrow(conf)
  mp <- mean(diag(conf) / rowSums(conf))
  100 * (mp - 1/C) / (1 - 1/C)
}

same_int_mat <- function(a, b) {
  if (is.null(a) || is.null(b)) return(FALSE)
  a <- unname(a); b <- unname(b)
  storage.mode(a) <- "integer"; storage.mode(b) <- "integer"
  identical(a, b)
}

harness_loo_refit_ordered_raw <- function(
    x, y, priors_on, degen, K_segments,
    boundary_mode = "right_closed"
) {
  n     <- length(y); y <- as.integer(y)
  y_pred <- rep(NA_integer_, n)
  for (i in seq_len(n)) {
    keep  <- seq_len(n) != i
    fit_i <- tryCatch(
      oda_multiclass_unioda_core(
        x = x[keep], y = y[keep], w = NULL,
        attr_type  = "ordered", priors_on = priors_on,
        K_segments = as.integer(K_segments), degen = degen,
        mcarlo = FALSE, loo = "off", boundary_mode = boundary_mode
      ),
      error = function(e) list(ok = FALSE, reason = conditionMessage(e))
    )
    if (!isTRUE(fit_i$ok))
      return(list(ok = FALSE, reason = paste0("fold_", i, "_", fit_i$reason)))
    y_pred[i] <- as.integer(
      oda_rule_predict_multiclass(x[i], fit_i$rule, boundary = boundary_mode))
  }
  if (anyNA(y_pred)) return(list(ok = FALSE, reason = "NA_predictions"))
  C    <- length(sort(unique(y)))
  conf <- confusion_raw(y, y_pred, C)
  list(ok = TRUE, confusion = conf, y_pred = y_pred)
}

.make_loo_gate_data <- function() {
  add <- function(cls, trap, stable, n)
    data.frame(Class = rep(cls, n), Trap = rep(trap, n), Stable = rep(stable, n))
  rbind(
    add(1L,  6, 0, 1), add(2L, 1, 0, 4), add(2L, 2, 0, 4),
    add(2L,  3, 0, 4), add(2L, 4, 0, 4), add(2L, 5, 0, 3),
    add(2L,  6, 0, 3),
    add(1L,  7, 0, 1), add(2L, 7, 0, 6),
    add(1L,  8, 1, 1),
    add(1L,  9, 1, 2), add(2L,  9, 0, 4),
    add(1L, 10, 0, 2), add(2L, 10, 0, 4)
  )
}

rule_is_nondegenerate <- function(rule, C) {
  if (is.null(rule) || is.null(rule$type)) return(FALSE)
  if (identical(rule$type, "multiclass_ordered"))
    return(length(unique(as.integer(rule$seg_classes))) == C)
  if (identical(rule$type, "multiclass_nominal"))
    return(length(unique(as.integer(rule$level_class))) == C)
  FALSE
}
