###############################################################################
# R/oda_s3.R  -  S3 methods and accessors for oda_fit objects
#
# Classes:
#   oda_fit             -  base class (all oda_fit() results)
#   oda_fit_binary      -  binary-class fit (engine = "binary")
#   oda_fit_multiclass  -  multiclass fit   (engine = "multiclass")
#   oda_fit_failed      -  failed/degenerate fit (ok = FALSE)
#   oda_fit_summary     -  output of summary.oda_fit()
#
# Public API added here:
#   predict.oda_fit()
#   print.oda_fit()
#   summary.oda_fit()
#   print.oda_fit_summary()
#   oda_predictions()
#   oda_confusion()
#   oda_metrics()
###############################################################################

# ---- Internal helpers ------------------------------------------------------- #

# Format a rule object as a human-readable string.
.oda_fmt_rule <- function(rule) {
  if (is.null(rule)) return("<no rule>")
  tp <- rule$type %||% ""

  if (tp == "ordered_cut") {
    cut <- rule$cut_value
    dir <- rule$direction %||% "?"
    l0  <- rule$label_0 %||% "0"
    l1  <- rule$label_1 %||% "1"
    if (dir == "0->1")
      sprintf("<= %.4g --> %s   |   > %.4g --> %s", cut, l0, cut, l1)
    else
      sprintf("<= %.4g --> %s   |   > %.4g --> %s", cut, l1, cut, l0)

  } else if (tp == "multiclass_ordered") {
    cuts <- rule$cut_values
    segs <- rule$seg_classes
    K    <- length(segs)
    parts <- character(K)
    for (i in seq_len(K)) {
      lo <- if (i == 1L) -Inf else cuts[i - 1L]
      hi <- if (i == K)   Inf else cuts[i]
      if (is.infinite(lo) && is.infinite(hi))
        parts[i] <- sprintf("all --> %d", segs[i])
      else if (is.infinite(lo))
        parts[i] <- sprintf("<= %.4g --> %d", hi, segs[i])
      else if (is.infinite(hi))
        parts[i] <- sprintf("> %.4g --> %d", lo, segs[i])
      else
        parts[i] <- sprintf("(%.4g, %.4g] --> %d", lo, hi, segs[i])
    }
    paste(parts, collapse = "   |   ")

  } else if (tp %in% c("binary_map", "nominal_cut")) {
    "<categorical/binary rule>"

  } else if (tp == "multiclass_nominal") {
    levs <- rule$levels %||% character(0)
    cls  <- rule$level_class %||% integer(0)
    if (length(levs) == 0L) return("<nominal multiclass rule>")
    parts <- sprintf("%s --> %d", levs, cls)
    paste(parts, collapse = "   |   ")

  } else {
    sprintf("<rule type=%s>", tp)
  }
}

# Format a p-value for display.
.fmt_p <- function(p) {
  if (is.null(p) || is.na(p)) return("NA")
  if (p < 0.001) return("< .001")
  sprintf("%.3f", p)
}

# Extract numeric value safely.
.num <- function(x) if (is.null(x) || !is.finite(x)) NA_real_ else as.numeric(x)

# Build the four-field LOO p-value info block.
#
# Canon (MPE p.34): "Hold-out p is one-tailed: the null hypothesis is that the
# training model will not replicate when it is used to classify observations in
# the hold-out sample."  LOO Fisher alternative is therefore always "greater",
# regardless of whether the training direction was specified or not.
#
# Note: for non-directional analyses, MC p is more conservative than LOO p
# because each MC permutation also searches both directions; LOO Fisher does not
# adjust for direction selection.  MC p and LOO p legitimately diverge.
#
# Rules:
#   Binary 2x2 LOO, p_value stored and not NA:
#     p_method="Fisher exact (2x2), one-tailed [MPE p.34]", p_status="computed".
#   Binary 2x2 LOO, p_value absent or NA:
#     p_status="not_computed", p_reason="...not available in fit object".
#   Multiclass/polychotomous:
#     p_value=NA, p_method="none", p_status="not_computed",
#     p_reason="LOO p-value is not reported for multicategorical/polychotomous ODA".
#
# lo must be the $loo list from fit (caller confirms allowed=TRUE before calling).
.loo_p_info <- function(fit, lo) {
  if (inherits(fit, "oda_fit_binary")) {
    pv <- lo$p_value
    if (!is.null(pv) && !is.na(pv)) {
      list(
        p_value  = pv,
        p_method = "Fisher exact (2x2), one-tailed; MPE p.34",
        p_status = "computed",
        p_reason = NA_character_
      )
    } else {
      list(
        p_value  = NA_real_,
        p_method = "Fisher exact (2x2), one-tailed; MPE p.34",
        p_status = "not_computed",
        p_reason = "2x2 Fisher exact LOO p-value not available in fit object"
      )
    }
  } else {
    list(
      p_value  = NA_real_,
      p_method = "none",
      p_status = "not_computed",
      p_reason = "LOO p-value is not reported for multicategorical/polychotomous ODA"
    )
  }
}

# ---- predict.oda_fit -------------------------------------------------------- #

#' Predict class labels from a fitted ODA model
#'
#' Applies the fitted ODA rule to new attribute values, returning predicted
#' class labels in the original label space.  Missing values and miss-coded
#' values return \code{NA_integer_}.  Failed fits return all \code{NA_integer_}
#' with a warning.
#'
#' @param object An \code{oda_fit} object from \code{\link{oda_fit}}.
#' @param newdata Numeric vector or single-column data frame of attribute values.
#' @param ... Unused.
#' @return Integer vector of predicted class labels, length \code{length(newdata)}
#'   or \code{nrow(newdata)}.
#' @export
predict.oda_fit <- function(object, newdata, ...) {
  # Coerce newdata to numeric vector
  if (is.data.frame(newdata)) {
    if (ncol(newdata) != 1L)
      stop("predict.oda_fit: newdata must be a numeric vector or single-column data frame",
           call. = FALSE)
    x_new <- newdata[[1L]]
  } else {
    x_new <- newdata
  }
  x_new <- as.numeric(x_new)
  n_new <- length(x_new)

  # Failed fit
  if (!isTRUE(object$ok)) {
    warning("predict.oda_fit: fit has ok=FALSE (reason: ",
            object$reason %||% "unknown", "); returning all NA_integer_",
            call. = FALSE)
    return(rep(NA_integer_, n_new))
  }

  # Apply miss_codes masking
  mc <- object$miss_codes
  if (!is.null(mc) && length(mc) > 0L)
    x_new[x_new %in% mc] <- NA_real_

  eng <- object$engine %||% "binary"

  if (eng == "binary") {
    coded <- oda_rule_predict(x_new, object$rule)   # 0/1/NA in coded space
    l0 <- object$rule$label_0
    l1 <- object$rule$label_1
    if (is.null(l0) || is.null(l1)) {
      # No original-label mapping (internal use); return coded directly
      return(as.integer(coded))
    }
    out <- coded
    out[!is.na(coded) & coded == 0L] <- l0
    out[!is.na(coded) & coded == 1L] <- l1
    return(as.integer(out))
  }

  if (eng == "multiclass") {
    bm <- object$boundary_mode %||% "megaoda_halfopen"
    return(as.integer(oda_rule_predict_multiclass(x_new, object$rule,
                                                   boundary = bm)))
  }

  warning("predict.oda_fit: unknown engine '", eng, "'; returning all NA_integer_",
          call. = FALSE)
  rep(NA_integer_, n_new)
}

# ---- print.oda_fit ---------------------------------------------------------- #

#' Print a fitted ODA model
#'
#' Compact display of rule, ESS/Mean PAC, and available MC/LOO metadata.
#' Does not recompute any quantities.
#'
#' @param x An \code{oda_fit} object.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.oda_fit <- function(x, ...) {
  eng   <- x$engine %||% "?"
  atyp  <- x$attr_type %||% "?"
  n_eff <- x$n_eff %||% NA_integer_
  pr    <- if (isTRUE(x$priors_on)) "TRUE" else "FALSE"
  wt    <- if (isTRUE(x$has_weights)) "  weights=TRUE" else ""

  cat(sprintf("\nODA (%s)  attr_type=%s  priors=%s  n=%s%s\n",
              eng, atyp, pr, n_eff, wt))

  if (!isTRUE(x$ok)) {
    if (is.null(x$rule)) {
      cat(sprintf("Status: FAILED  reason=%s\n\n", x$reason %||% "?"))
      return(invisible(x))
    }
    # ok=FALSE but rule present means LOO said not significant/stable.
    # Show the model in full; note LOO result below.
  }

  cat("\nRule: ", .oda_fmt_rule(x$rule), "\n\n", sep = "")

  if (eng == "binary") {
    l0   <- x$rule$label_0 %||% 0L
    l1   <- x$rule$label_1 %||% 1L
    conf <- x$confusion
    if (!is.null(conf)) {
      n0  <- conf$TN + conf$FP
      n1  <- conf$TP + conf$FN
      p0  <- if (!is.null(conf$specificity) && !is.na(conf$specificity))
               sprintf("%.1f%%", conf$specificity * 100) else "NA"
      p1  <- if (!is.null(conf$sensitivity) && !is.na(conf$sensitivity))
               sprintf("%.1f%%", conf$sensitivity * 100) else "NA"
      cat(sprintf("  CLASS  %6s  %6s\n", "n", "PAC"))
      cat(sprintf("  %5s  %6s  %6s\n", l0, n0, p0))
      cat(sprintf("  %5s  %6s  %6s\n", l1, n1, p1))
      cat("\n")
    }
    ess_val  <- .num(x$ess)
    pac_val  <- .num(x$pac)
    ess_str  <- if (!is.na(ess_val))  sprintf("ESS: %.2f%%", ess_val)  else "ESS: NA"
    pac_str  <- if (!is.na(pac_val))  sprintf("Mean PAC: %.2f%%", pac_val) else "Mean PAC: NA"
    p_str    <- if (!is.null(x$p_mc) && !is.na(x$p_mc))
                  sprintf("  p(MC): %s", .fmt_p(x$p_mc)) else ""
    cat(sprintf("  %s   %s%s\n", pac_str, ess_str, p_str))

    lo <- x$loo
    if (!is.null(lo) && isTRUE(lo$allowed)) {
      lc <- lo$confusion
      if (!is.null(lc)) {
        cat("\n  -- LOO --\n")
        # LOO binary confusion stores rates (0-1), not counts.
        # Use training n per class; PAC from LOO rates.
        ln0 <- if (!is.null(conf)) conf$TN + conf$FP else NA_integer_
        ln1 <- if (!is.null(conf)) conf$TP + conf$FN else NA_integer_
        lp0 <- if (!is.null(lc$specificity) && !is.na(lc$specificity))
                 sprintf("%.1f%%", lc$specificity * 100) else "NA"
        lp1 <- if (!is.null(lc$sensitivity) && !is.na(lc$sensitivity))
                 sprintf("%.1f%%", lc$sensitivity * 100) else "NA"
        cat(sprintf("  CLASS  %6s  %6s\n", "n", "PAC"))
        cat(sprintf("  %5s  %6s  %6s\n", l0, ln0, lp0))
        cat(sprintf("  %5s  %6s  %6s\n", l1, ln1, lp1))
        cat("\n")
      }
      ess_loo_str <- if (!is.null(lo$ess_loo) && !is.na(lo$ess_loo))
                       sprintf("LOO ESS: %.2f%%", lo$ess_loo) else ""
      p_loo <- lo$p_value
      p_loo_str <- if (!is.null(p_loo) && !is.na(p_loo))
                     sprintf("  p(LOO): %s", .fmt_p(p_loo)) else ""
      cat(sprintf("  %s%s\n", ess_loo_str, p_loo_str))
    }

  } else {
    # Multiclass
    if (!is.null(x$pac_by_class) && !is.null(x$classes)) {
      cat(sprintf("  CLASS  %6s\n", "PAC"))
      for (i in seq_along(x$classes))
        cat(sprintf("  %5s  %5.1f%%\n", x$classes[i],
                    x$pac_by_class[i] %||% NA_real_))
      cat("\n")
    }
    ess_val <- .num(x$ess)
    mp_val  <- .num(x$mean_pac)
    ess_str <- if (!is.na(ess_val)) sprintf("ESS: %.2f%%", ess_val) else "ESS: NA"
    pac_str <- if (!is.na(mp_val))  sprintf("Mean PAC: %.2f%%", mp_val) else "Mean PAC: NA"
    p_str   <- if (!is.null(x$p_mc) && !is.na(x$p_mc))
                 sprintf("  p(MC): %s", .fmt_p(x$p_mc)) else ""
    cat(sprintf("  %s   %s%s\n", pac_str, ess_str, p_str))

    lo <- x$loo
    if (!is.null(lo) && isTRUE(lo$allowed)) {
      # Categorical LOO stores $confusion; ordered LOO stores $confusion_raw.
      lc <- lo$confusion %||% lo$confusion_raw
      cat("\n  -- LOO --\n")
      if (!is.null(lc)) {
        pb <- lc$pac_by_class %||% NULL
        if (!is.null(pb) && !is.null(x$classes)) {
          cat(sprintf("  CLASS  %6s\n", "PAC"))
          for (i in seq_along(x$classes))
            cat(sprintf("  %5s  %5.1f%%\n", x$classes[i], pb[i] * 100))
          cat("\n")
        }
        mp      <- lc$mean_pac   # proportion scale from oda_confusion_multiclass
        C_loo   <- length(lc$classes %||% x$classes)
        ess_loo_val <- if (!is.null(mp) && !is.na(mp) && C_loo >= 2L)
                         oda_ess_from_mean(mp, C_loo) else NA_real_
        ess_loo_str <- if (!is.na(ess_loo_val))
                         sprintf("   LOO ESS: %.2f%%", ess_loo_val) else ""
        if (!is.null(mp) && !is.na(mp))
          cat(sprintf("  LOO Mean PAC: %.2f%%%s\n", mp * 100, ess_loo_str))
      }
      cat("  p(LOO): not reported for multicategorical ODA\n")
    }
  }

  cat("\n")
  invisible(x)
}

# ---- summary.oda_fit -------------------------------------------------------- #

#' Summarize a fitted ODA model
#'
#' Returns a structured list with class \code{"oda_fit_summary"} exposing train
#' and LOO sections.  Does not recompute any quantities; fields absent from the
#' fit appear as \code{NA} or \code{NULL}.
#'
#' @param object An \code{oda_fit} object.
#' @param ... Unused.
#' @return A list of class \code{"oda_fit_summary"}.
#' @export
summary.oda_fit <- function(object, ...) {
  eng    <- object$engine %||% "?"
  status <- if (isTRUE(object$ok)) "valid" else
              if (!is.null(object$rule)) object$reason %||% "loo_not_significant"
              else "failed"

  # ---- train section -------------------------------------------------------
  train <- NULL
  if (!is.null(object$rule)) {
    if (inherits(object, "oda_fit_binary")) {
      conf <- object$confusion    # oda_confusion_binary list
      train <- list(
        confusion_raw      = conf,
        confusion_weighted = object$confusion_wt,
        ess                = object$ess,
        pac                = object$pac,   # mean PAC (priors-weighted) * 100
        mean_pac_raw       = if (!is.null(conf)) conf$mean_pac * 100 else NA_real_,
        sensitivity        = if (!is.null(conf)) conf$sensitivity    else NA_real_,
        specificity        = if (!is.null(conf)) conf$specificity    else NA_real_,
        p_mc               = object$p_mc,
        mc_info            = object$mc_info
      )
    } else {
      train <- list(
        confusion_raw      = object$confusion,
        confusion_weighted = object$confusion_wt,
        ess                = object$ess,
        mean_pac           = object$mean_pac,
        pac_by_class       = object$pac_by_class,
        p_mc               = object$p_mc,
        mc_info            = object$mc_info
      )
    }
  }

  # ---- LOO section ---------------------------------------------------------
  loo_section <- NULL
  if (!is.null(object$loo)) {
    lo <- object$loo
    if (isTRUE(lo$allowed)) {
      p_info <- .loo_p_info(object, lo)
      if (inherits(object, "oda_fit_binary")) {
        lc <- lo$confusion   # oda_confusion_binary list
        loo_section <- list(
          allowed   = TRUE,
          confusion = lc,
          ess_loo   = lo$ess_loo,
          mean_pac  = if (!is.null(lc)) lc$mean_pac * 100 else NA_real_,
          p_value   = p_info$p_value,
          p_method  = p_info$p_method,
          p_status  = p_info$p_status,
          p_reason  = p_info$p_reason
        )
      } else {
        # Multiclass: categorical LOO stores $confusion; ordered stores $confusion_raw.
        lc <- lo$confusion %||% lo$confusion_raw
        mp      <- if (!is.null(lc)) lc$mean_pac else NULL   # proportion
        C_loo   <- length(if (!is.null(lc)) lc$classes else object$classes)
        ess_loo <- if (!is.null(mp) && !is.na(mp) && C_loo >= 2L)
                     oda_ess_from_mean(mp, C_loo) else NA_real_
        loo_section <- list(
          allowed      = TRUE,
          confusion    = if (!is.null(lc)) lc$confusion else NULL,
          pac_by_class = if (!is.null(lc)) lc$pac_by_class else NULL,
          mean_pac     = if (!is.null(mp)) mp * 100 else NA_real_,
          ess_loo      = ess_loo,
          p_value      = p_info$p_value,
          p_method     = p_info$p_method,
          p_status     = p_info$p_status,
          p_reason     = p_info$p_reason
        )
      }
    } else {
      loo_section <- list(allowed = FALSE, reason = lo$reason %||% "unavailable")
    }
  }

  s <- list(
    model_type     = eng,
    status         = status,
    reason         = object$reason,
    attr_type      = object$attr_type,
    n_eff          = object$n_eff,
    priors_on      = object$priors_on,
    has_weights    = object$has_weights,
    rule           = object$rule,
    rule_string    = if (!is.null(object$rule)) .oda_fmt_rule(object$rule) else NA_character_,
    objective      = object$ess,
    train          = train,
    loo            = loo_section,
    fields_present = names(object)
  )
  class(s) <- "oda_fit_summary"
  s
}

#' Print an ODA fit summary
#'
#' @param x An \code{oda_fit_summary} from \code{\link{summary.oda_fit}}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.oda_fit_summary <- function(x, ...) {
  cat(sprintf("\nODA Summary (%s)  status=%s  n=%s\n",
              x$model_type, x$status, x$n_eff %||% NA_integer_))
  if (x$status == "failed" && is.null(x$train)) {
    cat(sprintf("  reason: %s\n\n", x$reason %||% "?"))
    return(invisible(x))
  }
  cat(sprintf("  attr_type=%s  priors=%s  weights=%s\n",
              x$attr_type %||% "?",
              if (isTRUE(x$priors_on)) "TRUE" else "FALSE",
              if (isTRUE(x$has_weights)) "TRUE" else "FALSE"))
  cat(sprintf("  Rule: %s\n\n", x$rule_string %||% "<none>"))

  tr <- x$train
  if (!is.null(tr)) {
    cat("  -- Train --\n")
    if (x$model_type == "binary") {
      ess <- .num(tr$ess); pac <- .num(tr$pac)
      sens <- .num(tr$sensitivity); spec <- .num(tr$specificity)
      cat(sprintf("    Mean PAC (wt): %.2f%%   ESS: %.2f%%\n",
                  if (!is.na(pac)) pac else NA_real_,
                  if (!is.na(ess)) ess else NA_real_))
      cat(sprintf("    Sensitivity: %.3f   Specificity: %.3f\n",
                  if (!is.na(sens)) sens else NA_real_,
                  if (!is.na(spec)) spec else NA_real_))
    } else {
      ep <- .num(tr$ess); mp <- .num(tr$mean_pac)
      cat(sprintf("    Mean PAC: %.2f%%   ESS: %.2f%%\n",
                  if (!is.na(mp)) mp else NA_real_,
                  if (!is.na(ep)) ep else NA_real_))
    }
    p_mc <- tr$p_mc
    if (!is.null(p_mc) && !is.na(p_mc))
      cat(sprintf("    p(MC): %s  [MC permutation, one-tailed]\n", .fmt_p(p_mc)))
  }

  lo <- x$loo
  if (!is.null(lo)) {
    cat("  -- LOO --\n")
    if (!isTRUE(lo$allowed)) {
      cat(sprintf("    unavailable (%s)\n", lo$reason %||% "?"))
    } else {
      lc <- lo$confusion
      if (x$model_type == "binary") {
        # lc is an oda_confusion_binary list (stores rates, not counts).
        # Use training confusion for n per class.
        if (!is.null(lc)) {
          r  <- x$rule
          l0 <- if (!is.null(r)) r$label_0 %||% 0L else 0L
          l1 <- if (!is.null(r)) r$label_1 %||% 1L else 1L
          tc <- if (!is.null(tr)) tr$confusion_raw else NULL
          ln0 <- if (!is.null(tc)) (tc$TN %||% 0L) + (tc$FP %||% 0L) else NA_integer_
          ln1 <- if (!is.null(tc)) (tc$TP %||% 0L) + (tc$FN %||% 0L) else NA_integer_
          lp0 <- if (!is.null(lc$specificity) && !is.na(lc$specificity))
                   sprintf("%.1f%%", lc$specificity * 100) else "NA"
          lp1 <- if (!is.null(lc$sensitivity) && !is.na(lc$sensitivity))
                   sprintf("%.1f%%", lc$sensitivity * 100) else "NA"
          cat(sprintf("    CLASS  %6s  %6s\n", "n", "PAC"))
          cat(sprintf("    %5s  %6s  %6s\n", l0, ln0, lp0))
          cat(sprintf("    %5s  %6s  %6s\n", l1, ln1, lp1))
        }
        if (!is.null(lo$ess_loo) && !is.na(lo$ess_loo))
          cat(sprintf("    LOO ESS: %.2f%%\n", lo$ess_loo))
        if (!is.null(lo$mean_pac) && !is.na(lo$mean_pac))
          cat(sprintf("    LOO Mean PAC: %.2f%%\n", lo$mean_pac))
        pv      <- lo$p_value
        pstat   <- lo$p_status %||% "unknown"
        pmethod <- lo$p_method %||% ""
        if (pstat == "computed") {
          cat(sprintf("    p(LOO): %s  [%s]\n", .fmt_p(pv), pmethod))
        } else if (pstat == "not_computed") {
          reason <- lo$p_reason %||% "not computed"
          cat(sprintf("    p(LOO): NA  (%s)\n", reason))
        } else if (pstat == "not_applicable") {
          cat("    p(LOO): not applicable\n")
        }
      } else {
        # Multiclass: loo_section fields extracted by summary.oda_fit.
        # pac_by_class is proportion scale; mean_pac and ess_loo are percent scale.
        pb  <- lo$pac_by_class   # proportion
        cls <- x[["classes"]] %||% NULL   # classes not stored in summary — use NULL gracefully
        if (!is.null(pb) && !is.null(cls)) {
          cat(sprintf("    CLASS  %6s\n", "PAC"))
          for (i in seq_along(cls))
            cat(sprintf("    %5s  %5.1f%%\n", cls[i], pb[i] * 100))
        }
        if (!is.null(lo$ess_loo) && !is.na(lo$ess_loo))
          cat(sprintf("    LOO ESS: %.2f%%\n", lo$ess_loo))
        if (!is.null(lo$mean_pac) && !is.na(lo$mean_pac))
          cat(sprintf("    LOO Mean PAC: %.2f%%\n", lo$mean_pac))
        cat("    p(LOO): not reported for multicategorical ODA\n")
      }
    }
  }
  cat("\n")
  invisible(x)
}

# ---- Accessors -------------------------------------------------------------- #

#' Retrieve predictions from a fitted ODA model
#'
#' Returns stored LOO predictions when available, or calls
#' \code{predict.oda_fit()} on supplied \code{newdata}.  Training predictions
#' are not stored by the engine; supply \code{newdata} to obtain them.
#'
#' @param fit An \code{oda_fit} object.
#' @param split One of \code{"train"} or \code{"loo"}.
#' @param newdata For \code{split = "train"}: numeric vector or single-column
#'   data frame; if \code{NULL}, returns \code{NULL} with a message.
#' @param ... Passed to \code{predict.oda_fit()} when \code{newdata} is
#'   supplied.
#' @return Integer vector of predictions or \code{NULL}.
#' @export
oda_predictions <- function(fit, split = c("train", "loo"), newdata = NULL, ...) {
  split <- match.arg(split)
  if (split == "train") {
    if (!is.null(newdata)) return(predict(fit, newdata, ...))
    message("oda_predictions: training predictions are not stored; ",
            "supply newdata to obtain them via predict.oda_fit()")
    return(NULL)
  }
  # LOO
  lo <- fit$loo
  if (is.null(lo) || !isTRUE(lo$allowed)) {
    message("oda_predictions: LOO predictions not available")
    return(NULL)
  }
  # Multiclass LOO stores y_pred; binary LOO does not
  if (!is.null(lo$y_pred)) return(as.integer(lo$y_pred))
  message("oda_predictions: individual LOO predictions are not stored for this engine/mode")
  NULL
}

#' Retrieve a confusion matrix from a fitted ODA model
#'
#' @param fit An \code{oda_fit} object.
#' @param split One of \code{"train"} or \code{"loo"}.
#' @param weighted Logical; if \code{TRUE} return the priors-weighted confusion.
#'   Weighted LOO confusion is not stored separately.
#' @return The confusion object stored on the fit, or \code{NULL}.
#' @export
oda_confusion <- function(fit, split = c("train", "loo"), weighted = FALSE) {
  split <- match.arg(split)
  if (split == "train") {
    if (weighted) return(fit$confusion_wt)
    return(fit$confusion)
  }
  # LOO
  lo <- fit$loo
  if (is.null(lo) || !isTRUE(lo$allowed)) return(NULL)
  lo$confusion
}

#' Retrieve scalar performance metrics from a fitted ODA model
#'
#' Returns a list of scalar metrics present on the fit.  No quantities are
#' recomputed; absent fields appear as \code{NA_real_}.  LOO p-value uses 2x2
#' Fisher exact when stored and available (\code{p_status = "computed"}); if
#' the value is absent or \code{NA} the status is \code{"not_computed"} with an
#' explicit reason.  Multiclass/polychotomous LOO always returns
#' \code{p_status = "not_computed"}.
#'
#' @param fit An \code{oda_fit} object.
#' @param split One of \code{"train"} or \code{"loo"}.
#' @return Named list of scalar metrics.
#' @export
oda_metrics <- function(fit, split = c("train", "loo")) {
  split <- match.arg(split)

  if (split == "train") {
    if (inherits(fit, "oda_fit_binary")) {
      conf <- fit$confusion
      list(
        ess         = fit$ess %||% NA_real_,
        pac         = fit$pac %||% NA_real_,
        mean_pac_wt = fit$pac %||% NA_real_,
        sensitivity = if (!is.null(conf)) conf$sensitivity else NA_real_,
        specificity = if (!is.null(conf)) conf$specificity else NA_real_,
        mean_pac_raw = if (!is.null(conf)) conf$mean_pac * 100 else NA_real_,
        p_mc        = fit$p_mc %||% NA_real_
      )
    } else {
      list(
        ess          = fit$ess %||% NA_real_,
        mean_pac     = fit$mean_pac %||% NA_real_,
        pac_by_class = fit$pac_by_class,
        p_mc         = fit$p_mc %||% NA_real_
      )
    }
  } else {
    # LOO
    lo <- fit$loo
    if (is.null(lo) || !isTRUE(lo$allowed)) {
      message("oda_metrics: LOO metrics not available")
      return(NULL)
    }
    p_info <- .loo_p_info(fit, lo)
    if (inherits(fit, "oda_fit_binary")) {
      lc <- lo$confusion
      list(
        ess_loo     = lo$ess_loo %||% NA_real_,
        mean_pac    = if (!is.null(lc)) lc$mean_pac * 100 else NA_real_,
        sensitivity = if (!is.null(lc)) lc$sensitivity else NA_real_,
        specificity = if (!is.null(lc)) lc$specificity else NA_real_,
        p_value     = p_info$p_value,
        p_method    = p_info$p_method,
        p_status    = p_info$p_status,
        p_reason    = p_info$p_reason
      )
    } else {
      lc <- lo$confusion   # oda_confusion_multiclass result
      list(
        mean_pac = if (!is.null(lc) && !is.null(lc$mean_pac))
                     lc$mean_pac * 100 else NA_real_,
        p_value  = p_info$p_value,
        p_method = p_info$p_method,
        p_status = p_info$p_status,
        p_reason = p_info$p_reason
      )
    }
  }
}

# ---- oda_d_stat ------------------------------------------------------------- #

#' Compute the D statistic for a fitted ODA model
#'
#' D measures the distance between a model's classification accuracy (ESS) and
#' chance, expressed relative to the number of terminal prediction strata.
#' Formula: \eqn{D = \frac{100}{ESS / strata} - strata}, where \emph{strata}
#' counts terminal prediction endpoints only.
#'
#' Supported rule types and strata definitions:
#' \itemize{
#'   \item Binary (\code{oda_fit_binary}): strata = 2, ESS = \code{fit$ess}.
#'   \item Multiclass ordered (\code{multiclass_ordered} rule): strata =
#'     \code{length(fit$rule$seg_classes)}, ESS = \code{fit$ess}.
#'   \item Multiclass nominal/categorical: returns \code{NA_real_} (strata
#'     count is ambiguous without additional canon specification).
#'   \item Failed fit (\code{ok = FALSE}): returns \code{NA_real_}.
#' }
#'
#' @param fit An \code{oda_fit} object from \code{\link{oda_fit}}.
#' @return A scalar \code{numeric} D value, or \code{NA_real_} when the fit
#'   failed or the rule type does not have an unambiguous strata count.
#' @export
oda_d_stat <- function(fit) {
  # Failed / degenerate fits
  if (!isTRUE(fit$ok)) return(NA_real_)

  if (inherits(fit, "oda_fit_binary")) {
    ess <- fit$ess
    if (is.null(ess) || !is.finite(ess) || ess <= 0) return(NA_real_)
    strata <- 2L
    return(100 / (ess / strata) - strata)
  }

  if (inherits(fit, "oda_fit_multiclass")) {
    rule <- fit$rule
    if (is.null(rule)) return(NA_real_)

    if (identical(rule$type, "multiclass_ordered")) {
      strata <- length(rule$seg_classes)
      if (strata < 2L) return(NA_real_)
      ess <- fit$ess
      if (is.null(ess) || !is.finite(ess) || ess <= 0) return(NA_real_)
      return(100 / (ess / strata) - strata)
    }

    # multiclass_nominal / multiclass_categorical: strata count is ambiguous
    # without explicit canon  -  return NA rather than guess.
    return(NA_real_)
  }

  NA_real_
}
