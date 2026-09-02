## Input coercion, prediction-callback resolution and task detection.
## Internal; not exported.

.sv_as_df <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.matrix(x)) return(as.data.frame(x, stringsAsFactors = FALSE))
  stop("`X` / `bg` must be a data frame or matrix.", call. = FALSE)
}

.sv_default_features <- function(X, task, times) {
  nm <- names(X)
  if (identical(task, "survival") || !is.null(times)) {
    drop <- nm[grepl("^(time|status|event|Surv)$", nm, ignore.case = TRUE)]
    nm <- setdiff(nm, drop)
  }
  nm
}

.sv_check_features <- function(features, X) {
  miss <- setdiff(features, names(X))
  if (length(miss))
    stop("features not found in `X`: ", paste(miss, collapse = ", "), call. = FALSE)
  if (!length(features)) stop("no features to explain.", call. = FALSE)
  invisible(TRUE)
}

## make `bg` column set / types match `X` where possible
.sv_align_bg <- function(bg, X) {
  common <- intersect(names(X), names(bg))
  miss <- setdiff(names(X), names(bg))
  if (length(miss))
    stop("`bg` is missing columns present in `X`: ", paste(miss, collapse = ", "),
         call. = FALSE)
  bg <- bg[, names(X), drop = FALSE]
  for (nm in common) {
    if (is.factor(X[[nm]]) && !is.factor(bg[[nm]]))
      bg[[nm]] <- factor(bg[[nm]], levels = levels(X[[nm]]))
    if (is.factor(X[[nm]]) && is.factor(bg[[nm]]))
      bg[[nm]] <- factor(bg[[nm]], levels = levels(X[[nm]]))
  }
  bg
}

.sv_resolve_pred_fun <- function(object, pred_fun) {
  if (!is.null(pred_fun)) {
    if (!is.function(pred_fun)) stop("`pred_fun` must be a function.", call. = FALSE)
    return(pred_fun)
  }
  if (is.null(object))
    stop("supply `pred_fun` when `object` is not a fitted model.", call. = FALSE)
  cls <- class(object)
  if (inherits(object, "ranger")) {
    return(function(object, newdata, ...) {
      pr <- stats::predict(object, data = newdata)$predictions
      pr
    })
  }
  if (inherits(object, "xgb.Booster")) {
    return(function(object, newdata, ...) {
      m <- if (is.matrix(newdata)) newdata else stats::model.matrix(~ . - 1, newdata)
      stats::predict(object, m)
    })
  }
  if (inherits(object, "glm")) {
    return(function(object, newdata, ...) stats::predict(object, newdata, type = "response"))
  }
  if (inherits(object, "lm")) {
    return(function(object, newdata, ...) stats::predict(object, newdata))
  }
  ## generic fallback
  function(object, newdata, ...) stats::predict(object, newdata)
}

## call the callback, tolerating with/without a `times` argument
.sv_predict <- function(pf, object, newdata, times, task) {
  args <- names(formals(pf))
  use_times <- !is.null(times) &&
    ("times" %in% args || "..." %in% args || identical(task, "survival"))
  out <- tryCatch(
    if (use_times) pf(object, newdata, times = times) else pf(object, newdata),
    error = function(e) {
      if (use_times) pf(object, newdata) else stop(e)
    }
  )
  out
}

## normalise a prediction to a numeric matrix (rows = observations)
.sv_detect_task <- function(preds, task, times) {
  if (is.data.frame(preds)) preds <- as.matrix(preds)
  if (is.list(preds) && !is.data.frame(preds)) {
    stop("prediction callback returned a list; return a numeric vector or matrix.",
         call. = FALSE)
  }
  if (is.null(dim(preds))) preds <- matrix(as.numeric(preds), ncol = 1L)
  storage.mode(preds) <- "double"

  detected <- task
  if (identical(task, "auto")) {
    detected <- if (!is.null(times) && ncol(preds) > 1L) "survival"
    else if (ncol(preds) > 1L) "classification"
    else "regression"
  }
  if (identical(detected, "survival") && !is.null(times)) {
    if (ncol(preds) != length(times)) {
      stop(sprintf("survival prediction has %d columns but `times` has %d entries.",
                   ncol(preds), length(times)), call. = FALSE)
    }
  }
  list(preds = preds, task = detected)
}

.sv_output_names <- function(preds, task, times) {
  k <- ncol(preds)
  if (identical(task, "survival") && !is.null(times)) return(paste0("t=", times))
  if (k == 1L) return("value")
  if (!is.null(colnames(preds))) return(colnames(preds))
  paste0("class", seq_len(k))
}
