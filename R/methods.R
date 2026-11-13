## S3 methods and tidy tables for `shapvalues` objects.

#' @export
print.shapvalues <- function(x, ...) {
  arr <- x$phi_array
  d <- dim(arr)
  cat("<shapvalues>\n")
  cat(sprintf("  task     : %s\n", x$task))
  cat(sprintf("  method   : %s (nsim = %d)\n", x$method, x$nsim))
  cat(sprintf("  explained: %d obs x %d features\n", d[1], d[2]))
  if (d[3] > 1L)
    cat(sprintf("  outputs  : %d (%s)\n", d[3],
                paste(utils::head(x$outputs, 4), collapse = ", ")))
  imp <- sv_importance_table(x)
  cat("  top features by mean|SHAP|:\n")
  show <- utils::head(imp, 5)
  for (i in seq_len(nrow(show)))
    cat(sprintf("    %-20s %.4g\n", show$feature[i], show$mean_abs[i]))
  invisible(x)
}

#' @export
summary.shapvalues <- function(object, ...) {
  imp <- sv_importance_table(object)
  structure(list(importance = imp, task = object$task,
                 baseline = object$baseline), class = "summary.shapvalues")
}

#' @export
print.summary.shapvalues <- function(x, ...) {
  cat(sprintf("shapvalues summary (task: %s)\n\n", x$task))
  print(x$importance, row.names = FALSE)
  invisible(x)
}

#' Tidy SHAP attribution table
#'
#' @param x A `shapvalues` object.
#' @param ... Unused.
#' @return A data frame in long form with columns `observation`, `feature`,
#'   `phi`, `feature_value`, and `output` (present when there is more than one
#'   model output, e.g. classes or survival times).
#' @examples
#' fit <- lm(mpg ~ ., data = mtcars)
#' sv <- shapvalues(fit, mtcars[1:4, -1], mtcars[, -1], nsim = 10, nthreads = 2)
#' head(sv_table(sv))
#' @export
sv_table <- function(x, ...) {
  stopifnot(inherits(x, "shapvalues"))
  arr <- x$phi_array
  d <- dim(arr)
  obs <- rownames(x$feature_values)
  if (is.null(obs)) obs <- as.character(seq_len(d[1]))
  grid <- expand.grid(observation = obs, feature = x$features,
                      output = x$outputs, KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  grid$phi <- as.vector(arr)
  fv <- x$feature_values
  fv_long <- basetable::rbindfill(lapply(x$features, function(f) {
    data.frame(observation = obs, feature = f,
               feature_value = .sv_num(fv[[f]]),
               stringsAsFactors = FALSE)
  }))
  out <- basetable::merge(grid, fv_long, by = c("observation", "feature"))
  out <- as.data.frame(out)[c("observation", "feature", "phi", "feature_value", "output")]
  if (length(x$outputs) == 1L) out$output <- NULL
  rownames(out) <- NULL
  out
}

#' @export
as.data.frame.shapvalues <- function(x, ...) sv_table(x)

#' Feature importance table (mean absolute SHAP)
#'
#' @param x A `shapvalues` object.
#' @param by_output If `TRUE` and the model has multiple outputs, return one
#'   importance value per feature per output; otherwise average over outputs.
#' @return A data frame ordered by decreasing `mean_abs`.
#' @examples
#' fit <- lm(mpg ~ ., data = mtcars)
#' sv <- shapvalues(fit, mtcars[1:6, -1], mtcars[, -1], nsim = 10, nthreads = 2)
#' sv_importance_table(sv)
#' @export
sv_importance_table <- function(x, by_output = FALSE) {
  stopifnot(inherits(x, "shapvalues"))
  arr <- x$phi_array
  d <- dim(arr)
  if (by_output && d[3] > 1L) {
    res <- basetable::rbindfill(lapply(seq_len(d[3]), function(o) {
      data.frame(feature = x$features,
                 output = x$outputs[o],
                 mean_abs = as.numeric(colMeans(abs(arr[, , o, drop = FALSE]), dims = 1L)),
                 stringsAsFactors = FALSE)
    }))
  } else {
    ma <- apply(abs(arr), 2, mean)
    res <- data.frame(feature = x$features, mean_abs = as.numeric(ma),
                      stringsAsFactors = FALSE)
  }
  res <- as.data.frame(res)[order(-res$mean_abs), , drop = FALSE]
  rownames(res) <- NULL
  res
}

.sv_num <- function(v) {
  if (is.numeric(v)) return(v)
  if (is.factor(v)) return(as.numeric(v))
  suppressWarnings(as.numeric(as.character(v)))
}
