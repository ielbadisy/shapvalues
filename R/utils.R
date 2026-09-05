## small internal helper
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

## Silence R CMD check notes for bare ggplot2 aes() symbols.
utils::globalVariables(c(
  "feature", "shap", "mean_abs", "start", "end", "sign", "time",
  "y", "yend", "value01", "bx", "tx", "ty", "grp", "lx", "lcol", "lhj",
  "cx", "cy", "ux", "uy", "total_abs", "col", "axis", "mag", "direction",
  "panel", "node", "path_id", "nobs", "dom", "rank", "stratum", "count"
))
