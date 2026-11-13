#' Deprecated functions
#'
#' These underscore names are kept for one release and will be removed in the next.
#' Use the glued lowercase names instead.
#'
#' @param ... Passed to the new function.
#' @return The value of the new function.
#' @name shapvalues-deprecated
#' @keywords internal
NULL

#' @rdname shapvalues-deprecated
#' @export
sv_concentration <- function(...) {
  .Deprecated("svconcentration", package = "shapvalues")
  svconcentration(...)
}

#' @rdname shapvalues-deprecated
#' @export
sv_importance_table <- function(...) {
  .Deprecated("svimportancetable", package = "shapvalues")
  svimportancetable(...)
}

#' @rdname shapvalues-deprecated
#' @export
sv_palette <- function(...) {
  .Deprecated("svpalette", package = "shapvalues")
  svpalette(...)
}

#' @rdname shapvalues-deprecated
#' @export
sv_style <- function(...) {
  .Deprecated("svstyle", package = "shapvalues")
  svstyle(...)
}

#' @rdname shapvalues-deprecated
#' @export
sv_table <- function(...) {
  .Deprecated("svtable", package = "shapvalues")
  svtable(...)
}

#' @rdname shapvalues-deprecated
#' @export
theme_shapvalues <- function(...) {
  .Deprecated("themeshapvalues", package = "shapvalues")
  themeshapvalues(...)
}
