## The shapvalues colour system.
##
## One place decides colour for every plot. A *style* fixes three roles:
##   - a diverging pair for signed quantities (a contribution's sign, a feature
##     value from low to high),
##   - a sequential ramp for magnitudes and densities,
##   - a qualitative set for categories (clusters, feature domains, lines).
## Structural neutrals (axes, reference lines, text) are the same in every
## style, so the frame reads identically regardless of the data palette.
##
## Default style "shapvalues" is the package's own steel-blue / brick-rose
## pair. The other styles are the 'viridisLite' maps, addressed by their
## option letters: "viridis" (D), "magma" (A), "inferno" (B), "plasma" (C),
## "cividis" (E). Set globally with sv_style(), or per plot with `style =`.

.sv_styles <- c("shapvalues", "viridis", "magma", "inferno", "plasma", "cividis")

## structural neutrals (style-independent)
.sv_ink  <- "#333333"
.sv_axis <- "#5B5B5B"
.sv_line <- "#CBCBCB"
.sv_zero <- "#9A9A9A"
.sv_grid <- "#ECECEC"

.sv_vopt <- function(style)
  switch(style, viridis = "D", magma = "A", inferno = "B",
         plasma = "C", cividis = "E", "D")

#' Set the plotting style
#'
#' Chooses the colour style used by every `sv_*` plot. `"shapvalues"` (the
#' default) is the package's own diverging pair; the others are 'viridisLite'
#' maps. Individual plots also accept a `style` argument that overrides this
#' for one call.
#'
#' @param style One of `"shapvalues"`, `"viridis"`, `"magma"`, `"inferno"`,
#'   `"plasma"`, `"cividis"`.
#' @return The style, invisibly.
#' @examples
#' sv_style("plasma")
#' sv_style("shapvalues")
#' @export
sv_style <- function(style = c("shapvalues", "viridis", "magma", "inferno",
                               "plasma", "cividis")) {
  style <- match.arg(style)
  options(shapvalues.style = style)
  invisible(style)
}

.sv_style <- function() {
  s <- getOption("shapvalues.style", "shapvalues")
  if (length(s) != 1L || !s %in% .sv_styles) "shapvalues" else s
}

## resolve a style name to the full colour set used by the plots
.sv_colours <- function(style = NULL) {
  if (is.null(style)) style <- .sv_style()
  style <- match.arg(style, .sv_styles)
  neutral <- list(ink = .sv_ink, axis = .sv_axis, line = .sv_line,
                  zero = .sv_zero, grid = .sv_grid, style = style)
  if (style == "shapvalues") {
    lo <- "#3B6EA5"; hi <- "#C1495A"; mid <- "#EFEFEC"
    qual <- c("#3B6EA5", "#C1495A", "#5B8C5A", "#D6A24A",
              "#7E6699", "#C97B4A", "#4C4C4C", "#68B0AB")
    ramp <- function(n) grDevices::colorRampPalette(c(lo, mid, hi))(n)
  } else {
    opt <- .sv_vopt(style)
    v5 <- viridisLite::viridis(5L, option = opt)
    lo <- v5[1L]; hi <- v5[5L]; mid <- v5[3L]
    qual <- viridisLite::viridis(8L, option = opt, begin = 0.06, end = 0.94)
    ramp <- function(n) viridisLite::viridis(n, option = opt)
  }
  c(neutral, list(neg = lo, pos = hi, mid = mid, bar = hi,
                  grad_low = lo, grad_high = hi, ramp = ramp,
                  qualitative = qual))
}

#' The shapvalues colour palette
#'
#' The colours a given style resolves to. Inspection helper, not a user
#' setting; use [sv_style()] to change the active style.
#'
#' @param style Style name; defaults to the active style.
#' @return A named list: `neg`, `pos`, `mid`, `bar`, `grad_low`, `grad_high`,
#'   `ramp` (a function of `n`), `qualitative`, and the neutrals `ink`, `axis`,
#'   `line`, `zero`, `grid`.
#' @examples
#' sv_palette()$qualitative
#' sv_palette("inferno")$ramp(5)
#' @export
sv_palette <- function(style = NULL) .sv_colours(style)
