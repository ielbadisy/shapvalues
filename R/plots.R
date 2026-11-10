## Visualisations for `shapvalues` objects. Colour comes from the style system
## in palette.R: the default "shapvalues" style, or any 'viridisLite' map, set
## globally with sv_style() or per call with `style =`. Structural neutrals
## (axes, reference lines, text) are constant across styles. ggplot2 is a
## Suggests dependency; every entry point checks for it first.

.sv_need_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for shapvalues plots; install it to continue.",
         call. = FALSE)
}

#' The shapvalues plot theme
#'
#' The clean frame every `plot.shapvalues()` view uses by default: a white
#' panel with no grid or border and light-grey axis lines. Exported so it can
#' be re-applied after other tweaks, or used on your own plots. `plot()`
#' returns a plain `ggplot` object, so you may also drop this theme and add
#' any `ggplot2::theme_*()` afterwards with `+`.
#'
#' @param base_size Base font size, in points.
#' @return A `ggplot2` theme object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_shapvalues()
#'
#'   fit <- lm(mpg ~ ., data = mtcars)
#'   sv <- shapvalues(fit, mtcars[, -1], mtcars[, -1], nsim = 20, nthreads = 2)
#'   plot(sv) + theme_bw()   # swap in any ggplot2 theme
#' }
#' @export
theme_shapvalues <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid        = ggplot2::element_blank(),
      panel.border      = ggplot2::element_blank(),
      axis.line.x       = ggplot2::element_line(colour = .sv_line, linewidth = 0.5),
      axis.ticks.x      = ggplot2::element_line(colour = .sv_line, linewidth = 0.5),
      axis.ticks.length = ggplot2::unit(3, "pt"),
      axis.ticks.y      = ggplot2::element_blank(),
      axis.title        = ggplot2::element_text(colour = .sv_axis, size = base_size),
      axis.text         = ggplot2::element_text(colour = .sv_axis, size = base_size - 1),
      plot.title        = ggplot2::element_text(colour = .sv_ink, size = base_size + 1,
                                                face = "plain",
                                                margin = ggplot2::margin(b = 6)),
      plot.title.position = "plot",
      plot.subtitle     = ggplot2::element_text(colour = .sv_axis, size = base_size - 1),
      legend.title      = ggplot2::element_text(colour = .sv_axis, size = base_size - 1),
      legend.text       = ggplot2::element_text(colour = .sv_axis, size = base_size - 1),
      legend.key        = ggplot2::element_blank(),
      plot.margin       = ggplot2::margin(10, 14, 8, 10)
    )
}

## internal wrapper: the exported theme plus an optional left axis line
.sv_theme <- function(base_size = 12, y_axis = FALSE) {
  th <- theme_shapvalues(base_size)
  if (y_axis)
    th <- th + ggplot2::theme(
      axis.line.y  = ggplot2::element_line(colour = .sv_line, linewidth = 0.5),
      axis.ticks.y = ggplot2::element_line(colour = .sv_line, linewidth = 0.5))
  th
}

## slim, unframed vertical colour bar
.sv_colourbar <- function(title = "Feature value") {
  ggplot2::guide_colourbar(title = title, barwidth = ggplot2::unit(0.4, "lines"),
                           barheight = ggplot2::unit(7, "lines"), ticks = FALSE,
                           frame.colour = NA, draw.ulim = FALSE, draw.llim = FALSE)
}

## continuous low->high colour scale in the active style, on [0, 1]
.sv_grad <- function(pal, name) {
  ggplot2::scale_colour_gradientn(
    colours = pal$ramp(64), limits = c(0, 1), breaks = c(0, 1),
    labels = c("Low", "High"), name = name, guide = .sv_colourbar(name))
}

## robust 0-1 rescale of a feature value column, clipped at the 2nd/98th pct
.sv_scale01 <- function(v) {
  v <- .sv_num(v)
  q <- stats::quantile(v, c(0.02, 0.98), names = FALSE, na.rm = TRUE)
  if (diff(q) <= 0) q <- range(v, na.rm = TRUE)
  if (diff(q) <= 0) return(rep(0.5, length(v)))
  pmin(pmax((v - q[1]) / (q[2] - q[1]), 0), 1)
}

## deterministic beeswarm offset: bin the values, then spread points within a
## bin symmetrically about the row, centre-outwards, so the swarm reads as a
## violin rather than diagonal stripes. Denser bins spread wider.
.sv_swarm <- function(v, width = 0.4, nbins = 50L) {
  n <- length(v)
  if (n == 0L) return(numeric(0))
  if (diff(range(v)) == 0) return(rep(0, n))
  br  <- seq(min(v), max(v), length.out = nbins + 1L)
  bin <- findInterval(v, br, rightmost.closed = TRUE)
  off <- numeric(n)
  mx  <- max(tabulate(bin), 1L)
  for (b in unique(bin)) {
    ix <- which(bin == b)
    k  <- length(ix)
    if (k == 1L) { off[ix] <- 0; next }
    rank_out <- ceiling(seq_len(k) / 2) * rep_len(c(1, -1), k)
    rank_out <- rank_out - mean(rank_out)
    off[ix] <- rank_out / max(abs(rank_out)) * width * sqrt(k / mx)
  }
  off
}

.sv_pick_output <- function(x, output) {
  outs <- x$outputs
  if (length(outs) == 1L) return(1L)
  if (is.null(output)) {
    idx <- if (identical(x$task, "classification")) length(outs) else 1L
    message("multiple outputs; plotting '", outs[idx],
            "'. Set `output=` to choose another.")
    return(idx)
  }
  if (is.numeric(output)) return(as.integer(output))
  m <- match(output, outs)
  if (is.na(m)) stop("`output` not among: ", paste(outs, collapse = ", "),
                     call. = FALSE)
  m
}

.sv_slice <- function(x, o) {
  arr <- x$phi_array
  phi <- arr[, , o, drop = FALSE][, , 1]
  if (is.null(dim(phi))) phi <- matrix(phi, nrow = dim(arr)[1],
                                       dimnames = list(NULL, x$features))
  phi
}

.sv_order_features <- function(phi, max_features) {
  ma <- colMeans(abs(phi))
  names(sort(ma, decreasing = TRUE))[seq_len(min(max_features, length(ma)))]
}

#' SHAP beeswarm (summary) plot
#'
#' One row per feature, most important at the top. Each point is one
#' observation, positioned by its SHAP value and coloured by the feature value.
#' Points are spread within the row by a deterministic beeswarm offset.
#'
#' @param x A `shapvalues` object.
#' @param output Which model output to show (name or index) when there is more
#'   than one. Defaults to the last class for classification, the first time for
#'   survival.
#' @param max_features Maximum number of features to display, most important
#'   first.
#' @param title Optional plot title. `NULL` (default) draws none.
#' @param style Colour style; `NULL` uses the active style (see [sv_style()]).
#' @return A `ggplot` object.
#' @examples
#' fit <- lm(mpg ~ ., data = mtcars)
#' sv <- shapvalues(fit, mtcars[, -1], mtcars[, -1], nsim = 20, nthreads = 2)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(sv, type = "beeswarm")
#' @noRd
sv_beeswarm <- function(x, output = NULL, max_features = 15L, title = NULL,
                        style = NULL) {
  .sv_need_ggplot()
  stopifnot(inherits(x, "shapvalues"))
  pal <- .sv_colours(style)
  o <- .sv_pick_output(x, output)
  phi <- .sv_slice(x, o)
  feats <- rev(.sv_order_features(phi, max_features))
  fv <- x$feature_values
  df <- basetable::rbindfill(lapply(seq_along(feats), function(i) {
    f <- feats[i]; s <- phi[, f]
    data.frame(feature = f, row = i, shap = s, y = i + .sv_swarm(s),
               value01 = .sv_scale01(fv[[f]]), stringsAsFactors = FALSE)
  }))
  ggplot2::ggplot(df, ggplot2::aes(x = shap, y = y, colour = value01)) +
    ggplot2::geom_hline(yintercept = seq_along(feats), colour = .sv_line,
                        linewidth = 0.3, linetype = "22") +
    ggplot2::geom_vline(xintercept = 0, colour = .sv_zero, linewidth = 0.5) +
    ggplot2::geom_point(size = 1.1, alpha = 1, shape = 16) +
    ggplot2::scale_y_continuous(breaks = seq_along(feats), labels = feats,
                                expand = ggplot2::expansion(add = 0.7)) +
    .sv_grad(pal, "Feature value") +
    ggplot2::labs(x = "SHAP value (impact on model output)", y = NULL,
                  title = title) +
    .sv_theme()
}

#' SHAP importance bar plot
#'
#' Mean absolute SHAP value per feature.
#'
#' @inheritParams sv_beeswarm
#' @return A `ggplot` object.
#' @noRd
sv_importance <- function(x, output = NULL, max_features = 15L, title = NULL,
                          style = NULL) {
  .sv_need_ggplot()
  stopifnot(inherits(x, "shapvalues"))
  pal <- .sv_colours(style)
  o <- .sv_pick_output(x, output)
  phi <- .sv_slice(x, o)
  ma <- sort(colMeans(abs(phi)), decreasing = TRUE)
  ma <- utils::head(ma, max_features)
  df <- data.frame(feature = factor(names(ma), levels = rev(names(ma))),
                   mean_abs = as.numeric(ma))
  ggplot2::ggplot(df, ggplot2::aes(x = mean_abs, y = feature)) +
    ggplot2::geom_col(fill = pal$bar, width = 0.62) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", mean_abs)),
                       hjust = -0.2, size = 3.1, colour = .sv_axis) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.14))) +
    ggplot2::labs(x = "mean(|SHAP value|)", y = NULL, title = title) +
    .sv_theme()
}

#' SHAP waterfall plot for a single observation
#'
#' Builds the prediction from the baseline \eqn{E[f(X)]} by adding one feature
#' contribution at a time, largest first.
#'
#' @inheritParams sv_beeswarm
#' @param row Row index into the explained data.
#' @param max_features Features shown individually; the rest are pooled.
#' @return A `ggplot` object.
#' @noRd
sv_waterfall <- function(x, row = 1L, output = NULL, max_features = 10L,
                         title = NULL, style = NULL) {
  .sv_need_ggplot()
  stopifnot(inherits(x, "shapvalues"))
  pal <- .sv_colours(style)
  o <- .sv_pick_output(x, output)
  phi <- .sv_slice(x, o)[row, ]
  base <- x$baseline[o]
  pred <- x$prediction[row, o]
  fv <- x$feature_values[row, , drop = TRUE]
  ord <- order(abs(phi), decreasing = TRUE)
  phi <- phi[ord]
  if (length(phi) > max_features) {
    rest <- sum(phi[-seq_len(max_features)])
    n_rest <- length(phi) - max_features
    phi <- c(phi[seq_len(max_features)],
             stats::setNames(rest, sprintf("%d other features", n_rest)))
  }
  phi <- rev(phi)
  ny  <- length(phi)
  fmt_v <- function(nm) {
    if (!nm %in% names(fv)) return(nm)
    v <- fv[[nm]]
    if (is.numeric(v)) sprintf("%s = %s", formatC(v, format = "g", digits = 4), nm)
    else sprintf("%s = %s", as.character(v), nm)
  }
  ylabs <- vapply(names(phi), fmt_v, character(1))
  df <- data.frame(y = seq_len(ny), shap = as.numeric(phi))
  df$end   <- base + cumsum(df$shap)
  df$start <- df$end - df$shap
  df$sign  <- ifelse(df$shap >= 0, "pos", "neg")
  span <- diff(range(c(base, pred, df$end, df$start)))
  if (span == 0) span <- 1
  h <- 0.34; tipw <- 0.022 * span
  df$bx <- with(df, ifelse(abs(shap) > tipw, end - sign(shap) * tipw, start))
  tri <- basetable::rbindfill(lapply(seq_len(ny), function(i) data.frame(
    grp = i, sign = df$sign[i],
    tx = c(df$bx[i], df$end[i], df$bx[i]),
    ty = c(df$y[i] - h, df$y[i], df$y[i] + h))))
  lbl <- sprintf("%+.2f", df$shap)
  z <- abs(df$shap) < 5e-3
  lbl[z] <- ifelse(df$shap[z] >= 0, "+0", "-0")
  df$inside <- abs(df$shap) > 0.16 * span
  df$lx   <- ifelse(df$inside, (df$start + df$end) / 2, df$end)
  df$lcol <- ifelse(df$inside, "in", df$sign)
  df$lhj  <- ifelse(df$inside, 0.5, ifelse(df$shap >= 0, -0.25, 1.25))
  ggplot2::ggplot(df) +
    ggplot2::geom_vline(xintercept = base, colour = .sv_line, linewidth = 0.5) +
    ggplot2::geom_rect(ggplot2::aes(xmin = pmin(start, bx), xmax = pmax(start, bx),
                                    ymin = y - h, ymax = y + h, fill = sign)) +
    ggplot2::geom_polygon(data = tri,
                          ggplot2::aes(x = tx, y = ty, group = grp, fill = sign)) +
    ggplot2::geom_text(ggplot2::aes(x = lx, y = y, colour = lcol, label = lbl,
                                    hjust = lhj), size = 3, show.legend = FALSE) +
    ggplot2::annotate("text", x = base, y = 0.25,
                      label = sprintf("E[f(X)] = %.3g", base),
                      colour = .sv_axis, size = 3, hjust = 0.5) +
    ggplot2::annotate("text", x = pred, y = ny + 0.8,
                      label = sprintf("f(x) = %.3g", pred),
                      colour = .sv_axis, size = 3.1, hjust = 0.5) +
    ggplot2::scale_y_continuous(breaks = df$y, labels = ylabs,
                                expand = ggplot2::expansion(add = c(0.7, 1.1))) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.14)) +
    ggplot2::scale_fill_manual(values = c(pos = pal$pos, neg = pal$neg),
                               guide = "none") +
    ggplot2::scale_colour_manual(values = c(pos = pal$pos, neg = pal$neg,
                                            `in` = "white"), guide = "none") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    .sv_theme() +
    ggplot2::theme(axis.line.x = ggplot2::element_line(colour = .sv_line,
                                                       linewidth = 0.5))
}

#' SHAP dependence plot
#'
#' SHAP value of one feature against its observed value, optionally coloured by
#' a second feature.
#'
#' @inheritParams sv_beeswarm
#' @param feature Feature on the x axis.
#' @param colour_by Optional second feature used to colour the points.
#' @param smooth If `TRUE`, add a loess trend line.
#' @return A `ggplot` object.
#' @noRd
sv_dependence <- function(x, feature, colour_by = NULL, output = NULL,
                          smooth = FALSE, title = NULL, style = NULL) {
  .sv_need_ggplot()
  stopifnot(inherits(x, "shapvalues"))
  if (!feature %in% x$features) stop("unknown feature: ", feature, call. = FALSE)
  pal <- .sv_colours(style)
  o <- .sv_pick_output(x, output)
  phi <- .sv_slice(x, o)
  df <- data.frame(x = .sv_num(x$feature_values[[feature]]), shap = phi[, feature])
  coloured <- !is.null(colour_by)
  if (coloured) {
    if (!colour_by %in% x$features) stop("unknown feature: ", colour_by, call. = FALSE)
    df$value01 <- .sv_scale01(x$feature_values[[colour_by]])
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = shap, colour = value01))
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = shap))
  }
  p <- p + ggplot2::geom_hline(yintercept = 0, colour = .sv_zero, linewidth = 0.4)
  if (isTRUE(smooth))
    p <- p + ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                                  span = 0.9, colour = .sv_axis, linewidth = 0.6)
  p <- if (coloured)
    p + ggplot2::geom_point(size = 1.3, alpha = 0.9, shape = 16)
  else
    p + ggplot2::geom_point(size = 1.3, alpha = 0.6, shape = 16, colour = pal$neg)
  if (coloured) p <- p + .sv_grad(pal, colour_by)
  p +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.08)) +
    ggplot2::labs(x = feature, y = sprintf("SHAP value for %s", feature),
                  title = title) +
    .sv_theme(y_axis = TRUE)
}

#' SurvSHAP(t): SHAP values over time for survival models
#'
#' One line per feature across the evaluation times, for one observation.
#'
#' @param x A `shapvalues` object with `task = "survival"`.
#' @param row Row index to display.
#' @param max_features Number of features drawn, largest time-integrated
#'   absolute attribution first.
#' @param title Optional plot title.
#' @param style Colour style; `NULL` uses the active style.
#' @return A `ggplot` object.
#' @noRd
sv_surv <- function(x, row = 1L, max_features = 8L, title = NULL, style = NULL) {
  .sv_need_ggplot()
  stopifnot(inherits(x, "shapvalues"))
  if (!identical(x$task, "survival"))
    stop("sv_surv() needs a survival `shapvalues` object.", call. = FALSE)
  pal <- .sv_colours(style)
  arr <- x$phi_array
  times <- x$times
  mat <- arr[row, , , drop = FALSE][1, , ]
  score <- rowMeans(abs(mat))
  feats <- names(sort(score, decreasing = TRUE))[seq_len(min(max_features, length(score)))]
  df <- basetable::rbindfill(lapply(feats, function(f)
    data.frame(time = times, shap = mat[f, ], feature = f, stringsAsFactors = FALSE)))
  df$feature <- factor(df$feature, levels = feats)
  ggplot2::ggplot(df, ggplot2::aes(x = time, y = shap, colour = feature)) +
    ggplot2::geom_hline(yintercept = 0, colour = .sv_zero, linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_colour_manual(values = rep_len(pal$qualitative, length(feats)),
                                 name = NULL) +
    ggplot2::labs(x = "Time", y = "SHAP value", title = title) +
    .sv_theme(y_axis = TRUE)
}

#' Plot a shapvalues object
#'
#' Single entry point for every `shapvalues` visualisation; `type` picks the
#' view, and `...` is forwarded to it. All views take `output` (which model
#' output to show, when there is more than one) and `style` (colour style;
#' `NULL` uses the active style, see [sv_style()]); most take `title`.
#'
#' \describe{
#'   \item{`"beeswarm"`}{(default) one row per feature, points positioned by
#'     SHAP value and coloured by feature value. Args: `max_features`.}
#'   \item{`"importance"`}{bar plot of mean `|SHAP|` per feature. Args:
#'     `max_features`.}
#'   \item{`"waterfall"`}{builds one observation's prediction from the
#'     baseline, one feature at a time. Args: `row`, `max_features`.}
#'   \item{`"dependence"`}{SHAP value of one feature against its observed
#'     value. Args: `feature` (required), `colour_by`, `smooth`.}
#'   \item{`"interaction"`}{pairwise SHAP interaction scatter; needs
#'     `interactions = TRUE` in [shapvalues()]. Args: `feature` (required),
#'     `interact`.}
#'   \item{`"surv"`}{SurvSHAP(t): one line per feature over evaluation times,
#'     for one observation. Needs `task = "survival"`. Args: `row`,
#'     `max_features`.}
#'   \item{`"compass"`}{places each observation by the direction its
#'     attribution mass points; distance from centre is concentration (see
#'     [sv_concentration()]). Args: `groups`, `order`, `colour_by`.}
#'   \item{`"fingerprint"`}{radial `|SHAP|` signature, one facet per
#'     observation or per k-means cluster. Args: `rows`, `k`, `max_features`.}
#'   \item{`"pathway"`}{alluvial diagram of the top-ranked contributors per
#'     observation; needs the \pkg{ggalluvial} package. Args: `rank`,
#'     `min_n`, `signed`.}
#' }
#'
#' @param x A `shapvalues` object.
#' @param type One of `"beeswarm"`, `"importance"`, `"waterfall"`,
#'   `"dependence"`, `"interaction"`, `"surv"`, `"compass"`, `"fingerprint"`,
#'   `"pathway"`.
#' @param ... Arguments for the chosen `type`; see Details.
#' @return The `ggplot` object, invisibly (printed as a side effect). You may
#'   assign it and add further `ggplot2` layers with `+`.
#' @examples
#' fit <- lm(mpg ~ ., data = mtcars)
#' sv <- shapvalues(fit, mtcars[, -1], mtcars[, -1], nsim = 20, nthreads = 2)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot(sv)
#'   plot(sv, type = "importance")
#'   plot(sv, type = "dependence", feature = "wt", colour_by = "hp")
#'   plot(sv, type = "waterfall", row = 1) + ggplot2::theme_bw()
#' }
#' @export
#' @importFrom graphics plot
plot.shapvalues <- function(x, type = c("beeswarm", "importance", "waterfall",
                                        "dependence", "interaction", "surv",
                                        "compass", "fingerprint", "pathway"), ...) {
  type <- match.arg(type)
  fn <- switch(type,
               beeswarm = sv_beeswarm, importance = sv_importance,
               waterfall = sv_waterfall, dependence = sv_dependence,
               interaction = sv_interaction, surv = sv_surv,
               compass = sv_compass, fingerprint = sv_fingerprint,
               pathway = sv_pathway)
  p <- fn(x, ...)
  print(p)
  invisible(p)
}
