## Opinionated attribution views. These do not marginalise features one at a
## time; they summarise the multivariate attribution profile of each
## observation.
##
##  - sv_compass:      each observation is a point placed where its attribution
##                     mass points; distance from the centre is the explanation
##                     concentration EC in [0, 1].
##  - sv_concentration: the numbers behind the compass.
##  - sv_fingerprint:  a radial signature of |SHAP| per feature, per observation
##                     or per k-means cluster of SHAP profiles.
##  - sv_pathway:      an alluvial of the top-`rank` contributors: node =
##                     feature, flow = observations, colour = sign of the
##                     dominant contributor.

## shared core: attribution barycentre + concentration for every observation
.sv_compass_core <- function(x, groups = NULL, output = NULL,
                             order = c("cluster", "importance", "angle", "none")) {
  order <- match.arg(order)
  o <- .sv_pick_output(x, output)
  phi <- .sv_slice(x, o)                                   # n x p
  A <- abs(phi)
  if (!is.null(groups)) {
    if (is.list(groups))
      groups <- stats::setNames(rep(names(groups), lengths(groups)),
                                unlist(groups, use.names = FALSE))
    g <- groups[colnames(A)]
    g[is.na(g)] <- "other"
    A <- t(rowsum(t(A), group = g, reorder = FALSE))
  }
  axes <- colnames(A)
  m <- length(axes)
  if (m < 2L) stop("the compass needs at least two features or groups.",
                   call. = FALSE)
  ord <- switch(
    order,
    importance = order(colMeans(A), decreasing = TRUE),
    cluster = if (m >= 3L && nrow(A) >= 3L) {
      cc <- suppressWarnings(stats::cor(A))
      cc[!is.finite(cc)] <- 0
      stats::hclust(stats::as.dist(1 - cc), method = "average")$order
    } else seq_len(m),
    angle = seq_len(m),
    none = seq_len(m)
  )
  A <- A[, ord, drop = FALSE]
  axes <- axes[ord]
  ang <- seq(0, 2 * pi, length.out = m + 1L)[-(m + 1L)]
  ux <- cos(ang); uy <- sin(ang)
  tot <- rowSums(A)
  tot[tot == 0] <- .Machine$double.eps
  cx <- as.numeric(A %*% ux) / tot
  cy <- as.numeric(A %*% uy) / tot
  dom <- axes[max.col(A, ties.method = "first")]
  pts <- data.frame(
    observation = rownames(phi) %||% as.character(seq_len(nrow(phi))),
    cx = cx, cy = cy, ec = sqrt(cx^2 + cy^2),
    angle = atan2(cy, cx), dominant = dom, total_abs = rowSums(abs(phi)),
    stringsAsFactors = FALSE)
  geom <- data.frame(axis = axes, ux = ux, uy = uy, stringsAsFactors = FALSE)
  list(points = pts, geometry = geom, output = x$outputs[o])
}

#' Explanation concentration
#'
#' For each explained observation, the direction and concentration of its
#' attribution mass. With unit vectors \eqn{u_j} placed around a circle, one per
#' feature (or feature group), the barycentre is
#' \eqn{C_i = \sum_j |\phi_{ij}| u_j / \sum_j |\phi_{ij}|} and the concentration
#' index is \eqn{EC_i = \lVert C_i \rVert \in [0, 1]}: near 1 the prediction is
#' driven by one coherent feature direction, near 0 by competing contributions
#' spread across the feature space.
#'
#' @param x A `shapvalues` object.
#' @param groups Optional feature grouping: a named vector `feature -> group`,
#'   or a list `group -> features`. When given, the compass axes are the groups.
#' @param output Model output to use when there is more than one.
#' @param order Axis order: `"cluster"` (correlated attributions adjacent),
#'   `"importance"`, `"angle"`, or `"none"`.
#' @return A data frame with columns `observation`, `cx`, `cy`, `ec`, `angle`,
#'   `dominant`, `total_abs`.
#' @examples
#' fit <- lm(mpg ~ ., data = mtcars)
#' sv <- shapvalues(fit, mtcars[, -1], mtcars[, -1], nsim = 20, nthreads = 2)
#' head(sv_concentration(sv))
#' @export
sv_concentration <- function(x, groups = NULL, output = NULL, order = "cluster") {
  stopifnot(inherits(x, "shapvalues"))
  .sv_compass_core(x, groups, output, order)$points
}

#' SHAP compass
#'
#' Places every observation where its attribution mass points: features (or
#' feature groups) sit around a circle, and an observation is drawn at the
#' `|SHAP|`-weighted average of the feature directions. Distance from the centre
#' is the explanation concentration `EC` (see [sv_concentration()]); the size of
#' a point is the total attribution magnitude.
#'
#' @inheritParams sv_concentration
#' @param colour_by `"prediction"`, `"concentration"`, or `"none"`.
#' @param title Optional plot title.
#' @param style Colour style; `NULL` uses the active style.
#' @return A `ggplot` object.
#' @examples
#' fit <- lm(mpg ~ ., data = mtcars)
#' sv <- shapvalues(fit, mtcars[, -1], mtcars[, -1], nsim = 20, nthreads = 2)
#' if (requireNamespace("ggplot2", quietly = TRUE)) sv_compass(sv)
#' @export
sv_compass <- function(x, groups = NULL, output = NULL, order = "cluster",
                       colour_by = c("prediction", "concentration", "none"),
                       title = NULL, style = NULL) {
  .sv_need_ggplot()
  stopifnot(inherits(x, "shapvalues"))
  colour_by <- match.arg(colour_by)
  pal <- .sv_colours(style)
  cc <- .sv_compass_core(x, groups, output, order)
  pts <- cc$points; geom <- cc$geometry
  circ <- data.frame(a = seq(0, 2 * pi, length.out = 400))
  circ$cx <- cos(circ$a); circ$cy <- sin(circ$a)
  message(sprintf("mean explanation concentration EC = %.2f", mean(pts$ec)))
  has_col <- colour_by != "none"
  if (has_col)
    pts$col <- if (colour_by == "prediction")
      x$prediction[, .sv_pick_output(x, output)] else pts$ec
  p <- ggplot2::ggplot() +
    ggplot2::geom_path(data = circ, ggplot2::aes(cx, cy), colour = .sv_line,
                       linewidth = 0.5) +
    ggplot2::geom_segment(data = geom,
                          ggplot2::aes(x = 0, y = 0, xend = ux, yend = uy),
                          colour = .sv_grid, linewidth = 0.5)
  if (!has_col) {
    p <- p + ggplot2::geom_point(data = pts,
      ggplot2::aes(cx, cy, size = total_abs), colour = pal$neg, alpha = 0.55)
  } else {
    p <- p + ggplot2::geom_point(data = pts,
      ggplot2::aes(cx, cy, size = total_abs, colour = col), alpha = 0.65) +
      ggplot2::scale_colour_gradientn(colours = pal$ramp(64),
        name = if (colour_by == "prediction") "Prediction" else "EC")
  }
  p +
    ggplot2::geom_text(data = geom,
      ggplot2::aes(1.14 * ux, 1.14 * uy, label = axis),
      colour = .sv_ink, fontface = "bold", size = 3.3) +
    ggplot2::scale_size_continuous(range = c(0.6, 4.5), name = "|SHAP| total") +
    ggplot2::coord_equal(xlim = c(-1.32, 1.32), ylim = c(-1.32, 1.32)) +
    ggplot2::labs(title = title,
                  subtitle = "direction = where attribution mass points; radius = concentration (EC)") +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(colour = .sv_ink, size = 13),
      plot.subtitle = ggplot2::element_text(colour = .sv_axis, size = 10),
      legend.title = ggplot2::element_text(colour = .sv_axis, size = 10),
      legend.text = ggplot2::element_text(colour = .sv_axis, size = 9))
}

#' SHAP fingerprint
#'
#' A radial signature of the explanation: one spoke per feature, radius the mean
#' `|SHAP|`, points coloured by the sign of the contribution. Either one
#' fingerprint per observation (`rows`) or one per k-means cluster of the SHAP
#' profiles (`k`), which surfaces recurring attribution phenotypes.
#'
#' @inheritParams sv_beeswarm
#' @param rows Row indices to draw, one facet each. Ignored when `k` is set.
#'   Default: the three observations with the largest total `|SHAP|`.
#' @param k Number of k-means clusters of the (scaled) SHAP vectors; one facet
#'   per cluster.
#' @param max_features Features drawn, most important first.
#' @return A `ggplot` object.
#' @examples
#' fit <- lm(mpg ~ ., data = mtcars)
#' sv <- shapvalues(fit, mtcars[, -1], mtcars[, -1], nsim = 20, nthreads = 2)
#' if (requireNamespace("ggplot2", quietly = TRUE)) sv_fingerprint(sv, k = 3)
#' @export
sv_fingerprint <- function(x, rows = NULL, k = NULL, output = NULL,
                           max_features = 16L, title = NULL, style = NULL) {
  .sv_need_ggplot()
  stopifnot(inherits(x, "shapvalues"))
  pal <- .sv_colours(style)
  o <- .sv_pick_output(x, output)
  phi <- .sv_slice(x, o)
  feats <- .sv_order_features(phi, max_features)
  phi <- phi[, feats, drop = FALSE]
  if (!is.null(k)) {
    k <- as.integer(k)
    if (k < 2L || k >= nrow(phi)) stop("`k` must be between 2 and n - 1.",
                                       call. = FALSE)
    cl <- stats::kmeans(scale(phi), centers = k, nstart = 20L)$cluster
    df <- do.call(rbind, lapply(sort(unique(cl)), function(g) {
      sub <- phi[cl == g, , drop = FALSE]
      data.frame(panel = sprintf("cluster %d (n = %d)", g, nrow(sub)),
                 feature = feats, mag = colMeans(abs(sub)),
                 direction = ifelse(colMeans(sub) >= 0, "Positive", "Negative"),
                 stringsAsFactors = FALSE)
    }))
  } else {
    if (is.null(rows))
      rows <- utils::head(order(rowSums(abs(phi)), decreasing = TRUE), 3L)
    df <- do.call(rbind, lapply(rows, function(i) {
      data.frame(panel = sprintf("obs %s",
                   (rownames(phi) %||% as.character(seq_len(nrow(phi))))[i]),
                 feature = feats, mag = abs(phi[i, ]),
                 direction = ifelse(phi[i, ] >= 0, "Positive", "Negative"),
                 stringsAsFactors = FALSE)
    }))
  }
  df$feature <- factor(df$feature, levels = feats)
  ymax <- max(df$mag)
  ggplot2::ggplot(df, ggplot2::aes(x = feature, y = mag, group = panel)) +
    ggplot2::geom_polygon(fill = pal$grid, colour = .sv_axis, linewidth = 0.6,
                          alpha = 0.55) +
    ggplot2::geom_point(ggplot2::aes(colour = direction), size = 2.4) +
    ggplot2::coord_polar(clip = "off") +
    ggplot2::facet_wrap(~ panel) +
    ggplot2::scale_y_continuous(limits = c(0, ymax * 1.05)) +
    ggplot2::scale_colour_manual(
      values = c(Positive = pal$pos, Negative = pal$neg), name = NULL) +
    ggplot2::labs(x = NULL, y = "mean |SHAP|", title = title) +
    .sv_theme(base_size = 11) +
    ggplot2::theme(
      panel.spacing = ggplot2::unit(14, "pt"),
      axis.text.x   = ggplot2::element_text(colour = .sv_axis, size = 7),
      axis.text.y   = ggplot2::element_blank(),
      legend.text   = ggplot2::element_text(colour = .sv_axis, size = 9))
}

#' SHAP attribution pathways
#'
#' An alluvial diagram of the `rank` strongest contributors per observation.
#' Each node is a feature; the flow width is the number of observations; the
#' flow colour is the sign of the dominant (rank-1) contributor. Recurring
#' pathways stand out; rare ones (fewer than `min_n` observations) are dropped.
#'
#' @inheritParams sv_beeswarm
#' @param rank Number of ranked contributors per observation (2 to 5).
#' @param min_n Minimum observations for a pathway to be kept. Default about
#'   2 percent of the explained rows.
#' @param signed If `TRUE`, split each node into `feature (+)` / `feature (-)`.
#' @return A `ggplot` object. Requires the \pkg{ggalluvial} package.
#' @examples
#' fit <- lm(mpg ~ ., data = mtcars)
#' sv <- shapvalues(fit, mtcars[, -1], mtcars[, -1], nsim = 20, nthreads = 2)
#' if (requireNamespace("ggalluvial", quietly = TRUE)) sv_pathway(sv)
#' @export
sv_pathway <- function(x, rank = 3L, min_n = NULL, output = NULL,
                       signed = FALSE, title = NULL, style = NULL) {
  .sv_need_ggplot()
  if (!requireNamespace("ggalluvial", quietly = TRUE))
    stop("sv_pathway() needs the 'ggalluvial' package.", call. = FALSE)
  stopifnot(inherits(x, "shapvalues"))
  rank <- max(2L, min(5L, as.integer(rank)))
  pal <- .sv_colours(style)
  o <- .sv_pick_output(x, output)
  phi <- .sv_slice(x, o)
  n <- nrow(phi)
  if (is.null(min_n)) min_n <- max(3L, round(0.02 * n))

  top <- t(apply(phi, 1L, function(r) {
    ix <- order(abs(r), decreasing = TRUE)[seq_len(rank)]
    node <- if (signed)
      paste0(colnames(phi)[ix], ifelse(r[ix] >= 0, " (+)", " (-)"))
    else colnames(phi)[ix]
    c(node, if (r[ix[1L]] >= 0) "positive" else "negative")
  }))
  paths <- as.data.frame(top, stringsAsFactors = FALSE)
  names(paths) <- c(paste0("r", seq_len(rank)), "dom")
  agg <- stats::aggregate(list(nobs = seq_len(n)),
                          by = as.list(paths), FUN = length)
  keep <- agg$nobs >= min_n
  dropped <- sum(agg$nobs[!keep])
  agg <- agg[keep, , drop = FALSE]
  if (!nrow(agg)) stop("no pathway reaches min_n = ", min_n,
                       "; lower min_n or rank.", call. = FALSE)
  agg$path_id <- seq_len(nrow(agg))
  long <- do.call(rbind, lapply(seq_len(rank), function(r)
    data.frame(path_id = agg$path_id, nobs = agg$nobs, dom = agg$dom,
               rank = r, node = agg[[paste0("r", r)]],
               stringsAsFactors = FALSE)))
  rlabs <- c("Dominant", "Second", "Third", "Fourth", "Fifth")[seq_len(rank)]
  if (dropped > 0L)
    message(sprintf("%d observations in pathways below min_n = %d not shown",
                    dropped, min_n))
  lab_min <- 0.03 * sum(agg$nobs)                 # label only chunky strata
  ggplot2::ggplot(long, ggplot2::aes(x = factor(rank), y = nobs,
                                     stratum = node, alluvium = path_id,
                                     fill = dom)) +
    ggalluvial::geom_flow(width = 0.3, alpha = 0.62, curve_type = "sigmoid") +
    ggalluvial::geom_stratum(width = 0.3, fill = "white", colour = .sv_line,
                             linewidth = 0.4) +
    ggplot2::geom_text(stat = ggalluvial::StatStratum, size = 2.9,
                       colour = .sv_ink,
                       ggplot2::aes(label = ifelse(
                         ggplot2::after_stat(count) >= lab_min,
                         as.character(ggplot2::after_stat(stratum)), ""))) +
    ggplot2::scale_x_discrete(labels = rlabs, expand = c(0.08, 0.08)) +
    ggplot2::scale_fill_manual(
      values = c(positive = pal$pos, negative = pal$neg),
      name = "Dominant SHAP sign", na.translate = FALSE) +
    ggplot2::labs(x = "Attribution rank", y = "Observations", title = title) +
    .sv_theme() +
    ggplot2::theme(axis.line.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   legend.position = "top")
}
