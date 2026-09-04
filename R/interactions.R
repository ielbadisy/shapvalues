## Pairwise SHAP interaction values (single-output models).
##
## For an instance x and a pair of features (i, j),
##   Phi_ij = (1/2) E_{s, S}[ v(S u {i,j}) - v(S u {i}) - v(S u {j}) + v(S) ],
## where s is uniform on {0, ..., p-2} and S is a uniform subset of the other
## p-2 features of size s. This is the Grabisch-Roubens Shapley interaction
## index, halved and symmetrised as in Lundberg et al. (2020). The diagonal is
## the main effect, Phi_ii = phi_i - sum_{j != i} Phi_ij, so each row of the
## p x p interaction matrix sums to the SHAP value phi_i.

.sv_interactions <- function(pf, object, X, bg, features, phi_mat, nsim) {
  n <- nrow(X)
  p <- length(features)
  if (p < 2L) stop("interactions need at least two features.", call. = FALSE)
  B <- nrow(bg)
  pairs <- utils::combn(p, 2L)
  np <- ncol(pairs)

  ## one background row per sample, reused across observations
  pool <- sample.int(B)
  bg_idx <- if (nsim >= B) rep_len(pool, nsim)
            else pool[floor((seq_len(nsim) - 1) * B / nsim) + 1L]

  ## per (sample, pair): a subset S of the other p-2 features
  Smask <- matrix(FALSE, nsim * np, p)
  for (m in seq_len(nsim)) for (pp in seq_len(np)) {
    ij <- pairs[, pp]
    others <- seq_len(p)[-ij]
    s <- sample.int(p - 1L, 1L) - 1L            # uniform on {0, ..., p-2}
    if (s > 0L) Smask[(m - 1L) * np + pp, others[sample.int(p - 2L, s)]] <- TRUE
  }

  ## batch rows ordered (obs, sample, pair, coalition c = 1..4)
  ##   c: 1 = S, 2 = S u {i}, 3 = S u {j}, 4 = S u {i,j}
  L <- 4L
  N <- n * nsim * np * L
  cc  <- rep_len(1:4, N)
  rem <- (seq_len(N) - 1L) %/% L
  pp  <- (rem %% np) + 1L
  rem <- rem %/% np
  mm  <- (rem %% nsim) + 1L
  ii  <- (rem %/% nsim) + 1L
  sp  <- (mm - 1L) * np + pp                    # row into Smask

  batch <- bg[bg_idx[mm], , drop = FALSE]
  rownames(batch) <- NULL
  pi_ <- pairs[1, pp]; pj_ <- pairs[2, pp]
  for (f in seq_len(p)) {
    on <- Smask[cbind(sp, f)] |
      (pi_ == f & (cc == 2L | cc == 4L)) |
      (pj_ == f & (cc == 3L | cc == 4L))
    if (any(on)) batch[[features[f]]][on] <- X[[features[f]]][ii[on]]
  }

  preds <- as.numeric(pf(object, batch))
  if (length(preds) != N)
    stop("prediction callback did not return one value per row for interactions.",
         call. = FALSE)

  dim(preds) <- c(L, np, nsim, n)
  delta <- preds[4, , , , drop = FALSE] - preds[2, , , , drop = FALSE] -
           preds[3, , , , drop = FALSE] + preds[1, , , , drop = FALSE]
  dim(delta) <- c(np, nsim, n)
  Iij <- apply(delta, c(1, 3), mean)            # np x n

  Phi <- array(0, c(n, p, p))
  for (pp2 in seq_len(np)) {
    v <- Iij[pp2, ] / 2
    Phi[, pairs[1, pp2], pairs[2, pp2]] <- v
    Phi[, pairs[2, pp2], pairs[1, pp2]] <- v
  }
  for (f in seq_len(p)) Phi[, f, f] <- phi_mat[, f] - rowSums(Phi[, f, ])
  dimnames(Phi) <- list(rownames(X), features, features)
  Phi
}

#' SHAP interaction plot
#'
#' Scatter of the pairwise SHAP interaction value between `feature` and a second
#' feature against the value of `feature`, coloured by the second feature. The
#' `shapvalues` object must have been produced with `interactions = TRUE`.
#'
#' @param x A `shapvalues` object carrying an interaction array.
#' @param feature Feature on the x axis.
#' @param interact Second feature. `NULL` (default) picks the feature with the
#'   largest mean absolute interaction with `feature`.
#' @param title Optional plot title.
#' @param style Colour style; `NULL` uses the active style (see [sv_style()]).
#' @return A `ggplot` object.
#' @examples
#' fit <- lm(mpg ~ wt * hp + disp + drat, data = mtcars)
#' sv <- shapvalues(fit, mtcars[, c("wt", "hp", "disp", "drat")],
#'                  mtcars[, c("wt", "hp", "disp", "drat")], nsim = 20,
#'                  interactions = TRUE, nthreads = 2)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(sv, type = "interaction", feature = "wt")
#' @noRd
sv_interaction <- function(x, feature, interact = NULL, title = NULL,
                           style = NULL) {
  .sv_need_ggplot()
  stopifnot(inherits(x, "shapvalues"))
  if (is.null(x$interactions))
    stop("no interaction values; call shapvalues(..., interactions = TRUE).",
         call. = FALSE)
  if (!feature %in% x$features) stop("unknown feature: ", feature, call. = FALSE)
  pal <- .sv_colours(style)
  Phi <- x$interactions
  fi <- match(feature, x$features)
  if (is.null(interact)) {
    score <- colMeans(abs(Phi[, fi, ]))
    score[fi] <- -Inf
    interact <- x$features[which.max(score)]
    message("interacting feature: ", interact)
  }
  if (!interact %in% x$features) stop("unknown feature: ", interact, call. = FALSE)
  ki <- match(interact, x$features)
  df <- data.frame(x = .sv_num(x$feature_values[[feature]]),
                   shap = Phi[, fi, ki] * 2,       # full pairwise effect
                   value01 = .sv_scale01(x$feature_values[[interact]]))
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = shap, colour = value01)) +
    ggplot2::geom_hline(yintercept = 0, colour = .sv_zero, linewidth = 0.4) +
    ggplot2::geom_point(size = 1.3, alpha = 0.9, shape = 16) +
    .sv_grad(pal, interact) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.08)) +
    ggplot2::labs(x = feature,
                  y = sprintf("SHAP interaction value for %s and %s", feature, interact),
                  title = title) +
    .sv_theme(y_axis = TRUE)
}
