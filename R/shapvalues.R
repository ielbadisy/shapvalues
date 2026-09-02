#' Compute SHAP values
#'
#' Sampling (permutation) estimator of Shapley additive explanations for the
#' predictions of any supervised model. All coalition rows are scored in one
#' batched prediction call; the marginal-contribution reduction runs in
#' parallel C++.
#'
#' @param object A fitted model, or a function. If a function is supplied it is
#'   used directly as the prediction callback (see `pred_fun`) and `X` becomes
#'   its first data argument.
#' @param X Data frame (or matrix) of rows to explain.
#' @param bg Background data frame used to marginalise absent features. Defaults
#'   to `X` (a random subsample of size `bg_size` is taken).
#' @param pred_fun Prediction callback with signature `function(object, newdata)`
#'   returning either a numeric vector (regression / single score), or a numeric
#'   matrix with one column per class (classification) or per time point
#'   (survival). If `NULL`, a callback is chosen from the class of `object`
#'   (`lm`, `glm`, `ranger`, `xgb.Booster`) or falls back to `stats::predict`.
#' @param features Character vector of columns in `X` to attribute. Defaults to
#'   all columns of `X` (minus survival outcome columns when detected).
#' @param nsim Number of random feature permutations per explained row.
#' @param method Estimator. Only `"permutation"` is currently implemented.
#' @param task One of `"auto"`, `"regression"`, `"classification"`,
#'   `"survival"`. `"auto"` inspects a probe prediction and the `times`
#'   argument.
#' @param times Numeric vector of evaluation times. Required for, and only used
#'   by, `task = "survival"`; forwarded to `pred_fun` as `times = times`.
#' @param bg_size Number of background rows to sample from `bg`.
#' @param antithetic If `TRUE`, permutations are drawn in reversed pairs to
#'   reduce variance (`nsim` is rounded up to an even number).
#' @param interactions If `TRUE`, also estimate the pairwise SHAP interaction
#'   matrix (one `n_feature x n_feature` matrix per explained row, stored in
#'   `$interactions`). Single-output models only; costs `O(n_feature^2)` more
#'   predictions.
#' @param nthreads Number of C++ threads for the reduction. `NULL` uses the
#'   'RcppParallel' default.
#' @param seed Optional integer seed for permutation and background sampling.
#' @param verbose If `TRUE`, report batch size and timing.
#'
#' @return An object of class `"shapvalues"`: a list with
#'   \item{phi}{numeric array `n_obs x n_feature x n_output` of attributions
#'     (drop to a matrix when `n_output == 1`).}
#'   \item{baseline}{length-`n_output` vector `E[f(x)]` over the background.}
#'   \item{prediction}{`n_obs x n_output` matrix of model predictions for `X`.}
#'   \item{feature_values}{the explained rows, `X[, features]`.}
#'   \item{features, task, method, outputs, times}{metadata.}
#'
#' @references Strumbelj, E. and Kononenko, I. (2014). Explaining prediction
#'   models and individual predictions with feature contributions.
#'   \emph{Knowledge and Information Systems} 41(3), 647-665.
#'   \doi{10.1007/s10115-013-0679-x}
#'
#' @examples
#' fit <- lm(mpg ~ ., data = mtcars)
#' sv <- shapvalues(fit, X = mtcars[1:5, -1], bg = mtcars[, -1], nsim = 20,
#'                  nthreads = 2)
#' sv
#'
#' @export
shapvalues <- function(object, X, bg = NULL, pred_fun = NULL, features = NULL,
                       nsim = 50L, method = c("permutation"),
                       task = c("auto", "regression", "classification", "survival"),
                       times = NULL, bg_size = 100L, antithetic = TRUE,
                       interactions = FALSE, nthreads = NULL, seed = NULL,
                       verbose = FALSE) {
  method <- match.arg(method)
  task <- match.arg(task)
  if (is.function(object)) {
    if (is.null(pred_fun)) pred_fun <- object
    object <- NULL
  }
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
    set.seed(seed)
    if (!is.null(old_seed)) on.exit(assign(".Random.seed", old_seed, .GlobalEnv), add = TRUE)
  }
  if (!is.null(nthreads)) {
    RcppParallel::setThreadOptions(numThreads = max(1L, as.integer(nthreads)))
  }

  X <- .sv_as_df(X)
  pf <- .sv_resolve_pred_fun(object, pred_fun)

  if (is.null(features)) features <- .sv_default_features(X, task, times)
  .sv_check_features(features, X)

  if (is.null(bg)) bg <- X else bg <- .sv_as_df(bg)
  bg <- .sv_align_bg(bg, X)
  bg_size <- min(as.integer(bg_size), nrow(bg))
  bg <- bg[sample.int(nrow(bg), bg_size), , drop = FALSE]

  nsim <- as.integer(nsim)
  if (nsim < 1L) stop("`nsim` must be >= 1.", call. = FALSE)
  if (antithetic && nsim %% 2L == 1L) nsim <- nsim + 1L

  n_obs <- nrow(X)
  p <- length(features)
  L <- p + 1L
  nblk <- n_obs * nsim

  ## permutations: (n_obs * nsim) x p, 0-based feature indices into `features`
  perm_order <- .sv_permutations(n_obs, nsim, p, antithetic)
  ## One background row per permutation, and the SAME nsim rows for every
  ## explained observation, so E[f(empty)] is identical across observations and
  ## sum_f phi[i, f] + baseline == prediction[i] holds exactly (this also gives
  ## common random numbers across observations). When nsim >= bg_size the whole
  ## background is tiled so the empirical mean equals the population mean;
  ## otherwise a systematic sample of the shuffled pool keeps it close.
  pool <- sample.int(nrow(bg))                     # shuffle once
  if (nsim >= length(pool)) {
    bg_idx <- rep_len(pool, nsim)
  } else {
    bg_idx <- pool[floor((seq_len(nsim) - 1) * length(pool) / nsim) + 1L]
  }
  bg_row <- rep(bg_idx, times = n_obs)

  ## rank of each feature within its permutation (1-based position)
  rank_of <- matrix(0L, nrow = nblk, ncol = p)
  seqblk <- seq_len(nblk)
  for (k in seq_len(p)) rank_of[cbind(seqblk, perm_order[, k] + 1L)] <- k

  ## batched coalition rows: N = nblk * L
  within <- rep.int(0:p, nblk)                    # coalition step per row
  blk_of <- rep(seqblk, each = L)                 # block index per row
  obs_of <- ((blk_of - 1L) %/% nsim) + 1L

  batch <- bg[bg_row[blk_of], , drop = FALSE]     # start: all-background
  rownames(batch) <- NULL
  for (j in seq_len(p)) {
    fcol <- features[j]
    on <- within >= rank_of[blk_of, j]            # feature j switched on
    if (any(on)) batch[[fcol]][on] <- X[[fcol]][obs_of[on]]
  }

  if (verbose) message(sprintf("shapvalues: %d coalition rows (%d obs x %d perms x %d)",
                               nrow(batch), n_obs, nsim, L))

  t0 <- proc.time()[["elapsed"]]
  preds <- .sv_predict(pf, object, batch, times, task)
  det <- .sv_detect_task(preds, task, times)
  task <- det$task
  preds <- det$preds
  n_out <- ncol(preds)

  phi <- shap_reduce_cpp(preds, perm_order, n_obs = n_obs, nsim = nsim,
                         p = p, n_out = n_out)
  baseline <- shap_baseline_cpp(preds, n_obs = n_obs, nsim = nsim, p = p)
  if (verbose) message(sprintf("shapvalues: done in %.2fs", proc.time()[["elapsed"]] - t0))

  ## model prediction for each explained row = last (all-on) coalition row
  full_rows <- (seqblk - 1L) * L + L
  pred_full <- preds[full_rows, , drop = FALSE]
  pred_obs <- rowsum(pred_full, group = obs_of[full_rows], reorder = FALSE) / nsim
  rownames(pred_obs) <- NULL

  out_names <- .sv_output_names(preds, det$task, times)
  dimnames(phi) <- list(rownames(X), features, out_names)
  names(baseline) <- out_names
  colnames(pred_obs) <- out_names

  inter <- NULL
  if (isTRUE(interactions)) {
    if (n_out != 1L)
      stop("interactions currently support single-output models ",
           "(regression or binary).", call. = FALSE)
    if (verbose) message("shapvalues: computing pairwise interactions")
    phi_mat <- matrix(phi[, , 1L], nrow = dim(phi)[1L],
                      dimnames = list(NULL, features))
    inter <- .sv_interactions(pf, object, X, bg, features,
                              phi_mat = phi_mat, nsim = nsim)
  }

  structure(
    list(
      phi = if (n_out == 1L) phi[, , 1L, drop = TRUE] else phi,
      phi_array = phi,
      baseline = baseline,
      prediction = pred_obs,
      feature_values = X[, features, drop = FALSE],
      features = features,
      task = task,
      method = method,
      outputs = out_names,
      times = if (identical(task, "survival")) times else NULL,
      nsim = nsim,
      interactions = inter,
      call = match.call()
    ),
    class = "shapvalues"
  )
}

## permutation generator ------------------------------------------------------

.sv_permutations <- function(n_obs, nsim, p, antithetic) {
  nblk <- n_obs * nsim
  po <- matrix(0L, nrow = nblk, ncol = p)
  if (p == 1L) return(po)  # single feature: order is trivial
  if (antithetic) {
    half <- nsim %/% 2L
    for (i in seq_len(n_obs)) {
      for (h in seq_len(half)) {
        ord <- sample.int(p) - 1L
        r1 <- (i - 1L) * nsim + (2L * h - 1L)
        r2 <- (i - 1L) * nsim + (2L * h)
        po[r1, ] <- ord
        po[r2, ] <- rev(ord)
      }
    }
  } else {
    for (r in seq_len(nblk)) po[r, ] <- sample.int(p) - 1L
  }
  po
}
