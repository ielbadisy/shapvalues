# Benchmark shapvalues against other model-agnostic SHAP implementations.
#
#   Rscript inst/benchmarks/bench-shapvalues.R          # quick
#   Rscript inst/benchmarks/bench-shapvalues.R --full   # full grid (slow)
#
# Models explained: ranger, xgboost, fastgbm (a compiled histogram GBM).
# All three share the hyperparameters that can be shared, so that differences
# between them reflect the model family, not the tuning:
#   200 trees, max depth 4, min node size 10, all features per split,
#   no row/column subsampling, L2 = 1, 255 histogram bins.
# Model-intrinsic knobs that cannot be shared: the boosters use shrinkage
# lr = 0.1 (a random forest has none); the forest keeps bootstrap resampling
# (structurally required -- without it, depth-matched full-feature trees would
# be identical). ranger also splits on exact values rather than bins.
# All SHAP methods get the identical (model, rows, background) setup.
#
# Studies:
#   1. head-to-head  : median runtime + accuracy per model and feature count for
#                        shapvalues (1 and 4 threads),
#                        kernelshap (default: exact <=8 features, hybrid above),
#                        fastshap   (C++ permutation; CRAN-archived),
#                        iml        (Shapley, pure R, per row),
#                        treeshap   (path-dependent TreeSHAP; tree models only).
#   2. convergence   : RMSE and time vs nsim for the sampling estimators
#                        (shapvalues, fastshap).
#   3. thread scaling: shapvalues wall time at 1/2/4/8 threads, many-output model.
#
# Ground truth is INTERVENTIONAL SHAP, obtained from kernelshap(exact = TRUE)
# over the full background on a capped subset of rows. shapvalues, fastshap and
# hybrid kernelshap all target this same estimand; treeshap targets the
# path-dependent variant, so its small gap here is expected (features are
# independent by construction, which keeps that gap small).
#
# Writes inst/benchmarks/bench-shapvalues-results.rds.

suppressPackageStartupMessages(library(shapvalues))
stopifnot(requireNamespace("bench", quietly = TRUE))

FULL    <- "--full" %in% commandArgs(TRUE)
set.seed(20260901)
has     <- function(p) requireNamespace(p, quietly = TRUE)
NSIM    <- 50L
ITERS   <- if (FULL) 8L else 2L
IML_MAX <- if (FULL) 200L else 50L    # iml is per row
KS_MAX  <- if (FULL) 1500L else 150L  # kernelshap time ~ n_expl * n_bg * coalitions
N_TRUTH <- if (FULL) 80L  else 50L    # rows scored for the exact reference

## ------------------------------------------------------------------ data/models

make_data <- function(n, p) {
  X <- as.data.frame(matrix(rnorm(n * p), n, p))
  names(X) <- paste0("x", seq_len(p))
  beta <- rep_len(c(2, -1.5, 1, -0.5), p)
  X$y <- as.numeric(as.matrix(X) %*% beta + X$x1 * X$x2) + rnorm(n, 0, 0.3)
  X
}

## shared model hyperparameters
HP <- list(ntrees = 200L, max_depth = 4L, min_node = 10L, lr = 0.1,
           subsample = 1, colsample = 1, lambda = 1, nbins = 255L)

build_model <- function(model, BG, feats) {
  if (model == "ranger") {
    fit <- ranger::ranger(y ~ ., data = BG, num.trees = HP$ntrees,
                          max.depth = HP$max_depth, min.node.size = HP$min_node,
                          mtry = length(feats),            # all features per split
                          replace = TRUE, sample.fraction = 1,  # forest keeps bagging
                          num.threads = 1, seed = 1)
    list(fit = fit, tree = TRUE,
         pf = function(object, newdata) predict(object, newdata, num.threads = 1)$predictions,
         unified = treeshap::ranger.unify(fit, BG[, feats]))
  } else if (model == "xgboost") {
    dm <- xgboost::xgb.DMatrix(as.matrix(BG[, feats]), label = BG$y)
    fit <- xgboost::xgb.train(
      list(objective = "reg:squarederror", tree_method = "hist",
           max_depth = HP$max_depth, eta = HP$lr, min_child_weight = HP$min_node,
           subsample = HP$subsample, colsample_bytree = HP$colsample,
           colsample_bynode = HP$colsample, lambda = HP$lambda,
           max_bin = HP$nbins, nthread = 1),
      dm, nrounds = HP$ntrees, verbose = 0)
    list(fit = fit, tree = TRUE,
         pf = function(object, newdata) predict(object, as.matrix(newdata[, feats])),
         unified = treeshap::xgboost.unify(fit, as.matrix(BG[, feats])))
  } else {                                                          # fastgbm
    fit <- fastgbm::fastgbm(
      x = BG[, feats], y = BG$y, objective = "regression",
      ntrees = HP$ntrees, learning_rate = HP$lr, max_depth = HP$max_depth,
      min_node_size = HP$min_node, subsample = HP$subsample,
      colsample = HP$colsample, lambda = HP$lambda, min_child_weight = 1,
      max_bins = HP$nbins, threads = 1, seed = 1, verbose = FALSE)
    list(fit = fit, tree = FALSE, unified = NULL,
         pf = function(object, newdata)
           as.numeric(predict(object, newdata[, feats, drop = FALSE])))
  }
}

rmse  <- function(a, b) sqrt(mean((as.matrix(a) - as.matrix(b))^2))
maxae <- function(a, b) max(abs(as.matrix(a) - as.matrix(b)))
## `thunk` MUST be a zero-arg function so bench::mark re-runs it each iteration.
med_time <- function(thunk, iters) {
  gc(FALSE)
  as.numeric(bench::mark(thunk(), iterations = iters, check = FALSE,
                         filter_gc = FALSE, memory = FALSE)$median)
}

## ------------------------------------------------------------------ SHAP runners

sv_run <- function(m, Xe, BG, feats, nsim, th) {
  phi <- shapvalues(m$fit, Xe[, feats], BG[, feats], pred_fun = m$pf,
                    nsim = nsim, bg_size = nrow(BG), nthreads = th, seed = 1)$phi
  if (length(dim(phi)) == 2L) as.matrix(phi)[, feats, drop = FALSE] else phi
}

ks_run <- function(m, Xe, BG, feats, exact = FALSE) {
  S <- kernelshap::kernelshap(m$fit, X = Xe[, feats], bg_X = BG[, feats],
                              exact = exact, verbose = FALSE,
                              pred_fun = function(o, X) m$pf(o, X))$S
  `colnames<-`(as.matrix(S), feats)
}

fs_run <- function(m, Xe, BG, feats, nsim)
  `colnames<-`(as.matrix(fastshap::explain(
    m$fit, X = BG[, feats], nsim = nsim, newdata = Xe[, feats], shap_only = TRUE,
    pred_wrapper = function(o, newdata) m$pf(o, newdata))), feats)

iml_run <- function(m, Xe, BG, feats, nsim) {
  pr <- iml::Predictor$new(m$fit, data = BG[, feats],
                           predict.function = function(model, newdata) m$pf(model, newdata))
  t(vapply(seq_len(nrow(Xe)), function(i) {
    r <- iml::Shapley$new(pr, x.interest = Xe[i, feats, drop = FALSE],
                          sample.size = nsim)$results
    r$phi[match(feats, r$feature)]
  }, numeric(length(feats))))
}

ts_run <- function(m, Xe, feats)
  as.matrix(treeshap::treeshap(m$unified, Xe[, feats], verbose = FALSE)$shaps)[, feats, drop = FALSE]

## ================================================================= 1. head to head

hh_grid <- if (FULL) {
  expand.grid(model = c("ranger", "xgboost", "fastgbm"), p = c(5L, 10L, 16L),
              stringsAsFactors = FALSE)
} else {
  expand.grid(model = c("ranger", "xgboost", "fastgbm"), p = c(5L, 10L),
              stringsAsFactors = FALSE)
}
N_BG   <- if (FULL) 200L else 100L
N_EXPL <- if (FULL) 2000L else 400L

message(sprintf("[1/3] head-to-head: %d scenarios, n_bg=%d, n_expl=%d, nsim=%d",
                nrow(hh_grid), N_BG, N_EXPL, NSIM))

headtohead <- do.call(rbind, lapply(seq_len(nrow(hh_grid)), function(i) {
  model <- hh_grid$model[i]; p <- hh_grid$p[i]
  feats <- paste0("x", seq_len(p))
  dat <- make_data(N_BG + N_EXPL, p)
  BG <- dat[seq_len(N_BG), ]; Xe <- dat[N_BG + seq_len(N_EXPL), ]
  m <- build_model(model, BG, feats)
  message(sprintf("   %s p=%d", model, p))

  Xt <- Xe[seq_len(min(N_TRUTH, nrow(Xe))), , drop = FALSE]   # accuracy subset
  truth <- ks_run(m, Xt, BG, feats, exact = TRUE)             # interventional exact

  Xi <- Xe[seq_len(min(IML_MAX, nrow(Xe))), , drop = FALSE]
  Xk <- Xe[seq_len(min(KS_MAX,  nrow(Xe))), , drop = FALSE]

  runners <- list(
    `shapvalues (1c)` = list(val = function() sv_run(m, Xt, BG, feats, NSIM, 1L),
                             tim = function() sv_run(m, Xe, BG, feats, NSIM, 1L), it = ITERS),
    `shapvalues (4c)` = list(val = function() sv_run(m, Xt, BG, feats, NSIM, 4L),
                             tim = function() sv_run(m, Xe, BG, feats, NSIM, 4L), it = ITERS),
    kernelshap        = list(val = function() ks_run(m, Xt, BG, feats),
                             tim = function() ks_run(m, Xk, BG, feats), it = ITERS),
    fastshap          = list(val = function() fs_run(m, Xt, BG, feats, NSIM),
                             tim = function() fs_run(m, Xe, BG, feats, NSIM), it = ITERS),
    iml               = list(val = function() iml_run(m, Xt, BG, feats, NSIM),
                             tim = function() iml_run(m, Xi, BG, feats, NSIM), it = 1L)
  )
  if (isTRUE(m$tree))
    runners$treeshap <- list(val = function() ts_run(m, Xt, feats),
                             tim = function() ts_run(m, Xe, feats), it = ITERS)

  do.call(rbind, lapply(names(runners), function(nm) {
    r <- runners[[nm]]
    val <- tryCatch(r$val(), error = function(e) {
      message("     ", nm, " failed: ", conditionMessage(e)); NULL })
    if (is.null(val)) return(NULL)
    n_tim <- if (nm == "iml") nrow(Xi) else if (nm == "kernelshap") nrow(Xk) else N_EXPL
    tm <- tryCatch(med_time(r$tim, r$it), error = function(e) NA_real_)
    data.frame(model = model, p = p, method = nm,
               time_s = tm, ms_per_row = 1000 * tm / n_tim,
               rmse = rmse(val, truth), max_ae = maxae(val, truth),
               stringsAsFactors = FALSE)
  }))
}))
key <- with(headtohead, paste(model, p))
headtohead$rel_speed <- ave(headtohead$ms_per_row, key,
                            FUN = function(z) z / min(z, na.rm = TRUE))

## ================================================================= 2. convergence

message("[2/3] convergence vs nsim (xgboost, p=10)")
p <- 10L; feats <- paste0("x", seq_len(p))
NC <- if (FULL) 1000L else 400L
dat <- make_data(200L + NC, p)
BG <- dat[seq_len(200L), ]; Xe <- dat[200L + seq_len(NC), ]
m <- build_model("xgboost", BG, feats)
Xt <- Xe[seq_len(min(N_TRUTH, nrow(Xe))), , drop = FALSE]
truth <- ks_run(m, Xt, BG, feats, exact = TRUE)

nsims <- if (FULL) c(10L, 25L, 50L, 100L, 200L, 400L) else c(10L, 25L, 50L, 100L, 200L)
convergence <- do.call(rbind, lapply(nsims, function(M) {
  sv <- sv_run(m, Xt, BG, feats, M, 4L)
  fs <- fs_run(m, Xt, BG, feats, M)
  rbind(
    data.frame(method = "shapvalues", nsim = M, rmse = rmse(sv, truth),
               time_s = med_time(function() sv_run(m, Xe, BG, feats, M, 4L), ITERS)),
    data.frame(method = "fastshap", nsim = M, rmse = rmse(fs, truth),
               time_s = med_time(function() fs_run(m, Xe, BG, feats, M), ITERS))
  )
}))
attr(convergence, "kernelshap_time") <-
  med_time(function() ks_run(m, Xe[seq_len(min(KS_MAX, nrow(Xe))), ], BG, feats), ITERS)

saveRDS(list(headtohead = headtohead, convergence = convergence),
        "inst/benchmarks/bench-shapvalues-results.rds")

## ================================================================= 3. thread scaling
## The parallel C++ reduction only becomes the bottleneck (rather than the
## R-side predict call) when the model has many outputs -- multiclass, or
## survival over a fine time grid. A cheap analytic many-output "model"
## isolates that effect.

NT   <- if (FULL) 4000L else 2500L
NOUT <- 40L
p <- 18L; feats <- paste0("x", seq_len(p))
message(sprintf("[3/3] thread scaling (analytic %d-output model, p=%d, n_expl=%d, nsim=%d)",
                NOUT, p, NT, NSIM))
dat <- make_data(200L + NT, p)
BG <- dat[seq_len(200L), ]; Xe <- dat[200L + seq_len(NT), ]
W <- matrix(rnorm(p * 3), p, 3); tgrid <- seq(-2, 2, length.out = NOUT)
m <- list(fit = NULL, tree = FALSE, unified = NULL,
          pf = function(object, newdata) {
            z <- rowSums((as.matrix(newdata[, feats]) %*% W)^2)
            plogis(outer(z, tgrid, `-`))
          })
threads <- c(1L, 2L, 4L, 8L)
threads <- threads[threads <= max(1L, parallel::detectCores())]
thread_scaling <- do.call(rbind, lapply(threads, function(th)
  data.frame(threads = th, n_out = NOUT,
             time_s = med_time(function() sv_run(m, Xe, BG, feats, NSIM, th), ITERS))))
thread_scaling$speedup <- thread_scaling$time_s[1] / thread_scaling$time_s

## ------------------------------------------------------------------ save + print

saveRDS(list(
  headtohead = headtohead, convergence = convergence, thread_scaling = thread_scaling,
  meta = list(full = FULL, iters = ITERS, nsim = NSIM, n_bg = N_BG, hp = HP,
              date = Sys.Date(),
              versions = vapply(c("shapvalues", "kernelshap", "fastshap", "iml",
                                  "treeshap", "fastgbm"),
                                function(p) as.character(packageVersion(p)), ""))),
  "inst/benchmarks/bench-shapvalues-results.rds")

cat("\n===============  1. HEAD TO HEAD  ===============\n")
print(headtohead[order(headtohead$model, headtohead$p, headtohead$ms_per_row),
                 c("model", "p", "method", "time_s", "ms_per_row", "rel_speed",
                   "rmse", "max_ae")], row.names = FALSE, digits = 3)
cat("\n===============  2. CONVERGENCE (xgboost, p=10, vs interventional exact)  ===============\n")
print(convergence[order(convergence$method, convergence$nsim), ], row.names = FALSE, digits = 3)
cat(sprintf("kernelshap (default) time on %d rows: %.3fs\n",
            min(KS_MAX, NC), attr(convergence, "kernelshap_time")))
cat("\n===============  3. THREAD SCALING (analytic 40-output, p=18)  ===============\n")
print(thread_scaling, row.names = FALSE, digits = 3)
cat("\nsaved -> inst/benchmarks/bench-shapvalues-results.rds\n")

## ------------------------------------------------------------------ plots

if (has("ggplot2")) {
  suppressPackageStartupMessages(library(ggplot2))
  hh <- headtohead; hh$method <- reorder(hh$method, hh$ms_per_row)
  ggsave("inst/benchmarks/figures/runtime.png",
    ggplot(hh, aes(method, ms_per_row, fill = method)) +
      geom_col(width = .7, show.legend = FALSE) +
      facet_wrap(~ paste0(model, " (p=", p, ")"), scales = "free_x") +
      scale_y_log10() + coord_flip() +
      labs(title = "SHAP runtime per explained row (log scale)",
           x = NULL, y = "ms / row") + theme_minimal(base_size = 10),
    width = 9, height = 6, dpi = 120)

  cv <- convergence
  ggsave("inst/benchmarks/figures/convergence.png",
    ggplot(cv, aes(nsim, pmax(rmse, 1e-16), colour = method)) +
      geom_line(linewidth = .8) + geom_point() +
      scale_x_log10() + scale_y_log10() +
      labs(title = "Error vs interventional-exact SHAP as nsim grows (xgboost, p=10)",
           x = "nsim", y = "RMSE vs exact") + theme_minimal(base_size = 10),
    width = 7, height = 4.5, dpi = 120)

  ts <- thread_scaling
  ggsave("inst/benchmarks/figures/threads.png",
    ggplot(ts, aes(factor(threads), speedup)) +
      geom_col(width = .6, fill = "#2c7fb8") +
      geom_text(aes(label = sprintf("%.2fx", speedup)), vjust = -0.4, size = 3) +
      labs(title = "shapvalues thread scaling (analytic 40-output model)",
           x = "threads", y = "speedup vs 1 thread") + theme_minimal(base_size = 10),
    width = 6, height = 4, dpi = 120)
  cat("plots -> inst/benchmarks/figures/{runtime,convergence,threads}.png\n")
}
