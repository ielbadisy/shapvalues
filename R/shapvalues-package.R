#' shapvalues: fast model-agnostic Shapley additive explanations
#'
#' Estimates SHAP values for any supervised model via a single prediction
#' callback. Regression, binary and multiclass classification, and
#' right-censored survival (time-dependent SHAP) are supported. The sampling
#' (permutation) estimator is used: coalition predictions are issued in one
#' batched call and reduced to attributions in parallel C++.
#'
#' @section Entry points:
#' \itemize{
#'   \item [shapvalues()] - compute attributions (and, with
#'     `interactions = TRUE`, the pairwise interaction matrix).
#'   \item [sv_beeswarm()], [sv_importance()], [sv_waterfall()],
#'     [sv_dependence()], [sv_interaction()], [sv_surv()] - standard plots.
#'   \item [sv_compass()], [svconcentration()], [sv_fingerprint()],
#'     [sv_pathway()] - multivariate attribution views.
#'   \item [svtable()], [svimportancetable()] - tidy tables.
#'   \item [svstyle()], [svpalette()] - the shared colour system.
#' }
#'
#' @keywords internal
#' @name shapvalues-package
#' @useDynLib shapvalues, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs setThreadOptions defaultNumThreads
#' @importFrom stats predict model.frame terms delete.response na.pass runif
"_PACKAGE"
