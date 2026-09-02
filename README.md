
<!-- README.md is generated from README.Rmd. Please edit that file -->

# shapvalues

<!-- badges: start -->
<!-- badges: end -->

`shapvalues` computes Shapley additive explanations (SHAP) for the
predictions of any supervised model through a single prediction
callback. It covers regression, binary and multiclass classification,
and right-censored survival analysis with time-dependent attributions.

The sampling (permutation) estimator is used: every coalition prediction
is issued in one batched call, and the marginal-contribution reduction
runs in parallel C++ via `RcppParallel`.

## Installation

``` r
# install.packages("pak")
pak::pak("ielbadisy/shapvalues")
```

## Usage

``` r
library(shapvalues)

fit <- lm(mpg ~ ., data = mtcars)
sv  <- shapvalues(fit, X = mtcars[, -1], bg = mtcars[, -1], nsim = 50)

sv
sv_importance_table(sv)
sv_beeswarm(sv)
sv_waterfall(sv, row = 1)
```

For a model class `shapvalues` does not know, pass a `pred_fun`:

``` r
pf <- function(object, newdata) predict(object, newdata, type = "response")
sv <- shapvalues(fit, X = X, bg = bg, pred_fun = pf, task = "classification")
```

Survival models return one attribution per feature per time point:

``` r
pf <- function(object, newdata, times) {
  s <- summary(survival::survfit(object, newdata = newdata),
               times = times, extend = TRUE)$surv
  t(matrix(s, nrow = length(times)))
}
sv <- shapvalues(cox_fit, X = X, bg = bg, pred_fun = pf,
                 times = c(30, 90, 180, 365), task = "survival")
sv_surv(sv, row = 1)
```
