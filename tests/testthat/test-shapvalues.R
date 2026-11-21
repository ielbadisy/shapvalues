set.seed(1)

test_that("regression: SHAP satisfies local accuracy", {
  fit <- lm(mpg ~ ., data = mtcars)
  sv <- shapvalues(fit, X = mtcars[, -1], bg = mtcars[, -1],
                   nsim = 60, nthreads = 2, seed = 1)
  expect_s3_class(sv, "shapvalues")
  expect_equal(sv$task, "regression")
  expect_equal(dim(sv$phi), c(nrow(mtcars), ncol(mtcars) - 1L))
  # local accuracy: sum of phi + baseline == prediction, per observation
  recon <- rowSums(sv$phi) + sv$baseline[1]
  expect_equal(unname(recon), unname(sv$prediction[, 1]), tolerance = 1e-6)
})

test_that("function interface works without a model object", {
  fit <- lm(mpg ~ ., data = mtcars)
  pf <- function(object, newdata) predict(fit, newdata)
  sv <- shapvalues(pf, X = mtcars[1:5, -1], bg = mtcars[, -1],
                   nsim = 20, nthreads = 2, seed = 2)
  expect_equal(dim(sv$phi), c(5L, ncol(mtcars) - 1L))
})

test_that("classification: multi-output array and per-class local accuracy", {
  fit <- suppressWarnings(
    glm(vs ~ mpg + hp + wt + drat, data = mtcars, family = "binomial"))
  pf <- function(object, newdata) {
    p <- suppressWarnings(predict(fit, newdata, type = "response"))
    cbind(no = 1 - p, yes = p)
  }
  vars <- c("mpg", "hp", "wt", "drat")
  sv <- shapvalues(pf, X = mtcars[, vars], bg = mtcars[, vars],
                   nsim = 40, task = "classification", nthreads = 2, seed = 3)
  expect_equal(sv$task, "classification")
  expect_equal(dim(sv$phi), c(nrow(mtcars), 4L, 2L))
  recon <- rowSums(sv$phi[, , 2]) + sv$baseline[2]
  expect_equal(unname(recon), unname(sv$prediction[, 2]), tolerance = 1e-6)
  # attributions across the two classes are equal and opposite
  expect_equal(sv$phi[, , 1], -sv$phi[, , 2], tolerance = 1e-6)
})

test_that("survival: time-dependent SHAP shape and times metadata", {
  skip_if_not_installed("survival")
  data(cancer, package = "survival", envir = environment())
  fit <- survival::coxph(survival::Surv(time, status) ~ age + karno + diagtime,
                         data = veteran, x = TRUE)
  tt <- c(30, 90, 180, 365)
  pf <- function(object, newdata, times) {
    sf <- survival::survfit(object, newdata = newdata)
    s <- summary(sf, times = times, extend = TRUE)$surv
    t(matrix(s, nrow = length(times)))
  }
  sv <- shapvalues(fit, X = veteran[1:6, c("age", "karno", "diagtime")],
                   bg = veteran[, c("age", "karno", "diagtime")],
                   pred_fun = pf, times = tt, task = "survival",
                   nsim = 20, nthreads = 2, seed = 4)
  expect_equal(sv$task, "survival")
  expect_equal(dim(sv$phi), c(6L, 3L, 4L))
  expect_equal(sv$times, tt)
})

test_that("tables and importance are consistent", {
  fit <- lm(mpg ~ ., data = mtcars)
  sv <- shapvalues(fit, mtcars[1:8, -1], mtcars[, -1], nsim = 20,
                   nthreads = 2, seed = 5)
  tab <- svtable(sv)
  expect_equal(nrow(tab), 8L * (ncol(mtcars) - 1L))
  expect_true(all(c("observation", "feature", "phi", "feature_value") %in% names(tab)))
  imp <- svimportancetable(sv)
  expect_equal(nrow(imp), ncol(mtcars) - 1L)
  expect_true(!is.unsorted(rev(imp$mean_abs)))
})

test_that("plots build when ggplot2 is available", {
  skip_if_not_installed("ggplot2")
  fit <- lm(mpg ~ ., data = mtcars)
  sv <- shapvalues(fit, mtcars[, -1], mtcars[, -1], nsim = 20,
                   nthreads = 2, seed = 6)
  expect_s3_class(sv_beeswarm(sv), "ggplot")
  expect_s3_class(sv_importance(sv), "ggplot")
  expect_s3_class(sv_waterfall(sv, row = 1), "ggplot")
  expect_s3_class(sv_dependence(sv, "wt", colour_by = "hp"), "ggplot")
})

test_that("interaction matrix rows sum to the SHAP values", {
  vars <- c("wt", "hp", "disp", "drat", "qsec")
  fit <- lm(mpg ~ wt * hp + disp + drat + qsec, data = mtcars)
  sv <- shapvalues(fit, mtcars[, vars], mtcars[, vars], nsim = 80,
                   interactions = TRUE, nthreads = 2, seed = 7)
  Phi <- sv$interactions
  expect_equal(dim(Phi), c(nrow(mtcars), 5L, 5L))
  # symmetry
  expect_equal(Phi, aperm(Phi, c(1, 3, 2)), tolerance = 1e-8)
  # each row of the interaction matrix sums to the SHAP value
  rowsum_phi <- apply(Phi, c(1, 2), sum)
  expect_equal(unname(rowsum_phi), unname(as.matrix(sv$phi)), tolerance = 1e-6)
  # everything sums to prediction - baseline
  total <- apply(Phi, 1, sum) + sv$baseline[1]
  expect_equal(unname(total), unname(sv$prediction[, 1]), tolerance = 1e-6)
  skip_if_not_installed("ggplot2")
  expect_s3_class(sv_interaction(sv, "wt", "hp"), "ggplot")
})

test_that("colour styles resolve", {
  expect_named(svpalette(), c("ink", "axis", "line", "zero", "grid", "style",
                               "neg", "pos", "mid", "bar", "grad_low",
                               "grad_high", "ramp", "qualitative"),
               ignore.order = TRUE)
  expect_length(svpalette("plasma")$ramp(7), 7L)
  old <- svstyle("inferno")
  expect_identical(getOption("shapvalues.style"), "inferno")
  svstyle("shapvalues")
})

test_that("concentration and the multivariate views build", {
  skip_if_not_installed("ggplot2")
  fit <- lm(mpg ~ ., data = mtcars)
  sv <- shapvalues(fit, mtcars[, -1], mtcars[, -1], nsim = 20,
                   nthreads = 2, seed = 8)
  ec <- svconcentration(sv)
  expect_equal(nrow(ec), nrow(mtcars))
  expect_true(all(ec$ec >= 0 & ec$ec <= 1 + 1e-8))
  expect_true(all(c("observation", "cx", "cy", "ec", "dominant") %in% names(ec)))
  expect_s3_class(sv_compass(sv), "ggplot")
  expect_s3_class(sv_compass(sv, groups = list(a = c("hp", "disp", "cyl"),
                                               b = c("wt", "drat", "qsec"))),
                  "ggplot")
  expect_s3_class(sv_fingerprint(sv, k = 3), "ggplot")
  expect_s3_class(sv_fingerprint(sv, rows = c(1, 5)), "ggplot")
  skip_if_not_installed("ggalluvial")
  expect_s3_class(sv_pathway(sv, rank = 3), "ggplot")
  expect_s3_class(sv_pathway(sv, rank = 2, signed = TRUE), "ggplot")
})

test_that("single explained row keeps phi as a 1-row matrix", {
  fit <- lm(mpg ~ ., data = mtcars)
  sv <- shapvalues(fit, X = mtcars[1, -1], bg = mtcars[, -1],
                   nsim = 20, nthreads = 2, seed = 1)
  expect_true(is.matrix(sv$phi))
  expect_equal(dim(sv$phi), c(1L, ncol(mtcars) - 1L))
  expect_equal(colnames(sv$phi), colnames(mtcars)[-1])
  expect_output(print(sv), "1 obs")
  expect_s3_class(sv_waterfall(sv), "ggplot")
  expect_s3_class(sv_importance(sv), "ggplot")
  expect_equal(nrow(svimportancetable(sv)), ncol(mtcars) - 1L)
})

test_that("plot() returns an unprinted ggplot with editable text and shared colours", {
  skip_if_not_installed("ggplot2")
  fit <- lm(mpg ~ ., data = mtcars)
  sv <- shapvalues(fit, mtcars[, -1], mtcars[, -1], nsim = 20, nthreads = 1)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(g <- plot(sv, type = "waterfall", row = 1))
  expect_s3_class(g, "ggplot")
  expect_null(g$labels$title)
  expect_null(g$labels$subtitle)
  g2 <- plot(sv, type = "importance", title = "T", subtitle = "S", caption = "C") +
    ggplot2::labs(x = "Mean |SHAP| (mpg)")
  expect_identical(unlist(g2$labels[c("title", "subtitle", "caption", "x")], use.names = FALSE),
                   c("T", "S", "C", "Mean |SHAP| (mpg)"))
  expect_null(suppressMessages(plot(sv, type = "compass"))$labels$subtitle)
  pal <- svpalette("shapvalues")
  expect_identical(c(pal$neg, pal$pos), c("#0072B2", "#D55E00"))
  b <- ggplot2::ggplot_build(g)
  expect_setequal(unique(b$data[[2]]$fill), c("#0072B2", "#D55E00"))
})

test_that("waterfall labels read 'feature = value' and tiny contributions are not '+0' / '-0'", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  d <- transform(mtcars, z = rnorm(32))
  d$y <- 30 - 3 * d$wt - 0.03 * d$hp + 0.002 * d$z
  X <- d[c("wt", "hp", "z")]
  sv <- shapvalues(lm(y ~ wt + hp + z, d), X, X, nsim = 20, nthreads = 1, seed = 1)
  g <- plot(sv, type = "waterfall", row = 3)
  yl <- ggplot2::ggplot_build(g)$layout$panel_params[[1]]$y$get_labels()
  expect_true(all(c("wt = 2.32", "hp = 93") %in% yl))
  lb <- ggplot2::ggplot_build(g)$data[[4]]$label
  expect_false(any(lb %in% c("+0", "-0", "+0.00", "-0.00")))
  expect_equal(.sv_signed(c(-0.0017, 3.191, 0), 5), c("-0.002", "+3.19", "0"))
})
