#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
if (!length(file_arg)) stop("This script must be run with Rscript.", call. = FALSE)
repo_root <- normalizePath(file.path(dirname(file_arg), ".."))

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required to run this harness.", call. = FALSE)
}

devtools::load_all(repo_root, quiet = TRUE)

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

expect_error <- function(expr, pattern = NULL) {
  err <- NULL
  tryCatch(force(expr), error = function(e) err <<- e)
  if (is.null(err)) stop("Expected an error, but the expression succeeded.", call. = FALSE)
  if (!is.null(pattern) && !grepl(pattern, conditionMessage(err), fixed = TRUE)) {
    stop(
      sprintf("Error message did not match '%s': %s", pattern, conditionMessage(err)),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

is_symmetric <- function(x, tol = 1e-8) {
  max(abs(x - t(x))) < tol
}

is_spd <- function(x, tol = 1e-8) {
  x <- (x + t(x)) / 2
  vals <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  all(vals > tol)
}

offdiag_mean <- function(x) {
  idx <- row(x) != col(x)
  mean(abs(x[idx]))
}

fit_corrshrink <- function(formula, data, tree) {
  mvgls(
    formula,
    data = data,
    tree = tree,
    model = "BMM",
    method = "LL",
    REML = FALSE,
    bmm.structure = "corrshrink",
    echo = FALSE
  )
}

fit_proportional <- function(formula, data, tree) {
  mvgls(
    formula,
    data = data,
    tree = tree,
    model = "BMM",
    method = "LL",
    REML = FALSE,
    echo = FALSE
  )
}

set.seed(20260319)
n_tips <- 18
tree <- phytools::pbtree(n = n_tips, scale = 1)
states <- setNames(rep(c("A", "B"), each = n_tips / 2), tree$tip.label)
simmap <- suppressMessages(phytools::make.simmap(tree, states, model = "ER", nsim = 1))

sigma_base <- matrix(c(1.30, 0.55, 0.55, 1.00), 2, 2)
sigma_prop <- 1.8 * sigma_base
sigma_shrink <- matrix(c(1.70, 0.10, 0.10, 1.20), 2, 2)
theta <- c(0, 0)

make_response <- function(sigmas) {
  y <- mvSIM(
    simmap,
    nsim = 1,
    model = "BMM",
    param = list(ntraits = 2, sigma = sigmas, theta = theta)
  )
  if (is.list(y)) y <- y[[1]]
  y
}

Y_prop <- make_response(list(A = sigma_base, B = sigma_prop))
Y_shrink <- make_response(list(A = sigma_base, B = sigma_shrink))
X_shrink <- as.numeric(scale(rnorm(nrow(Y_shrink))))
names(X_shrink) <- rownames(Y_shrink)

fit_prop <- fit_proportional(Y ~ 1, data = list(Y = Y_prop), tree = simmap)
fit_corr_prop <- fit_corrshrink(Y ~ 1, data = list(Y = Y_prop), tree = simmap)

assert_true(identical(fit_corr_prop$bmm.structure, "corrshrink"), "corr-shrink fit did not mark the structure")
assert_true(all(c("df.free", "df.free_cov", "df.free_beta", "df.free_model") %in% names(fit_corr_prop)),
  "corr-shrink fit is missing one or more df.free fields"
)
assert_true(is.finite(fit_corr_prop$logLik), "corr-shrink fit produced a non-finite logLik")
assert_true(abs(as.numeric(fit_corr_prop$logLik) - as.numeric(fit_prop$logLik)) < 5,
  "corr-shrink logLik is not comparable to the proportional fit"
)

rho_names <- grep("\\.rho$", names(fit_corr_prop$param), value = TRUE)
assert_true(length(rho_names) >= 2, "expected at least one non-reference rho parameter")
assert_true(max(abs(fit_corr_prop$param[rho_names[-1]] - 1)) < 0.20,
  "proportional-case rho estimates did not remain close to 1"
)

regime_mats <- fit_corr_prop$sigma$regime
assert_true(length(regime_mats) >= 2, "expected regime covariance matrices in the fitted object")
for (nm in names(regime_mats)) {
  mat <- regime_mats[[nm]]
  assert_true(is_symmetric(mat), sprintf("regime covariance '%s' is not symmetric", nm))
  assert_true(is_spd(mat), sprintf("regime covariance '%s' is not positive definite", nm))
}

base_cor <- cov2cor(regime_mats[[1]])
for (nm in names(regime_mats)[-1]) {
  reg_cor <- cov2cor(regime_mats[[nm]])
  assert_true(max(abs(reg_cor - base_cor)) < 0.20,
    sprintf("regime '%s' is not close to the proportional correlation structure", nm)
  )
}

Y_pred <- Y_shrink
fit_corr_shrink <- fit_corrshrink(Y ~ 1, data = list(Y = Y_pred), tree = simmap)
fit_corr_reg <- fit_corrshrink(Y ~ x, data = list(Y = Y_shrink, x = X_shrink), tree = simmap)

assert_true(identical(fit_corr_reg$bmm.structure, "corrshrink"), "regression fit did not preserve corr-shrink mode")
assert_true(nrow(coef(fit_corr_reg)) == 2L, "unexpected coefficient row count for regression fit")
assert_true(ncol(coef(fit_corr_reg)) == ncol(Y_shrink), "unexpected coefficient column count for regression fit")
assert_true(identical(dim(fitted(fit_corr_reg)), dim(Y_shrink)), "fitted values have the wrong dimensions")
assert_true(identical(dim(residuals(fit_corr_reg)), dim(Y_shrink)), "response residuals have the wrong dimensions")

gic_reg <- GIC(fit_corr_reg)
aic_reg <- AIC(fit_corr_reg)
assert_true(is.finite(gic_reg$GIC), "GIC returned a non-finite value")
assert_true(is.finite(aic_reg$AIC), "AIC returned a non-finite value")
assert_true(abs(gic_reg$GIC - aic_reg$AIC) < 1e-8, "GIC and AIC should coincide for corr-shrink ML fits")
assert_true(gic_reg$bias_cov == fit_corr_reg$df.free_cov, "GIC bias_cov does not match df.free_cov")
assert_true(gic_reg$bias == fit_corr_reg$df.free, "GIC bias does not match the total free parameter count")
assert_true(
  fit_corr_reg$df.free == fit_corr_reg$df.free_cov + fit_corr_reg$df.free_beta + fit_corr_reg$df.free_model,
  "df.free fields do not add up correctly"
)

summary_reg <- summary(fit_corr_reg)
assert_true(inherits(summary_reg, "summary.mvgls"), "summary() did not return a summary.mvgls object")
assert_true(all(c("AIC", "GIC", "logLik") %in% colnames(summary_reg$results.fit)),
  "summary output is missing the expected fit statistics"
)

regime_mats_shrink <- fit_corr_shrink$sigma$regime
assert_true(length(regime_mats_shrink) >= 2, "expected regime covariance matrices in the shrink fit")
base_cor_shrink <- cov2cor(regime_mats_shrink[[1]])
for (nm in names(regime_mats_shrink)) {
  mat <- regime_mats_shrink[[nm]]
  assert_true(is_symmetric(mat), sprintf("shrink regime covariance '%s' is not symmetric", nm))
  assert_true(is_spd(mat), sprintf("shrink regime covariance '%s' is not positive definite", nm))
}
rho_shrink <- fit_corr_shrink$param[grep("\\.rho$", names(fit_corr_shrink$param))]
assert_true(any(rho_shrink[-1] < 0.95), "shrink case did not move rho away from 1")
assert_true(offdiag_mean(cov2cor(regime_mats_shrink[[names(regime_mats_shrink)[2]]])) < offdiag_mean(base_cor_shrink),
  "shrink case did not reduce the off-diagonal correlation strength"
)

expect_error(
  mvgls(Y ~ 1, data = list(Y = Y_shrink), tree = simmap, model = "BMM", method = "PL-LOOCV",
        REML = FALSE, bmm.structure = "corrshrink", echo = FALSE),
  "corr-shrink"
)
expect_error(
  mvgls(Y ~ 1, data = list(Y = Y_shrink), tree = simmap, model = "BMM", method = "H&L",
        REML = FALSE, bmm.structure = "corrshrink", echo = FALSE),
  "corr-shrink"
)
expect_error(
  mvgls(Y ~ 1, data = list(Y = Y_shrink), tree = simmap, model = "BMM", method = "Mahalanobis",
        REML = FALSE, bmm.structure = "corrshrink", echo = FALSE),
  "corr-shrink"
)
expect_error(
  mvgls(Y ~ 1, data = list(Y = Y_shrink), tree = simmap, model = "BMM", method = "EmpBayes",
        REML = FALSE, bmm.structure = "corrshrink", echo = FALSE),
  "corr-shrink"
)
expect_error(
  mvgls(Y ~ 1, data = list(Y = Y_shrink), tree = simmap, model = "BMM", method = "LL",
        REML = TRUE, bmm.structure = "corrshrink", echo = FALSE),
  "corr-shrink"
)
expect_error(
  mvgls(Y ~ 1, data = list(Y = Y_shrink), tree = simmap, model = "BMM", method = "LL",
        REML = FALSE, error = TRUE, bmm.structure = "corrshrink", echo = FALSE),
  "corr-shrink"
)

expect_error(predict(fit_corr_reg), "corr-shrink")
expect_error(ancestral(fit_corr_reg), "corr-shrink")
expect_error(residuals(fit_corr_reg, type = "normalized"), "corr-shrink")
expect_error(vcov(fit_corr_reg, type = "coef"), "corr-shrink")

sim_one <- simulate(fit_corr_reg, nsim = 1)
assert_true(is.matrix(sim_one), "simulate(..., nsim=1) should return a matrix")
assert_true(identical(dim(sim_one), dim(Y_shrink)), "simulate(..., nsim=1) returned the wrong dimensions")
assert_true(identical(rownames(sim_one), rownames(Y_shrink)), "simulate(..., nsim=1) did not preserve row names")
assert_true(all(is.finite(sim_one)), "simulate(..., nsim=1) returned non-finite values")

sim_three <- simulate(fit_corr_reg, nsim = 3)
assert_true(is.list(sim_three), "simulate(..., nsim=3) should return a list")
assert_true(length(sim_three) == 3L, "simulate(..., nsim=3) returned the wrong number of replicates")
for (i in seq_along(sim_three)) {
  sim_i <- sim_three[[i]]
  assert_true(is.matrix(sim_i), sprintf("simulate replicate %d is not a matrix", i))
  assert_true(identical(dim(sim_i), dim(Y_shrink)), sprintf("simulate replicate %d has the wrong dimensions", i))
  assert_true(identical(rownames(sim_i), rownames(Y_shrink)), sprintf("simulate replicate %d did not preserve row names", i))
  assert_true(all(is.finite(sim_i)), sprintf("simulate replicate %d returned non-finite values", i))
}

eic_base <- EIC(fit_corr_reg, nboot = 4L, nbcores = 1L)
assert_true(is.finite(eic_base$EIC), "EIC returned a non-finite criterion")
assert_true(is.finite(eic_base$LogLikelihood), "EIC returned a non-finite log-likelihood")
assert_true(is.finite(eic_base$se), "EIC returned a non-finite standard error")
assert_true(length(eic_base$bias) >= 1L, "EIC did not return bootstrap bias samples")

eic_restricted <- EIC(fit_corr_reg, nboot = 4L, nbcores = 1L, restricted = TRUE)
assert_true(is.finite(eic_restricted$EIC), "restricted EIC returned a non-finite criterion")
assert_true(is.finite(eic_restricted$LogLikelihood), "restricted EIC returned a non-finite log-likelihood")
assert_true(is.finite(eic_restricted$se), "restricted EIC returned a non-finite standard error")

eic_warn <- NULL
eic_sqmf <- withCallingHandlers(
  EIC(fit_corr_reg, nboot = 4L, nbcores = 1L, eigSqm = FALSE),
  warning = function(w) {
    eic_warn <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  }
)
assert_true(!is.null(eic_warn), "EIC(..., eigSqm = FALSE) should warn")
assert_true(is.finite(eic_sqmf$EIC), "EIC(..., eigSqm = FALSE) returned a non-finite criterion")
assert_true(is.finite(eic_sqmf$LogLikelihood), "EIC(..., eigSqm = FALSE) returned a non-finite log-likelihood")
assert_true(is.finite(eic_sqmf$se), "EIC(..., eigSqm = FALSE) returned a non-finite standard error")

cat("corr-shrink mvgls docs/test harness checks passed\n")
