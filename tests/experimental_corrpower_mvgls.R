#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
if (!length(file_arg)) stop("This script must be run with Rscript.", call. = FALSE)
repo_root <- normalizePath(file.path(dirname(file_arg), ".."))

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required to run this harness.", call. = FALSE)
}

pkgload::load_all(repo_root, quiet = TRUE)

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

fit_corrpower <- function(formula, data, tree, bmm.reference = NULL) {
  mvgls(
    formula,
    data = data,
    tree = tree,
    model = "BMM",
    method = "LL",
    REML = FALSE,
    bmm.structure = "corrpower",
    bmm.reference = bmm.reference,
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

make_response <- function(tree, sigma_by_regime, beta = NULL, x = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  y <- mvSIM(
    tree,
    nsim = 1,
    model = "BMM",
    param = list(
      ntraits = nrow(sigma_by_regime[[1]]),
      sigma = sigma_by_regime,
      theta = rep(0, nrow(sigma_by_regime[[1]]))
    )
  )
  if (is.list(y)) y <- y[[1]]
  y <- as.matrix(y)
  if (!is.null(beta)) {
    X <- cbind("(Intercept)" = 1, x = as.numeric(x))
    y <- y + X %*% beta
  }
  y
}

build_omega <- function(fit) {
  kron_sum <- get("kroneckerSum", envir = asNamespace("mvMORPH"))
  regime_cov <- vcvSplit(fit$variables$tree)
  .Call(
    kron_sum,
    R = fit$sigma$regime,
    C = regime_cov,
    Rrows = as.integer(fit$dims$p),
    Crows = as.integer(fit$dims$n),
    dimlist = as.integer(length(regime_cov))
  )
}

whiten_residuals <- function(fit) {
  resid <- residuals(fit, type = "response")
  omega <- build_omega(fit)
  chol_omega <- chol(omega)
  whitened <- backsolve(chol_omega, as.numeric(resid), transpose = TRUE)
  matrix(whitened, nrow = nrow(resid), ncol = ncol(resid), dimnames = dimnames(resid))
}

set.seed(20260319)
n_tips <- 20L
tree <- phytools::pbtree(n = n_tips, scale = 1)
states <- setNames(sample(rep(c("A", "B"), each = n_tips / 2L), n_tips), tree$tip.label)
simmap <- suppressMessages(phytools::make.simmap(tree, states, model = "ER", nsim = 1))

set.seed(20260319)
weak_n_tips <- 30L
weak_tree <- phytools::pbtree(n = weak_n_tips, scale = 1)
weak_states <- setNames(sample(rep(c("A", "B"), each = weak_n_tips / 2L), weak_n_tips), weak_tree$tip.label)
weak_simmap <- suppressMessages(suppressWarnings(phytools::make.simmap(weak_tree, weak_states, model = "ER", nsim = 1, message = FALSE)))

sigma_base <- matrix(c(1.30, 0.55, 0.55, 1.00), 2, 2)
sigma_prop <- 1.8 * sigma_base
sigma_weak <- matrix(c(1.70, 0.00, 0.00, 1.20), 2, 2)
sigma_strong <- matrix(c(1.70, 0.85, 0.85, 1.20), 2, 2)

Y_prop <- make_response(simmap, list(A = sigma_base, B = sigma_prop), seed = 20260320)
Y_weak <- make_response(weak_simmap, list(A = sigma_base, B = sigma_weak), seed = 20260328)
Y_strong <- make_response(simmap, list(A = sigma_base, B = sigma_strong), seed = 20260322)
X_prop <- as.numeric(scale(rnorm(nrow(Y_prop))))
names(X_prop) <- rownames(Y_prop)

fit_prop <- fit_proportional(Y ~ 1, data = list(Y = Y_prop), tree = simmap)
fit_corr_prop <- fit_corrpower(Y ~ 1, data = list(Y = Y_prop), tree = simmap)
fit_corr_weak <- fit_corrpower(Y ~ 1, data = list(Y = Y_weak), tree = weak_simmap)
fit_corr_strong <- fit_corrpower(Y ~ 1, data = list(Y = Y_strong), tree = simmap)

assert_true(identical(fit_corr_prop$bmm.structure, "corrpower"), "corrpower fit did not mark the structure")
assert_true(all(c("df.free", "df.free_cov", "df.free_beta", "df.free_model") %in% names(fit_corr_prop)),
  "corrpower fit is missing one or more df.free fields"
)
assert_true(is.finite(fit_corr_prop$logLik), "corrpower fit produced a non-finite logLik")

corr_power_names <- grep("\\.corr_power$", names(fit_corr_prop$param), value = TRUE)
assert_true(length(corr_power_names) >= 2, "expected at least one non-reference corr_power parameter")
assert_true(all(is.finite(fit_corr_prop$param[corr_power_names])),
  "proportional-case corr_power estimates should remain finite"
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

assert_true("corr_power" %in% colnames(fit_corr_prop$regime.summary),
  "regime.summary does not expose the corr_power column"
)
assert_true(identical(rownames(fit_corr_prop$regime.summary), names(fit_corr_prop$sigma$regime)),
  "regime.summary row names do not match the fitted regimes"
)

corr_power_weak <- fit_corr_weak$param[grep("\\.corr_power$", names(fit_corr_weak$param))]
assert_true(corr_power_weak[[2]] < corr_power_weak[[1]],
  "weaker-than-reference case did not estimate a smaller corr_power for the derived regime"
)
weak_regime_mats <- fit_corr_weak$sigma$regime
assert_true(offdiag_mean(cov2cor(weak_regime_mats[[names(weak_regime_mats)[2]]])) < offdiag_mean(cov2cor(weak_regime_mats[[1]])),
  "weaker-than-reference case did not reduce the off-diagonal correlation power"
)

corr_power_strong <- fit_corr_strong$param[grep("\\.corr_power$", names(fit_corr_strong$param))]
assert_true(corr_power_strong[[2]] > corr_power_strong[[1]],
  "stronger-than-reference case did not estimate a larger corr_power for the derived regime"
)
strong_regime_mats <- fit_corr_strong$sigma$regime
assert_true(offdiag_mean(cov2cor(strong_regime_mats[[names(strong_regime_mats)[2]]])) > offdiag_mean(cov2cor(strong_regime_mats[[1]])),
  "stronger-than-reference case did not increase the off-diagonal correlation power"
)

fit_corr_reg <- fit_corrpower(Y ~ x, data = list(Y = Y_prop, x = X_prop), tree = simmap)
assert_true(identical(fit_corr_reg$bmm.structure, "corrpower"), "regression fit did not preserve corr-power mode")
assert_true(nrow(coef(fit_corr_reg)) == 2L, "unexpected coefficient row count for regression fit")
assert_true(ncol(coef(fit_corr_reg)) == ncol(Y_prop), "unexpected coefficient column count for regression fit")
assert_true(identical(dim(fitted(fit_corr_reg)), dim(Y_prop)), "fitted values have the wrong dimensions")
assert_true(identical(dim(residuals(fit_corr_reg)), dim(Y_prop)), "response residuals have the wrong dimensions")

summary_reg <- summary(fit_corr_reg)
assert_true(inherits(summary_reg, "summary.mvgls"), "summary() did not return a summary.mvgls object")
assert_true(all(c("AIC", "GIC", "logLik") %in% colnames(summary_reg$results.fit)),
  "summary output is missing the expected fit statistics"
)
assert_true(is.data.frame(fit_corr_reg$regime.summary), "corr-power fit did not store regime.summary")
assert_true(is.data.frame(summary_reg$regime.summary), "summary() did not preserve regime.summary")
assert_true(all(c("reference", "scale", "corr_power", "mean_rate", "mean_variance", "mean_covariance",
                  "mean_correlation", "mean_abs_correlation") %in% colnames(summary_reg$regime.summary)),
  "regime.summary is missing one or more expected columns"
)
assert_true(all(is.finite(summary_reg$regime.summary$corr_power[summary_reg$regime.summary$reference])),
  "reference regime should retain a finite corr_power estimate"
)

gic_reg <- GIC(fit_corr_reg)
aic_reg <- AIC(fit_corr_reg)
assert_true(is.finite(gic_reg$GIC), "GIC returned a non-finite value")
assert_true(is.finite(aic_reg$AIC), "AIC returned a non-finite value")
assert_true(abs(gic_reg$GIC - aic_reg$AIC) < 1e-8, "GIC and AIC should coincide for corr-power ML fits")
assert_true(gic_reg$bias_cov == fit_corr_reg$df.free_cov, "GIC bias_cov does not match df.free_cov")
assert_true(gic_reg$bias == fit_corr_reg$df.free, "GIC bias does not match the total free parameter count")
assert_true(
  fit_corr_reg$df.free == fit_corr_reg$df.free_cov + fit_corr_reg$df.free_beta + fit_corr_reg$df.free_model,
  "df.free fields do not add up correctly"
)

diag_api <- corrpower_diagnostics(fit_corr_reg, nboot = 4L, nbcores = 1L, profile_points = 5L)
assert_true(inherits(diag_api, "corrpower_diagnostics"), "corrpower_diagnostics() did not return the expected class")
assert_true(all(c("parameter_summary", "regime_summary", "anchor_summary", "anchor_pairwise_cov", "acceptance") %in% names(diag_api)),
  "corrpower_diagnostics() is missing one or more expected top-level components"
)
assert_true(is.data.frame(diag_api$parameter_summary), "parameter_summary should be a data.frame")
assert_true(is.data.frame(diag_api$regime_summary), "regime_summary diagnostics should be a data.frame")
assert_true(is.data.frame(diag_api$anchor_summary), "anchor_summary should be a data.frame")
assert_true(any(diag_api$parameter_summary$label == "B.scale"), "parameter_summary is missing the non-reference scale parameter")
assert_true(any(diag_api$parameter_summary$label == "B.corr_power"), "parameter_summary is missing the non-reference corr_power parameter")
assert_true(any(diag_api$regime_summary$label == "A.mean_rate"), "regime_summary diagnostics are missing the derived mean-rate metric")
assert_true(all(c("selected", "boundary_corr_power", "pathological_scale") %in% names(diag_api$acceptance)),
  "acceptance checks are missing one or more expected fields"
)

diag_print <- capture.output(print(diag_api))
assert_true(any(grepl("Corr-power diagnostics", diag_print, fixed = TRUE)), "print.corrpower_diagnostics() did not print the expected header")
assert_true(any(grepl("parameter_summary", diag_print, fixed = TRUE)), "print.corrpower_diagnostics() did not print the parameter summary")

printed <- capture.output(print(fit_corr_reg))
assert_true(any(grepl("Corr-power", printed, fixed = TRUE)), "print.mvgls() did not mention corr-power mode")
summary_printed <- capture.output(print(summary(fit_corr_reg)))
assert_true(any(grepl("Corr-power", summary_printed, fixed = TRUE)), "summary.mvgls() did not mention corr-power mode")

norm_resid <- whiten_residuals(fit_corr_reg)
assert_true(identical(dim(norm_resid), dim(residuals(fit_corr_reg))), "normalized residuals have the wrong dimensions")
assert_true(all(is.finite(norm_resid)), "normalized residuals contain non-finite values")

pred_no_tree <- predict(fit_corr_reg)
assert_true(identical(dim(pred_no_tree), dim(Y_prop)), "predict() without newdata returned the wrong dimensions")

cat("corr-power mvgls harness checks passed\n")
