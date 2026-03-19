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

fit_corrstrength <- function(formula, data, tree, bmm.reference = NULL) {
  mvgls(
    formula,
    data = data,
    tree = tree,
    model = "BMM",
    method = "LL",
    REML = FALSE,
    bmm.structure = "corrstrength",
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
    param = list(ntraits = nrow(sigma_by_regime[[1]]), sigma = sigma_by_regime, theta = rep(0, nrow(sigma_by_regime[[1]])))
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
fit_corr_prop <- fit_corrstrength(Y ~ 1, data = list(Y = Y_prop), tree = simmap)
fit_corr_weak <- fit_corrstrength(Y ~ 1, data = list(Y = Y_weak), tree = weak_simmap)
fit_corr_strong <- fit_corrstrength(Y ~ 1, data = list(Y = Y_strong), tree = simmap)

assert_true(identical(fit_corr_prop$bmm.structure, "corrstrength"), "corr-strength fit did not mark the structure")
assert_true(all(c("df.free", "df.free_cov", "df.free_beta", "df.free_model") %in% names(fit_corr_prop)),
  "corr-strength fit is missing one or more df.free fields"
)
assert_true(is.finite(fit_corr_prop$logLik), "corr-strength fit produced a non-finite logLik")
assert_true(abs(as.numeric(fit_corr_prop$logLik) - as.numeric(fit_prop$logLik)) < 5,
  "corr-strength logLik is not comparable to the proportional fit"
)

corr_strength_names <- grep("\\.corr_strength$", names(fit_corr_prop$param), value = TRUE)
assert_true(length(corr_strength_names) >= 2, "expected at least one non-reference corr_strength parameter")
assert_true(all(is.finite(fit_corr_prop$param[corr_strength_names])),
  "proportional-case corr_strength estimates should remain finite"
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

assert_true("corr_strength" %in% colnames(fit_corr_prop$regime.summary),
  "regime.summary does not expose the corr_strength column"
)
assert_true(identical(rownames(fit_corr_prop$regime.summary), names(fit_corr_prop$sigma$regime)),
  "regime.summary row names do not match the fitted regimes"
)

corr_strength_weak <- fit_corr_weak$param[grep("\\.corr_strength$", names(fit_corr_weak$param))]
assert_true(corr_strength_weak[[2]] < corr_strength_weak[[1]],
  "weaker-than-reference case did not estimate a smaller corr_strength for the derived regime"
)
weak_regime_mats <- fit_corr_weak$sigma$regime
assert_true(offdiag_mean(cov2cor(weak_regime_mats[[names(weak_regime_mats)[2]]])) < offdiag_mean(cov2cor(weak_regime_mats[[1]])),
  "weaker-than-reference case did not reduce the off-diagonal correlation strength"
)

corr_strength_strong <- fit_corr_strong$param[grep("\\.corr_strength$", names(fit_corr_strong$param))]
assert_true(corr_strength_strong[[2]] > corr_strength_strong[[1]],
  "stronger-than-reference case did not estimate a larger corr_strength for the derived regime"
)
strong_regime_mats <- fit_corr_strong$sigma$regime
assert_true(offdiag_mean(cov2cor(strong_regime_mats[[names(strong_regime_mats)[2]]])) > offdiag_mean(cov2cor(strong_regime_mats[[1]])),
  "stronger-than-reference case did not increase the off-diagonal correlation strength"
)

fit_corr_reg <- fit_corrstrength(Y ~ x, data = list(Y = Y_prop, x = X_prop), tree = simmap)
assert_true(identical(fit_corr_reg$bmm.structure, "corrstrength"), "regression fit did not preserve corr-strength mode")
assert_true(nrow(coef(fit_corr_reg)) == 2L, "unexpected coefficient row count for regression fit")
assert_true(ncol(coef(fit_corr_reg)) == ncol(Y_prop), "unexpected coefficient column count for regression fit")
assert_true(identical(dim(fitted(fit_corr_reg)), dim(Y_prop)), "fitted values have the wrong dimensions")
assert_true(identical(dim(residuals(fit_corr_reg)), dim(Y_prop)), "response residuals have the wrong dimensions")

summary_reg <- summary(fit_corr_reg)
assert_true(inherits(summary_reg, "summary.mvgls"), "summary() did not return a summary.mvgls object")
assert_true(all(c("AIC", "GIC", "logLik") %in% colnames(summary_reg$results.fit)),
  "summary output is missing the expected fit statistics"
)
assert_true(is.data.frame(fit_corr_reg$regime.summary), "corr-strength fit did not store regime.summary")
assert_true(is.data.frame(summary_reg$regime.summary), "summary() did not preserve regime.summary")
assert_true(all(c("reference", "scale", "corr_strength", "mean_rate", "mean_variance", "mean_covariance",
                  "mean_correlation", "mean_abs_correlation") %in% colnames(summary_reg$regime.summary)),
  "regime.summary is missing one or more expected columns"
)
assert_true(all(is.finite(summary_reg$regime.summary$corr_strength[summary_reg$regime.summary$reference])),
  "reference regime should retain a finite corr_strength estimate"
)

gic_reg <- GIC(fit_corr_reg)
aic_reg <- AIC(fit_corr_reg)
assert_true(is.finite(gic_reg$GIC), "GIC returned a non-finite value")
assert_true(is.finite(aic_reg$AIC), "AIC returned a non-finite value")
assert_true(abs(gic_reg$GIC - aic_reg$AIC) < 1e-8, "GIC and AIC should coincide for corr-strength ML fits")
assert_true(gic_reg$bias_cov == fit_corr_reg$df.free_cov, "GIC bias_cov does not match df.free_cov")
assert_true(gic_reg$bias == fit_corr_reg$df.free, "GIC bias does not match the total free parameter count")
assert_true(
  fit_corr_reg$df.free == fit_corr_reg$df.free_cov + fit_corr_reg$df.free_beta + fit_corr_reg$df.free_model,
  "df.free fields do not add up correctly"
)

diag_api <- corrstrength_diagnostics(fit_corr_reg, nboot = 4L, nbcores = 1L, profile_points = 5L)
assert_true(inherits(diag_api, "corrstrength_diagnostics"), "corrstrength_diagnostics() did not return the expected class")
assert_true(all(c("parameter_summary", "regime_summary", "anchor_summary", "anchor_pairwise_cov", "acceptance") %in% names(diag_api)),
  "corrstrength_diagnostics() is missing one or more expected top-level components"
)
assert_true(is.data.frame(diag_api$parameter_summary), "parameter_summary should be a data.frame")
assert_true(is.data.frame(diag_api$regime_summary), "regime_summary diagnostics should be a data.frame")
assert_true(is.data.frame(diag_api$anchor_summary), "anchor_summary should be a data.frame")
assert_true(any(diag_api$parameter_summary$label == "B.scale"), "parameter_summary is missing the non-reference scale parameter")
assert_true(any(diag_api$parameter_summary$label == "B.corr_strength"), "parameter_summary is missing the non-reference corr_strength parameter")
assert_true(any(diag_api$regime_summary$label == "A.mean_rate"), "regime_summary diagnostics are missing the derived mean-rate metric")
assert_true(all(c("logLik_spread", "logLik_spread_ok", "no_pathological_scale") %in% names(diag_api$acceptance)),
  "corrstrength_diagnostics() did not return the expected acceptance checks"
)
diag_print <- capture.output(print(diag_api))
assert_true(any(grepl("Corr-strength diagnostics", diag_print, fixed = TRUE)), "print.corrstrength_diagnostics() did not print the expected header")

ci_profile <- confint(fit_corr_reg, method = "profile", profile_points = 5L)
ci_boot <- confint(fit_corr_reg, method = "bootstrap", nboot = 4L, nbcores = 1L)
ci_both <- confint(fit_corr_reg, method = "both", nboot = 4L, nbcores = 1L, profile_points = 5L)
assert_true(all(c("B.scale", "B.corr_strength") %in% rownames(ci_profile)), "profile confint() is missing the non-reference parameters")
assert_true(all(c("lower", "upper") %in% colnames(ci_profile)), "profile confint() returned unexpected columns")
assert_true(all(c("B.scale", "B.corr_strength") %in% rownames(ci_boot)), "bootstrap confint() is missing the non-reference parameters")
assert_true(all(c("lower", "upper") %in% colnames(ci_boot)), "bootstrap confint() returned unexpected columns")
assert_true(all(c("estimate", "profile_low", "profile_high", "bootstrap_low", "bootstrap_high") %in% colnames(ci_both)),
  "confint(..., method = \"both\") returned unexpected columns"
)

predict_no_tree <- predict(fit_corr_reg, newdata = data.frame(x = X_prop, row.names = rownames(Y_prop)))
assert_true(is.matrix(predict_no_tree), "predict() without a tree should return a matrix")
assert_true(identical(dim(predict_no_tree), dim(fitted(fit_corr_reg))),
  "predict() without a tree returned the wrong dimensions"
)
assert_true(max(abs(predict_no_tree - fitted(fit_corr_reg))) < 1e-8,
  "predict() without a tree should agree with fitted values"
)

full_tree <- phytools::pbtree(n = 24L, scale = 1)
full_states <- setNames(sample(rep(c("A", "B"), each = 12L), 24L), full_tree$tip.label)
full_simmap <- suppressMessages(phytools::make.simmap(full_tree, full_states, model = "ER", nsim = 1))
x_full <- as.numeric(scale(rnorm(Ntip(full_simmap))))
names(x_full) <- full_simmap$tip.label
beta_full <- rbind(
  "(Intercept)" = c(0.4, -0.2),
  x = c(0.6, -0.3)
)
Y_full <- make_response(full_simmap, list(A = sigma_base, B = sigma_strong), beta = beta_full, x = x_full, seed = 20260323)
train_tips <- full_simmap$tip.label[1:18]
hold_tips <- setdiff(full_simmap$tip.label, train_tips)
train_tree <- drop.tip(full_simmap, hold_tips)
fit_predict <- fit_corrstrength(
  Y ~ x,
  data = list(Y = Y_full[train_tips, , drop = FALSE], x = x_full[train_tips]),
  tree = train_tree
)
predict_with_tree <- predict(
  fit_predict,
  newdata = data.frame(x = x_full[hold_tips], row.names = hold_tips),
  tree = full_simmap
)
assert_true(is.matrix(predict_with_tree), "predict() with a tree should return a matrix")
assert_true(identical(rownames(predict_with_tree), hold_tips), "predict() with a tree did not preserve row names")
assert_true(ncol(predict_with_tree) == ncol(Y_full), "predict() with a tree returned the wrong number of traits")
assert_true(all(is.finite(predict_with_tree)), "predict() with a tree returned non-finite values")

anc_int <- ancestral(fit_corr_weak)
assert_true(is.matrix(anc_int), "ancestral() for an intercept-only fit should return a matrix")
assert_true(ncol(anc_int) == ncol(Y_weak), "ancestral() for an intercept-only fit returned the wrong number of traits")
assert_true(nrow(anc_int) == Nnode(weak_simmap), "ancestral() for an intercept-only fit returned the wrong number of nodes")
assert_true(all(is.finite(anc_int)), "ancestral() for an intercept-only fit returned non-finite values")

node_labels <- paste0("node_", Ntip(train_tree) + seq_len(Nnode(train_tree)))
node_newdata <- data.frame(
  x = seq(-0.5, 0.5, length.out = length(node_labels)),
  row.names = node_labels
)
anc_reg <- ancestral(fit_predict, newdata = node_newdata)
assert_true(is.matrix(anc_reg), "ancestral() for a regression fit should return a matrix")
assert_true(ncol(anc_reg) == ncol(Y_full), "ancestral() for a regression fit returned the wrong number of traits")
assert_true(all(is.finite(anc_reg)), "ancestral() for a regression fit returned non-finite values")

norm_resid <- residuals(fit_predict, type = "normalized")
direct_norm <- whiten_residuals(fit_predict)
assert_true(identical(dim(norm_resid), dim(direct_norm)), "normalized residuals have the wrong dimensions")
assert_true(max(abs(norm_resid - direct_norm)) < 1e-6,
  "normalized residuals do not match direct whitening by the reconstructed Omega"
)

expect_error(vcov(fit_corr_reg, type = "coef"), "corr-strength")

sim_one <- simulate(fit_corr_reg, nsim = 1)
assert_true(is.matrix(sim_one), "simulate(..., nsim=1) should return a matrix")
assert_true(identical(dim(sim_one), dim(Y_prop)), "simulate(..., nsim=1) returned the wrong dimensions")
assert_true(identical(rownames(sim_one), rownames(Y_prop)), "simulate(..., nsim=1) did not preserve row names")
assert_true(all(is.finite(sim_one)), "simulate(..., nsim=1) returned non-finite values")

sim_three <- simulate(fit_corr_reg, nsim = 3)
assert_true(is.list(sim_three), "simulate(..., nsim=3) should return a list")
assert_true(length(sim_three) == 3L, "simulate(..., nsim=3) returned the wrong number of replicates")
for (i in seq_along(sim_three)) {
  sim_i <- sim_three[[i]]
  assert_true(is.matrix(sim_i), sprintf("simulate replicate %d is not a matrix", i))
  assert_true(identical(dim(sim_i), dim(Y_prop)), sprintf("simulate replicate %d has the wrong dimensions", i))
  assert_true(identical(rownames(sim_i), rownames(Y_prop)), sprintf("simulate replicate %d did not preserve row names", i))
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

cat("corr-strength mvgls docs/test harness checks passed\n")
