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
    stop(sprintf("Error message did not match '%s': %s", pattern, conditionMessage(err)), call. = FALSE)
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

make_simmap <- function(seed, n_tips, states_per_regime) {
  set.seed(seed)
  tree <- phytools::pbtree(n = n_tips, scale = 1)
  sampled_states <- sample(rep(names(states_per_regime), times = states_per_regime), n_tips)
  states <- setNames(sampled_states, tree$tip.label)
  suppressMessages(phytools::make.simmap(tree, states, model = "ER", nsim = 1))
}

apply_corrpower <- function(base_sigma, corr_power) {
  eig <- eigen(stats::cov2cor(base_sigma), symmetric = TRUE)
  vals <- pmax(eig$values, .Machine$double.eps)
  corr_raw <- eig$vectors %*% diag(vals ^ corr_power, nrow = length(vals)) %*% t(eig$vectors)
  corr_mat <- stats::cov2cor(0.5 * (corr_raw + t(corr_raw)))
  D <- diag(sqrt(diag(base_sigma)))
  D %*% corr_mat %*% D
}

simulate_response <- function(tree, sigma_by_regime, seed) {
  set.seed(seed)
  y <- mvSIM(
    tree,
    nsim = 1,
    model = "BMM",
    param = list(ntraits = nrow(sigma_by_regime[[1]]), sigma = sigma_by_regime, theta = rep(0, nrow(sigma_by_regime[[1]])))
  )
  if (is.list(y)) y <- y[[1]]
  as.matrix(y)
}

fit_coronly <- function(formula, data, tree, bmm.reference = NULL) {
  mvgls(
    formula,
    data = data,
    tree = tree,
    model = "BMM",
    method = "LL",
    REML = FALSE,
    bmm.structure = "corrpower_coronly",
    bmm.reference = bmm.reference,
    echo = FALSE
  )
}

set.seed(20260326)
tree <- make_simmap(seed = 20260326, n_tips = 20L, states_per_regime = c(A = 10L, B = 10L))
base_sigma <- matrix(c(1.30, 0.55, 0.55, 1.00), 2, 2)
derived_sigma <- apply_corrpower(base_sigma, 0.55)
Y <- simulate_response(tree, list(A = base_sigma, B = derived_sigma), seed = 20260327)

fit <- fit_coronly(Y ~ 1, data = list(Y = Y), tree = tree)

assert_true(identical(fit$bmm.structure, "corrpower_coronly"), "fit did not record corrpower_coronly as the structure")
assert_true(all(!grepl("\\.scale$", names(fit$param))), "corrpower_coronly should not expose raw scale parameters")
assert_true(all(grepl("\\.corr_power$", names(fit$param))), "corrpower_coronly should expose corr_power parameters")
assert_true(is.data.frame(fit$regime.summary), "corrpower_coronly fit did not store a regime.summary table")
assert_true(!"scale" %in% colnames(fit$regime.summary), "corrpower_coronly regime.summary should not include a scale column")
assert_true(all(c("reference", "corr_power", "mean_rate", "mean_variance", "mean_covariance",
                  "mean_correlation", "mean_abs_correlation") %in% colnames(fit$regime.summary)),
            "corrpower_coronly regime.summary is missing expected columns")

for (nm in names(fit$sigma$regime)) {
  mat <- fit$sigma$regime[[nm]]
  assert_true(is_symmetric(mat), sprintf("regime covariance '%s' is not symmetric", nm))
  assert_true(is_spd(mat), sprintf("regime covariance '%s' is not positive definite", nm))
}

gic <- GIC(fit)
aic <- AIC(fit)
eic <- EIC(fit, nboot = 3L, nbcores = 1L)
assert_true(is.finite(gic$GIC), "GIC returned a non-finite value")
assert_true(is.finite(aic$AIC), "AIC returned a non-finite value")
assert_true(is.finite(eic$EIC), "EIC returned a non-finite value")
assert_true(abs(gic$GIC - aic$AIC) < 1e-8, "GIC and AIC should coincide for corrpower_coronly ML fits")

sim <- simulate(fit, nsim = 2, seed = 1)
assert_true(is.list(sim) && length(sim) == 2L, "simulate() should return a list of length 2 for nsim = 2")
assert_true(identical(dim(sim[[1]]), dim(Y)), "simulate() returned the wrong dimensions")
assert_true(identical(dim(sim[[2]]), dim(Y)), "simulate() returned the wrong dimensions for the second draw")

pred <- predict(fit)
assert_true(identical(dim(pred), dim(Y)), "predict() returned the wrong dimensions")

anc <- ancestral(fit)
assert_true(is.matrix(anc), "ancestral() should return a matrix")
assert_true(ncol(anc) == ncol(Y), "ancestral() returned the wrong number of trait columns")
assert_true(nrow(anc) == Nnode(tree), "ancestral() returned the wrong number of node rows")

norm_resid <- residuals(fit, type = "normalized")
assert_true(identical(dim(norm_resid), dim(Y)), "normalized residuals have the wrong dimensions")
assert_true(all(is.finite(norm_resid)), "normalized residuals contain non-finite values")

diag_api <- corrpower_coronly_diagnostics(fit, nboot = 4L, nbcores = 1L, profile_points = 5L)
assert_true(inherits(diag_api, "corrpower_coronly_diagnostics"), "corrpower_coronly_diagnostics() returned the wrong class")
assert_true(all(c("parameter_summary", "regime_summary", "anchor_summary", "anchor_pairwise_cov", "acceptance") %in% names(diag_api)),
            "corrpower_coronly_diagnostics() is missing expected top-level components")
assert_true(any(diag_api$parameter_summary$label == "B.corr_power"), "diagnostics missing the derived corr_power parameter")
assert_true(!any(diag_api$parameter_summary$label == "B.scale"), "corrpower_coronly diagnostics should not report a scale parameter")

conf <- confint(fit, method = "both", nboot = 4L, profile_points = 5L)
assert_true(is.data.frame(conf), "confint() should return a data.frame when method = 'both'")
assert_true(all(grepl("corr_power", rownames(conf), fixed = TRUE)), "confint() should only report corr_power parameters")
assert_true(all(c("estimate", "profile_low", "profile_high", "bootstrap_low", "bootstrap_high") %in% colnames(conf)),
            "confint() is missing expected summary columns")

cat("corrpower_coronly mvgls harness checks passed\n")
