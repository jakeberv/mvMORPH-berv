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

simulate_response <- function(tree, sigma_by_regime, seed, x = NULL, beta = NULL) {
  set.seed(seed)
  y <- mvSIM(
    tree,
    nsim = 1,
    model = "BMM",
    param = list(ntraits = nrow(sigma_by_regime[[1]]), sigma = sigma_by_regime, theta = rep(0, nrow(sigma_by_regime[[1]])))
  )
  if (is.list(y)) y <- y[[1]]
  y <- as.matrix(y)
  if (!is.null(x) && !is.null(beta)) {
    X <- cbind("(Intercept)" = 1, x = as.numeric(x))
    beta_mat <- if (is.matrix(beta)) beta else matrix(rep(beta, ncol(y)), nrow = ncol(X), ncol = ncol(y))
    y <- y + X %*% beta_mat
  }
  y
}

fit_coronly <- function(formula, data, tree, bmm.reference = NULL) {
  mvgls(
    formula,
    data = data,
    tree = tree,
    model = "BMM",
    method = "LL",
    REML = FALSE,
    bmm.structure = "corrpower",
    bmm.scale = FALSE,
    bmm.reference = bmm.reference,
    echo = FALSE
  )
}

scenario_fit <- function(label, tree, Y, formula, data, expected_direction) {
  fit <- fit_coronly(formula, data = data, tree = tree)
  corr_power_names <- grep("\\.corr_power$", names(fit$param), value = TRUE)
  corr_power_vals <- unname(fit$param[corr_power_names])
  assert_true(length(corr_power_vals) >= 2L, sprintf("%s did not produce a derived corr_power estimate", label))
  ref_value <- corr_power_vals[[1]]
  derived_value <- corr_power_vals[[2]]
  if (expected_direction == "lower") {
    assert_true(derived_value < ref_value, sprintf("%s did not estimate a lower corr_power for the derived regime", label))
  } else if (expected_direction == "higher") {
    assert_true(derived_value > ref_value, sprintf("%s did not estimate a higher corr_power for the derived regime", label))
  } else {
    stop("unexpected expected_direction", call. = FALSE)
  }
  assert_true(!"scale" %in% colnames(fit$regime.summary), sprintf("%s should not report a scale column", label))
  assert_true(all(c("reference", "corr_power", "mean_rate", "mean_variance", "mean_covariance",
                    "mean_correlation", "mean_abs_correlation") %in% colnames(fit$regime.summary)),
              sprintf("%s regime.summary is missing expected columns", label))
  assert_true(is.finite(as.numeric(fit$logLik)), sprintf("%s produced a non-finite logLik", label))
  diag <- corrpower_diagnostics(fit, nboot = 4L, nbcores = 1L, profile_points = 5L)
  assert_true(any(diag$parameter_summary$label == "B.corr_power"), sprintf("%s diagnostics missing B.corr_power", label))
  assert_true(!any(diag$parameter_summary$label == "B.scale"), sprintf("%s diagnostics should not include scale", label))
  conf <- confint(fit, method = "profile", profile_points = 5L)
  assert_true(is.matrix(conf) || is.data.frame(conf), sprintf("%s confint() returned an unexpected type", label))
  invisible(fit)
}

balanced_tree <- make_simmap(seed = 20260334, n_tips = 24L, states_per_regime = c(A = 12L, B = 12L))
balanced_base <- matrix(c(1.3, 0.60, 0.60, 1.0), 2, 2)

weaker_Y <- simulate_response(
  balanced_tree,
  list(A = balanced_base, B = apply_corrpower(balanced_base, 0.45)),
  seed = 20260335
)
scenario_fit(
  "balanced_weaker",
  tree = balanced_tree,
  Y = weaker_Y,
  formula = Y ~ 1,
  data = list(Y = weaker_Y),
  expected_direction = "lower"
)

stronger_Y <- simulate_response(
  balanced_tree,
  list(A = balanced_base, B = apply_corrpower(balanced_base, 1.45)),
  seed = 20260336
)
scenario_fit(
  "balanced_stronger",
  tree = balanced_tree,
  Y = stronger_Y,
  formula = Y ~ 1,
  data = list(Y = stronger_Y),
  expected_direction = "higher"
)

covariate_tree <- make_simmap(seed = 20260337, n_tips = 30L, states_per_regime = c(A = 15L, B = 15L))
covariate_base <- matrix(c(
  1.2, 0.20, 0.10, 0.05,
  0.20, 1.1, 0.15, 0.10,
  0.10, 0.15, 1.0, 0.20,
  0.05, 0.10, 0.20, 0.9
), 4, 4)
x <- as.numeric(scale(seq_len(30L)))
names(x) <- covariate_tree$tip.label
covariate_Y <- simulate_response(
  covariate_tree,
  list(A = covariate_base, B = apply_corrpower(covariate_base, 1.20)),
  seed = 20260338,
  x = x,
  beta = c(0.1, -0.05)
)
covariate_fit <- fit_coronly(Y ~ x, data = list(Y = covariate_Y, x = x), tree = covariate_tree)
pred <- predict(covariate_fit)
assert_true(identical(dim(pred), dim(covariate_Y)), "predict() returned the wrong dimensions for the covariate scenario")
node_labels <- paste0("node_", Ntip(covariate_tree) + seq_len(Nnode(covariate_tree)))
node_newdata <- data.frame(
  x = seq(-0.5, 0.5, length.out = length(node_labels)),
  row.names = node_labels
)
anc <- ancestral(covariate_fit, newdata = node_newdata)
assert_true(is.matrix(anc), "ancestral() should return a matrix for the covariate scenario")
assert_true(ncol(anc) == ncol(covariate_Y), "ancestral() returned the wrong number of trait columns for the covariate scenario")

cat("corrpower (bmm.scale=FALSE) identifiability harness checks passed\n")
