#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
if (!length(file_arg)) stop("This script must be run with Rscript.", call. = FALSE)
repo_root <- normalizePath(file.path(dirname(file_arg), ".."))

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required to run this study.", call. = FALSE)
}

devtools::load_all(repo_root, quiet = TRUE)

`%||%` <- function(x, y) if (is.null(x)) y else x

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed <= 0L) default else parsed
}

reps <- env_int("CORRSHRINK_IDENT_REPS", 6L)
seed_base <- env_int("CORRSHRINK_IDENT_SEED", 20260319L)

make_sigma <- function(sds, corr) {
  diag(sds) %*% corr %*% diag(sds)
}

apply_corrshrink <- function(base_sigma, scale, rho) {
  V <- diag(diag(base_sigma))
  W <- base_sigma - V
  scale * (V + rho * W)
}

make_simmap <- function(n_tips, state_counts, seed) {
  set.seed(seed)
  tree <- phytools::pbtree(n = n_tips, scale = 1)
  labels <- tree$tip.label
  sampled_states <- sample(rep(names(state_counts), times = state_counts), length(labels))
  states <- setNames(sampled_states, labels)
  invisible(capture.output({
    simmap <- suppressMessages(phytools::make.simmap(tree, states, model = "ER", nsim = 1))
  }))
  simmap
}

simulate_response <- function(tree, sigma_by_regime, beta = NULL, x = NULL) {
  p <- nrow(sigma_by_regime[[1]])
  Y <- mvSIM(
    tree,
    nsim = 1,
    model = "BMM",
    param = list(ntraits = p, sigma = sigma_by_regime, theta = rep(0, p))
  )
  if (is.list(Y)) Y <- Y[[1]]
  Y <- as.matrix(Y)
  if (!is.null(beta)) {
    X <- cbind("(Intercept)" = 1, x = as.numeric(x))
    Y <- Y + X %*% beta
  }
  Y
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

try_fit <- function(expr) {
  tryCatch(expr, error = function(e) e)
}

run_replicate <- function(scenario, rep_id) {
  tree <- make_simmap(scenario$n_tips, scenario$state_counts, scenario$seed_offset + rep_id)
  derived_regime <- scenario$derived_regime %||% names(scenario$state_counts)[2]
  sigma_by_regime <- list(
    A = scenario$base_sigma,
    B = apply_corrshrink(scenario$base_sigma, scenario$true_scale, scenario$true_rho)
  )

  x <- NULL
  beta <- NULL
  data <- list()
  formula <- if (isTRUE(scenario$with_covariate)) Y ~ x else Y ~ 1
  if (isTRUE(scenario$with_covariate)) {
    set.seed(scenario$seed_offset + 1000L + rep_id)
    x <- as.numeric(scale(rnorm(Ntip(tree))))
    names(x) <- tree$tip.label
    beta <- scenario$beta
    data$x <- x
  }

  Y <- simulate_response(tree, sigma_by_regime, beta = beta, x = x)
  data$Y <- Y
  data <- data[c("Y", setdiff(names(data), "Y"))]

  fit_cs <- try_fit(fit_corrshrink(formula, data = data, tree = tree))
  fit_prop <- try_fit(fit_proportional(formula, data = data, tree = tree))

  success_cs <- !inherits(fit_cs, "error")
  success_prop <- !inherits(fit_prop, "error")

  est_scale <- NA_real_
  est_rho <- NA_real_
  aic_delta <- NA_real_
  loglik_gain <- NA_real_
  rho_boundary <- NA_real_

  if (success_cs) {
    est_scale <- unname(fit_cs$param[paste0(derived_regime, ".scale")])
    est_rho <- unname(fit_cs$param[paste0(derived_regime, ".rho")])
    rho_boundary <- as.numeric(est_rho < 0.05 || est_rho > 0.95)
  }
  if (success_cs && success_prop) {
    aic_delta <- AIC(fit_prop)$AIC - AIC(fit_cs)$AIC
    loglik_gain <- as.numeric(fit_cs$logLik - fit_prop$logLik)
  }

  mapped_fraction <- sum(tree$mapped.edge[, derived_regime]) / sum(tree$edge.length)

  data.frame(
    scenario = scenario$name,
    description = scenario$description,
    rep = rep_id,
    success_corrshrink = success_cs,
    success_proportional = success_prop,
    n_tips = scenario$n_tips,
    p = nrow(scenario$base_sigma),
    with_covariate = isTRUE(scenario$with_covariate),
    true_scale = scenario$true_scale,
    true_rho = scenario$true_rho,
    est_scale = est_scale,
    est_rho = est_rho,
    scale_error = est_scale - scenario$true_scale,
    rho_error = est_rho - scenario$true_rho,
    abs_scale_error = abs(est_scale - scenario$true_scale),
    abs_rho_error = abs(est_rho - scenario$true_rho),
    aic_delta = aic_delta,
    loglik_gain = loglik_gain,
    rho_boundary = rho_boundary,
    mapped_fraction = mapped_fraction,
    stringsAsFactors = FALSE
  )
}

summarise_scenario <- function(df) {
  success_df <- df[df$success_corrshrink, , drop = FALSE]
  data.frame(
    scenario = df$scenario[1],
    description = df$description[1],
    n_tips = df$n_tips[1],
    p = df$p[1],
    covariate = df$with_covariate[1],
    true_scale = df$true_scale[1],
    true_rho = df$true_rho[1],
    mapped_fraction_mean = mean(df$mapped_fraction),
    success_rate = mean(df$success_corrshrink),
    est_scale_mean = mean(success_df$est_scale),
    est_scale_median = stats::median(success_df$est_scale),
    est_scale_sd = stats::sd(success_df$est_scale),
    est_rho_mean = mean(success_df$est_rho),
    est_rho_median = stats::median(success_df$est_rho),
    est_rho_sd = stats::sd(success_df$est_rho),
    scale_rmse = sqrt(mean(success_df$scale_error^2)),
    rho_rmse = sqrt(mean(success_df$rho_error^2)),
    rho_mae = mean(success_df$abs_rho_error),
    rho_boundary_rate = mean(success_df$rho_boundary),
    mean_delta_aic_prop_minus_corr = mean(success_df$aic_delta),
    mean_loglik_gain = mean(success_df$loglik_gain),
    stringsAsFactors = FALSE
  )
}

base_sigma_2 <- matrix(c(1.30, 0.55, 0.55, 1.00), 2, 2)
base_sigma_2_weak <- matrix(c(1.30, 0.08, 0.08, 1.00), 2, 2)
base_sigma_4 <- make_sigma(
  sds = c(1.3, 1.1, 1.0, 0.9),
  corr = matrix(
    c(
      1.00, 0.50, 0.35, 0.25,
      0.50, 1.00, 0.45, 0.30,
      0.35, 0.45, 1.00, 0.40,
      0.25, 0.30, 0.40, 1.00
    ),
    4, 4, byrow = TRUE
  )
)

scenarios <- list(
  list(
    name = "balanced_rho025_p2",
    description = "Balanced regimes, 2 traits, strong off-diagonal signal",
    n_tips = 18L,
    state_counts = c(A = 9L, B = 9L),
    base_sigma = base_sigma_2,
    true_scale = 1.8,
    true_rho = 0.25,
    with_covariate = FALSE,
    seed_offset = seed_base + 0L
  ),
  list(
    name = "balanced_rho095_p2",
    description = "Balanced regimes, 2 traits, near-proportional boundary",
    n_tips = 18L,
    state_counts = c(A = 9L, B = 9L),
    base_sigma = base_sigma_2,
    true_scale = 1.8,
    true_rho = 0.95,
    with_covariate = FALSE,
    seed_offset = seed_base + 100L
  ),
  list(
    name = "sparse_regime_rho025_p2",
    description = "Sparse derived regime, 2 traits, strong off-diagonal signal",
    n_tips = 18L,
    state_counts = c(A = 15L, B = 3L),
    base_sigma = base_sigma_2,
    true_scale = 1.8,
    true_rho = 0.25,
    with_covariate = FALSE,
    seed_offset = seed_base + 200L
  ),
  list(
    name = "weak_offdiag_rho025_p2",
    description = "Balanced regimes, 2 traits, weak base off-diagonal signal",
    n_tips = 18L,
    state_counts = c(A = 9L, B = 9L),
    base_sigma = base_sigma_2_weak,
    true_scale = 1.8,
    true_rho = 0.25,
    with_covariate = FALSE,
    seed_offset = seed_base + 300L
  ),
  list(
    name = "balanced_rho025_p4_cov",
    description = "Balanced regimes, 4 traits, one predictor",
    n_tips = 20L,
    state_counts = c(A = 10L, B = 10L),
    base_sigma = base_sigma_4,
    true_scale = 1.6,
    true_rho = 0.25,
    with_covariate = TRUE,
    beta = rbind(
      "(Intercept)" = c(0.4, -0.2, 0.1, 0.3),
      "x" = c(0.6, -0.3, 0.2, 0.4)
    ),
    seed_offset = seed_base + 400L
  )
)

cat("Running corr-shrink identifiability study with", reps, "replicates per scenario\n")

results <- do.call(
  rbind,
  lapply(scenarios, function(scenario) {
    do.call(
      rbind,
      lapply(seq_len(reps), function(rep_id) run_replicate(scenario, rep_id))
    )
  })
)

summary_table <- do.call(
  rbind,
  lapply(split(results, results$scenario), summarise_scenario)
)
summary_table <- summary_table[match(vapply(scenarios, `[[`, character(1), "name"), summary_table$scenario), ]
rownames(summary_table) <- NULL

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

print_table <- summary_table[, c(
  "scenario", "success_rate", "true_rho", "est_rho_mean", "est_rho_median", "rho_rmse",
  "true_scale", "est_scale_mean", "est_scale_median", "scale_rmse",
  "mapped_fraction_mean", "rho_boundary_rate", "mean_delta_aic_prop_minus_corr"
)]
print_table$success_rate <- fmt(print_table$success_rate, 2)
print_table$true_rho <- fmt(print_table$true_rho, 2)
print_table$est_rho_mean <- fmt(print_table$est_rho_mean, 2)
print_table$est_rho_median <- fmt(print_table$est_rho_median, 2)
print_table$rho_rmse <- fmt(print_table$rho_rmse, 2)
print_table$true_scale <- fmt(print_table$true_scale, 2)
print_table$est_scale_mean <- fmt(print_table$est_scale_mean, 2)
print_table$est_scale_median <- fmt(print_table$est_scale_median, 2)
print_table$scale_rmse <- fmt(print_table$scale_rmse, 2)
print_table$mapped_fraction_mean <- fmt(print_table$mapped_fraction_mean, 2)
print_table$rho_boundary_rate <- fmt(print_table$rho_boundary_rate, 2)
print_table$mean_delta_aic_prop_minus_corr <- fmt(print_table$mean_delta_aic_prop_minus_corr, 2)

cat("\nScenario summary\n")
print(print_table, row.names = FALSE)

best_rho <- summary_table[which.min(summary_table$rho_rmse), , drop = FALSE]
worst_rho <- summary_table[which.max(summary_table$rho_rmse), , drop = FALSE]
best_scale <- summary_table[which.min(summary_table$scale_rmse), , drop = FALSE]

cat("\nTakeaways\n")
cat(
  "- Best rho recovery:", best_rho$scenario,
  "with mean rho", fmt(best_rho$est_rho_mean, 2),
  "for true rho", fmt(best_rho$true_rho, 2),
  "and rho RMSE", fmt(best_rho$rho_rmse, 2), "\n"
)
cat(
  "- Weakest rho recovery:", worst_rho$scenario,
  "with mean rho", fmt(worst_rho$est_rho_mean, 2),
  "for true rho", fmt(worst_rho$true_rho, 2),
  "and rho RMSE", fmt(worst_rho$rho_rmse, 2), "\n"
)
cat(
  "- Best scale recovery:", best_scale$scenario,
  "with mean scale", fmt(best_scale$est_scale_mean, 2),
  "for true scale", fmt(best_scale$true_scale, 2),
  "and scale RMSE", fmt(best_scale$scale_rmse, 2), "\n"
)

low_info <- summary_table[
  summary_table$rho_rmse == max(summary_table$rho_rmse) |
    summary_table$mapped_fraction_mean == min(summary_table$mapped_fraction_mean),
  c("scenario", "description", "mapped_fraction_mean", "rho_rmse", "rho_boundary_rate")
]
rownames(low_info) <- NULL

cat("\nPotential weak-identifiability flags\n")
print(low_info, row.names = FALSE)

if (!any(results$success_corrshrink)) {
  stop("The identifiability study produced no successful corr-shrink fits.", call. = FALSE)
}

cat("\ncorr-shrink identifiability study completed\n")
