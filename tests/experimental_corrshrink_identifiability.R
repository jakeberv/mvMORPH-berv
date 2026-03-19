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

apply_corrstrength <- function(base_sigma, scale, kappa) {
  V <- diag(diag(base_sigma))
  W <- base_sigma - V
  scale * (V + kappa * W)
}

start_helper <- getFromNamespace(".mvgls_bmm_corrshrink_start", "mvMORPH")

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

fit_corrshrink <- function(formula, data, tree, start = NULL) {
  mvgls(
    formula,
    data = data,
    tree = tree,
    model = "BMM",
    method = "LL",
    REML = FALSE,
    bmm.structure = "corrshrink",
    start = start,
    echo = FALSE
  )
}

try_fit <- function(expr) {
  tryCatch(expr, error = function(e) e)
}

prepare_design <- function(formula, data) {
  model_fr <- model.frame(formula = formula, data = data)
  X <- model.matrix(attr(model_fr, "terms"), data = model_fr)
  Y <- model.response(model_fr)
  list(model_frame = model_fr, X = X, Y = as.matrix(Y))
}

fit_modes <- function(formula, data, tree, X, Y) {
  single_start <- start_helper(tree, Y, X)
  list(
    `single-start` = list(fit = try_fit(fit_corrshrink(formula, data = data, tree = tree, start = single_start)),
      start = single_start, mode = "single-start"),
    `multi-start` = list(fit = try_fit(fit_corrshrink(formula, data = data, tree = tree)),
      start = NULL, mode = "multi-start")
  )
}

run_replicate <- function(scenario, rep_id) {
  tree <- make_simmap(scenario$n_tips, scenario$state_counts, scenario$seed_offset + rep_id)
  derived_regime <- scenario$derived_regime %||% names(scenario$state_counts)[2]
  sigma_by_regime <- list(
    A = scenario$base_sigma,
    B = apply_corrstrength(scenario$base_sigma, scenario$true_scale, scenario$true_kappa)
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

  design <- prepare_design(formula, data)
  mode_fits <- fit_modes(formula, data, tree, design$X, design$Y)
  mapped_fraction <- sum(tree$mapped.edge[, derived_regime]) / sum(tree$edge.length)

  rows <- lapply(mode_fits, function(entry) {
    fit <- entry$fit
    success <- !inherits(fit, "error")
    est_scale <- NA_real_
    est_kappa <- NA_real_
    if (success) {
      est_scale <- unname(fit$param[paste0(derived_regime, ".scale")])
      est_kappa <- unname(fit$param[paste0(derived_regime, ".kappa")])
    }
    data.frame(
      scenario = scenario$name,
      description = scenario$description,
      rep = rep_id,
      mode = entry$mode,
      success_corrshrink = success,
      n_tips = scenario$n_tips,
      p = nrow(scenario$base_sigma),
      with_covariate = isTRUE(scenario$with_covariate),
      true_scale = scenario$true_scale,
      true_kappa = scenario$true_kappa,
      est_scale = est_scale,
      est_kappa = est_kappa,
      scale_error = est_scale - scenario$true_scale,
      kappa_error = est_kappa - scenario$true_kappa,
      abs_scale_error = abs(est_scale - scenario$true_scale),
      abs_kappa_error = abs(est_kappa - scenario$true_kappa),
      runaway_scale = as.numeric(is.finite(est_scale) && est_scale > 20),
      kappa_boundary = as.numeric(is.finite(est_kappa) && (est_kappa < 0.05 || est_kappa > 1.95)),
      mapped_fraction = mapped_fraction,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

summarise_mode <- function(df) {
  success_df <- df[df$success_corrshrink, , drop = FALSE]
  data.frame(
    scenario = df$scenario[1],
    description = df$description[1],
    mode = df$mode[1],
    n_tips = df$n_tips[1],
    p = df$p[1],
    covariate = df$with_covariate[1],
    true_scale = df$true_scale[1],
    true_kappa = df$true_kappa[1],
    mapped_fraction_mean = mean(df$mapped_fraction),
    success_rate = mean(df$success_corrshrink),
    est_scale_mean = mean(success_df$est_scale),
    est_scale_median = stats::median(success_df$est_scale),
    est_scale_sd = stats::sd(success_df$est_scale),
    est_kappa_mean = mean(success_df$est_kappa),
    est_kappa_median = stats::median(success_df$est_kappa),
    est_kappa_sd = stats::sd(success_df$est_kappa),
    scale_rmse = sqrt(mean(success_df$scale_error^2)),
    kappa_rmse = sqrt(mean(success_df$kappa_error^2)),
    kappa_mae = mean(success_df$abs_kappa_error),
    kappa_boundary_rate = mean(success_df$kappa_boundary),
    runaway_scale_rate = mean(success_df$runaway_scale),
    stringsAsFactors = FALSE
  )
}

base_sigma_2 <- matrix(c(1.30, 0.55, 0.55, 1.00), 2, 2)
base_sigma_2_weak <- matrix(c(1.30, 0.08, 0.08, 1.00), 2, 2)
base_sigma_2_strong <- matrix(c(1.30, 0.80, 0.80, 1.00), 2, 2)
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
    name = "balanced_kappa100_p2",
    description = "Balanced regimes, 2 traits, proportional correlation-strength",
    n_tips = 18L,
    state_counts = c(A = 9L, B = 9L),
    base_sigma = base_sigma_2,
    true_scale = 1.8,
    true_kappa = 1.0,
    with_covariate = FALSE,
    seed_offset = seed_base + 0L
  ),
  list(
    name = "balanced_kappa025_p2",
    description = "Balanced regimes, 2 traits, weaker-than-reference correlation-strength",
    n_tips = 18L,
    state_counts = c(A = 9L, B = 9L),
    base_sigma = base_sigma_2,
    true_scale = 1.8,
    true_kappa = 0.25,
    with_covariate = FALSE,
    seed_offset = seed_base + 100L
  ),
  list(
    name = "stronger_kappa140_p2",
    description = "Balanced regimes, 2 traits, stronger-than-reference correlation-strength",
    n_tips = 18L,
    state_counts = c(A = 9L, B = 9L),
    base_sigma = base_sigma_2_strong,
    true_scale = 1.8,
    true_kappa = 1.40,
    with_covariate = FALSE,
    seed_offset = seed_base + 200L
  ),
  list(
    name = "sparse_regime_kappa025_p2",
    description = "Sparse derived regime, 2 traits, weaker-than-reference correlation-strength",
    n_tips = 18L,
    state_counts = c(A = 15L, B = 3L),
    base_sigma = base_sigma_2,
    true_scale = 1.8,
    true_kappa = 0.25,
    with_covariate = FALSE,
    seed_offset = seed_base + 300L
  ),
  list(
    name = "weak_offdiag_kappa025_p2",
    description = "Balanced regimes, 2 traits, weak base off-diagonal signal",
    n_tips = 18L,
    state_counts = c(A = 9L, B = 9L),
    base_sigma = base_sigma_2_weak,
    true_scale = 1.8,
    true_kappa = 0.25,
    with_covariate = FALSE,
    seed_offset = seed_base + 400L
  ),
  list(
    name = "balanced_kappa100_p4_cov",
    description = "Balanced regimes, 4 traits, one predictor",
    n_tips = 20L,
    state_counts = c(A = 10L, B = 10L),
    base_sigma = base_sigma_4,
    true_scale = 1.6,
    true_kappa = 1.0,
    with_covariate = TRUE,
    beta = rbind(
      "(Intercept)" = c(0.4, -0.2, 0.1, 0.3),
      "x" = c(0.6, -0.3, 0.2, 0.4)
    ),
    seed_offset = seed_base + 500L
  )
)

cat("Running corr-strength identifiability study with", reps, "replicates per scenario\n")

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
  lapply(split(results, list(results$scenario, results$mode), drop = TRUE), summarise_mode)
)
summary_table <- summary_table[order(match(summary_table$scenario, vapply(scenarios, `[[`, character(1), "name")),
                                    match(summary_table$mode, c("single-start", "multi-start"))), ]
rownames(summary_table) <- NULL

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

print_table <- summary_table[, c(
  "scenario", "mode", "success_rate", "true_kappa", "est_kappa_mean", "est_kappa_median", "kappa_rmse",
  "true_scale", "est_scale_mean", "est_scale_median", "scale_rmse",
  "mapped_fraction_mean", "kappa_boundary_rate", "runaway_scale_rate"
)]
print_table$success_rate <- fmt(print_table$success_rate, 2)
print_table$true_kappa <- fmt(print_table$true_kappa, 2)
print_table$est_kappa_mean <- fmt(print_table$est_kappa_mean, 2)
print_table$est_kappa_median <- fmt(print_table$est_kappa_median, 2)
print_table$kappa_rmse <- fmt(print_table$kappa_rmse, 2)
print_table$true_scale <- fmt(print_table$true_scale, 2)
print_table$est_scale_mean <- fmt(print_table$est_scale_mean, 2)
print_table$est_scale_median <- fmt(print_table$est_scale_median, 2)
print_table$scale_rmse <- fmt(print_table$scale_rmse, 2)
print_table$mapped_fraction_mean <- fmt(print_table$mapped_fraction_mean, 2)
print_table$kappa_boundary_rate <- fmt(print_table$kappa_boundary_rate, 2)
print_table$runaway_scale_rate <- fmt(print_table$runaway_scale_rate, 2)

cat("\nScenario summary\n")
print(print_table, row.names = FALSE)

comparison_table <- do.call(
  rbind,
  lapply(split(summary_table, summary_table$scenario), function(df) {
    single <- df[df$mode == "single-start", , drop = FALSE]
    multi <- df[df$mode == "multi-start", , drop = FALSE]
    data.frame(
      scenario = single$scenario[1],
      kappa_rmse_single = single$kappa_rmse,
      kappa_rmse_multi = multi$kappa_rmse,
      scale_rmse_single = single$scale_rmse,
      scale_rmse_multi = multi$scale_rmse,
      kappa_boundary_single = single$kappa_boundary_rate,
      kappa_boundary_multi = multi$kappa_boundary_rate,
      runaway_scale_single = single$runaway_scale_rate,
      runaway_scale_multi = multi$runaway_scale_rate,
      stringsAsFactors = FALSE
    )
  })
)

cat("\nSingle-start vs multi-start comparison\n")
cmp_print <- comparison_table
cmp_print$kappa_rmse_single <- fmt(cmp_print$kappa_rmse_single, 2)
cmp_print$kappa_rmse_multi <- fmt(cmp_print$kappa_rmse_multi, 2)
cmp_print$scale_rmse_single <- fmt(cmp_print$scale_rmse_single, 2)
cmp_print$scale_rmse_multi <- fmt(cmp_print$scale_rmse_multi, 2)
cmp_print$kappa_boundary_single <- fmt(cmp_print$kappa_boundary_single, 2)
cmp_print$kappa_boundary_multi <- fmt(cmp_print$kappa_boundary_multi, 2)
cmp_print$runaway_scale_single <- fmt(cmp_print$runaway_scale_single, 2)
cmp_print$runaway_scale_multi <- fmt(cmp_print$runaway_scale_multi, 2)
print(cmp_print, row.names = FALSE)

improvement_table <- do.call(
  rbind,
  lapply(split(comparison_table, comparison_table$scenario), function(df) {
    data.frame(
      scenario = df$scenario[1],
      delta_kappa_rmse = df$kappa_rmse_multi - df$kappa_rmse_single,
      delta_scale_rmse = df$scale_rmse_multi - df$scale_rmse_single,
      delta_kappa_boundary = df$kappa_boundary_multi - df$kappa_boundary_single,
      delta_runaway_scale = df$runaway_scale_multi - df$runaway_scale_single,
      stringsAsFactors = FALSE
    )
  })
)

cat("\nMulti-start deltas (multi minus single)\n")
imp_print <- improvement_table
imp_print$delta_kappa_rmse <- fmt(imp_print$delta_kappa_rmse, 2)
imp_print$delta_scale_rmse <- fmt(imp_print$delta_scale_rmse, 2)
imp_print$delta_kappa_boundary <- fmt(imp_print$delta_kappa_boundary, 2)
imp_print$delta_runaway_scale <- fmt(imp_print$delta_runaway_scale, 2)
print(imp_print, row.names = FALSE)

low_info <- summary_table[
  summary_table$kappa_rmse == max(summary_table$kappa_rmse) |
    summary_table$runaway_scale_rate == max(summary_table$runaway_scale_rate),
  c("scenario", "mode", "mapped_fraction_mean", "kappa_rmse", "kappa_boundary_rate", "runaway_scale_rate")
]
rownames(low_info) <- NULL

cat("\nPotential weak-identifiability flags\n")
print(low_info, row.names = FALSE)

if (!any(results$success_corrshrink)) {
  stop("The identifiability study produced no successful corr-strength fits.", call. = FALSE)
}

cat("\ncorr-strength identifiability study completed\n")
