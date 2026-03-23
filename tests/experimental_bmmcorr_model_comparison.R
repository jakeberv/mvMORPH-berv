#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
if (!length(file_arg)) stop("This script must be run with Rscript.", call. = FALSE)
repo_root <- normalizePath(file.path(dirname(file_arg), ".."))

use_installed <- identical(toupper(Sys.getenv("MV_MORPH_USE_INSTALLED", "FALSE")), "TRUE")

if (use_installed) {
  suppressPackageStartupMessages(library(mvMORPH))
} else {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required to run this harness.", call. = FALSE)
  }
  pkgload::load_all(repo_root, quiet = TRUE)
}

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

env_flag <- function(name, default = FALSE) {
  val <- toupper(Sys.getenv(name, if (default) "TRUE" else "FALSE"))
  identical(val, "TRUE")
}

env_int <- function(name, default = NA_integer_) {
  val <- Sys.getenv(name, "")
  if (!nzchar(val)) return(default)
  as.integer(val)
}

write_outputs <- function(results_df) {
  save_csv <- Sys.getenv("CORRPOWER_COMPARISON_SAVE_CSV", "")
  save_rds <- Sys.getenv("CORRPOWER_COMPARISON_SAVE_RDS", "")
  if (nzchar(save_csv)) {
    dir.create(dirname(save_csv), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(results_df, save_csv, row.names = FALSE)
  }
  if (nzchar(save_rds)) {
    dir.create(dirname(save_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(results_df, save_rds)
  }
}

regime_summary_from_sigma <- function(sigma_by_regime) {
  regime_names <- names(sigma_by_regime)
  if (is.null(regime_names)) regime_names <- seq_along(sigma_by_regime)
  p <- nrow(sigma_by_regime[[1]])
  rows <- lapply(regime_names, function(regime_name) {
    Sigma <- sigma_by_regime[[regime_name]]
    diag_vals <- diag(Sigma)
    cov_vals <- Sigma[upper.tri(Sigma)]
    cor_vals <- try({
      corr <- stats::cov2cor(Sigma)
      corr[upper.tri(corr)]
    }, silent = TRUE)
    if (inherits(cor_vals, "try-error")) cor_vals <- rep(NA_real_, length(cov_vals))
    data.frame(
      regime = regime_name,
      mean_rate = sum(diag_vals) / p,
      mean_variance = mean(diag_vals),
      mean_covariance = if (length(cov_vals)) mean(cov_vals) else 0,
      mean_correlation = if (length(cor_vals)) mean(cor_vals) else 0,
      mean_abs_correlation = if (length(cor_vals)) mean(abs(cor_vals)) else 0,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- out$regime
  out
}

cov_error_metrics <- function(fitted_sigma, true_sigma) {
  regime_names <- names(true_sigma)
  errs <- vapply(regime_names, function(regime_name) {
    sqrt(sum((fitted_sigma[[regime_name]] - true_sigma[[regime_name]])^2))
  }, numeric(1))
  list(
    err_A = unname(errs[[1]]),
    err_B = unname(errs[[2]]),
    err_mean = mean(errs)
  )
}

summary_error_metric <- function(fit_summary, true_summary) {
  metrics <- c(
    "mean_rate",
    "mean_variance",
    "mean_covariance",
    "mean_correlation",
    "mean_abs_correlation"
  )
  common_regimes <- intersect(rownames(true_summary), rownames(fit_summary))
  if (!length(common_regimes)) return(NA_real_)
  fit_vals <- as.matrix(fit_summary[common_regimes, metrics, drop = FALSE])
  true_vals <- as.matrix(true_summary[common_regimes, metrics, drop = FALSE])
  sqrt(mean((fit_vals - true_vals)^2))
}

apply_corrpower <- function(base_sigma, scale, corr_power) {
  eig <- eigen(stats::cov2cor(base_sigma), symmetric = TRUE)
  vals <- pmax(eig$values, .Machine$double.eps)
  corr_raw <- eig$vectors %*% diag(vals ^ corr_power, nrow = length(vals)) %*% t(eig$vectors)
  corr_mat <- stats::cov2cor(0.5 * (corr_raw + t(corr_raw)))
  D <- diag(sqrt(diag(base_sigma)))
  scale * (D %*% corr_mat %*% D)
}

make_fossil_tree <- function(seed, n_tips, fossil_fraction) {
  set.seed(seed)
  tree <- phytools::pbtree(n = n_tips, scale = 1)
  if (fossil_fraction <= 0) return(tree)
  n_fossil <- max(1L, round(n_tips * fossil_fraction))
  fossil_tips <- sample(tree$tip.label, n_fossil)
  fossil_edge_ids <- match(fossil_tips, tree$tip.label)
  tip_rows <- match(fossil_edge_ids, tree$edge[, 2])
  tree$edge.length[tip_rows] <- tree$edge.length[tip_rows] * 0.25
  tree
}

make_simmap <- function(seed, n_tips, states_per_regime) {
  set.seed(seed)
  tree <- phytools::pbtree(n = n_tips, scale = 1)
  sampled_states <- sample(rep(names(states_per_regime), times = states_per_regime), n_tips)
  states <- setNames(sampled_states, tree$tip.label)
  suppressMessages(phytools::make.simmap(tree, states, model = "ER", nsim = 1))
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

fit_corrpower <- function(formula, data, tree) {
  mvgls(
    formula,
    data = data,
    tree = tree,
    model = "BMM",
    method = "LL",
    REML = FALSE,
    bmm.structure = "corrpower",
    echo = FALSE
  )
}

fit_corrstrength <- function(formula, data, tree) {
  mvgls(
    formula,
    data = data,
    tree = tree,
    model = "BMM",
    method = "LL",
    REML = FALSE,
    bmm.structure = "corrstrength",
    echo = FALSE
  )
}

compare_case <- function(n_tips, p, fossil_fraction, relation, seed, replicate) {
  tree <- make_fossil_tree(seed, n_tips = n_tips, fossil_fraction = fossil_fraction)
  if (p == 2L) {
    base_sigma <- if (relation == "weaker") matrix(c(1.3, 0.08, 0.08, 1.0), 2, 2) else matrix(c(1.3, 0.55, 0.55, 1.0), 2, 2)
    derived_sigma <- if (relation == "weaker") apply_corrpower(base_sigma, 1.4, 0.45) else apply_corrpower(base_sigma, 1.8, 1.4)
  } else {
    base_sigma <- matrix(c(
      1.2, 0.20, 0.10, 0.05,
      0.20, 1.1, 0.15, 0.10,
      0.10, 0.15, 1.0, 0.20,
      0.05, 0.10, 0.20, 0.9
    ), 4, 4)
    derived_sigma <- if (relation == "weaker") apply_corrpower(base_sigma, 1.4, 0.45) else apply_corrpower(base_sigma, 1.9, 1.3)
  }
  sampled_states <- c(rep("A", round(length(tree$tip.label) * 0.6)), rep("B", length(tree$tip.label) - round(length(tree$tip.label) * 0.6)))
  sampled_states <- sample(sampled_states, length(sampled_states), replace = FALSE)
  tree <- phytools::make.simmap(tree, setNames(sampled_states, tree$tip.label), model = "ER", nsim = 1)
  true_sigma <- list(A = base_sigma, B = derived_sigma)
  Y <- simulate_response(tree, true_sigma, seed = seed + 50L)
  fit_power <- fit_corrpower(Y ~ 1, data = list(Y = Y), tree = tree)
  fit_strength <- fit_corrstrength(Y ~ 1, data = list(Y = Y), tree = tree)
  power_cov_metrics <- cov_error_metrics(fit_power$sigma$regime[names(true_sigma)], true_sigma)
  strength_cov_metrics <- cov_error_metrics(fit_strength$sigma$regime[names(true_sigma)], true_sigma)
  true_summary <- regime_summary_from_sigma(true_sigma)
  data.frame(
    replicate = replicate,
    n_tips = n_tips,
    p = p,
    fossil_fraction = fossil_fraction,
    relation = relation,
    corrpower_cov_err_A = power_cov_metrics$err_A,
    corrpower_cov_err_B = power_cov_metrics$err_B,
    corrpower_cov_err = power_cov_metrics$err_mean,
    corrstrength_cov_err_A = strength_cov_metrics$err_A,
    corrstrength_cov_err_B = strength_cov_metrics$err_B,
    corrstrength_cov_err = strength_cov_metrics$err_mean,
    corrpower_summary_err = summary_error_metric(fit_power$regime.summary, true_summary),
    corrstrength_summary_err = summary_error_metric(fit_strength$regime.summary, true_summary),
    corrpower_pathological = isTRUE(getElement(fit_power$diagnostics$corrpower, "pathological_scale")),
    corrstrength_pathological = isTRUE(getElement(fit_strength$diagnostics$corrstrength, "pathological_scale")),
    corrpower_boundary = isTRUE(getElement(fit_power$diagnostics$corrpower, "boundary_corr_power")),
    corrstrength_boundary = isTRUE(getElement(fit_strength$diagnostics$corrstrength, "boundary_corr_strength")),
    logLik_power = as.numeric(fit_power$logLik),
    logLik_strength = as.numeric(fit_strength$logLik),
    stringsAsFactors = FALSE
  )
}

reps <- env_int("CORRPOWER_COMPARISON_REPS", 1L)
if (is.na(reps) || reps < 1L) stop("CORRPOWER_COMPARISON_REPS must be >= 1", call. = FALSE)

full_grid <- expand.grid(
  n_tips = c(50L, 100L),
  p = c(4L),
  fossil_fraction = c(0, 0.5),
  relation = c("weaker", "stronger"),
  replicate = seq_len(reps),
  stringsAsFactors = FALSE
)
full_grid$scenario_id <- seq_len(nrow(full_grid))

grid <- if (env_flag("CORRPOWER_COMPARISON_FULL", FALSE)) {
  full_grid
} else {
  subset(full_grid, n_tips == 50L)
}

chunk_total <- env_int("CORRPOWER_COMPARISON_CHUNK_TOTAL", 1L)
chunk_index <- env_int("CORRPOWER_COMPARISON_CHUNK_INDEX", 1L)
if (is.na(chunk_total) || chunk_total < 1L) stop("CORRPOWER_COMPARISON_CHUNK_TOTAL must be >= 1", call. = FALSE)
if (is.na(chunk_index) || chunk_index < 1L || chunk_index > chunk_total) {
  stop("CORRPOWER_COMPARISON_CHUNK_INDEX must be between 1 and CORRPOWER_COMPARISON_CHUNK_TOTAL", call. = FALSE)
}
if (chunk_total > 1L) {
  keep <- ((seq_len(nrow(grid)) - 1L) %% chunk_total) == (chunk_index - 1L)
  grid <- grid[keep, , drop = FALSE]
}
if (!nrow(grid)) {
  res <- full_grid[0, , drop = FALSE]
  res$corrpower_cov_err_A <- numeric(0)
  res$corrpower_cov_err_B <- numeric(0)
  res$corrpower_cov_err <- numeric(0)
  res$corrstrength_cov_err_A <- numeric(0)
  res$corrstrength_cov_err_B <- numeric(0)
  res$corrstrength_cov_err <- numeric(0)
  res$corrpower_summary_err <- numeric(0)
  res$corrstrength_summary_err <- numeric(0)
  res$corrpower_pathological <- logical(0)
  res$corrstrength_pathological <- logical(0)
  res$corrpower_boundary <- logical(0)
  res$corrstrength_boundary <- logical(0)
  res$logLik_power <- numeric(0)
  res$logLik_strength <- numeric(0)
  write_outputs(res)
  cat("No scenarios assigned to this corr-power comparison chunk\n")
  quit(save = "no", status = 0L)
}

res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  compare_case(
    g$n_tips,
    g$p,
    g$fossil_fraction,
    g$relation,
    seed = 20260700 + g$scenario_id * 100L,
    replicate = g$replicate
  )
}))

write_outputs(res)
print(res)

cat("\nComparison summary\n")
cat(sprintf("corrpower median covariance error: %.4f\n", median(res$corrpower_cov_err, na.rm = TRUE)))
cat(sprintf("corrstrength median covariance error: %.4f\n", median(res$corrstrength_cov_err, na.rm = TRUE)))
cat(sprintf("corrpower median regime-summary error: %.4f\n", median(res$corrpower_summary_err, na.rm = TRUE)))
cat(sprintf("corrstrength median regime-summary error: %.4f\n", median(res$corrstrength_summary_err, na.rm = TRUE)))
cat(sprintf("corrpower pathology rate: %.3f\n", mean(res$corrpower_pathological, na.rm = TRUE)))
cat(sprintf("corrstrength pathology rate: %.3f\n", mean(res$corrstrength_pathological, na.rm = TRUE)))

assert_true(any(is.finite(res$corrpower_cov_err)), "corrpower comparison produced no finite covariance errors")
assert_true(any(is.finite(res$corrstrength_cov_err)), "corrstrength comparison produced no finite covariance errors")

cat("bmm correlation model comparison harness checks passed\n")
