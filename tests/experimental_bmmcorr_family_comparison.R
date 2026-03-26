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
`%||%` <- function(x, y) if (is.null(x)) y else x

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
  save_csv <- Sys.getenv("BMMCORR_FAMILY_COMPARISON_SAVE_CSV", "")
  save_rds <- Sys.getenv("BMMCORR_FAMILY_COMPARISON_SAVE_RDS", "")
  if (nzchar(save_csv)) {
    dir.create(dirname(save_csv), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(results_df, save_csv, row.names = FALSE)
  }
  if (nzchar(save_rds)) {
    dir.create(dirname(save_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(results_df, save_rds)
  }
}

make_simmap <- function(seed, n_tips, states_per_regime) {
  set.seed(seed)
  tree <- phytools::pbtree(n = n_tips, scale = 1)
  sampled_states <- sample(rep(names(states_per_regime), times = states_per_regime), n_tips)
  states <- setNames(sampled_states, tree$tip.label)
  suppressMessages(phytools::make.simmap(tree, states, model = "ER", nsim = 1))
}

base_sigma_from_signal <- function(signal) {
  matrix(c(
    1.20, 0.25 * signal, 0.14 * signal, 0.08 * signal,
    0.25 * signal, 1.10, 0.18 * signal, 0.12 * signal,
    0.14 * signal, 0.18 * signal, 1.00, 0.22 * signal,
    0.08 * signal, 0.12 * signal, 0.22 * signal, 0.90
  ), 4, 4)
}

apply_corrpower <- function(base_sigma, corr_power, scale = 1) {
  eig <- eigen(stats::cov2cor(base_sigma), symmetric = TRUE)
  vals <- pmax(eig$values, .Machine$double.eps)
  corr_raw <- eig$vectors %*% diag(vals ^ corr_power, nrow = length(vals)) %*% t(eig$vectors)
  corr_mat <- stats::cov2cor(0.5 * (corr_raw + t(corr_raw)))
  D <- diag(sqrt(diag(base_sigma)))
  scale * (D %*% corr_mat %*% D)
}

simulate_response <- function(tree, sigma_by_regime, seed) {
  set.seed(seed)
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
  as.matrix(y)
}

fit_bmm <- function(formula, data, tree, structure) {
  mvgls(
    formula,
    data = data,
    tree = tree,
    model = "BMM",
    method = "LL",
    REML = FALSE,
    bmm.structure = structure,
    echo = FALSE
  )
}

cov_error <- function(fitted_sigma, true_sigma) {
  regime_names <- names(true_sigma)
  errs <- vapply(regime_names, function(regime_name) {
    sqrt(sum((fitted_sigma[[regime_name]] - true_sigma[[regime_name]])^2))
  }, numeric(1))
  mean(errs)
}

diag_flag <- function(fit, field, default = NA) {
  if (identical(fit$bmm.structure, "corrpower")) {
    return(fit$diagnostics$corrpower[[field]] %||% default)
  }
  if (identical(fit$bmm.structure, "corrpower_coronly")) {
    return(fit$diagnostics$corrpower_coronly[[field]] %||% default)
  }
  default
}

compare_case <- function(case_type, relation, n_tips, signal, replicate, seed) {
  tree <- make_simmap(seed, n_tips = n_tips, states_per_regime = c(A = n_tips %/% 2, B = n_tips - n_tips %/% 2))
  base_sigma <- base_sigma_from_signal(signal)
  if (relation == "weaker") {
    corr_power_B <- 0.45
    scale_B <- if (case_type == "full") 1.55 else 1.00
  } else {
    corr_power_B <- 1.40
    scale_B <- if (case_type == "full") 1.70 else 1.00
  }
  true_sigma <- list(
    A = base_sigma,
    B = apply_corrpower(base_sigma, corr_power = corr_power_B, scale = scale_B)
  )
  Y <- simulate_response(tree, true_sigma, seed = seed + 100L)

  fit_prop <- fit_bmm(Y ~ 1, data = list(Y = Y), tree = tree, structure = "proportional")
  fit_coronly <- fit_bmm(Y ~ 1, data = list(Y = Y), tree = tree, structure = "corrpower_coronly")
  fit_corrpower <- fit_bmm(Y ~ 1, data = list(Y = Y), tree = tree, structure = "corrpower")

  data.frame(
    case_type = case_type,
    relation = relation,
    n_tips = n_tips,
    signal = signal,
    replicate = replicate,
    proportional_cov_err = cov_error(fit_prop$sigma$regime[names(true_sigma)], true_sigma),
    corrpower_coronly_cov_err = cov_error(fit_coronly$sigma$regime[names(true_sigma)], true_sigma),
    corrpower_cov_err = cov_error(fit_corrpower$sigma$regime[names(true_sigma)], true_sigma),
    proportional_logLik = as.numeric(fit_prop$logLik),
    corrpower_coronly_logLik = as.numeric(fit_coronly$logLik),
    corrpower_logLik = as.numeric(fit_corrpower$logLik),
    corrpower_coronly_boundary = isTRUE(diag_flag(fit_coronly, "boundary_corr_power", FALSE)),
    corrpower_boundary = isTRUE(diag_flag(fit_corrpower, "boundary_corr_power", FALSE)),
    corrpower_coronly_pathological = isTRUE(diag_flag(fit_coronly, "pathological_corr_power", FALSE)),
    corrpower_pathological = isTRUE(diag_flag(fit_corrpower, "pathological_scale", FALSE)) ||
      isTRUE(diag_flag(fit_corrpower, "pathological_corr_power", FALSE)),
    stringsAsFactors = FALSE
  )
}

full_run <- env_flag("BMMCORR_FAMILY_COMPARISON_FULL", FALSE)
reps <- env_int("BMMCORR_FAMILY_COMPARISON_REPS", if (full_run) 5L else 2L)
if (is.na(reps) || reps < 1L) stop("BMMCORR_FAMILY_COMPARISON_REPS must be >= 1", call. = FALSE)

grid <- expand.grid(
  case_type = c("coronly", "full"),
  relation = c("weaker", "stronger"),
  n_tips = if (full_run) c(50L, 100L) else 50L,
  signal = if (full_run) c(0.20, 0.50) else 0.20,
  replicate = seq_len(reps),
  stringsAsFactors = FALSE
)

results <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  row <- grid[i, , drop = FALSE]
  compare_case(
    case_type = row$case_type[[1]],
    relation = row$relation[[1]],
    n_tips = row$n_tips[[1]],
    signal = row$signal[[1]],
    replicate = row$replicate[[1]],
    seed = 20260326L + i
  )
}))

write_outputs(results)

coronly_cases <- results[results$case_type == "coronly", , drop = FALSE]
full_cases <- results[results$case_type == "full", , drop = FALSE]

assert_true(
  stats::median(coronly_cases$corrpower_coronly_cov_err) <= stats::median(coronly_cases$corrpower_cov_err) + 0.05,
  "corrpower_coronly should match or beat corrpower in correlation-only scenarios"
)
assert_true(
  mean(coronly_cases$corrpower_coronly_boundary) <= mean(coronly_cases$corrpower_boundary) + 0.10,
  "corrpower_coronly should not have materially worse boundary behavior in correlation-only scenarios"
)
assert_true(
  mean(coronly_cases$corrpower_coronly_pathological) <= mean(coronly_cases$corrpower_pathological) + 0.05,
  "corrpower_coronly should not have materially worse pathology in correlation-only scenarios"
)
assert_true(
  stats::median(full_cases$corrpower_cov_err) <= stats::median(full_cases$corrpower_coronly_cov_err) + 0.05,
  "full corrpower should retain the advantage in scale-plus-correlation scenarios"
)

cat("bmm correlation-family comparison harness checks passed\n")
