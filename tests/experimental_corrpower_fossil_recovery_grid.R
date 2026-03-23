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

env_flag <- function(name, default = FALSE) {
  val <- toupper(Sys.getenv(name, if (default) "TRUE" else "FALSE"))
  identical(val, "TRUE")
}

env_int <- function(name, default = NA_integer_) {
  val <- Sys.getenv(name, "")
  if (!nzchar(val)) return(default)
  as.integer(val)
}

env_num <- function(name, default = NA_real_) {
  val <- Sys.getenv(name, "")
  if (!nzchar(val)) return(default)
  as.numeric(val)
}

lambda_scale <- env_num("CORRPOWER_FOSSIL_RECOVERY_LAMBDA_SCALE", 0.05)
lambda_corr_power <- env_num("CORRPOWER_FOSSIL_RECOVERY_LAMBDA_CORR_POWER", 0.05)
options(
  mvMORPH.corrpower.lambda_scale = lambda_scale,
  mvMORPH.corrpower.lambda_corr_power = lambda_corr_power
)

write_outputs <- function(results_df) {
  save_csv <- Sys.getenv("CORRPOWER_FOSSIL_RECOVERY_SAVE_CSV", "")
  save_rds <- Sys.getenv("CORRPOWER_FOSSIL_RECOVERY_SAVE_RDS", "")
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

cov_rel_frob_metrics <- function(fitted_sigma, true_sigma) {
  regime_names <- names(true_sigma)
  errs <- vapply(regime_names, function(regime_name) {
    num <- sqrt(sum((fitted_sigma[[regime_name]] - true_sigma[[regime_name]])^2))
    den <- sqrt(sum(true_sigma[[regime_name]]^2))
    num / den
  }, numeric(1))
  names(errs) <- regime_names
  list(
    err_A = unname(errs[["A"]]),
    err_B = unname(errs[["B"]]),
    err_mean = mean(errs)
  )
}

summary_error_metrics <- function(fit_summary, true_summary) {
  metrics <- c(
    "mean_rate",
    "mean_variance",
    "mean_covariance",
    "mean_correlation",
    "mean_abs_correlation"
  )
  common_regimes <- intersect(rownames(true_summary), rownames(fit_summary))
  if (!length(common_regimes)) {
    return(list(
      rmse = NA_real_,
      B_mean_rate_abs_err = NA_real_,
      B_mean_abs_correlation_abs_err = NA_real_
    ))
  }
  fit_vals <- as.matrix(fit_summary[common_regimes, metrics, drop = FALSE])
  true_vals <- as.matrix(true_summary[common_regimes, metrics, drop = FALSE])
  list(
    rmse = sqrt(mean((fit_vals - true_vals)^2)),
    B_mean_rate_abs_err = abs(fit_summary["B", "mean_rate"] - true_summary["B", "mean_rate"]),
    B_mean_abs_correlation_abs_err = abs(fit_summary["B", "mean_abs_correlation"] - true_summary["B", "mean_abs_correlation"])
  )
}

apply_corrpower <- function(base_sigma, scale, corr_power) {
  eig <- eigen(stats::cov2cor(base_sigma), symmetric = TRUE)
  vals <- pmax(eig$values, .Machine$double.eps)
  corr_raw <- eig$vectors %*% diag(vals ^ corr_power, nrow = length(vals)) %*% t(eig$vectors)
  corr_mat <- stats::cov2cor(0.5 * (corr_raw + t(corr_raw)))
  D <- diag(sqrt(diag(base_sigma)))
  scale * (D %*% corr_mat %*% D)
}

make_extant_tree <- function(seed, n_extant) {
  set.seed(seed)
  phytools::pbtree(n = n_extant, scale = 1)
}

choose_attachment_height <- function(tree_height, fossil_depth) {
  frac <- switch(
    fossil_depth,
    deep = stats::runif(1, 0.20, 0.45),
    shallow = stats::runif(1, 0.70, 0.95),
    stats::runif(1, 0.20, 0.95)
  )
  tree_height * frac
}

add_fossil_tips <- function(extant_tree, n_fossil, fossil_depth, seed) {
  if (n_fossil <= 0L) {
    return(list(full = extant_tree, fossil_labels = character(0), extant_labels = extant_tree$tip.label))
  }

  set.seed(seed)
  tree <- extant_tree
  extant_labels <- extant_tree$tip.label
  fossil_labels <- character(n_fossil)

  for (i in seq_len(n_fossil)) {
    heights <- phytools::nodeHeights(tree)
    tree_height <- max(heights)
    target_height <- choose_attachment_height(tree_height, fossil_depth)
    eligible <- which(heights[, 1] < target_height & heights[, 2] > target_height)
    if (!length(eligible)) {
      eligible <- which(heights[, 2] > target_height)
    }
    if (!length(eligible)) {
      eligible <- seq_len(nrow(tree$edge))
    }
    edge_id <- sample(eligible, 1L)
    where <- tree$edge[edge_id, 2]
    position <- max(heights[edge_id, 2] - target_height, 0)
    fossil_label <- sprintf("fossil_%03d", i)
    fossil_labels[[i]] <- fossil_label
    tree <- phytools::bind.tip(
      tree,
      tip.label = fossil_label,
      edge.length = 0,
      where = where,
      position = position
    )
  }

  list(full = tree, fossil_labels = fossil_labels, extant_labels = extant_labels)
}

choose_shift_node <- function(tree, target_fraction) {
  total_length <- sum(tree$edge.length)
  internal_nodes <- seq_len(tree$Nnode) + ape::Ntip(tree)
  root_node <- ape::Ntip(tree) + 1L
  candidates <- setdiff(internal_nodes, root_node)
  if (!length(candidates)) return(root_node)

  score <- vapply(candidates, function(node) {
    clade <- ape::extract.clade(tree, node)
    clade_len <- sum(clade$edge.length)
    stem_row <- which(tree$edge[, 2] == node)
    stem_len <- if (length(stem_row)) tree$edge.length[stem_row[[1]]] else 0
    abs((clade_len + stem_len) / total_length - target_fraction)
  }, numeric(1))
  candidates[[which.min(score)]]
}

paint_two_regime_tree <- function(tree, mapped_fraction) {
  node <- choose_shift_node(tree, mapped_fraction)
  phytools::paintSubTree(tree, node = node, state = "B", anc.state = "A", stem = TRUE)
}

mapped_fraction_for_regime <- function(tree, regime_name) {
  mapped <- tree$mapped.edge
  if (is.null(mapped) || is.null(colnames(mapped)) || !regime_name %in% colnames(mapped)) return(NA_real_)
  sum(mapped[, regime_name], na.rm = TRUE) / sum(tree$edge.length)
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

make_base_sigma <- function(p, base_signal) {
  if (p == 2L) {
    if (base_signal < 0.3) {
      matrix(c(1.3, 0.08, 0.08, 1.0), 2, 2)
    } else {
      matrix(c(1.3, 0.35, 0.35, 1.0), 2, 2)
    }
  } else {
    if (base_signal < 0.3) {
      matrix(c(
        1.2, 0.08, 0.05, 0.03,
        0.08, 1.1, 0.06, 0.04,
        0.05, 0.06, 1.0, 0.08,
        0.03, 0.04, 0.08, 0.9
      ), 4, 4)
    } else {
      matrix(c(
        1.2, 0.22, 0.12, 0.08,
        0.22, 1.1, 0.18, 0.12,
        0.12, 0.18, 1.0, 0.20,
        0.08, 0.12, 0.20, 0.9
      ), 4, 4)
    }
  }
}

make_true_sigma <- function(p, base_signal, relation) {
  base_sigma <- make_base_sigma(p, base_signal)
  if (relation == "weaker") {
    scale_b <- if (p == 2L) 1.4 else 1.5
    corr_b <- if (base_signal < 0.3) 0.45 else 0.70
  } else {
    scale_b <- if (p == 2L) 1.8 else 1.9
    corr_b <- if (base_signal < 0.3) 1.30 else 1.50
  }
  list(
    sigma = list(A = base_sigma, B = apply_corrpower(base_sigma, scale = scale_b, corr_power = corr_b)),
    scale_b = scale_b,
    corr_power_b = corr_b
  )
}

compare_case <- function(sc) {
  base_seed <- 20261000L + sc$scenario_id * 1000L + sc$replicate
  n_fossil <- if (sc$fossil_fraction <= 0) 0L else max(1L, round(sc$n_extant * sc$fossil_fraction))

  extant_plain <- make_extant_tree(base_seed, sc$n_extant)
  fossilized <- add_fossil_tips(extant_plain, n_fossil = n_fossil, fossil_depth = sc$fossil_depth, seed = base_seed + 1L)
  full_plain <- fossilized$full
  full_simmap <- paint_two_regime_tree(full_plain, sc$mapped_fraction)
  extant_simmap <- if (length(fossilized$fossil_labels)) {
    ape::drop.tip(full_simmap, fossilized$fossil_labels)
  } else {
    full_simmap
  }

  true_par <- make_true_sigma(sc$p, sc$base_signal, sc$relation)
  true_sigma <- true_par$sigma
  true_summary <- regime_summary_from_sigma(true_sigma)

  Y_full <- simulate_response(full_simmap, true_sigma, seed = base_seed + 2L)
  Y_extant <- Y_full[extant_simmap$tip.label, , drop = FALSE]

  fit_full <- fit_corrpower(Y ~ 1, data = list(Y = Y_full), tree = full_simmap)
  fit_extant <- fit_corrpower(Y ~ 1, data = list(Y = Y_extant), tree = extant_simmap)

  full_cov <- cov_rel_frob_metrics(fit_full$sigma$regime[names(true_sigma)], true_sigma)
  extant_cov <- cov_rel_frob_metrics(fit_extant$sigma$regime[names(true_sigma)], true_sigma)
  full_summary <- summary_error_metrics(fit_full$regime.summary, true_summary)
  extant_summary <- summary_error_metrics(fit_extant$regime.summary, true_summary)

  data.frame(
    scenario_id = sc$scenario_id,
    replicate = sc$replicate,
    n_extant = sc$n_extant,
    n_fossil = n_fossil,
    p = sc$p,
    fossil_fraction = sc$fossil_fraction,
    fossil_depth = sc$fossil_depth,
    mapped_fraction_target = sc$mapped_fraction,
    mapped_fraction_full_B = mapped_fraction_for_regime(full_simmap, "B"),
    mapped_fraction_extant_B = mapped_fraction_for_regime(extant_simmap, "B"),
    base_signal = sc$base_signal,
    relation = sc$relation,
    lambda_scale = lambda_scale,
    lambda_corr_power = lambda_corr_power,
    true_scale_B = true_par$scale_b,
    true_corr_power_B = true_par$corr_power_b,
    full_cov_rel_frob_A = full_cov$err_A,
    full_cov_rel_frob_B = full_cov$err_B,
    full_cov_rel_frob = full_cov$err_mean,
    extant_cov_rel_frob_A = extant_cov$err_A,
    extant_cov_rel_frob_B = extant_cov$err_B,
    extant_cov_rel_frob = extant_cov$err_mean,
    delta_cov_rel_frob_B = extant_cov$err_B - full_cov$err_B,
    full_summary_rmse = full_summary$rmse,
    extant_summary_rmse = extant_summary$rmse,
    delta_summary_rmse = extant_summary$rmse - full_summary$rmse,
    full_B_mean_rate_abs_err = full_summary$B_mean_rate_abs_err,
    extant_B_mean_rate_abs_err = extant_summary$B_mean_rate_abs_err,
    full_B_mean_abs_corr_abs_err = full_summary$B_mean_abs_correlation_abs_err,
    extant_B_mean_abs_corr_abs_err = extant_summary$B_mean_abs_correlation_abs_err,
    fossil_better_cov_B = full_cov$err_B < extant_cov$err_B,
    fossil_better_summary = full_summary$rmse < extant_summary$rmse,
    logLik_full = as.numeric(fit_full$logLik),
    logLik_extant = as.numeric(fit_extant$logLik),
    full_pathological_scale = isTRUE(getElement(fit_full$diagnostics$corrpower, "pathological_scale")),
    extant_pathological_scale = isTRUE(getElement(fit_extant$diagnostics$corrpower, "pathological_scale")),
    full_boundary = isTRUE(getElement(fit_full$diagnostics$corrpower, "boundary_corr_power")),
    extant_boundary = isTRUE(getElement(fit_extant$diagnostics$corrpower, "boundary_corr_power")),
    stringsAsFactors = FALSE
  )
}

reps <- env_int("CORRPOWER_FOSSIL_RECOVERY_REPS", 20L)
if (is.na(reps) || reps < 1L) stop("CORRPOWER_FOSSIL_RECOVERY_REPS must be >= 1", call. = FALSE)

base_grid <- expand.grid(
  n_extant = c(50L, 100L, 200L),
  p = c(2L, 4L),
  fossil_fraction = c(0, 0.25, 0.50, 0.75),
  mapped_fraction = c(0.20, 0.80),
  base_signal = c(0.20, 0.50),
  relation = c("weaker", "stronger"),
  replicate = seq_len(reps),
  stringsAsFactors = FALSE
)

expanded_grid <- do.call(rbind, lapply(seq_len(nrow(base_grid)), function(i) {
  sc <- base_grid[i, , drop = FALSE]
  depth_values <- if (sc$fossil_fraction[[1]] <= 0) "none" else c("deep", "shallow")
  cbind(sc[rep(1L, length(depth_values)), , drop = FALSE], fossil_depth = depth_values, stringsAsFactors = FALSE)
}))
expanded_grid$scenario_id <- seq_len(nrow(expanded_grid))

scenario_grid <- if (env_flag("CORRPOWER_FOSSIL_RECOVERY_FULL", FALSE)) {
  expanded_grid
} else {
  subset(
    expanded_grid,
    n_extant == 50L &
      p %in% c(2L, 4L) &
      fossil_fraction %in% c(0, 0.50) &
      fossil_depth %in% c("none", "deep") &
      mapped_fraction == 0.20 &
      base_signal == 0.20 &
      relation == "weaker" &
      replicate == 1L
  )
}

chunk_total <- env_int("CORRPOWER_FOSSIL_RECOVERY_CHUNK_TOTAL", 1L)
chunk_index <- env_int("CORRPOWER_FOSSIL_RECOVERY_CHUNK_INDEX", 1L)
if (is.na(chunk_total) || chunk_total < 1L) stop("CORRPOWER_FOSSIL_RECOVERY_CHUNK_TOTAL must be >= 1", call. = FALSE)
if (is.na(chunk_index) || chunk_index < 1L || chunk_index > chunk_total) {
  stop("CORRPOWER_FOSSIL_RECOVERY_CHUNK_INDEX must be between 1 and CORRPOWER_FOSSIL_RECOVERY_CHUNK_TOTAL", call. = FALSE)
}
if (chunk_total > 1L) {
  keep <- ((seq_len(nrow(scenario_grid)) - 1L) %% chunk_total) == (chunk_index - 1L)
  scenario_grid <- scenario_grid[keep, , drop = FALSE]
}

empty_results <- function() {
  data.frame(
    scenario_id = integer(0),
    replicate = integer(0),
    n_extant = integer(0),
    n_fossil = integer(0),
    p = integer(0),
    fossil_fraction = numeric(0),
    fossil_depth = character(0),
    mapped_fraction_target = numeric(0),
    mapped_fraction_full_B = numeric(0),
    mapped_fraction_extant_B = numeric(0),
    base_signal = numeric(0),
    relation = character(0),
    lambda_scale = numeric(0),
    lambda_corr_power = numeric(0),
    true_scale_B = numeric(0),
    true_corr_power_B = numeric(0),
    full_cov_rel_frob_A = numeric(0),
    full_cov_rel_frob_B = numeric(0),
    full_cov_rel_frob = numeric(0),
    extant_cov_rel_frob_A = numeric(0),
    extant_cov_rel_frob_B = numeric(0),
    extant_cov_rel_frob = numeric(0),
    delta_cov_rel_frob_B = numeric(0),
    full_summary_rmse = numeric(0),
    extant_summary_rmse = numeric(0),
    delta_summary_rmse = numeric(0),
    full_B_mean_rate_abs_err = numeric(0),
    extant_B_mean_rate_abs_err = numeric(0),
    full_B_mean_abs_corr_abs_err = numeric(0),
    extant_B_mean_abs_corr_abs_err = numeric(0),
    fossil_better_cov_B = logical(0),
    fossil_better_summary = logical(0),
    logLik_full = numeric(0),
    logLik_extant = numeric(0),
    full_pathological_scale = logical(0),
    extant_pathological_scale = logical(0),
    full_boundary = logical(0),
    extant_boundary = logical(0),
    stringsAsFactors = FALSE
  )
}

if (!nrow(scenario_grid)) {
  results_df <- empty_results()
  write_outputs(results_df)
  cat("No scenarios assigned to this corrpower fossil-recovery chunk\n")
  quit(save = "no", status = 0L)
}

results_df <- do.call(rbind, lapply(seq_len(nrow(scenario_grid)), function(i) {
  compare_case(scenario_grid[i, , drop = FALSE])
}))

write_outputs(results_df)
print(results_df)

cat("\nFossil recovery summary\n")
cat(sprintf("rows: %d\n", nrow(results_df)))
cat(sprintf("fossil better on regime-B covariance: %.3f\n", mean(results_df$fossil_better_cov_B, na.rm = TRUE)))
cat(sprintf("fossil better on regime summary: %.3f\n", mean(results_df$fossil_better_summary, na.rm = TRUE)))
cat(sprintf("median delta cov rel frob B (extant - fossil): %.4f\n", stats::median(results_df$delta_cov_rel_frob_B, na.rm = TRUE)))
cat(sprintf("median delta summary rmse (extant - fossil): %.4f\n", stats::median(results_df$delta_summary_rmse, na.rm = TRUE)))
cat(sprintf("full pathological rate: %.3f\n", mean(results_df$full_pathological_scale, na.rm = TRUE)))
cat(sprintf("extant pathological rate: %.3f\n", mean(results_df$extant_pathological_scale, na.rm = TRUE)))
cat(sprintf("full boundary rate: %.3f\n", mean(results_df$full_boundary, na.rm = TRUE)))
cat(sprintf("extant boundary rate: %.3f\n", mean(results_df$extant_boundary, na.rm = TRUE)))

if (!any(is.finite(results_df$full_cov_rel_frob_B))) {
  stop("No corrpower fossil-recovery fits succeeded.", call. = FALSE)
}

cat("corrpower fossil-recovery grid harness checks passed\n")
