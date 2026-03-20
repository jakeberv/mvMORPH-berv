#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
if (!length(file_arg)) stop("This script must be run with Rscript.", call. = FALSE)
script_path <- normalizePath(file_arg)
repo_root <- normalizePath(file.path(dirname(script_path), ".."))

desc_packages <- function(value) {
  if (is.null(value) || is.na(value) || !nzchar(value)) return(character())
  packages <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  packages <- sub("\\s*\\(.*\\)$", "", packages)
  packages <- trimws(packages)
  packages[nzchar(packages) & packages != "R"]
}

bootstrap_runtime <- function(repo_root) {
  desc_path <- file.path(repo_root, "DESCRIPTION")
  if (!file.exists(desc_path)) {
    stop("Could not find DESCRIPTION at ", desc_path, call. = FALSE)
  }

  desc <- read.dcf(desc_path)
  base_pkgs <- rownames(installed.packages(priority = c("base", "recommended")))
  required_pkgs <- unique(c(
    "pkgload",
    desc_packages(if ("Depends" %in% colnames(desc)) desc[1, "Depends"] else ""),
    desc_packages(if ("Imports" %in% colnames(desc)) desc[1, "Imports"] else "")
  ))
  required_pkgs <- setdiff(required_pkgs, base_pkgs)

  install_lib <- Sys.getenv("OUM_GRID_R_LIB", unset = "")
  if (!nzchar(install_lib)) {
    install_lib <- Sys.getenv("R_LIBS_USER", unset = "")
  }
  if (!nzchar(install_lib)) {
    install_lib <- .libPaths()[1]
  }
  dir.create(install_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(unique(c(install_lib, .libPaths())))

  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    message(
      "Installing missing R packages into ",
      install_lib,
      ": ",
      paste(missing_pkgs, collapse = ", ")
    )
    utils::install.packages(
      missing_pkgs,
      repos = Sys.getenv("OUM_GRID_CRAN_REPO", unset = "https://cloud.r-project.org"),
      lib = install_lib,
      quiet = TRUE
    )
  }

  still_missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing)) {
    stop(
      "Missing required packages after attempted install: ",
      paste(still_missing, collapse = ", "),
      call. = FALSE
    )
  }
}

# Remove staged local build products so pkgload compiles a host-native shared object.
clean_staged_artifacts <- function(repo_root) {
  src_dir <- file.path(repo_root, "src")
  if (dir.exists(src_dir)) {
    stale_files <- list.files(
      src_dir,
      pattern = "\\.(o|so|dll|dylib)$",
      full.names = TRUE
    )
    if (length(stale_files)) {
      unlink(stale_files, force = TRUE)
    }
  }
  if (requireNamespace("pkgbuild", quietly = TRUE)) {
    pkgbuild::clean_dll(repo_root)
  }
}

# Bootstrap runtime dependencies, then force a clean host-native package build.
bootstrap_runtime(repo_root)
clean_staged_artifacts(repo_root)
pkgload::load_all(repo_root, quiet = TRUE, compile = TRUE, recompile = TRUE)

`%||%` <- function(x, y) if (is.null(x)) y else x

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed <= 0L) return(default)
  parsed
}

env_num_subset <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1]]))
  parsed <- parsed[is.finite(parsed)]
  if (!length(parsed)) return(default)
  parsed
}

env_chr_subset <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  parsed <- parsed[nzchar(parsed)]
  if (!length(parsed)) return(default)
  parsed
}

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

safe_mean <- function(x) {
  if (!length(x) || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (sum(is.finite(x)) <= 1L) return(NA_real_)
  sd(x, na.rm = TRUE)
}

safe_rate <- function(x) {
  if (!length(x) || all(is.na(x))) return(NA_real_)
  mean(as.numeric(x), na.rm = TRUE)
}

has_text <- function(x) {
  !is.na(x) & nzchar(x)
}

signal_scale <- function(level) {
  switch(
    level,
    null = 0,
    weak = 0.10,
    medium = 0.22,
    strong = 0.35,
    stop("Unknown signal level", call. = FALSE)
  )
}

make_simmap <- function(n_tips, n_regimes, seed) {
  set.seed(seed)
  tree <- phytools::pbtree(n = n_tips, scale = 1)
  regime_names <- LETTERS[seq_len(n_regimes)]
  counts <- rep(floor(n_tips / n_regimes), n_regimes)
  counts[seq_len(n_tips %% n_regimes)] <- counts[seq_len(n_tips %% n_regimes)] + 1L
  states <- setNames(sample(rep(regime_names, counts)), tree$tip.label)
  invisible(capture.output({
    simmap <- suppressMessages(phytools::make.simmap(tree, states, model = "ER", nsim = 1))
  }))
  simmap
}

make_predictor <- function(tips, type, seed) {
  set.seed(seed)
  n <- length(tips)
  if (identical(type, "numeric")) {
    x <- as.numeric(scale(stats::rnorm(n)))
    names(x) <- tips
    frame <- data.frame(x = x, row.names = tips)
    X <- stats::model.matrix(~ x, data = frame)
    list(
      full_formula = Y ~ x,
      null_formula = Y ~ 1,
      data = as.list(frame),
      training_frame = frame,
      node_frame = function(nnode) data.frame(x = seq(-1, 1, length.out = nnode)),
      X = X
    )
  } else if (identical(type, "factor3")) {
    grp <- factor(sample(rep(c("g1", "g2", "g3"), length.out = n)), levels = c("g1", "g2", "g3"))
    names(grp) <- tips
    frame <- data.frame(grp = grp, row.names = tips)
    X <- stats::model.matrix(~ grp, data = frame)
    list(
      full_formula = Y ~ grp,
      null_formula = Y ~ 1,
      data = as.list(frame),
      training_frame = frame,
      node_frame = function(nnode) data.frame(grp = factor(rep("g1", nnode), levels = levels(grp))),
      X = X
    )
  } else {
    stop("Unknown predictor type", call. = FALSE)
  }
}

build_theta <- function(regimes, p) {
  offsets <- seq(-(length(regimes) - 1) / 2, (length(regimes) - 1) / 2, length.out = length(regimes))
  loadings <- seq(0.30, 0.08, length.out = p)
  theta <- outer(offsets, loadings, "*")
  theta <- theta + outer(seq_along(regimes) - 1, (-1) ^ seq_len(p), "*") * 0.03
  rownames(theta) <- regimes
  colnames(theta) <- paste0("trait", seq_len(p))
  theta
}

build_beta <- function(term_names, p, signal_level) {
  scale <- signal_scale(signal_level)
  if (length(term_names) == 0L) {
    return(matrix(numeric(0), nrow = 0, ncol = p, dimnames = list(character(0), paste0("trait", seq_len(p)))))
  }
  base <- seq(0.45, 0.15, length.out = p)
  beta <- vapply(seq_along(term_names), function(i) {
    sign_pattern <- (-1) ^ (seq_len(p) + i)
    scale * sign_pattern * base * (1 + 0.25 * (i - 1L))
  }, numeric(p))
  beta <- t(beta)
  rownames(beta) <- term_names
  colnames(beta) <- paste0("trait", seq_len(p))
  beta
}

simulate_oum_response <- function(tree, alpha, Sigma, theta, X_full, beta, seed) {
  set.seed(seed)
  tips <- tree$tip.label
  p <- ncol(theta)
  W <- mvMORPH:::.mvgls_oum_weight_matrix(tree, param = alpha, root = "stationary", std = 1L)
  W <- W[tips, , drop = FALSE]

  mean_mat <- W %*% theta
  if (nrow(beta) > 0L) {
    mean_mat <- mean_mat + X_full[, rownames(beta), drop = FALSE] %*% beta
  }
  rownames(mean_mat) <- tips

  Vphy <- .Call(
    "mvmorph_covar_ou_fixed",
    A = ape::vcv.phylo(tree),
    alpha = as.double(alpha),
    sigma = 1,
    PACKAGE = "mvMORPH"
  )
  chol_full <- t(chol(kronecker(Sigma, Vphy)))
  noise <- matrix(chol_full %*% stats::rnorm(length(tips) * p), nrow = length(tips), ncol = p)
  rownames(noise) <- tips
  colnames(noise) <- colnames(theta)
  mean_mat + noise
}

fit_mvgls <- function(formula, data, tree, method) {
  common_args <- list(
    formula = formula,
    data = data,
    tree = tree,
    model = "OUM",
    method = method,
    REML = FALSE,
    echo = FALSE
  )
  if (identical(method, "LOOCV")) common_args$penalty <- "RidgeArch"
  tryCatch(
    suppressWarnings(do.call(mvgls, common_args)),
    error = function(e) e
  )
}

evaluate_method <- function(full_fit, null_fit, beta_true, predictor_frame, tree) {
  if (inherits(full_fit, "error") || inherits(null_fit, "error")) {
    return(list(
      success = FALSE,
      delta_gic = NA_real_,
      coef_present = NA_real_,
      beta_rmse = NA_real_,
      predict_ok = FALSE,
      ancestral_ok = FALSE,
      manova_ok = FALSE,
      full_error = if (inherits(full_fit, "error")) conditionMessage(full_fit) else NA_character_,
      null_error = if (inherits(null_fit, "error")) conditionMessage(null_fit) else NA_character_
    ))
  }

  delta_gic <- GIC(null_fit)$GIC - GIC(full_fit)$GIC
  coef_est <- coef(full_fit)
  coef_present <- all(rownames(beta_true) %in% rownames(coef_est))
  beta_rmse <- if (!coef_present || nrow(beta_true) == 0L) {
    0
  } else {
    sqrt(mean((coef_est[rownames(beta_true), , drop = FALSE] - beta_true)^2))
  }

  predict_ok <- FALSE
  ancestral_ok <- FALSE
  manova_ok <- FALSE

  pred <- tryCatch(predict(full_fit, newdata = predictor_frame), error = function(e) e)
  if (!inherits(pred, "error")) {
    predict_ok <- identical(dim(pred), dim(full_fit$variables$Y))
  }

  node_data <- predictor_frame[rep(1L, ape::Nnode(tree)), , drop = FALSE]
  rownames(node_data) <- NULL
  anc <- tryCatch(ancestral(full_fit, newdata = node_data), error = function(e) e)
  if (!inherits(anc, "error")) {
    ancestral_ok <- identical(dim(anc), c(ape::Nnode(tree), full_fit$dims$p))
  }

  mv_test <- tryCatch(manova.gls(full_fit, type = "II", test = "Pillai"), error = function(e) e)
  if (!inherits(mv_test, "error")) {
    manova_ok <- length(mv_test$terms) >= 1L
  }

  list(
    success = TRUE,
    delta_gic = delta_gic,
    coef_present = as.numeric(coef_present),
    beta_rmse = beta_rmse,
    predict_ok = predict_ok,
    ancestral_ok = ancestral_ok,
    manova_ok = manova_ok,
    full_error = NA_character_,
    null_error = NA_character_
  )
}

cores <- env_int("OUM_GRID_CORES", max(1L, parallel::detectCores(logical = FALSE) %||% 1L))
reps <- env_int("OUM_GRID_REPS", 3L)
seed_base <- env_int("OUM_GRID_SEED", 20260319L)
scenario_grid <- expand.grid(
  n_tips = env_num_subset("OUM_GRID_NTIPS", c(30L, 60L, 120L)),
  p = env_num_subset("OUM_GRID_P", c(4L, 8L, 12L)),
  n_regimes = env_num_subset("OUM_GRID_REGIMES", c(2L, 3L)),
  predictor = env_chr_subset("OUM_GRID_PREDICTOR", c("numeric", "factor3")),
  signal = env_chr_subset("OUM_GRID_SIGNAL", c("null", "weak", "medium", "strong")),
  stringsAsFactors = FALSE
)

scenario_grid$scenario_id <- sprintf(
  "n%d_p%d_r%d_%s_%s",
  scenario_grid$n_tips,
  scenario_grid$p,
  scenario_grid$n_regimes,
  scenario_grid$predictor,
  scenario_grid$signal
)

methods <- c("LL", "LOOCV", "EmpBayes")
tasks <- do.call(
  rbind,
  lapply(seq_len(nrow(scenario_grid)), function(i) {
    data.frame(
      scenario_index = i,
      rep = seq_len(reps),
      stringsAsFactors = FALSE
    )
  })
)

run_replicate <- function(scenario, rep_id) {
  seed_offset <- seed_base + 1000L * rep_id + 10000L * match(scenario$scenario_id, scenario_grid$scenario_id)
  tree <- make_simmap(scenario$n_tips, scenario$n_regimes, seed = seed_offset)
  predictor <- make_predictor(tree$tip.label, scenario$predictor, seed = seed_offset + 1L)
  p <- scenario$p
  Sigma <- 0.04 * stats::toeplitz(0.65 ^ (0:(p - 1L)))
  theta <- build_theta(colnames(mvMORPH:::.mvgls_oum_weight_matrix(tree, param = 0.8, root = "stationary", std = 1L)), p)
  beta_true <- build_beta(colnames(predictor$X)[colnames(predictor$X) != "(Intercept)"], p, scenario$signal)
  Y <- simulate_oum_response(tree, alpha = 0.8, Sigma = Sigma, theta = theta, X_full = predictor$X, beta = beta_true, seed = seed_offset + 2L)

  data_full <- c(list(Y = Y), predictor$data)
  data_null <- list(Y = Y)

  rows <- lapply(methods, function(method) {
    full_fit <- fit_mvgls(predictor$full_formula, data_full, tree, method)
    null_fit <- fit_mvgls(predictor$null_formula, data_null, tree, method)
    metrics <- evaluate_method(full_fit, null_fit, beta_true, predictor$training_frame, tree)
    data.frame(
      scenario_id = scenario$scenario_id,
      n_tips = scenario$n_tips,
      p = scenario$p,
      n_regimes = scenario$n_regimes,
      predictor = scenario$predictor,
      signal = scenario$signal,
      rep = rep_id,
      method = method,
      success = metrics$success,
      delta_gic = metrics$delta_gic,
      select_full = ifelse(is.na(metrics$delta_gic), NA, metrics$delta_gic > 0),
      coef_present = metrics$coef_present,
      beta_rmse = metrics$beta_rmse,
      predict_ok = metrics$predict_ok,
      ancestral_ok = metrics$ancestral_ok,
      manova_ok = metrics$manova_ok,
      full_error = metrics$full_error,
      null_error = metrics$null_error,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

cat(
  "Running OUM covariate grid with",
  nrow(scenario_grid), "scenario combinations x", reps, "replicate(s)",
  "across", nrow(tasks), "task(s) on", cores, "worker(s)\n"
)

task_runner <- function(task_row) {
  scenario <- scenario_grid[task_row$scenario_index, , drop = FALSE]
  cat(
    sprintf(
      "[worker %s] scenario %d/%d rep %d: %s\n",
      Sys.getpid(),
      task_row$scenario_index,
      nrow(scenario_grid),
      task_row$rep,
      scenario$scenario_id[[1]]
    )
  )
  run_replicate(scenario, task_row$rep)
}

task_rows <- split(tasks, seq_len(nrow(tasks)))
if (cores > 1L) {
  result_list <- parallel::mclapply(task_rows, task_runner, mc.cores = cores, mc.preschedule = FALSE)
} else {
  result_list <- lapply(task_rows, task_runner)
}
results <- do.call(rbind, result_list)

summary_table <- do.call(
  rbind,
  lapply(split(results, list(results$scenario_id, results$method), drop = TRUE), function(df) {
    data.frame(
      scenario_id = df$scenario_id[[1]],
      method = df$method[[1]],
      n_tips = df$n_tips[[1]],
      p = df$p[[1]],
      n_regimes = df$n_regimes[[1]],
      predictor = df$predictor[[1]],
      signal = df$signal[[1]],
      reps = nrow(df),
      success_rate = safe_mean(df$success),
      select_full_rate = safe_mean(df$select_full),
      mean_delta_gic = safe_mean(df$delta_gic),
      sd_delta_gic = safe_sd(df$delta_gic),
      coef_present_rate = safe_mean(df$coef_present),
      beta_rmse_mean = safe_mean(df$beta_rmse),
      predict_rate = safe_mean(df$predict_ok),
      ancestral_rate = safe_mean(df$ancestral_ok),
      manova_rate = safe_mean(df$manova_ok),
      stringsAsFactors = FALSE
    )
  })
)

selection_summary <- do.call(
  rbind,
  lapply(split(results, list(results$signal, results$method, results$predictor), drop = TRUE), function(df) {
    data.frame(
      signal = df$signal[[1]],
      method = df$method[[1]],
      predictor = df$predictor[[1]],
      scenarios = length(unique(df$scenario_id)),
      reps = nrow(df),
      success_rate = safe_mean(df$success),
      select_full_rate = safe_mean(df$select_full),
      mean_delta_gic = safe_mean(df$delta_gic),
      coef_present_rate = safe_mean(df$coef_present),
      stringsAsFactors = FALSE
    )
  })
)

agreement_summary <- do.call(
  rbind,
  lapply(split(results, list(results$scenario_id, results$rep), drop = TRUE), function(df) {
    picks <- df$select_full[!is.na(df$select_full)]
    data.frame(
      scenario_id = df$scenario_id[[1]],
      rep = df$rep[[1]],
      signal = df$signal[[1]],
      predictor = df$predictor[[1]],
      all_methods_agree = length(unique(picks)) <= 1L,
      stringsAsFactors = FALSE
    )
  })
)

agreement_by_group <- do.call(
  rbind,
  lapply(split(agreement_summary, list(agreement_summary$signal, agreement_summary$predictor), drop = TRUE), function(df) {
    data.frame(
      signal = df$signal[[1]],
      predictor = df$predictor[[1]],
      scenario_reps = nrow(df),
      all_methods_agree_rate = safe_rate(df$all_methods_agree),
      stringsAsFactors = FALSE
    )
  })
)

ambiguous_summary <- summary_table
ambiguous_summary$abs_mean_delta_gic <- abs(ambiguous_summary$mean_delta_gic)
ambiguous_summary <- ambiguous_summary[order(ambiguous_summary$abs_mean_delta_gic, ambiguous_summary$beta_rmse_mean), ]

failure_rows <- results[!results$success | has_text(results$full_error) | has_text(results$null_error), , drop = FALSE]

summary_table <- summary_table[order(summary_table$signal, summary_table$predictor, summary_table$method, summary_table$n_tips, summary_table$p, summary_table$n_regimes), ]
selection_summary <- selection_summary[order(selection_summary$signal, selection_summary$predictor, selection_summary$method), ]
agreement_by_group <- agreement_by_group[order(agreement_by_group$signal, agreement_by_group$predictor), ]

cat("\nScenario-by-method summary\n")
print(summary_table, row.names = FALSE)

cat("\nCollapsed selection summary\n")
print(selection_summary, row.names = FALSE)

cat("\nMethod-agreement summary\n")
print(agreement_by_group, row.names = FALSE)

cat("\nMost ambiguous scenario-method cells by |mean delta GIC|\n")
print(utils::head(ambiguous_summary[, c(
  "scenario_id", "method", "signal", "predictor", "success_rate",
  "select_full_rate", "mean_delta_gic", "sd_delta_gic", "beta_rmse_mean"
)], 12L), row.names = FALSE)

if (nrow(failure_rows) > 0L) {
  cat("\nFailure rows\n")
  print(failure_rows[, c("scenario_id", "rep", "method", "full_error", "null_error")], row.names = FALSE)
} else {
  cat("\nNo failed fits were recorded.\n")
}

output_csv <- Sys.getenv("OUM_GRID_OUT", unset = "")
if (nzchar(output_csv)) {
  utils::write.csv(results, output_csv, row.names = FALSE)
  cat("\nSaved replicate-level results to", output_csv, "\n")
}

summary_csv <- Sys.getenv("OUM_GRID_SUMMARY_OUT", unset = "")
if (nzchar(summary_csv)) {
  utils::write.csv(summary_table, summary_csv, row.names = FALSE)
  cat("Saved scenario summary to", summary_csv, "\n")
}

selection_csv <- Sys.getenv("OUM_GRID_SELECTION_OUT", unset = "")
if (nzchar(selection_csv)) {
  utils::write.csv(selection_summary, selection_csv, row.names = FALSE)
  cat("Saved collapsed selection summary to", selection_csv, "\n")
}

agreement_csv <- Sys.getenv("OUM_GRID_AGREEMENT_OUT", unset = "")
if (nzchar(agreement_csv)) {
  utils::write.csv(agreement_by_group, agreement_csv, row.names = FALSE)
  cat("Saved method-agreement summary to", agreement_csv, "\n")
}
