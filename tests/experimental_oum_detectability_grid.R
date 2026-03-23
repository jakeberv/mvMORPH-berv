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

  install_lib <- Sys.getenv("OUM_DET_R_LIB", unset = "")
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
      repos = Sys.getenv("OUM_DET_CRAN_REPO", unset = "https://cloud.r-project.org"),
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

clean_staged_artifacts <- function(repo_root) {
  src_dir <- file.path(repo_root, "src")
  if (dir.exists(src_dir)) {
    stale_files <- list.files(src_dir, pattern = "\\.(o|so|dll|dylib)$", full.names = TRUE)
    if (length(stale_files)) {
      unlink(stale_files, force = TRUE)
    }
  }
  if (requireNamespace("pkgbuild", quietly = TRUE)) {
    pkgbuild::clean_dll(repo_root)
  }
}

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

safe_mean <- function(x) {
  if (!length(x) || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_rate <- function(x) {
  if (!length(x) || all(is.na(x))) return(NA_real_)
  mean(as.numeric(x), na.rm = TRUE)
}

has_text <- function(x) {
  !is.na(x) & nzchar(x)
}

format_level_id <- function(x, scale = 1000) {
  sprintf("%04d", as.integer(round(scale * x)))
}

collapse_messages <- function(x) {
  x <- trimws(x)
  x <- unique(x[nzchar(x)])
  if (!length(x)) return(NA_character_)
  paste(x, collapse = " | ")
}

normalize_counts <- function(n, weights) {
  weights <- weights / sum(weights)
  raw_counts <- weights * n
  counts <- floor(raw_counts)
  remainder <- n - sum(counts)
  if (remainder > 0L) {
    frac_order <- order(raw_counts - counts, decreasing = TRUE)
    counts[frac_order[seq_len(remainder)]] <- counts[frac_order[seq_len(remainder)]] + 1L
  }
  counts
}

make_regime_counts <- function(n_tips, n_regimes, regime_balance) {
  if (identical(regime_balance, "balanced")) {
    return(normalize_counts(n_tips, rep(1, n_regimes)))
  }

  if (identical(regime_balance, "imbalanced")) {
    if (n_regimes == 2L) {
      return(normalize_counts(n_tips, c(0.8, 0.2)))
    }
    if (n_regimes == 3L) {
      return(normalize_counts(n_tips, c(0.7, 0.2, 0.1)))
    }
    return(normalize_counts(n_tips, rev(seq_len(n_regimes))))
  }

  stop("Unknown regime_balance value: ", regime_balance, call. = FALSE)
}

make_simmap <- function(n_tips, n_regimes, regime_balance, seed) {
  set.seed(seed)
  tree <- phytools::pbtree(n = n_tips, scale = 1)
  regime_names <- LETTERS[seq_len(n_regimes)]
  counts <- make_regime_counts(n_tips, n_regimes, regime_balance)
  states <- setNames(sample(rep(regime_names, counts)), tree$tip.label)
  invisible(capture.output({
    simmap <- suppressMessages(phytools::make.simmap(tree, states, model = "ER", nsim = 1))
  }))
  list(tree = simmap, tip_states = states, regime_names = regime_names)
}

make_predictor <- function(tips, type, tip_states, seed) {
  set.seed(seed)
  n <- length(tips)
  regime_factor <- factor(tip_states[tips], levels = unique(tip_states))
  regime_index <- as.numeric(regime_factor)

  if (identical(type, "cont_indep")) {
    x <- as.numeric(scale(stats::rnorm(n)))
    names(x) <- tips
    frame <- data.frame(x = x, row.names = tips)
    X <- stats::model.matrix(~ x, data = frame)
    return(list(formula = Y ~ x, data = as.list(frame), training_frame = frame, X = X))
  }

  if (identical(type, "cont_confounded")) {
    x <- regime_index + stats::rnorm(n, sd = 0.35)
    x <- as.numeric(scale(x))
    names(x) <- tips
    frame <- data.frame(x = x, row.names = tips)
    X <- stats::model.matrix(~ x, data = frame)
    return(list(formula = Y ~ x, data = as.list(frame), training_frame = frame, X = X))
  }

  if (identical(type, "disc_binary")) {
    grp <- factor(sample(rep(c("g1", "g2"), length.out = n)), levels = c("g1", "g2"))
    names(grp) <- tips
    frame <- data.frame(grp = grp, row.names = tips)
    X <- stats::model.matrix(~ grp, data = frame)
    return(list(formula = Y ~ grp, data = as.list(frame), training_frame = frame, X = X))
  }

  if (identical(type, "disc3_sparse_confounded")) {
    level_names <- c("g1", "g2", "g3")
    probs_by_regime <- lapply(seq_len(length(levels(regime_factor))), function(i) {
      if (i == 1L) return(c(0.72, 0.23, 0.05))
      if (i == 2L) return(c(0.18, 0.72, 0.10))
      c(0.08, 0.34, 0.58)
    })
    grp <- character(n)
    for (attempt in seq_len(100L)) {
      grp <- vapply(seq_len(n), function(i) {
        sample(level_names, size = 1L, prob = probs_by_regime[[regime_index[i]]])
      }, character(1))
      if (length(unique(grp)) == length(level_names)) break
    }
    grp <- factor(grp, levels = level_names)
    names(grp) <- tips
    frame <- data.frame(grp = grp, row.names = tips)
    X <- stats::model.matrix(~ grp, data = frame)
    return(list(formula = Y ~ grp, data = as.list(frame), training_frame = frame, X = X))
  }

  stop("Unknown predictor type: ", type, call. = FALSE)
}

build_sigma_base <- function(p) {
  0.04 * stats::toeplitz(0.65 ^ (0:(p - 1L)))
}

build_rate_scales <- function(n_regimes, rate_ratio) {
  if (!is.finite(rate_ratio) || rate_ratio < 1) {
    stop("rate_ratio must be a finite value >= 1", call. = FALSE)
  }
  if (n_regimes == 1L || identical(rate_ratio, 1)) {
    return(rep(1, n_regimes))
  }
  exp(seq(log(1), log(rate_ratio), length.out = n_regimes))
}

build_bmm_sigmas <- function(base_sigma, regimes, rate_ratio) {
  scales <- build_rate_scales(length(regimes), rate_ratio)
  sigmas <- lapply(scales, function(scale) scale * base_sigma)
  names(sigmas) <- regimes
  sigmas
}

build_theta <- function(regimes, p, theta_sep) {
  theta <- matrix(0, nrow = length(regimes), ncol = p)
  rownames(theta) <- regimes
  colnames(theta) <- paste0("trait", seq_len(p))

  if (!is.finite(theta_sep) || theta_sep < 0) {
    stop("theta_sep must be a non-negative finite number", call. = FALSE)
  }
  if (identical(theta_sep, 0)) return(theta)

  offsets <- seq(-0.5, 0.5, length.out = length(regimes))
  loadings <- seq(1, 0.4, length.out = p) * ((-1) ^ seq_len(p))
  loadings <- loadings / sqrt(sum(loadings^2))
  theta[,] <- theta_sep * outer(offsets, loadings)
  theta
}

build_beta <- function(X_full, p, covariate_effect) {
  term_names <- colnames(X_full)[colnames(X_full) != "(Intercept)"]
  if (length(term_names) == 0L) {
    return(matrix(
      numeric(0),
      nrow = 0,
      ncol = p,
      dimnames = list(character(0), paste0("trait", seq_len(p)))
    ))
  }

  if (!is.finite(covariate_effect) || covariate_effect < 0) {
    stop("covariate_effect must be a non-negative finite number", call. = FALSE)
  }

  if (identical(covariate_effect, 0)) {
    beta <- matrix(0, nrow = length(term_names), ncol = p)
    rownames(beta) <- term_names
    colnames(beta) <- paste0("trait", seq_len(p))
    return(beta)
  }

  X_terms <- X_full[, term_names, drop = FALSE]
  term_scales <- apply(X_terms, 2, stats::sd)
  term_scales[!is.finite(term_scales) | term_scales <= 0] <- 1
  base <- seq(0.45, 0.15, length.out = p)
  beta <- vapply(seq_along(term_names), function(i) {
    sign_pattern <- (-1) ^ (seq_len(p) + i)
    (covariate_effect / term_scales[[i]]) * sign_pattern * base * (1 + 0.25 * (i - 1L))
  }, numeric(p))
  beta <- t(beta)
  rownames(beta) <- term_names
  colnames(beta) <- paste0("trait", seq_len(p))
  beta
}

get_oum_weight_matrix <- function(tree, alpha) {
  W <- mvMORPH:::.mvgls_oum_weight_matrix(tree, param = alpha, root = "stationary", std = 1L)
  W[tree$tip.label, , drop = FALSE]
}

simulate_gaussian_process <- function(Vphy, sigma, mean_mat, seed) {
  set.seed(seed)
  n <- nrow(Vphy)
  p <- ncol(sigma)
  chol_full <- t(chol(kronecker(sigma, Vphy)))
  noise <- matrix(chol_full %*% stats::rnorm(n * p), nrow = n, ncol = p)
  rownames(noise) <- rownames(Vphy)
  colnames(noise) <- colnames(mean_mat)
  mean_mat + noise
}

simulate_scalar_ou <- function(tree, alpha, sigma, mean_mat, seed) {
  Vphy <- .Call(
    "mvmorph_covar_ou_fixed",
    A = ape::vcv.phylo(tree),
    alpha = as.double(alpha),
    sigma = 1,
    PACKAGE = "mvMORPH"
  )
  simulate_gaussian_process(Vphy = Vphy, sigma = sigma, mean_mat = mean_mat, seed = seed)
}

simulate_truth <- function(truth, tree, alpha, base_sigma, theta_sep, rate_ratio, X_full, beta, seed) {
  p <- ncol(base_sigma)
  tip_order <- tree$tip.label
  trait_names <- paste0("trait", seq_len(p))
  zero_mean <- matrix(0, nrow = length(tip_order), ncol = p, dimnames = list(tip_order, trait_names))

  W <- get_oum_weight_matrix(tree, alpha = max(alpha, 1e-8))
  regimes <- colnames(W)
  theta <- build_theta(regimes, p, theta_sep)
  regime_signal <- W %*% theta
  rownames(regime_signal) <- tip_order
  colnames(regime_signal) <- trait_names

  if (truth == "shared_OU") {
    Y <- simulate_scalar_ou(tree = tree, alpha = alpha, sigma = base_sigma, mean_mat = zero_mean, seed = seed)
  } else if (truth == "pure_OUM") {
    Y <- simulate_scalar_ou(tree = tree, alpha = alpha, sigma = base_sigma, mean_mat = regime_signal, seed = seed)
  } else if (truth == "mixed") {
    bmm_sigmas <- build_bmm_sigmas(base_sigma, regimes, rate_ratio)
    Y <- mvSIM(
      tree,
      model = "BMM",
      nsim = 1,
      param = list(sigma = bmm_sigmas, theta = rep(0, p))
    )
    Y <- Y + regime_signal
  } else {
    stop("Unknown truth class: ", truth, call. = FALSE)
  }

  if (nrow(beta) > 0L) {
    Y <- Y + X_full[, rownames(beta), drop = FALSE] %*% beta
  }

  rownames(Y) <- tip_order
  colnames(Y) <- trait_names
  Y
}

fit_candidate <- function(formula, data, tree, model, method) {
  warnings <- character()
  fit <- withCallingHandlers(
    tryCatch(
      mvgls(
        formula = formula,
        data = data,
        tree = tree,
        model = model,
        method = method,
        REML = FALSE,
        echo = FALSE
      ),
      error = function(e) e
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = warnings)
}

compute_gic_safe <- function(fit_object) {
  warnings <- character()
  gic <- withCallingHandlers(
    tryCatch(GIC(fit_object), error = function(e) e),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (!grepl("GIC criterion with multiple predictors has not been fully tested", msg, fixed = TRUE)) {
        warnings <<- c(warnings, msg)
      }
      invokeRestart("muffleWarning")
    }
  )
  list(value = gic, warnings = warnings)
}

extract_beta_rmse <- function(fit_object, beta_true) {
  if (inherits(fit_object, "error")) return(NA_real_)
  if (!nrow(beta_true)) return(0)
  coef_est <- coef(fit_object)
  if (!all(rownames(beta_true) %in% rownames(coef_est))) return(NA_real_)
  sqrt(mean((coef_est[rownames(beta_true), , drop = FALSE] - beta_true)^2))
}

extract_alpha_hat <- function(fit_object) {
  if (inherits(fit_object, "error")) return(NA_real_)
  param <- fit_object$param %||% NULL
  if (is.null(param) || !length(param)) return(NA_real_)
  if (!is.null(names(param)) && "alpha" %in% names(param)) {
    return(as.numeric(param[["alpha"]]))
  }
  if (length(param) == 1L && is.finite(param[[1]])) {
    return(as.numeric(param[[1]]))
  }
  NA_real_
}

choose_best_model <- function(values, tol = 1e-8) {
  finite <- values[is.finite(values)]
  if (!length(finite)) return(list(winner = NA_character_, best = NA_real_))
  best <- min(finite)
  winners <- names(finite)[abs(finite - best) <= tol]
  list(winner = if (length(winners) == 1L) winners[[1]] else "tie", best = best)
}

expected_winner <- function(theta_sep) {
  if (!is.finite(theta_sep) || theta_sep <= 0) "OU" else "OUM"
}

write_checkpoint <- function(df, checkpoint_dir, scenario_id, rep_id) {
  if (!nzchar(checkpoint_dir)) return(invisible(NULL))
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(checkpoint_dir, sprintf("result_%s_rep%03d.csv", scenario_id, as.integer(rep_id)))
  utils::write.csv(df, out_path, row.names = FALSE)
  invisible(out_path)
}

summarize_results <- function(results) {
  summary_table <- do.call(
    rbind,
    lapply(split(results, results$scenario_id, drop = TRUE), function(df) {
      data.frame(
        phase = df$phase[[1]],
        truth = df$truth[[1]],
        scenario_id = df$scenario_id[[1]],
        method = df$method[[1]],
        n_tips = df$n_tips[[1]],
        p = df$p[[1]],
        n_regimes = df$n_regimes[[1]],
        regime_balance = df$regime_balance[[1]],
        predictor = df$predictor[[1]],
        covariate_effect = df$covariate_effect[[1]],
        alpha = df$alpha[[1]],
        theta_sep = df$theta_sep[[1]],
        rate_ratio = df$rate_ratio[[1]],
        reps = nrow(df),
        ou_fit_rate = safe_rate(df$ou_fit_success),
        oum_fit_rate = safe_rate(df$oum_fit_success),
        both_gic_rate = safe_rate(df$both_gic_available),
        select_ou_rate = safe_rate(df$winner == "OU"),
        select_oum_rate = safe_rate(df$winner == "OUM"),
        tie_rate = safe_rate(df$winner == "tie"),
        correct_rate = safe_rate(df$correct),
        mean_delta_gic_ou_minus_oum = safe_mean(df$delta_gic_ou_minus_oum),
        beta_rmse_ou_mean = safe_mean(df$beta_rmse_ou),
        beta_rmse_oum_mean = safe_mean(df$beta_rmse_oum),
        alpha_hat_ou_mean = safe_mean(df$alpha_hat_ou),
        alpha_hat_oum_mean = safe_mean(df$alpha_hat_oum),
        warning_rate = safe_rate(df$has_issue),
        stringsAsFactors = FALSE
      )
    })
  )

  selection_summary <- do.call(
    rbind,
    lapply(split(
      results,
      list(results$phase, results$truth, results$method, results$predictor, results$regime_balance),
      drop = TRUE
    ), function(df) {
      data.frame(
        phase = df$phase[[1]],
        truth = df$truth[[1]],
        method = df$method[[1]],
        predictor = df$predictor[[1]],
        regime_balance = df$regime_balance[[1]],
        scenario_reps = nrow(df),
        select_ou_rate = safe_rate(df$winner == "OU"),
        select_oum_rate = safe_rate(df$winner == "OUM"),
        correct_rate = safe_rate(df$correct),
        mean_delta_gic_ou_minus_oum = safe_mean(df$delta_gic_ou_minus_oum),
        warning_rate = safe_rate(df$has_issue),
        stringsAsFactors = FALSE
      )
    })
  )

  boundary_rows <- list()

  for (phase_name in c("oum_frontier", "oum_misspec")) {
    phase_df <- summary_table[summary_table$phase == phase_name, , drop = FALSE]
    if (!nrow(phase_df)) next

    split_args <- list(
      phase_df$method,
      phase_df$n_tips,
      phase_df$p,
      phase_df$n_regimes,
      phase_df$regime_balance,
      phase_df$predictor,
      phase_df$covariate_effect,
      phase_df$alpha
    )
    if (phase_name == "oum_misspec") {
      split_args <- c(split_args, list(phase_df$rate_ratio))
    }

    boundary_rows[[phase_name]] <- do.call(
      rbind,
      lapply(split(phase_df, split_args, drop = TRUE), function(df) {
        df <- df[order(df$theta_sep), , drop = FALSE]
        signal <- df$theta_sep
        select_rate <- df$select_oum_rate
        empirical50 <- if (any(select_rate >= 0.5, na.rm = TRUE)) min(signal[select_rate >= 0.5], na.rm = TRUE) else NA_real_
        empirical80 <- if (any(select_rate >= 0.8, na.rm = TRUE)) min(signal[select_rate >= 0.8], na.rm = TRUE) else NA_real_
        fit <- tryCatch(
          stats::glm(
            cbind(round(select_rate * df$reps), df$reps - round(select_rate * df$reps)) ~ signal,
            family = stats::binomial()
          ),
          error = function(e) e
        )
        slope <- NA_real_
        model50 <- NA_real_
        model80 <- NA_real_
        if (!inherits(fit, "error")) {
          coef_hat <- stats::coef(fit)
          slope <- unname(coef_hat[[2]])
          if (
            is.finite(slope) &&
            abs(slope) > 1e-6 &&
            length(unique(select_rate[is.finite(select_rate)])) > 1L
          ) {
            model50 <- unname(-coef_hat[[1]] / slope)
            model80 <- unname((stats::qlogis(0.8) - coef_hat[[1]]) / slope)
          }
        }
        data.frame(
          phase = phase_name,
          truth = df$truth[[1]],
          method = df$method[[1]],
          n_tips = df$n_tips[[1]],
          p = df$p[[1]],
          n_regimes = df$n_regimes[[1]],
          regime_balance = df$regime_balance[[1]],
          predictor = df$predictor[[1]],
          covariate_effect = df$covariate_effect[[1]],
          alpha = df$alpha[[1]],
          rate_ratio = df$rate_ratio[[1]],
          signal_parameter = "theta_sep",
          select_rate_at_null_signal = select_rate[match(0, signal)],
          empirical_signal50 = empirical50,
          empirical_signal80 = empirical80,
          logistic_slope = slope,
          model_signal50 = model50,
          model_signal80 = model80,
          stringsAsFactors = FALSE
        )
      })
    )
  }

  boundary_summary <- if (length(boundary_rows)) do.call(rbind, boundary_rows) else data.frame()
  issue_rows <- results[results$has_issue | !results$ou_fit_success | !results$oum_fit_success, , drop = FALSE]

  summary_table <- summary_table[order(
    summary_table$phase,
    summary_table$truth,
    summary_table$predictor,
    summary_table$regime_balance,
    summary_table$n_tips,
    summary_table$p,
    summary_table$n_regimes,
    summary_table$covariate_effect,
    summary_table$alpha,
    summary_table$rate_ratio,
    summary_table$theta_sep
  ), ]

  selection_summary <- selection_summary[order(
    selection_summary$phase,
    selection_summary$truth,
    selection_summary$predictor,
    selection_summary$regime_balance,
    selection_summary$method
  ), ]

  if (nrow(boundary_summary)) {
    boundary_summary <- boundary_summary[order(
      boundary_summary$phase,
      boundary_summary$predictor,
      boundary_summary$regime_balance,
      boundary_summary$n_tips,
      boundary_summary$p,
      boundary_summary$n_regimes,
      boundary_summary$covariate_effect,
      boundary_summary$alpha,
      boundary_summary$rate_ratio
    ), ]
  }

  list(
    summary_table = summary_table,
    selection_summary = selection_summary,
    boundary_summary = boundary_summary,
    issue_rows = issue_rows
  )
}

cores <- env_int("OUM_DET_CORES", max(1L, parallel::detectCores(logical = FALSE) %||% 1L))
seed_base <- env_int("OUM_DET_SEED", 20260323L)
checkpoint_dir <- Sys.getenv("OUM_DET_CHECKPOINT_DIR", unset = "")
phase_filter <- env_chr_subset("OUM_DET_PHASES", c("oum_null", "oum_frontier", "oum_misspec"))
methods <- env_chr_subset("OUM_DET_METHODS", c("LL"))
if (!all(methods %in% "LL")) {
  stop("This detectability study currently supports method='LL' only.", call. = FALSE)
}

regime_balances <- env_chr_subset("OUM_DET_BALANCE", c("balanced", "imbalanced"))
null_predictors <- env_chr_subset(
  "OUM_DET_NULL_PREDICTOR",
  c("cont_indep", "cont_confounded", "disc_binary", "disc3_sparse_confounded")
)
frontier_predictors <- env_chr_subset(
  "OUM_DET_FRONTIER_PREDICTOR",
  c("cont_confounded", "disc_binary", "disc3_sparse_confounded")
)
misspec_predictors <- env_chr_subset(
  "OUM_DET_MISSPEC_PREDICTOR",
  c("cont_confounded", "disc_binary", "disc3_sparse_confounded")
)

covariate_levels <- env_num_subset("OUM_DET_COVARIATE_EFFECT", c(0.00, 0.05))
regime_levels <- env_num_subset("OUM_DET_REGIMES", c(2L, 3L))

null_reps <- env_int("OUM_DET_NULL_REPS", 8L)
frontier_reps <- env_int("OUM_DET_FRONTIER_REPS", 8L)
misspec_reps <- env_int("OUM_DET_MISSPEC_REPS", 6L)

null_grid <- expand.grid(
  phase = "oum_null",
  truth = "shared_OU",
  n_tips = env_num_subset("OUM_DET_NULL_NTIPS", c(60L, 120L)),
  p = env_num_subset("OUM_DET_NULL_P", c(8L, 24L)),
  n_regimes = regime_levels,
  regime_balance = regime_balances,
  predictor = null_predictors,
  covariate_effect = covariate_levels,
  alpha = env_num_subset("OUM_DET_NULL_ALPHA", c(0.10, 0.25, 0.50, 1.00, 1.50)),
  theta_sep = 0,
  rate_ratio = 1,
  reps = null_reps,
  stringsAsFactors = FALSE
)

frontier_grid <- expand.grid(
  phase = "oum_frontier",
  truth = "pure_OUM",
  n_tips = env_num_subset("OUM_DET_FRONTIER_NTIPS", c(60L, 120L)),
  p = env_num_subset("OUM_DET_FRONTIER_P", c(24L)),
  n_regimes = regime_levels,
  regime_balance = regime_balances,
  predictor = frontier_predictors,
  covariate_effect = covariate_levels,
  alpha = env_num_subset("OUM_DET_FRONTIER_ALPHA", c(0.10, 0.25, 0.50, 1.00, 1.50)),
  theta_sep = env_num_subset("OUM_DET_FRONTIER_THETA", c(0.0000, 0.0025, 0.0050, 0.0100, 0.0200, 0.0400, 0.0800, 0.1200, 0.1800, 0.2600, 0.3600)),
  rate_ratio = 1,
  reps = frontier_reps,
  stringsAsFactors = FALSE
)

misspec_grid <- expand.grid(
  phase = "oum_misspec",
  truth = "mixed",
  n_tips = env_num_subset("OUM_DET_MISSPEC_NTIPS", c(120L)),
  p = env_num_subset("OUM_DET_MISSPEC_P", c(24L)),
  n_regimes = regime_levels,
  regime_balance = regime_balances,
  predictor = misspec_predictors,
  covariate_effect = covariate_levels,
  alpha = env_num_subset("OUM_DET_MISSPEC_ALPHA", c(0.25, 0.50, 1.00, 1.50)),
  theta_sep = env_num_subset("OUM_DET_MISSPEC_THETA", c(0.0000, 0.0100, 0.0200, 0.0400, 0.0800, 0.1200, 0.1800)),
  rate_ratio = env_num_subset("OUM_DET_MISSPEC_RATE", c(1.00, 1.05, 1.10, 1.20, 1.35, 1.50, 2.00)),
  reps = misspec_reps,
  stringsAsFactors = FALSE
)

scenario_grid <- rbind(
  null_grid[null_grid$phase %in% phase_filter, , drop = FALSE],
  frontier_grid[frontier_grid$phase %in% phase_filter, , drop = FALSE],
  misspec_grid[misspec_grid$phase %in% phase_filter, , drop = FALSE]
)

balance_id <- ifelse(scenario_grid$regime_balance == "balanced", "bal", "imb")
scenario_grid$scenario_id <- sprintf(
  "%s_%s_n%d_p%d_r%d_%s_%s_c%s_a%s_t%s_rr%s",
  scenario_grid$phase,
  scenario_grid$truth,
  scenario_grid$n_tips,
  scenario_grid$p,
  scenario_grid$n_regimes,
  balance_id,
  scenario_grid$predictor,
  format_level_id(scenario_grid$covariate_effect),
  format_level_id(scenario_grid$alpha),
  format_level_id(scenario_grid$theta_sep, scale = 10000),
  format_level_id(scenario_grid$rate_ratio)
)

tasks <- do.call(
  rbind,
  lapply(seq_len(nrow(scenario_grid)), function(i) {
    data.frame(scenario_index = i, rep = seq_len(scenario_grid$reps[[i]]), stringsAsFactors = FALSE)
  })
)
set.seed(seed_base)
tasks <- tasks[sample.int(nrow(tasks)), , drop = FALSE]

run_replicate <- function(scenario, rep_id, scenario_index) {
  seed_offset <- seed_base + 1000L * rep_id + 10000L * scenario_index
  simmap_info <- make_simmap(
    n_tips = scenario$n_tips,
    n_regimes = scenario$n_regimes,
    regime_balance = scenario$regime_balance,
    seed = seed_offset
  )
  tree <- simmap_info$tree
  predictor <- make_predictor(tree$tip.label, scenario$predictor, simmap_info$tip_states, seed = seed_offset + 1L)
  base_sigma <- build_sigma_base(scenario$p)
  beta_true <- build_beta(predictor$X, scenario$p, scenario$covariate_effect)
  Y <- simulate_truth(
    truth = scenario$truth,
    tree = tree,
    alpha = scenario$alpha,
    base_sigma = base_sigma,
    theta_sep = scenario$theta_sep,
    rate_ratio = scenario$rate_ratio,
    X_full = predictor$X,
    beta = beta_true,
    seed = seed_offset + 2L
  )

  data_full <- c(list(Y = Y), predictor$data)

  rows <- lapply(methods, function(method) {
    fit_ou <- fit_candidate(predictor$formula, data_full, tree, "OU", method)
    fit_oum <- fit_candidate(predictor$formula, data_full, tree, "OUM", method)

    gic_ou <- if (inherits(fit_ou$fit, "error")) list(value = fit_ou$fit, warnings = character()) else compute_gic_safe(fit_ou$fit)
    gic_oum <- if (inherits(fit_oum$fit, "error")) list(value = fit_oum$fit, warnings = character()) else compute_gic_safe(fit_oum$fit)

    gic_values <- c(
      OU = if (inherits(gic_ou$value, "error")) NA_real_ else as.numeric(gic_ou$value$GIC),
      OUM = if (inherits(gic_oum$value, "error")) NA_real_ else as.numeric(gic_oum$value$GIC)
    )
    choice <- choose_best_model(gic_values)
    warnings_by_model <- list(
      OU = collapse_messages(c(
        fit_ou$warnings,
        gic_ou$warnings,
        if (inherits(fit_ou$fit, "error")) conditionMessage(fit_ou$fit) else character()
      )),
      OUM = collapse_messages(c(
        fit_oum$warnings,
        gic_oum$warnings,
        if (inherits(fit_oum$fit, "error")) conditionMessage(fit_oum$fit) else character()
      ))
    )
    expected <- expected_winner(scenario$theta_sep)

    data.frame(
      phase = scenario$phase,
      truth = scenario$truth,
      scenario_id = scenario$scenario_id,
      n_tips = scenario$n_tips,
      p = scenario$p,
      n_regimes = scenario$n_regimes,
      regime_balance = scenario$regime_balance,
      predictor = scenario$predictor,
      covariate_effect = scenario$covariate_effect,
      alpha = scenario$alpha,
      theta_sep = scenario$theta_sep,
      rate_ratio = scenario$rate_ratio,
      rep = rep_id,
      method = method,
      ou_fit_success = !inherits(fit_ou$fit, "error"),
      oum_fit_success = !inherits(fit_oum$fit, "error"),
      both_gic_available = all(is.finite(gic_values)),
      gic_ou = gic_values[["OU"]],
      gic_oum = gic_values[["OUM"]],
      winner = choice$winner,
      correct = if (is.na(choice$winner) || identical(choice$winner, "tie")) NA else identical(choice$winner, expected),
      delta_gic_ou_minus_oum = if (is.finite(gic_values[["OU"]]) && is.finite(gic_values[["OUM"]])) gic_values[["OU"]] - gic_values[["OUM"]] else NA_real_,
      beta_rmse_ou = extract_beta_rmse(fit_ou$fit, beta_true),
      beta_rmse_oum = extract_beta_rmse(fit_oum$fit, beta_true),
      alpha_hat_ou = extract_alpha_hat(fit_ou$fit),
      alpha_hat_oum = extract_alpha_hat(fit_oum$fit),
      ou_messages = warnings_by_model[["OU"]],
      oum_messages = warnings_by_model[["OUM"]],
      has_issue = any(vapply(warnings_by_model, has_text, logical(1))),
      stringsAsFactors = FALSE
    )
  })

  replicate_results <- do.call(rbind, rows)
  write_checkpoint(replicate_results, checkpoint_dir, scenario$scenario_id, rep_id)
  replicate_results
}

cat(
  "Running OU-vs-OUM detectability grid with",
  nrow(scenario_grid), "scenario combinations across", nrow(tasks), "task(s)",
  "on", cores, "worker(s)\n"
)

task_runner <- function(task_row) {
  scenario <- scenario_grid[task_row$scenario_index, , drop = FALSE]
  if (task_row$rep == 1L) {
    cat(sprintf(
      "[worker %s] scenario %d/%d: %s\n",
      Sys.getpid(),
      task_row$scenario_index,
      nrow(scenario_grid),
      scenario$scenario_id[[1]]
    ))
  }
  run_replicate(scenario, task_row$rep, task_row$scenario_index)
}

task_rows <- split(tasks, seq_len(nrow(tasks)))
if (cores > 1L) {
  result_list <- parallel::mclapply(task_rows, task_runner, mc.cores = cores, mc.preschedule = FALSE)
} else {
  result_list <- lapply(task_rows, task_runner)
}

results <- do.call(rbind, result_list)
summary_objects <- summarize_results(results)
summary_table <- summary_objects$summary_table
selection_summary <- summary_objects$selection_summary
boundary_summary <- summary_objects$boundary_summary
issue_rows <- summary_objects$issue_rows

cat("\nCompleted OU-vs-OUM detectability grid.\n")
cat("Scenario summary rows:", nrow(summary_table), "\n")
cat("Selection summary rows:", nrow(selection_summary), "\n")
cat("Boundary summary rows:", nrow(boundary_summary), "\n")
cat("Issue rows:", nrow(issue_rows), "\n")

cat("\nCollapsed selection summary\n")
print(selection_summary, row.names = FALSE)

if (nrow(boundary_summary) > 0L) {
  cat("\nBoundary summary\n")
  print(boundary_summary, row.names = FALSE)
}

output_csv <- Sys.getenv("OUM_DET_OUT", unset = "")
if (nzchar(output_csv)) {
  utils::write.csv(results, output_csv, row.names = FALSE)
  cat("\nSaved replicate-level results to", output_csv, "\n")
}

summary_csv <- Sys.getenv("OUM_DET_SUMMARY_OUT", unset = "")
if (nzchar(summary_csv)) {
  utils::write.csv(summary_table, summary_csv, row.names = FALSE)
  cat("Saved scenario summary to", summary_csv, "\n")
}

selection_csv <- Sys.getenv("OUM_DET_SELECTION_OUT", unset = "")
if (nzchar(selection_csv)) {
  utils::write.csv(selection_summary, selection_csv, row.names = FALSE)
  cat("Saved collapsed selection summary to", selection_csv, "\n")
}

boundary_csv <- Sys.getenv("OUM_DET_BOUNDARY_OUT", unset = "")
if (nzchar(boundary_csv)) {
  utils::write.csv(boundary_summary, boundary_csv, row.names = FALSE)
  cat("Saved boundary summary to", boundary_csv, "\n")
}

issue_csv <- Sys.getenv("OUM_DET_ISSUE_OUT", unset = "")
if (nzchar(issue_csv)) {
  utils::write.csv(issue_rows, issue_csv, row.names = FALSE)
  cat("Saved warning and issue rows to", issue_csv, "\n")
}
