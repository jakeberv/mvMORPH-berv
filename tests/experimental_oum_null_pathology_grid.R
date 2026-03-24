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

  install_lib <- Sys.getenv("OUM_NULL_R_LIB", unset = "")
  if (!nzchar(install_lib)) install_lib <- Sys.getenv("R_LIBS_USER", unset = "")
  if (!nzchar(install_lib)) install_lib <- .libPaths()[1]
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
      repos = Sys.getenv("OUM_NULL_CRAN_REPO", unset = "https://cloud.r-project.org"),
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
    if (length(stale_files)) unlink(stale_files, force = TRUE)
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
    if (n_regimes == 2L) return(normalize_counts(n_tips, c(0.8, 0.2)))
    if (n_regimes == 3L) return(normalize_counts(n_tips, c(0.7, 0.2, 0.1)))
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

simulate_shared_ou <- function(tree, alpha, sigma, mean_mat, seed) {
  Vphy <- .Call(
    "mvmorph_covar_ou_fixed",
    A = ape::vcv.phylo(tree),
    alpha = as.double(alpha),
    sigma = 1,
    PACKAGE = "mvMORPH"
  )
  simulate_gaussian_process(Vphy = Vphy, sigma = sigma, mean_mat = mean_mat, seed = seed)
}

simulate_truth <- function(tree, alpha, base_sigma, X_full, beta, seed) {
  p <- ncol(base_sigma)
  tip_order <- tree$tip.label
  trait_names <- paste0("trait", seq_len(p))
  zero_mean <- matrix(0, nrow = length(tip_order), ncol = p, dimnames = list(tip_order, trait_names))
  Y <- simulate_shared_ou(tree = tree, alpha = alpha, sigma = base_sigma, mean_mat = zero_mean, seed = seed)
  if (nrow(beta) > 0L) {
    Y <- Y + X_full[, rownames(beta), drop = FALSE] %*% beta
  }
  rownames(Y) <- tip_order
  colnames(Y) <- trait_names
  list(Y = Y)
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

extract_alpha_hat <- function(fit_object) {
  if (inherits(fit_object, "error")) return(NA_real_)
  param <- fit_object$param %||% NULL
  if (is.null(param) || !length(param)) return(NA_real_)
  if (!is.null(names(param)) && "alpha" %in% names(param)) return(as.numeric(param[["alpha"]]))
  if (length(param) == 1L && is.finite(param[[1]])) return(as.numeric(param[[1]]))
  NA_real_
}

extract_tuning <- function(fit_object) {
  if (inherits(fit_object, "error")) return(NA_real_)
  tuning <- fit_object$tuning %||% NULL
  if (is.null(tuning) || !length(tuning)) return(NA_real_)
  tuning_vec <- suppressWarnings(as.numeric(unlist(tuning, use.names = FALSE)))
  tuning_vec <- tuning_vec[is.finite(tuning_vec)]
  if (!length(tuning_vec)) return(NA_real_)
  tuning_vec[[1]]
}

extract_penalty <- function(fit_object) {
  if (inherits(fit_object, "error")) return(NA_character_)
  penalty <- fit_object$penalty %||% NULL
  if (is.null(penalty) || !length(penalty)) return(NA_character_)
  as.character(unlist(penalty, use.names = FALSE)[[1]])
}

extract_sigma_metrics <- function(fit_object) {
  default <- list(
    sigma_kappa = NA_real_,
    sigma_min_eig = NA_real_,
    sigma_max_eig = NA_real_
  )
  if (inherits(fit_object, "error")) return(default)
  sigma_hat <- fit_object$sigma$Pinv %||% NULL
  if (is.null(sigma_hat)) return(default)
  sigma_hat <- tryCatch(as.matrix(sigma_hat), error = function(e) NULL)
  if (is.null(sigma_hat) || !nrow(sigma_hat)) return(default)
  eig <- tryCatch(eigen((sigma_hat + t(sigma_hat)) / 2, symmetric = TRUE, only.values = TRUE)$values, error = function(e) rep(NA_real_, nrow(sigma_hat)))
  eig <- Re(eig)
  eig <- eig[is.finite(eig)]
  if (!length(eig)) return(default)
  min_eig <- min(eig)
  max_eig <- max(eig)
  list(
    sigma_kappa = if (is.finite(min_eig) && abs(min_eig) > .Machine$double.eps) max_eig / min_eig else NA_real_,
    sigma_min_eig = min_eig,
    sigma_max_eig = max_eig
  )
}

theta_spread <- function(theta_mat) {
  theta_mat <- as.matrix(theta_mat)
  if (nrow(theta_mat) < 2L) return(0)
  d <- stats::dist(theta_mat)
  mean(as.numeric(d))
}

extract_spurious_theta <- function(fit_object) {
  default <- list(theta_spread_hat = NA_real_, theta_mean_abs = NA_real_)
  if (inherits(fit_object, "error")) return(default)
  coef_est <- coef(fit_object)
  if (!identical(fit_object$model %||% "", "OUM")) return(default)
  beta_rows <- tryCatch({
    mf <- fit_object$model.frame
    form_terms <- stats::delete.response(stats::terms(fit_object))
    X_formula <- stats::model.matrix(form_terms, data = mf)
    colnames(X_formula)[colnames(X_formula) != "(Intercept)"]
  }, error = function(e) character())
  theta_rows <- setdiff(rownames(coef_est), c(beta_rows, "root"))
  theta_rows <- theta_rows[theta_rows %in% rownames(coef_est)]
  if (!length(theta_rows)) return(default)
  theta_est <- as.matrix(coef_est[theta_rows, , drop = FALSE])
  list(
    theta_spread_hat = theta_spread(theta_est),
    theta_mean_abs = mean(abs(theta_est))
  )
}

choose_best_model <- function(values, tol = 1e-8) {
  finite <- values[is.finite(values)]
  if (!length(finite)) return(list(winner = NA_character_, best = NA_real_))
  best <- min(finite)
  winners <- names(finite)[abs(finite - best) <= tol]
  list(winner = if (length(winners) == 1L) winners[[1]] else "tie", best = best)
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
    lapply(split(
      results,
      list(results$scenario_id, results$method),
      drop = TRUE
    ), function(df) {
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
        reps = nrow(df),
        ou_fit_rate = safe_rate(df$ou_fit_success),
        oum_fit_rate = safe_rate(df$oum_fit_success),
        both_gic_rate = safe_rate(df$both_gic_available),
        select_ou_rate = safe_rate(df$winner == "OU"),
        select_oum_rate = safe_rate(df$winner == "OUM"),
        tie_rate = safe_rate(df$winner == "tie"),
        mean_delta_gic_ou_minus_oum = safe_mean(df$delta_gic_ou_minus_oum),
        alpha_hat_ou_mean = safe_mean(df$alpha_hat_ou),
        alpha_hat_oum_mean = safe_mean(df$alpha_hat_oum),
        alpha_hat_oum_selected_mean = safe_mean(df$alpha_hat_oum[df$winner == "OUM"]),
        tuning_ou_mean = safe_mean(df$tuning_ou),
        tuning_oum_mean = safe_mean(df$tuning_oum),
        sigma_kappa_ou_mean = safe_mean(df$sigma_kappa_ou),
        sigma_kappa_oum_mean = safe_mean(df$sigma_kappa_oum),
        sigma_min_eig_ou_mean = safe_mean(df$sigma_min_eig_ou),
        sigma_min_eig_oum_mean = safe_mean(df$sigma_min_eig_oum),
        theta_spread_hat_oum_mean = safe_mean(df$theta_spread_hat_oum),
        theta_spread_hat_oum_selected_mean = safe_mean(df$theta_spread_hat_oum[df$winner == "OUM"]),
        theta_mean_abs_oum_mean = safe_mean(df$theta_mean_abs_oum),
        theta_mean_abs_oum_selected_mean = safe_mean(df$theta_mean_abs_oum[df$winner == "OUM"]),
        theta_selected_count = sum(df$winner == "OUM", na.rm = TRUE),
        warning_rate = safe_rate(df$has_issue),
        stringsAsFactors = FALSE
      )
    })
  )

  design_summary <- do.call(
    rbind,
    lapply(split(
      results,
      list(
        results$method,
        results$p,
        results$n_regimes,
        results$regime_balance,
        results$predictor,
        results$covariate_effect,
        results$alpha
      ),
      drop = TRUE
    ), function(df) {
      data.frame(
        method = df$method[[1]],
        p = df$p[[1]],
        n_regimes = df$n_regimes[[1]],
        regime_balance = df$regime_balance[[1]],
        predictor = df$predictor[[1]],
        covariate_effect = df$covariate_effect[[1]],
        alpha = df$alpha[[1]],
        reps = nrow(df),
        select_oum_rate = safe_rate(df$winner == "OUM"),
        mean_delta_gic_ou_minus_oum = safe_mean(df$delta_gic_ou_minus_oum),
        alpha_hat_oum_mean = safe_mean(df$alpha_hat_oum),
        tuning_oum_mean = safe_mean(df$tuning_oum),
        sigma_kappa_oum_mean = safe_mean(df$sigma_kappa_oum),
        theta_spread_hat_oum_mean = safe_mean(df$theta_spread_hat_oum),
        theta_spread_hat_oum_selected_mean = safe_mean(df$theta_spread_hat_oum[df$winner == "OUM"]),
        warning_rate = safe_rate(df$has_issue),
        stringsAsFactors = FALSE
      )
    })
  )

  compact_summary <- do.call(
    rbind,
    lapply(split(
      results,
      list(results$method, results$p, results$n_regimes, results$regime_balance),
      drop = TRUE
    ), function(df) {
      data.frame(
        method = df$method[[1]],
        p = df$p[[1]],
        n_regimes = df$n_regimes[[1]],
        regime_balance = df$regime_balance[[1]],
        reps = nrow(df),
        select_oum_rate = safe_rate(df$winner == "OUM"),
        mean_delta_gic_ou_minus_oum = safe_mean(df$delta_gic_ou_minus_oum),
        alpha_hat_oum_mean = safe_mean(df$alpha_hat_oum),
        tuning_oum_mean = safe_mean(df$tuning_oum),
        sigma_kappa_oum_mean = safe_mean(df$sigma_kappa_oum),
        theta_spread_hat_oum_selected_mean = safe_mean(df$theta_spread_hat_oum[df$winner == "OUM"]),
        warning_rate = safe_rate(df$has_issue),
        stringsAsFactors = FALSE
      )
    })
  )

  issue_rows <- results[results$has_issue | !results$ou_fit_success | !results$oum_fit_success, , drop = FALSE]

  summary_table <- summary_table[order(
    summary_table$method,
    summary_table$p,
    summary_table$n_regimes,
    summary_table$regime_balance,
    summary_table$predictor,
    summary_table$covariate_effect,
    summary_table$alpha
  ), ]

  design_summary <- design_summary[order(
    design_summary$method,
    design_summary$p,
    design_summary$n_regimes,
    design_summary$regime_balance,
    design_summary$predictor,
    design_summary$covariate_effect,
    design_summary$alpha
  ), ]

  compact_summary <- compact_summary[order(
    compact_summary$method,
    compact_summary$p,
    compact_summary$n_regimes,
    compact_summary$regime_balance
  ), ]

  list(
    summary_table = summary_table,
    design_summary = design_summary,
    compact_summary = compact_summary,
    issue_rows = issue_rows
  )
}

cores <- env_int("OUM_NULL_CORES", max(1L, parallel::detectCores(logical = FALSE) %||% 1L))
seed_base <- env_int("OUM_NULL_SEED", 20260324L)
checkpoint_dir <- Sys.getenv("OUM_NULL_CHECKPOINT_DIR", unset = "")
methods <- env_chr_subset("OUM_NULL_METHODS", c("LL", "LOOCV", "EmpBayes"))

scenario_grid <- expand.grid(
  phase = "null_pathology",
  truth = "shared_OU",
  n_tips = env_num_subset("OUM_NULL_NTIPS", c(60L)),
  p = env_num_subset("OUM_NULL_P", c(24L, 36L, 48L, 64L)),
  n_regimes = env_num_subset("OUM_NULL_REGIMES", c(2L, 3L)),
  regime_balance = env_chr_subset("OUM_NULL_BALANCE", c("balanced", "imbalanced")),
  predictor = env_chr_subset("OUM_NULL_PREDICTOR", c("cont_indep", "cont_confounded", "disc_binary", "disc3_sparse_confounded")),
  covariate_effect = env_num_subset("OUM_NULL_COVARIATE_EFFECT", c(0.00, 0.05)),
  alpha = env_num_subset("OUM_NULL_ALPHA", c(0.25, 0.50, 1.00, 1.50)),
  reps = env_int("OUM_NULL_REPS", 30L),
  stringsAsFactors = FALSE
)

balance_id <- ifelse(scenario_grid$regime_balance == "balanced", "bal", "imb")
scenario_grid$scenario_id <- sprintf(
  "%s_n%d_p%d_r%d_%s_%s_c%s_a%s",
  scenario_grid$truth,
  scenario_grid$n_tips,
  scenario_grid$p,
  scenario_grid$n_regimes,
  balance_id,
  scenario_grid$predictor,
  format_level_id(scenario_grid$covariate_effect),
  format_level_id(scenario_grid$alpha)
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
  truth <- simulate_truth(
    tree = tree,
    alpha = scenario$alpha,
    base_sigma = base_sigma,
    X_full = predictor$X,
    beta = beta_true,
    seed = seed_offset + 2L
  )

  data_full <- c(list(Y = truth$Y), predictor$data)

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
    sigma_ou <- extract_sigma_metrics(fit_ou$fit)
    sigma_oum <- extract_sigma_metrics(fit_oum$fit)
    theta_oum <- extract_spurious_theta(fit_oum$fit)

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
      rep = rep_id,
      method = method,
      ou_fit_success = !inherits(fit_ou$fit, "error"),
      oum_fit_success = !inherits(fit_oum$fit, "error"),
      both_gic_available = all(is.finite(gic_values)),
      gic_ou = gic_values[["OU"]],
      gic_oum = gic_values[["OUM"]],
      winner = choice$winner,
      delta_gic_ou_minus_oum = if (is.finite(gic_values[["OU"]]) && is.finite(gic_values[["OUM"]])) gic_values[["OU"]] - gic_values[["OUM"]] else NA_real_,
      alpha_hat_ou = extract_alpha_hat(fit_ou$fit),
      alpha_hat_oum = extract_alpha_hat(fit_oum$fit),
      tuning_ou = extract_tuning(fit_ou$fit),
      tuning_oum = extract_tuning(fit_oum$fit),
      penalty_ou = extract_penalty(fit_ou$fit),
      penalty_oum = extract_penalty(fit_oum$fit),
      sigma_kappa_ou = sigma_ou$sigma_kappa,
      sigma_kappa_oum = sigma_oum$sigma_kappa,
      sigma_min_eig_ou = sigma_ou$sigma_min_eig,
      sigma_min_eig_oum = sigma_oum$sigma_min_eig,
      sigma_max_eig_ou = sigma_ou$sigma_max_eig,
      sigma_max_eig_oum = sigma_oum$sigma_max_eig,
      theta_spread_hat_oum = theta_oum$theta_spread_hat,
      theta_mean_abs_oum = theta_oum$theta_mean_abs,
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
  "Running focused OUM null pathology grid with",
  nrow(scenario_grid), "scenario combinations across", nrow(tasks), "task(s)",
  "and", length(methods), "method(s) on", cores, "worker(s)\n"
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
design_summary <- summary_objects$design_summary
compact_summary <- summary_objects$compact_summary
issue_rows <- summary_objects$issue_rows

cat("\nCompleted focused OUM null pathology grid.\n")
cat("Scenario summary rows:", nrow(summary_table), "\n")
cat("Design summary rows:", nrow(design_summary), "\n")
cat("Compact summary rows:", nrow(compact_summary), "\n")
cat("Issue rows:", nrow(issue_rows), "\n")

cat("\nCompact summary\n")
print(compact_summary, row.names = FALSE)

output_csv <- Sys.getenv("OUM_NULL_OUT", unset = "")
if (nzchar(output_csv)) {
  utils::write.csv(results, output_csv, row.names = FALSE)
  cat("\nSaved replicate-level results to", output_csv, "\n")
}

summary_csv <- Sys.getenv("OUM_NULL_SUMMARY_OUT", unset = "")
if (nzchar(summary_csv)) {
  utils::write.csv(summary_table, summary_csv, row.names = FALSE)
  cat("Saved scenario summary to", summary_csv, "\n")
}

design_csv <- Sys.getenv("OUM_NULL_DESIGN_OUT", unset = "")
if (nzchar(design_csv)) {
  utils::write.csv(design_summary, design_csv, row.names = FALSE)
  cat("Saved design summary to", design_csv, "\n")
}

compact_csv <- Sys.getenv("OUM_NULL_COMPACT_OUT", unset = "")
if (nzchar(compact_csv)) {
  utils::write.csv(compact_summary, compact_csv, row.names = FALSE)
  cat("Saved compact summary to", compact_csv, "\n")
}

issue_csv <- Sys.getenv("OUM_NULL_ISSUE_OUT", unset = "")
if (nzchar(issue_csv)) {
  utils::write.csv(issue_rows, issue_csv, row.names = FALSE)
  cat("Saved warning and issue rows to", issue_csv, "\n")
}
