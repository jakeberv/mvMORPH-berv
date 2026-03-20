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
    stop(
      sprintf("Error message did not match '%s': %s", pattern, conditionMessage(err)),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

start_helper <- NULL
for (nm in c(".mvgls_bmm_corrpower_start", ".mvgls_bmm_corrstrength_start", ".mvgls_bmm_corrshrink_start")) {
  if (exists(nm, envir = asNamespace("mvMORPH"), inherits = FALSE)) {
    start_helper <- getFromNamespace(nm, "mvMORPH")
    break
  }
}
if (is.null(start_helper)) {
  stop("corr-power start helper was not found in the mvMORPH namespace", call. = FALSE)
}

make_simmap <- function(seed, n_tips, states_per_regime) {
  set.seed(seed)
  tree <- phytools::pbtree(n = n_tips, scale = 1)
  sampled_states <- sample(rep(names(states_per_regime), times = states_per_regime), n_tips)
  states <- setNames(sampled_states, tree$tip.label)
  invisible(capture.output({
    simmap <- suppressMessages(phytools::make.simmap(tree, states, model = "ER", nsim = 1))
  }))
  simmap
}

apply_corrpower <- function(base_sigma, scale, corr_power) {
  eig <- eigen(stats::cov2cor(base_sigma), symmetric = TRUE)
  vals <- pmax(eig$values, .Machine$double.eps)
  corr_raw <- eig$vectors %*% diag(vals ^ corr_power, nrow = length(vals)) %*% t(eig$vectors)
  corr_mat <- stats::cov2cor(0.5 * (corr_raw + t(corr_raw)))
  D <- diag(sqrt(diag(base_sigma)))
  scale * (D %*% corr_mat %*% D)
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

fit_corrpower <- function(formula, data, tree, start = NULL) {
  mvgls(
    formula,
    data = data,
    tree = tree,
    model = "BMM",
    method = "LL",
    REML = FALSE,
    bmm.structure = "corrpower",
    start = start,
    echo = FALSE
  )
}

multi_fit <- function(label, formula, data, tree, single_start) {
  fits <- list(
    `single-start` = tryCatch(fit_corrpower(formula, data = data, tree = tree, start = single_start), error = identity),
    `multi-start` = tryCatch(fit_corrpower(formula, data = data, tree = tree), error = identity)
  )
  out <- lapply(fits, function(fit) {
    if (inherits(fit, "error")) {
      return(list(
        success = FALSE,
        logLik = NA_real_,
        corr_power = NA_real_,
        scale = NA_real_,
        fit = fit
      ))
    }
    param <- fit$param
    corr_power <- param[grep("\\.corr_power$", names(param))]
    scale <- param[grep("\\.scale$", names(param))]
    list(
      success = TRUE,
      logLik = as.numeric(fit$logLik),
      corr_power = unname(corr_power),
      scale = unname(scale),
      fit = fit
    )
  })
  data.frame(
    scenario = label,
    mode = names(out),
    success = vapply(out, `[[`, logical(1), "success"),
    logLik = vapply(out, `[[`, numeric(1), "logLik"),
    corr_power = I(lapply(out, `[[`, "corr_power")),
    scale = I(lapply(out, `[[`, "scale")),
    stringsAsFactors = FALSE
  )
}

scenarios <- list(
  list(
    label = "balanced_2_trait_proportional",
    n_tips = 24L, p = 2L, relation = "equal",
    states = c(A = 12L, B = 12L),
    base_sigma = matrix(c(1.3, 0.55, 0.55, 1.0), 2, 2),
    derived_scale = 1.8,
    derived_corr_power = 1.0,
    x = NULL,
    beta = NULL
  ),
  list(
    label = "balanced_2_trait_weaker",
    n_tips = 24L, p = 2L, relation = "weaker",
    states = c(A = 12L, B = 12L),
    base_sigma = matrix(c(1.3, 0.60, 0.60, 1.0), 2, 2),
    derived_scale = 1.5,
    derived_corr_power = 0.45,
    x = NULL,
    beta = NULL
  ),
  list(
    label = "balanced_2_trait_stronger",
    n_tips = 24L, p = 2L, relation = "stronger",
    states = c(A = 12L, B = 12L),
    base_sigma = matrix(c(1.3, 0.40, 0.40, 1.0), 2, 2),
    derived_scale = 1.5,
    derived_corr_power = 1.45,
    x = NULL,
    beta = NULL
  ),
  list(
    label = "sparse_derived_regime",
    n_tips = 24L, p = 2L, relation = "weaker",
    states = c(A = 20L, B = 4L),
    base_sigma = matrix(c(1.3, 0.50, 0.50, 1.0), 2, 2),
    derived_scale = 1.9,
    derived_corr_power = 0.4,
    x = NULL,
    beta = NULL
  ),
  list(
    label = "weak_offdiag_signal",
    n_tips = 24L, p = 2L, relation = "weak_offdiag",
    states = c(A = 12L, B = 12L),
    base_sigma = matrix(c(1.3, 0.05, 0.05, 1.0), 2, 2),
    derived_scale = 1.8,
    derived_corr_power = 0.8,
    x = NULL,
    beta = NULL
  ),
  list(
    label = "4_trait_covariate",
    n_tips = 30L, p = 4L, relation = "covariate",
    states = c(A = 15L, B = 15L),
    base_sigma = diag(c(1.2, 1.0, 0.9, 0.8)) + 0.25,
    derived_scale = 1.6,
    derived_corr_power = 1.2,
    x = as.numeric(scale(seq_len(30L))),
    beta = c(0.1, -0.05)
  )
)

scenario_results <- list()
for (i in seq_along(scenarios)) {
  scenario <- scenarios[[i]]
  tree <- make_simmap(seed = 20260320 + i, n_tips = scenario$n_tips, states_per_regime = scenario$states)
  if (is.null(scenario$x)) {
    Y <- simulate_response(tree, list(
      A = scenario$base_sigma,
      B = apply_corrpower(scenario$base_sigma, scenario$derived_scale, scenario$derived_corr_power)
    ), seed = 20260350 + i)
    dat <- list(Y = Y)
    formula <- Y ~ 1
  } else {
    x <- scenario$x
    names(x) <- tree$tip.label
    Y <- simulate_response(tree, list(
      A = scenario$base_sigma,
      B = apply_corrpower(scenario$base_sigma, scenario$derived_scale, scenario$derived_corr_power)
    ), seed = 20260350 + i, x = x, beta = scenario$beta)
    dat <- list(Y = Y, x = x)
    formula <- Y ~ x
  }
  single_start <- start_helper(tree, Y, if (is.null(scenario$x)) matrix(1, nrow = nrow(Y), ncol = 1L, dimnames = list(rownames(Y), "(Intercept)")) else cbind("(Intercept)" = 1, x = dat$x))
  scenario_results[[i]] <- multi_fit(scenario$label, formula, dat, tree, single_start)
}

results <- do.call(rbind, scenario_results)
print(results)

summarise_mode <- function(df) {
  ok <- df$success
  data.frame(
    n = sum(ok),
    success_rate = mean(df$success),
    mean_logLik = mean(df$logLik[ok], na.rm = TRUE),
    median_logLik = median(df$logLik[ok], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

cat("\nScenario summary\n")
print(do.call(rbind, lapply(split(results, results$scenario), summarise_mode)))

if (!any(results$success)) {
  stop("No corr-power fits succeeded in the identifiability study.", call. = FALSE)
}

cat("corr-power identifiability harness checks passed\n")
