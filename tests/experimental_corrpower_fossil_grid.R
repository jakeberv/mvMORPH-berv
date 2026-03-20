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

apply_corrpower <- function(base_sigma, scale, corr_power) {
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

make_fossil_tree <- function(seed, n_extant = 16L, n_fossil = 8L, fossil_depth = c("deep", "mixed", "shallow")) {
  set.seed(seed)
  full <- phytools::pbtree(n = n_extant, scale = 1)
  n_fossil <- max(1L, min(n_fossil, n_extant - 1L))
  fossil_depth <- match.arg(fossil_depth)
  if (fossil_depth == "deep") {
    fossil_branches <- head(order(full$edge.length, decreasing = TRUE), n_fossil)
  } else if (fossil_depth == "shallow") {
    fossil_branches <- head(order(full$edge.length, decreasing = FALSE), n_fossil)
  } else {
    fossil_branches <- sample(seq_len(nrow(full$edge)), n_fossil)
  }
  fossil_tips <- unique(full$edge[fossil_branches, 2])
  fossil_tips <- fossil_tips[fossil_tips <= length(full$tip.label)]
  fossil_labels <- full$tip.label[fossil_tips]
  extant <- full
  keep <- setdiff(full$tip.label, fossil_labels)
  extant <- ape::drop.tip(full, fossil_labels)
  list(full = full, extant = extant, fossil_labels = fossil_labels)
}

build_simmap <- function(tree, seed, states_per_regime = c(A = 10L, B = 6L)) {
  set.seed(seed)
  sampled <- sample(rep(names(states_per_regime), times = states_per_regime), length(tree$tip.label))
  states <- setNames(sampled, tree$tip.label)
  suppressMessages(phytools::make.simmap(tree, states, model = "ER", nsim = 1))
}

full_grid <- expand.grid(
  n_extant = c(16L, 24L),
  p = c(2L, 4L),
  fossil_fraction = c(0.25, 0.50),
  fossil_depth = c("deep", "mixed"),
  mapped_fraction = c(0.20, 0.80),
  base_signal = c(0.20, 0.80),
  relation = c("weaker", "stronger"),
  stringsAsFactors = FALSE
)

scenario_grid <- if (identical(toupper(Sys.getenv("CORRPOWER_FOSSIL_GRID_FULL", "FALSE")), "TRUE")) {
  full_grid
} else {
  subset(
    full_grid,
    n_extant == 16L &
      p == 2L &
      fossil_depth == "deep" &
      mapped_fraction == 0.20 &
      base_signal == 0.20
  )
}

results <- list()
counter <- 1L
for (i in seq_len(nrow(scenario_grid))) {
  sc <- scenario_grid[i, ]
  n_fossil <- max(1L, round(sc$n_extant * sc$fossil_fraction))
  tr <- make_fossil_tree(20260380 + i, n_extant = sc$n_extant, n_fossil = n_fossil, fossil_depth = sc$fossil_depth)
  tree_full <- build_simmap(tr$full, seed = 20260400 + i, states_per_regime = c(A = round(sc$n_extant * 0.6), B = round(sc$n_extant * 0.4)))
  tree_extant <- build_simmap(tr$extant, seed = 20260420 + i, states_per_regime = c(A = round(length(tr$extant$tip.label) * 0.6), B = max(2L, round(length(tr$extant$tip.label) * 0.4))))
  base_sigma <- if (sc$p == 2L) {
    if (sc$base_signal < 0.5) matrix(c(1.3, 0.08, 0.08, 1.0), 2, 2) else matrix(c(1.3, 0.55, 0.55, 1.0), 2, 2)
  } else {
    matrix(c(
      1.2, 0.20, 0.10, 0.05,
      0.20, 1.1, 0.15, 0.10,
      0.10, 0.15, 1.0, 0.20,
      0.05, 0.10, 0.20, 0.9
    ), 4, 4)
  }
  if (sc$relation == "weaker") {
    scale_b <- 1.5
    corr_b <- if (sc$base_signal < 0.5) 0.45 else 0.70
  } else {
    scale_b <- 1.8
    corr_b <- if (sc$base_signal < 0.5) 1.30 else 1.50
  }
  sigma_fossil <- list(A = base_sigma, B = apply_corrpower(base_sigma, scale = scale_b, corr_power = corr_b))
  Y_full <- simulate_response(tree_full, sigma_fossil, seed = 20260500 + i)
  Y_extant <- Y_full[tree_extant$tip.label, , drop = FALSE]
  fit_full <- fit_corrpower(Y ~ 1, data = list(Y = Y_full), tree = tree_full)
  fit_extant <- fit_corrpower(Y ~ 1, data = list(Y = Y_extant), tree = tree_extant)
  results[[counter]] <- data.frame(
    scenario = paste(sc, collapse = "|"),
    n_extant = sc$n_extant,
    p = sc$p,
    fossil_fraction = sc$fossil_fraction,
    fossil_depth = sc$fossil_depth,
    mapped_fraction = sc$mapped_fraction,
    base_signal = sc$base_signal,
    relation = sc$relation,
    logLik_full = as.numeric(fit_full$logLik),
    logLik_extant = as.numeric(fit_extant$logLik),
    full_pathological_scale = isTRUE(getElement(fit_full$diagnostics$corrpower, "pathological_scale")),
    extant_pathological_scale = isTRUE(getElement(fit_extant$diagnostics$corrpower, "pathological_scale")),
    full_boundary = isTRUE(getElement(fit_full$diagnostics$corrpower, "boundary_corr_power")),
    extant_boundary = isTRUE(getElement(fit_extant$diagnostics$corrpower, "boundary_corr_power")),
    stringsAsFactors = FALSE
  )
  counter <- counter + 1L
}

results_df <- do.call(rbind, results)
print(results_df)

cat("\nFossil vs extant summary\n")
cat(sprintf("mean logLik delta: %.4f\n", mean(results_df$logLik_full - results_df$logLik_extant, na.rm = TRUE)))
cat(sprintf("fossil pathological rate: %.3f\n", mean(results_df$full_pathological_scale, na.rm = TRUE)))
cat(sprintf("extant pathological rate: %.3f\n", mean(results_df$extant_pathological_scale, na.rm = TRUE)))

if (!any(is.finite(results_df$logLik_full))) {
  stop("No corr-power fossil fits succeeded.", call. = FALSE)
}

cat("corr-power fossil grid harness checks passed\n")
