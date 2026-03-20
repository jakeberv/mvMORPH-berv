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

compare_case <- function(n_tips, p, fossil_fraction, relation, seed) {
  tree <- make_fossil_tree(seed, n_tips = n_tips, fossil_fraction = fossil_fraction)
  simmap <- make_simmap(seed + 1L, n_tips = n_tips, states_per_regime = c(A = round(n_tips * 0.6), B = round(n_tips * 0.4)))
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
  tree <- phytools::make.simmap(tree, setNames(sample(c("A", "B"), size = length(tree$tip.label), replace = TRUE), tree$tip.label), model = "ER", nsim = 1)
  Y <- simulate_response(tree, list(A = base_sigma, B = derived_sigma), seed = seed + 50L)
  fit_power <- fit_corrpower(Y ~ 1, data = list(Y = Y), tree = tree)
  fit_strength <- fit_corrstrength(Y ~ 1, data = list(Y = Y), tree = tree)
  cov_power <- fit_power$sigma$regime[[2]]
  cov_strength <- fit_strength$sigma$regime[[2]]
  cov_true <- derived_sigma
  frob <- function(a, b) sqrt(sum((a - b)^2))
  data.frame(
    n_tips = n_tips,
    p = p,
    fossil_fraction = fossil_fraction,
    relation = relation,
    corrpower_cov_err = frob(cov_power, cov_true),
    corrstrength_cov_err = frob(cov_strength, cov_true),
    corrpower_pathological = isTRUE(getElement(fit_power$diagnostics$corrpower, "pathological_scale")),
    corrstrength_pathological = isTRUE(getElement(fit_strength$diagnostics$corrstrength, "pathological_scale")),
    corrpower_boundary = isTRUE(getElement(fit_power$diagnostics$corrpower, "boundary_corr_power")),
    corrstrength_boundary = isTRUE(getElement(fit_strength$diagnostics$corrstrength, "boundary_corr_strength")),
    logLik_power = as.numeric(fit_power$logLik),
    logLik_strength = as.numeric(fit_strength$logLik),
    stringsAsFactors = FALSE
  )
}

full_grid <- expand.grid(
  n_tips = c(50L, 100L),
  p = c(4L),
  fossil_fraction = c(0, 0.5),
  relation = c("weaker", "stronger"),
  stringsAsFactors = FALSE
)

grid <- if (identical(toupper(Sys.getenv("CORRPOWER_COMPARISON_FULL", "FALSE")), "TRUE")) {
  full_grid
} else {
  subset(full_grid, n_tips == 50L)
}

res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  compare_case(g$n_tips, g$p, g$fossil_fraction, g$relation, seed = 20260700 + i)
}))

print(res)

cat("\nComparison summary\n")
cat(sprintf("corrpower median covariance error: %.4f\n", median(res$corrpower_cov_err, na.rm = TRUE)))
cat(sprintf("corrstrength median covariance error: %.4f\n", median(res$corrstrength_cov_err, na.rm = TRUE)))
cat(sprintf("corrpower pathology rate: %.3f\n", mean(res$corrpower_pathological, na.rm = TRUE)))
cat(sprintf("corrstrength pathology rate: %.3f\n", mean(res$corrstrength_pathological, na.rm = TRUE)))

assert_true(any(is.finite(res$corrpower_cov_err)), "corrpower comparison produced no finite covariance errors")
assert_true(any(is.finite(res$corrstrength_cov_err)), "corrstrength comparison produced no finite covariance errors")

cat("bmm correlation model comparison harness checks passed\n")
