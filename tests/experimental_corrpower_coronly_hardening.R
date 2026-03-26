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
    stop(sprintf("Error message did not match '%s': %s", pattern, conditionMessage(err)), call. = FALSE)
  }
  invisible(TRUE)
}

start_helper <- getFromNamespace(".mvgls_bmm_corrpower_start", "mvMORPH")
profile_helper <- getFromNamespace(".mvgls_bmm_corrpower_profile", "mvMORPH")

make_simmap <- function(seed, n_tips, states_per_regime) {
  set.seed(seed)
  tree <- phytools::pbtree(n = n_tips, scale = 1)
  sampled_states <- sample(rep(names(states_per_regime), times = states_per_regime), n_tips)
  states <- setNames(sampled_states, tree$tip.label)
  suppressMessages(phytools::make.simmap(tree, states, model = "ER", nsim = 1))
}

apply_corrpower <- function(base_sigma, corr_power) {
  eig <- eigen(stats::cov2cor(base_sigma), symmetric = TRUE)
  vals <- pmax(eig$values, .Machine$double.eps)
  corr_raw <- eig$vectors %*% diag(vals ^ corr_power, nrow = length(vals)) %*% t(eig$vectors)
  corr_mat <- stats::cov2cor(0.5 * (corr_raw + t(corr_raw)))
  D <- diag(sqrt(diag(base_sigma)))
  D %*% corr_mat %*% D
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

fit_coronly <- function(formula, data, tree, bmm.reference = NULL) {
  mvgls(
    formula,
    data = data,
    tree = tree,
    model = "BMM",
    method = "LL",
    REML = FALSE,
    bmm.structure = "corrpower",
    bmm.scale = FALSE,
    bmm.reference = bmm.reference,
    echo = FALSE
  )
}

get_diag <- function(fit) {
  diag <- fit$diagnostics[["corrpower"]]
  if (is.null(diag)) stop("corrpower fit with bmm.scale=FALSE did not populate diagnostics$corrpower", call. = FALSE)
  diag
}

set.seed(20260326)
base_sigma <- matrix(c(1.30, 0.55, 0.55, 1.00), 2, 2)

easy_tree <- make_simmap(seed = 20260330, n_tips = 18L, states_per_regime = c(A = 9L, B = 9L))
hard_tree <- make_simmap(seed = 20260331, n_tips = 18L, states_per_regime = c(A = 15L, B = 3L))

easy_Y <- simulate_response(easy_tree, list(A = base_sigma, B = apply_corrpower(base_sigma, 1.35)), seed = 20260332)
hard_Y <- simulate_response(hard_tree, list(A = base_sigma, B = apply_corrpower(base_sigma, 2.50)), seed = 20262133)

easy_fit <- fit_coronly(Y ~ 1, data = list(Y = easy_Y), tree = easy_tree)
easy_diag <- get_diag(easy_fit)

assert_true(is.data.frame(easy_diag$start_table), "diagnostics$corrpower$start_table must be a data.frame")
assert_true(nrow(easy_diag$start_table) >= 4L, "corrpower with bmm.scale=FALSE should populate deterministic multistart candidates")
assert_true(!"scale_multiplier" %in% names(easy_diag$start_table), "corrpower with bmm.scale=FALSE should not expose a scale_multiplier column")
assert_true("corr_power_seed" %in% names(easy_diag$start_table), "corrpower with bmm.scale=FALSE start_table should expose corr_power_seed")
assert_true(sum(as.logical(easy_diag$start_table$selected)) == 1L, "exactly one corrpower start should be marked selected")

best_objective <- min(easy_diag$start_table$objective_value[is.finite(easy_diag$start_table$objective_value)])
assert_true(
  abs(easy_diag$start_table$objective_value[easy_diag$start_table$selected] - best_objective) < 1e-6,
  "selected corrpower fit with bmm.scale=FALSE does not match the best finite candidate"
)
assert_true(is.finite(as.numeric(easy_fit$logLik)), "selected corrpower fit with bmm.scale=FALSE should report a finite log-likelihood")

alt_ref_fit_name <- fit_coronly(Y ~ 1, data = list(Y = easy_Y), tree = easy_tree, bmm.reference = "B")
alt_ref_fit_index <- fit_coronly(Y ~ 1, data = list(Y = easy_Y), tree = easy_tree, bmm.reference = 2L)
assert_true(identical(alt_ref_fit_name$reference_regime, "B"), "named bmm.reference did not set the expected anchor")
assert_true(identical(alt_ref_fit_index$reference_regime, "B"), "indexed bmm.reference did not set the expected anchor")
assert_true(abs(as.numeric(alt_ref_fit_name$logLik) - as.numeric(easy_fit$logLik)) < 2e-2,
            "switching the corrpower reference by name should preserve fit quality up to optimizer tolerance")
assert_true(abs(as.numeric(alt_ref_fit_index$logLik) - as.numeric(alt_ref_fit_name$logLik)) < 1e-5,
            "switching the corrpower reference by index should match the named reference fit")

user_X <- matrix(1, nrow = nrow(easy_Y), ncol = 1L, dimnames = list(rownames(easy_Y), "(Intercept)"))
user_start <- start_helper(easy_tree, easy_Y, user_X, include_scale = FALSE)
user_fit <- mvgls(
  Y ~ 1,
  data = list(Y = easy_Y),
  tree = easy_tree,
  model = "BMM",
  method = "LL",
  REML = FALSE,
  bmm.structure = "corrpower",
  bmm.scale = FALSE,
  start = user_start,
  echo = FALSE
)
user_diag <- get_diag(user_fit)
assert_true(nrow(user_diag$start_table) == 1L, "user-provided corrpower start should disable multistart")
assert_true(all(as.character(user_diag$start_table$source) == "user-provided"),
            "user-provided corrpower start was not marked correctly")

summary_out <- capture.output(summary(easy_fit))
assert_true(any(grepl("Corr-power", summary_out, fixed = TRUE)),
            "summary output did not mention corr-power terminology")
assert_true(any(grepl("selected", summary_out, fixed = TRUE)), "summary output did not report the selected start")
assert_true(any(grepl("corr_power", summary_out, fixed = TRUE)), "summary output did not report corr_power terminology")

conf <- confint(easy_fit, method = "profile", profile_points = 5L)
assert_true(is.matrix(conf) || is.data.frame(conf), "confint() should return a matrix or data.frame")

corr_profile <- profile_helper(
  object = easy_fit,
  parameter = "corr_power",
  regime = "B",
  grid = seq(0.25, 1.5, length.out = 5L)
)
assert_true(is.data.frame(corr_profile), "corrpower corr_power profile helper should return a data.frame")
expect_error(
  profile_helper(
    object = easy_fit,
    parameter = "scale",
    regime = "B",
    grid = c(0.5, 1, 2)
  ),
  "arg"
)

hard_fit <- fit_coronly(Y ~ 1, data = list(Y = hard_Y), tree = hard_tree)
hard_diag <- get_diag(hard_fit)
assert_true(
  isTRUE(hard_diag$boundary_corr_power) || isTRUE(hard_diag$pathological_corr_power),
  "difficult corrpower case did not trigger a corr_power stability flag"
)

expect_error(
  mvgls(
    Y ~ 1,
    data = list(Y = easy_Y),
    tree = easy_tree,
    model = "BMM",
    method = "LL",
    REML = FALSE,
    bmm.structure = "corrstrength",
    echo = FALSE
  ),
  "retired on this branch"
)

cat("corrpower (bmm.scale=FALSE) hardening harness checks passed\n")
