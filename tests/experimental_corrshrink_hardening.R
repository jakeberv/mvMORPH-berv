#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
if (!length(file_arg)) stop("This script must be run with Rscript.", call. = FALSE)
repo_root <- normalizePath(file.path(dirname(file_arg), ".."))

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required to run this harness.", call. = FALSE)
}

devtools::load_all(repo_root, quiet = TRUE)

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

`%||%` <- function(x, y) if (is.null(x)) y else x

start_helper <- getFromNamespace(".mvgls_bmm_corrshrink_start", "mvMORPH")
profile_helper <- NULL
for (nm in c(".mvgls_bmm_corrshrink_profile", ".mvgls_corrshrink_profile")) {
  if (exists(nm, envir = asNamespace("mvMORPH"), inherits = FALSE)) {
    profile_helper <- getFromNamespace(nm, "mvMORPH")
    break
  }
}
if (is.null(profile_helper)) {
  stop("corr-strength profile helper was not found in the mvMORPH namespace", call. = FALSE)
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

simulate_response <- function(tree, sigma_by_regime, seed) {
  set.seed(seed)
  y <- mvSIM(
    tree,
    nsim = 1,
    model = "BMM",
    param = list(ntraits = 2, sigma = sigma_by_regime, theta = c(0, 0))
  )
  if (is.list(y)) y <- y[[1]]
  as.matrix(y)
}

fit_corrshrink <- function(formula, data, tree, start = NULL, bmm.reference = NULL) {
  mvgls(
    formula,
    data = data,
    tree = tree,
    model = "BMM",
    method = "LL",
    REML = FALSE,
    bmm.structure = "corrshrink",
    bmm.reference = bmm.reference,
    start = start,
    echo = FALSE
  )
}

get_diag <- function(fit) {
  diag <- fit$diagnostics[["corrshrink"]]
  if (is.null(diag)) stop("corr-strength fit did not populate diagnostics$corrshrink", call. = FALSE)
  diag
}

check_start_table <- function(fit, expected_nstarts, expected_source = NULL) {
  diag <- get_diag(fit)
  st <- diag$start_table
  assert_true(is.data.frame(st), "diagnostics$corrshrink$start_table must be a data.frame")
  assert_true(nrow(st) == expected_nstarts, sprintf("expected %d starts, got %d", expected_nstarts, nrow(st)))
  required <- c("start_id", "source", "scale_multiplier", "kappa_seed", "convergence", "nloglik", "max_scale", "selected")
  assert_true(all(required %in% names(st)), "start_table is missing one or more required columns")
  assert_true(sum(as.logical(st$selected)) == 1L, "exactly one start should be marked selected")
  if (!is.null(expected_source)) {
    assert_true(all(as.character(st$source) == expected_source), "user-supplied fit was not marked as user-provided")
  }
  invisible(diag)
}

base_sigma <- matrix(c(1.30, 0.55, 0.55, 1.00), 2, 2)
hard_sigma <- base_sigma
easy_sigma_derived <- 1.8 * base_sigma
hard_sigma_derived <- matrix(c(2.40, 0.00, 0.00, 1.92), 2, 2)

easy_tree <- make_simmap(seed = 20260321, n_tips = 18L, states_per_regime = c(A = 9L, B = 9L))
hard_tree <- make_simmap(seed = 20260323, n_tips = 18L, states_per_regime = c(A = 15L, B = 3L))

easy_Y <- simulate_response(easy_tree, list(A = base_sigma, B = easy_sigma_derived), seed = 20260323)
hard_Y <- simulate_response(hard_tree, list(A = hard_sigma, B = hard_sigma_derived), seed = 20260324)

easy_fit <- fit_corrshrink(Y ~ 1, data = list(Y = easy_Y), tree = easy_tree)
easy_diag <- check_start_table(easy_fit, expected_nstarts = 12L)

summary_out <- capture.output(summary(easy_fit))
assert_true(any(grepl("Regime summary", summary_out, fixed = TRUE)), "summary output did not report the corr-strength regime summary")
assert_true(any(grepl("nstarts", summary_out, fixed = TRUE)), "summary output did not report corr-strength start diagnostics")
assert_true(any(grepl("selected", summary_out, fixed = TRUE)), "summary output did not report the selected start")
assert_true(any(grepl("kappa", summary_out, fixed = TRUE)), "summary output did not report corr-strength kappa terminology")

best_nloglik <- min(easy_diag$start_table$nloglik[is.finite(easy_diag$start_table$nloglik)])
assert_true(
  abs(as.numeric(easy_fit$logLik) - (-best_nloglik)) < 1e-6,
  "selected fit logLik does not match the best finite candidate in start_table"
)
assert_true(
  identical(easy_diag$selected_start_id, easy_diag$start_table$start_id[easy_diag$start_table$selected]),
  "selected_start_id does not match the selected row in start_table"
)

alt_ref_fit_name <- fit_corrshrink(Y ~ 1, data = list(Y = easy_Y), tree = easy_tree, bmm.reference = "B")
alt_ref_fit_index <- fit_corrshrink(Y ~ 1, data = list(Y = easy_Y), tree = easy_tree, bmm.reference = 2L)
alt_ref_diag <- get_diag(alt_ref_fit_name)
assert_true(identical(alt_ref_fit_name$reference_regime, "B"), "named bmm.reference did not set the expected anchor")
assert_true(identical(alt_ref_fit_index$reference_regime, "B"), "indexed bmm.reference did not set the expected anchor")
assert_true(identical(alt_ref_diag$reference_regime, "B"), "diagnostics did not record the selected reference regime")
assert_true(abs(as.numeric(alt_ref_fit_name$logLik) - as.numeric(easy_fit$logLik)) < 1e-2,
  "switching the corr-strength reference by name should preserve fit quality up to optimizer tolerance"
)
assert_true(abs(as.numeric(alt_ref_fit_index$logLik) - as.numeric(alt_ref_fit_name$logLik)) < 1e-5,
  "switching the corr-strength reference by index should match the named reference fit"
)
assert_true(abs(unname(alt_ref_fit_name$param["B.scale"]) - 1) < 1e-8,
  "the selected reference regime must have scale fixed at 1"
)
assert_true(abs(unname(alt_ref_fit_name$param["B.kappa"]) - 1) < 1e-8,
  "the selected reference regime must have kappa fixed at 1"
)
assert_true(all(names(alt_ref_fit_name$sigma$regime) == names(alt_ref_fit_index$sigma$regime)),
  "named and indexed reference fits should preserve the same regime ordering"
)
for (nm in names(alt_ref_fit_name$sigma$regime)) {
  assert_true(
    max(abs(alt_ref_fit_name$sigma$regime[[nm]] - alt_ref_fit_index$sigma$regime[[nm]])) < 1e-6,
    sprintf("named and indexed reference fits disagreed for regime %s", nm)
  )
}

user_X <- matrix(1, nrow = nrow(easy_Y), ncol = 1L, dimnames = list(rownames(easy_Y), "(Intercept)"))
user_start <- start_helper(easy_tree, easy_Y, user_X)
user_fit <- fit_corrshrink(Y ~ 1, data = list(Y = easy_Y), tree = easy_tree, start = user_start)
user_diag <- check_start_table(user_fit, expected_nstarts = 1L, expected_source = "user-provided")
assert_true(isTRUE(user_diag$start_table$selected[1]), "user-provided start should be marked selected")

hard_fit <- fit_corrshrink(Y ~ 1, data = list(Y = hard_Y), tree = hard_tree)
hard_diag <- get_diag(hard_fit)
assert_true(
  isTRUE(hard_diag$boundary_kappa) || isTRUE(hard_diag$pathological_scale),
  "difficult seeded case did not trigger boundary/pathology diagnostics"
)

profile_grid_kappa <- seq(0.05, 0.95, length.out = 5L)
profile_grid_scale <- c(0.5, 1, 2, 4)

kappa_profile <- tryCatch(
  do.call(profile_helper, list(
    object = easy_fit,
    parameter = "kappa",
    regime = "B",
    grid = profile_grid_kappa
  )),
  error = function(e) e
)
if (inherits(kappa_profile, "error")) stop(conditionMessage(kappa_profile), call. = FALSE)
assert_true(is.data.frame(kappa_profile), "kappa profile helper did not return a data.frame")
assert_true(all(c("fixed_value", "logLik", "convergence", "max_scale") %in% names(kappa_profile)),
  "kappa profile helper returned unexpected columns"
)
assert_true(nrow(kappa_profile) == length(profile_grid_kappa), "kappa profile helper returned the wrong number of rows")
assert_true(sum(is.finite(kappa_profile$logLik)) >= 4L, "kappa profile helper should return mostly finite log-likelihoods")
expect_error(
  do.call(profile_helper, list(
    object = alt_ref_fit_name,
    parameter = "kappa",
    regime = "B",
    grid = profile_grid_kappa
  )),
  "reference regime"
)

scale_profile <- tryCatch(
  do.call(profile_helper, list(
    object = easy_fit,
    parameter = "scale",
    regime = "B",
    grid = profile_grid_scale
  )),
  error = function(e) e
)
if (inherits(scale_profile, "error")) stop(conditionMessage(scale_profile), call. = FALSE)
assert_true(is.data.frame(scale_profile), "scale profile helper did not return a data.frame")
assert_true(all(c("fixed_value", "logLik", "convergence", "max_scale") %in% names(scale_profile)),
  "scale profile helper returned unexpected columns"
)
assert_true(nrow(scale_profile) == length(profile_grid_scale), "scale profile helper returned the wrong number of rows")
assert_true(sum(is.finite(scale_profile$logLik)) >= 3L, "scale profile helper should return mostly finite log-likelihoods")

cat("corr-strength hardening harness checks passed\n")
