#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
if (!length(file_arg)) stop("This script must be run with Rscript.", call. = FALSE)
repo_root <- normalizePath(file.path(dirname(file_arg), ".."))

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required to run this script.", call. = FALSE)
}

devtools::load_all(repo_root, quiet = TRUE)

`%||%` <- function(x, y) if (is.null(x)) y else x

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed <= 0L) default else parsed
}

get_internal <- function(candidates) {
  for (nm in candidates) {
    if (exists(nm, envir = asNamespace("mvMORPH"), inherits = FALSE)) {
      return(getFromNamespace(nm, "mvMORPH"))
    }
  }
  stop(sprintf("Could not find any of: %s", paste(candidates, collapse = ", ")), call. = FALSE)
}

profile_helper <- get_internal(c(".mvgls_bmm_corrshrink_profile", ".mvgls_corrshrink_profile"))
refit_helper <- get_internal(c(".mvgls_corrshrink_refit"))

nboot <- env_int("CORRSHRINK_ANCHOR_BOOT", 4L)
profile_points <- env_int("CORRSHRINK_ANCHOR_PROFILE_POINTS", 9L)
seed_base <- env_int("CORRSHRINK_ANCHOR_SEED", 20260320L)

make_sigma <- function(sds, corr) {
  diag(sds) %*% corr %*% diag(sds)
}

apply_corrstrength <- function(base_sigma, scale, kappa) {
  V <- diag(diag(base_sigma))
  W <- base_sigma - V
  scale * (V + kappa * W)
}

make_simmap <- function(n_tips, state_counts, seed) {
  set.seed(seed)
  tree <- phytools::pbtree(n = n_tips, scale = 1)
  labels <- tree$tip.label
  sampled_states <- sample(rep(names(state_counts), times = state_counts), length(labels))
  states <- setNames(sampled_states, labels)
  invisible(capture.output({
    simmap <- suppressMessages(phytools::make.simmap(tree, states, model = "ER", nsim = 1))
  }))
  simmap
}

simulate_response <- function(tree, sigma_by_regime, seed) {
  set.seed(seed)
  Y <- mvSIM(
    tree,
    nsim = 1,
    model = "BMM",
    param = list(ntraits = nrow(sigma_by_regime[[1]]), sigma = sigma_by_regime, theta = rep(0, nrow(sigma_by_regime[[1]])))
  )
  if (is.list(Y)) Y <- Y[[1]]
  as.matrix(Y)
}

fit_corrshrink <- function(formula, data, tree, bmm.reference = NULL) {
  mvgls(
    formula,
    data = data,
    tree = tree,
    model = "BMM",
    method = "LL",
    REML = FALSE,
    bmm.structure = "corrshrink",
    bmm.reference = bmm.reference,
    echo = FALSE
  )
}

summarise_regime_cov <- function(Sigma) {
  offdiag <- Sigma[upper.tri(Sigma)]
  data.frame(
    mean_var = mean(diag(Sigma)),
    mean_cov = if (length(offdiag)) mean(offdiag) else 0,
    frobenius = sqrt(sum(Sigma^2)),
    stringsAsFactors = FALSE
  )
}

make_scale_grid <- function(est, n = profile_points) {
  center <- max(as.numeric(est), 1e-4)
  lo <- max(center / 4, 1e-4)
  hi <- max(center * 4, lo * 1.1)
  exp(seq(log(lo), log(hi), length.out = n))
}

make_kappa_grid <- function(est, n = profile_points) {
  center <- max(as.numeric(est), 0.02)
  lo <- max(0.01, center - 0.6)
  hi <- max(center + 0.6, lo * 1.5)
  unique(round(seq(lo, hi, length.out = n), 6))
}

profile_interval <- function(fit, regime, parameter, level = 0.95) {
  est <- unname(fit$param[paste0(regime, ".", parameter)])
  grid <- if (identical(parameter, "scale")) make_scale_grid(est) else make_kappa_grid(est)
  prof <- tryCatch(
    do.call(profile_helper, list(object = fit, parameter = parameter, regime = regime, grid = grid)),
    error = function(e) e
  )
  if (inherits(prof, "error")) {
    return(list(low = NA_real_, high = NA_real_, table = NULL, error = conditionMessage(prof)))
  }
  finite <- prof[is.finite(prof$logLik), , drop = FALSE]
  if (!nrow(finite)) {
    return(list(low = NA_real_, high = NA_real_, table = prof, error = NULL))
  }
  cutoff <- max(finite$logLik) - stats::qchisq(level, df = 1) / 2
  inside <- finite$fixed_value[finite$logLik >= cutoff]
  list(
    low = if (length(inside)) min(inside) else NA_real_,
    high = if (length(inside)) max(inside) else NA_real_,
    table = prof,
    error = NULL
  )
}

bootstrap_free_params <- function(fit, nboot) {
  sims <- simulate(fit, nsim = nboot)
  if (!is.list(sims)) sims <- list(sims)
  free_regimes <- setdiff(names(fit$sigma$regime), fit$reference_regime)
  rows <- vector("list", length(sims) * max(length(free_regimes), 1L))
  idx <- 1L
  for (draw_id in seq_along(sims)) {
    refit <- tryCatch(
      refit_helper(fit, as.matrix(sims[[draw_id]]), start = fit$opt$par, echo = FALSE),
      error = function(e) e
    )
    if (inherits(refit, "error")) next
    for (regime in free_regimes) {
      scale_est <- unname(refit$param[paste0(regime, ".scale")])
      kappa_est <- unname(refit$param[paste0(regime, ".kappa")])
      rows[[idx]] <- data.frame(
        draw = draw_id,
        anchor = fit$reference_regime,
        regime = regime,
        scale = scale_est,
        kappa = kappa_est,
        kappa_boundary = as.numeric(is.finite(kappa_est) && (kappa_est < 0.05 || kappa_est > 1.95)),
        runaway_scale = as.numeric(is.finite(scale_est) && scale_est > 20),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  rows <- rows[seq_len(max(idx - 1L, 0L))]
  if (!length(rows)) {
    data.frame(
      draw = integer(0),
      anchor = character(0),
      regime = character(0),
      scale = numeric(0),
      kappa = numeric(0),
      kappa_boundary = numeric(0),
      runaway_scale = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, rows)
  }
}

bootstrap_summary <- function(boot_df, regime, parameter) {
  sub <- boot_df[boot_df$regime == regime, , drop = FALSE]
  values <- sub[[parameter]]
  values <- values[is.finite(values)]
  if (!length(values)) {
    return(list(low = NA_real_, high = NA_real_, mean = NA_real_, n = 0L))
  }
  list(
    low = unname(stats::quantile(values, probs = 0.025, names = FALSE)),
    high = unname(stats::quantile(values, probs = 0.975, names = FALSE)),
    mean = mean(values),
    n = length(values)
  )
}

base_sigma <- make_sigma(
  sds = c(1.3, 1.0),
  corr = matrix(c(1, 0.45, 0.45, 1), 2, 2)
)
sigma_by_regime <- list(
  A = base_sigma,
  B = apply_corrstrength(base_sigma, scale = 1.7, kappa = 1.25),
  C = apply_corrstrength(base_sigma, scale = 1.25, kappa = 0.35)
)

tree <- make_simmap(n_tips = 24L, state_counts = c(A = 8L, B = 8L, C = 8L), seed = seed_base)
Y <- simulate_response(tree, sigma_by_regime, seed = seed_base + 1L)
data <- list(Y = Y)

anchors <- names(sigma_by_regime)
fits <- lapply(anchors, function(anchor) fit_corrshrink(Y ~ 1, data = data, tree = tree, bmm.reference = anchor))
names(fits) <- anchors

anchor_rows <- lapply(fits, function(fit) {
  diag <- fit$diagnostics$corrshrink
  data.frame(
    anchor = fit$reference_regime,
    logLik = as.numeric(fit$logLik),
    AIC = AIC(fit)$AIC,
    GIC = GIC(fit)$GIC,
    selected_start = diag$selected_start_id %||% NA_integer_,
    max_scale = diag$max_scale %||% NA_real_,
    boundary_kappa = diag$boundary_kappa %||% NA,
    pathological_scale = diag$pathological_scale %||% NA,
    stringsAsFactors = FALSE
  )
})
anchor_summary <- do.call(rbind, anchor_rows)
rownames(anchor_summary) <- NULL

param_rows <- list()
cov_rows <- list()
boot_cache <- vector("list", length(fits))
names(boot_cache) <- names(fits)

for (anchor in names(fits)) {
  fit <- fits[[anchor]]
  free_regimes <- setdiff(names(fit$sigma$regime), fit$reference_regime)
  boot_df <- bootstrap_free_params(fit, nboot = nboot)
  boot_cache[[anchor]] <- boot_df
  for (regime in names(fit$sigma$regime)) {
    cov_stats <- summarise_regime_cov(fit$sigma$regime[[regime]])
    cov_rows[[length(cov_rows) + 1L]] <- cbind(anchor = anchor, regime = regime, cov_stats, stringsAsFactors = FALSE)
  }
  for (regime in free_regimes) {
    scale_profile <- profile_interval(fit, regime = regime, parameter = "scale")
    kappa_profile <- profile_interval(fit, regime = regime, parameter = "kappa")
    scale_boot <- bootstrap_summary(boot_df, regime = regime, parameter = "scale")
    kappa_boot <- bootstrap_summary(boot_df, regime = regime, parameter = "kappa")
    boot_sub <- boot_df[boot_df$regime == regime, , drop = FALSE]
    param_rows[[length(param_rows) + 1L]] <- data.frame(
      anchor = anchor,
      regime = regime,
      scale_est = unname(fit$param[paste0(regime, ".scale")]),
      scale_profile_low = scale_profile$low,
      scale_profile_high = scale_profile$high,
      scale_boot_low = scale_boot$low,
      scale_boot_high = scale_boot$high,
      kappa_est = unname(fit$param[paste0(regime, ".kappa")]),
      kappa_profile_low = kappa_profile$low,
      kappa_profile_high = kappa_profile$high,
      kappa_boot_low = kappa_boot$low,
      kappa_boot_high = kappa_boot$high,
      kappa_boot_boundary_rate = if (nrow(boot_sub)) mean(boot_sub$kappa_boundary) else NA_real_,
      scale_boot_runaway_rate = if (nrow(boot_sub)) mean(boot_sub$runaway_scale) else NA_real_,
      boot_success = scale_boot$n,
      stringsAsFactors = FALSE
    )
  }
}

param_summary <- do.call(rbind, param_rows)
cov_summary <- do.call(rbind, cov_rows)
rownames(param_summary) <- NULL
rownames(cov_summary) <- NULL

pairwise_rows <- list()
anchor_pairs <- utils::combn(names(fits), 2L, simplify = FALSE)
for (pair in anchor_pairs) {
  left <- fits[[pair[1]]]
  right <- fits[[pair[2]]]
  for (regime in names(left$sigma$regime)) {
    delta <- left$sigma$regime[[regime]] - right$sigma$regime[[regime]]
    pairwise_rows[[length(pairwise_rows) + 1L]] <- data.frame(
      anchor_left = pair[1],
      anchor_right = pair[2],
      regime = regime,
      max_abs_diff = max(abs(delta)),
      frobenius_diff = sqrt(sum(delta^2)),
      stringsAsFactors = FALSE
    )
  }
}
pairwise_cov <- do.call(rbind, pairwise_rows)
rownames(pairwise_cov) <- NULL

fmt <- function(x, digits = 3) {
  if (is.logical(x)) {
    return(ifelse(is.na(x), "NA", ifelse(x, "TRUE", "FALSE")))
  }
  if (is.character(x)) {
    return(ifelse(is.na(x), "NA", x))
  }
  x <- as.numeric(x)
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

cat("Running corr-strength anchor-sensitivity workflow\n")
cat("Bootstrap replicates per anchor:", nboot, "\n")
cat("Profile points per parameter:", profile_points, "\n\n")

anchor_print <- anchor_summary
anchor_print$logLik <- fmt(anchor_print$logLik, 3)
anchor_print$AIC <- fmt(anchor_print$AIC, 3)
anchor_print$GIC <- fmt(anchor_print$GIC, 3)
anchor_print$max_scale <- fmt(anchor_print$max_scale, 2)
if ("boundary_kappa" %in% names(anchor_print)) anchor_print$boundary_kappa <- fmt(anchor_print$boundary_kappa, 2)
cat("Anchor-level fit summary\n")
print(anchor_print, row.names = FALSE)

param_print <- param_summary
for (nm in setdiff(names(param_print), c("anchor", "regime", "boot_success"))) {
  param_print[[nm]] <- fmt(param_print[[nm]], 3)
}
cat("\nFree-regime parameter summary with profile and bootstrap intervals\n")
print(param_print, row.names = FALSE)

cov_print <- cov_summary
cov_print$mean_var <- fmt(cov_print$mean_var, 3)
cov_print$mean_cov <- fmt(cov_print$mean_cov, 3)
cov_print$frobenius <- fmt(cov_print$frobenius, 3)
cat("\nRegime covariance summaries by anchor\n")
print(cov_print, row.names = FALSE)

pairwise_print <- pairwise_cov
pairwise_print$max_abs_diff <- fmt(pairwise_print$max_abs_diff, 4)
pairwise_print$frobenius_diff <- fmt(pairwise_print$frobenius_diff, 4)
cat("\nPairwise covariance sensitivity across anchors\n")
print(pairwise_print, row.names = FALSE)

worst_pair <- pairwise_cov[which.max(pairwise_cov$frobenius_diff), , drop = FALSE]
cat("\nLargest anchor-to-anchor covariance deviation\n")
print(worst_pair, row.names = FALSE)

loglik_spread <- diff(range(anchor_summary$logLik))
accept_loglik <- is.finite(loglik_spread) && loglik_spread < 5
accept_pathology <- !any(as.logical(anchor_summary$pathological_scale), na.rm = TRUE)

cat("\nAcceptance checks\n")
cat("logLik spread < 5:", accept_loglik, "(spread =", fmt(loglik_spread, 3), ")\n")
cat("no pathological scale flags:", accept_pathology, "\n")
if (!(accept_loglik && accept_pathology)) {
  warning("corr-strength anchor-sensitivity acceptance checks were not met for this seeded scenario")
}

cat("\ncorr-strength anchor-sensitivity workflow completed\n")
