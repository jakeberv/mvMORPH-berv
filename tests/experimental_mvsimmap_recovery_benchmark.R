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

  install_lib <- Sys.getenv("MVSIMMAP_BENCH_R_LIB", unset = "")
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
      repos = Sys.getenv("MVSIMMAP_BENCH_CRAN_REPO", unset = "https://cloud.r-project.org"),
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

`%||%` <- function(x, y) if (is.null(x)) y else x

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed <= 0L) return(default)
  parsed
}

env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  value
}

normalize_trait_names <- function(ntraits) {
  paste0("y", seq_len(ntraits))
}

drop_simmap_mapping <- function(tree) {
  base_tree <- ape::reorder.phylo(tree, order = "cladewise")
  class(base_tree) <- setdiff(class(base_tree), "simmap")
  base_tree$maps <- NULL
  base_tree$mapped.edge <- NULL
  base_tree
}

make_benchmark_simmap <- function(seed, n_tips, min_derived_tips, target_derived_tips) {
  set.seed(seed)
  tree <- phytools::pbtree(n = n_tips, scale = 1)
  internal_nodes <- seq_len(tree$Nnode) + ape::Ntip(tree)
  root_node <- ape::Ntip(tree) + 1L
  candidates <- setdiff(internal_nodes, root_node)
  sizes <- vapply(candidates, function(node) ape::Ntip(ape::extract.clade(tree, node)), integer(1))
  eligible <- candidates[sizes >= min_derived_tips & sizes < n_tips]
  if (!length(eligible)) {
    stop("No eligible derived subtree found for the requested benchmark setup", call. = FALSE)
  }
  eligible_sizes <- sizes[sizes >= min_derived_tips & sizes < n_tips]
  node <- eligible[which.min(abs(eligible_sizes - target_derived_tips))]
  simmap <- phytools::paintSubTree(tree, node = node, state = "B", anc.state = "A", stem = TRUE)
  list(
    tree = simmap,
    shift_node = node,
    derived_tips = sum(phytools::getStates(simmap, "tips") == "B")
  )
}

make_true_param <- function(ntraits) {
  if (ntraits < 1L || ntraits > 3L) {
    stop("This benchmark currently supports between 1 and 3 traits", call. = FALSE)
  }
  trait_names <- normalize_trait_names(ntraits)
  theta_B_base <- c(1.5, -1.0, 0.75)
  sigma_base <- matrix(
    c(0.040, 0.012, 0.008,
      0.012, 0.030, 0.009,
      0.008, 0.009, 0.025),
    3, 3, byrow = TRUE
  )
  alpha_base <- diag(c(3.0, 2.5, 2.0), 3)
  sigma <- sigma_base[seq_len(ntraits), seq_len(ntraits), drop = FALSE]
  alpha <- alpha_base[seq_len(ntraits), seq_len(ntraits), drop = FALSE]
  dimnames(sigma) <- list(trait_names, trait_names)
  dimnames(alpha) <- list(trait_names, trait_names)
  theta <- rbind(root = rep(0, ntraits), "theta.B" = theta_B_base[seq_len(ntraits)])
  colnames(theta) <- trait_names
  list(
    theta = theta,
    sigma = list(A = sigma, B = sigma),
    alpha = list(B = alpha)
  )
}

fit_bm1_baseline <- function(tree, data, method) {
  fit_primary <- mvBM(
    tree, data,
    model = "BM1",
    method = method,
    optimization = "L-BFGS-B",
    control = list(maxit = 20000),
    diagnostic = FALSE,
    echo = FALSE
  )
  if (isTRUE(fit_primary$convergence == 0)) {
    fit_primary$benchmark.optimizer <- "L-BFGS-B"
    return(fit_primary)
  }

  fit_retry <- mvBM(
    tree, data,
    model = "BM1",
    method = method,
    optimization = "subplex",
    control = list(maxit = 20000),
    diagnostic = FALSE,
    echo = FALSE
  )
  fit_retry$benchmark.optimizer <- "subplex"
  if (fit_retry$AIC <= fit_primary$AIC) {
    return(fit_retry)
  }
  fit_primary$benchmark.optimizer <- "L-BFGS-B"
  fit_primary
}

flatten_fit <- function(fit, fit_bm1, seeds, i, ntraits) {
  trait_names <- normalize_trait_names(ntraits)
  out <- data.frame(rep = i, seed = seeds[i], stringsAsFactors = FALSE)

  for (trait in trait_names) {
    out[[paste0("root_", trait)]] <- fit$theta["root", trait]
    out[[paste0("thetaB_", trait)]] <- fit$theta["theta.B", trait]
  }

  for (group in c("A", "B")) {
    sigma <- fit$sigma[[group]]
    for (row in seq_len(ntraits)) {
      for (col in row:ntraits) {
        out[[sprintf("sigma%s_%d%d", group, row, col)]] <- sigma[row, col]
      }
    }
  }

  alpha <- fit$alpha$B
  for (idx in seq_len(ntraits)) {
    out[[sprintf("alphaB_%d%d", idx, idx)]] <- alpha[idx, idx]
  }

  out$logLik <- as.numeric(logLik(fit))
  out$convergence <- fit$convergence
  out$hessian <- fit$diagnostics$hessian$status
  out$AIC_mixed <- fit$AIC
  out$AICc_mixed <- fit$AICc
  out$nparam_mixed <- fit$param$nparam
  out$AIC_BM1 <- fit_bm1$AIC
  out$AICc_BM1 <- fit_bm1$AICc
  out$nparam_BM1 <- fit_bm1$param$nparam
  out$logLik_BM1 <- as.numeric(logLik(fit_bm1))
  out$bm1_convergence <- fit_bm1$convergence
  out$bm1_optimizer <- fit_bm1$benchmark.optimizer %||% "L-BFGS-B"
  out$delta_AIC_BM1_minus_mixed <- fit_bm1$AIC - fit$AIC
  out$delta_AICc_BM1_minus_mixed <- fit_bm1$AICc - fit$AICc
  delta <- c(mixed = 0, BM1 = out$delta_AIC_BM1_minus_mixed)
  weights <- exp(-0.5 * delta)
  weights <- weights / sum(weights)
  out$weight_mixed <- weights[["mixed"]]
  out$weight_BM1 <- weights[["BM1"]]
  out
}

make_truth_vector <- function(true_param, ntraits) {
  out <- c()
  for (idx in seq_len(ntraits)) out[[paste0("root_y", idx)]] <- 0
  for (idx in seq_len(ntraits)) out[[paste0("thetaB_y", idx)]] <- true_param$theta["theta.B", idx]
  for (group in c("A", "B")) {
    sigma <- true_param$sigma[[group]]
    for (row in seq_len(ntraits)) {
      for (col in row:ntraits) {
        out[[sprintf("sigma%s_%d%d", group, row, col)]] <- sigma[row, col]
      }
    }
  }
  for (idx in seq_len(ntraits)) out[[sprintf("alphaB_%d%d", idx, idx)]] <- true_param$alpha$B[idx, idx]
  out
}

summarize_recovery <- function(results, truth) {
  metric_names <- names(truth)
  do.call(rbind, lapply(metric_names, function(metric) {
    x <- results[[metric]]
    t <- truth[[metric]]
    data.frame(
      parameter = metric,
      truth = t,
      mean_est = mean(x),
      median_est = stats::median(x),
      sd_est = stats::sd(x),
      bias = mean(x - t),
      rmse = sqrt(mean((x - t)^2)),
      min_est = min(x),
      max_est = max(x),
      stringsAsFactors = FALSE
    )
  }))
}

summarize_aic_comparison <- function(results) {
  data.frame(
    metric = c(
      "mean_AIC_mixed",
      "mean_AIC_BM1",
      "mean_delta_AIC_BM1_minus_mixed",
      "median_delta_AIC_BM1_minus_mixed",
      "min_delta_AIC_BM1_minus_mixed",
      "max_delta_AIC_BM1_minus_mixed",
      "mixed_wins",
      "mixed_wins_delta_gt_2",
      "mixed_wins_delta_gt_10",
      "mean_weight_mixed",
      "mean_weight_BM1",
      "mean_logLik_gain_mixed_minus_BM1"
    ),
    value = c(
      mean(results$AIC_mixed),
      mean(results$AIC_BM1),
      mean(results$delta_AIC_BM1_minus_mixed),
      stats::median(results$delta_AIC_BM1_minus_mixed),
      min(results$delta_AIC_BM1_minus_mixed),
      max(results$delta_AIC_BM1_minus_mixed),
      sum(results$delta_AIC_BM1_minus_mixed > 0),
      sum(results$delta_AIC_BM1_minus_mixed > 2),
      sum(results$delta_AIC_BM1_minus_mixed > 10),
      mean(results$weight_mixed),
      mean(results$weight_BM1),
      mean(results$logLik - results$logLik_BM1)
    ),
    stringsAsFactors = FALSE
  )
}

mean_theta_matrix <- function(results, ntraits) {
  out <- rbind(
    root = vapply(seq_len(ntraits), function(i) mean(results[[paste0("root_y", i)]]), numeric(1)),
    "theta.B" = vapply(seq_len(ntraits), function(i) mean(results[[paste0("thetaB_y", i)]]), numeric(1))
  )
  colnames(out) <- normalize_trait_names(ntraits)
  out
}

mean_symmetric_matrix <- function(results, prefix, ntraits) {
  out <- matrix(0, ntraits, ntraits)
  for (row in seq_len(ntraits)) {
    for (col in row:ntraits) {
      value <- mean(results[[sprintf("%s_%d%d", prefix, row, col)]])
      out[row, col] <- value
      out[col, row] <- value
    }
  }
  dimnames(out) <- list(normalize_trait_names(ntraits), normalize_trait_names(ntraits))
  out
}

mean_alpha_matrix <- function(results, ntraits) {
  out <- matrix(0, ntraits, ntraits)
  for (idx in seq_len(ntraits)) {
    out[idx, idx] <- mean(results[[sprintf("alphaB_%d%d", idx, idx)]])
  }
  dimnames(out) <- list(normalize_trait_names(ntraits), normalize_trait_names(ntraits))
  out
}

save_outputs <- function(output_dir, results, recovery, aic_summary, metadata) {
  if (!nzchar(output_dir)) return(invisible(NULL))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(results, file.path(output_dir, "replicate_estimates.csv"), row.names = FALSE)
  utils::write.csv(recovery, file.path(output_dir, "recovery_summary.csv"), row.names = FALSE)
  utils::write.csv(aic_summary, file.path(output_dir, "aic_comparison_summary.csv"), row.names = FALSE)
  writeLines(metadata, file.path(output_dir, "run_metadata.txt"))
}

bootstrap_runtime(repo_root)
clean_staged_artifacts(repo_root)
pkgload::load_all(repo_root, quiet = TRUE, compile = TRUE, recompile = TRUE)

tree_seed <- env_int("MVSIMMAP_BENCH_TREE_SEED", 20260330L)
n_tips <- env_int("MVSIMMAP_BENCH_NTIPS", 50L)
min_derived_tips <- env_int("MVSIMMAP_BENCH_MIN_DERIVED_TIPS", 15L)
target_derived_tips <- env_int("MVSIMMAP_BENCH_TARGET_DERIVED_TIPS", 28L)
ntraits <- env_int("MVSIMMAP_BENCH_NTRAITS", 3L)
nsim <- env_int("MVSIMMAP_BENCH_REPS", 10L)
method <- env_chr("MVSIMMAP_BENCH_METHOD", "rpf")
output_dir <- env_chr("MVSIMMAP_BENCH_OUTPUT_DIR", "")

if (!method %in% c("rpf", "inverse", "pseudoinverse")) {
  stop("MVSIMMAP_BENCH_METHOD must be one of rpf, inverse, pseudoinverse", call. = FALSE)
}

benchmark <- make_benchmark_simmap(
  seed = tree_seed,
  n_tips = n_tips,
  min_derived_tips = min_derived_tips,
  target_derived_tips = target_derived_tips
)

process <- c(A = "BM", B = "OU")
true_param <- make_true_param(ntraits)
seeds <- seq.int(1001L, length.out = nsim)
ncores <- min(env_int("MVSIMMAP_BENCH_CORES", 4L), parallel::detectCores(logical = FALSE), nsim)
base_tree <- drop_simmap_mapping(benchmark$tree)

scaffold <- mvSIMMAP(
  benchmark$tree,
  data = NULL,
  process = process,
  optimization = "fixed",
  param = list(
    ntraits = ntraits,
    names_traits = normalize_trait_names(ntraits),
    decomp = "diagonalPositive",
    decompSigma = "cholesky"
  ),
  method = method,
  echo = FALSE
)

Ylist <- lapply(seeds, function(seed) simulate(scaffold, nsim = 1, seed = seed, param = true_param))

run_one <- function(i) {
  fit <- mvSIMMAP(
    benchmark$tree, Ylist[[i]],
    process = process,
    method = method,
    optimization = "L-BFGS-B",
    param = list(decomp = "diagonalPositive", decompSigma = "cholesky"),
    control = list(
      maxit = 500,
      retry.unreliable = TRUE,
      retry.max = 2,
      retry.jitter = 0.05,
      retry.seed = 1
    ),
    diagnostic = FALSE,
    echo = FALSE
  )
  fit_bm1 <- fit_bm1_baseline(base_tree, Ylist[[i]], method = method)
  flatten_fit(fit, fit_bm1, seeds = seeds, i = i, ntraits = ntraits)
}

timing <- system.time({
  results <- do.call(rbind, parallel::mclapply(seq_along(seeds), run_one, mc.cores = ncores, mc.set.seed = FALSE))
})

truth <- make_truth_vector(true_param, ntraits)
recovery <- summarize_recovery(results, truth)
aic_summary <- summarize_aic_comparison(results)
mean_theta <- mean_theta_matrix(results, ntraits)
mean_sigma_A <- mean_symmetric_matrix(results, "sigmaA", ntraits)
mean_sigma_B <- mean_symmetric_matrix(results, "sigmaB", ntraits)
mean_alpha_B <- mean_alpha_matrix(results, ntraits)
depths <- ape::node.depth.edgelength(benchmark$tree)

metadata <- c(
  paste("tree_seed:", tree_seed),
  paste("n_tips:", n_tips),
  paste("ntraits:", ntraits),
  paste("nsim:", nsim),
  paste("method:", method),
  paste("shift_node:", benchmark$shift_node),
  paste("derived_tips:", benchmark$derived_tips),
  paste("ncores:", ncores),
  paste("elapsed_sec:", unname(timing[["elapsed"]])),
  paste("mean_delta_AIC_BM1_minus_mixed:", mean(results$delta_AIC_BM1_minus_mixed))
)

save_outputs(output_dir, results, recovery, aic_summary, metadata)

cat("TREE_SUMMARY\n")
cat("tree_seed", tree_seed, "\n")
cat("total_tips", ape::Ntip(benchmark$tree), "\n")
cat("tree_height", max(depths), "\n")
cat("shift_node", benchmark$shift_node, "\n")
cat("derived_tips", benchmark$derived_tips, "\n")
cat("regime_B_root_depth", depths[benchmark$shift_node], "\n")
cat("regime_B_remaining_height", max(depths) - depths[benchmark$shift_node], "\n")
cat("method", method, "\n")
cat("ncores", ncores, "\n")
cat("elapsed_sec", unname(timing[["elapsed"]]), "\n")

cat("\nTRUE_THETA\n")
print(true_param$theta)

cat("\nMEAN_EST_THETA\n")
print(mean_theta)

cat("\nTRUE_SIGMA_A\n")
print(true_param$sigma$A)

cat("\nMEAN_EST_SIGMA_A\n")
print(mean_sigma_A)

cat("\nTRUE_SIGMA_B\n")
print(true_param$sigma$B)

cat("\nMEAN_EST_SIGMA_B\n")
print(mean_sigma_B)

cat("\nTRUE_ALPHA_B\n")
print(true_param$alpha$B)

cat("\nMEAN_EST_ALPHA_B\n")
print(mean_alpha_B)

cat("\nRECOVERY_SUMMARY\n")
print(recovery, row.names = FALSE)

cat("\nCONVERGENCE_TABLE\n")
print(table(results$convergence))

cat("\nHESSIAN_TABLE\n")
print(table(results$hessian))

cat("\nBM1_CONVERGENCE_TABLE\n")
print(table(results$bm1_convergence))

cat("\nBM1_OPTIMIZER_TABLE\n")
print(table(results$bm1_optimizer))

cat("\nAIC_COMPARISON\n")
print(
  results[, c(
    "rep",
    "AIC_mixed",
    "AIC_BM1",
    "delta_AIC_BM1_minus_mixed",
    "weight_mixed",
    "convergence",
    "bm1_convergence",
    "hessian",
    "bm1_optimizer"
  )],
  row.names = FALSE
)

cat("\nAIC_COMPARISON_SUMMARY\n")
print(aic_summary, row.names = FALSE)

if (nzchar(output_dir)) {
  cat("\nOUTPUT_DIR\n")
  cat(normalizePath(output_dir), "\n")
}
