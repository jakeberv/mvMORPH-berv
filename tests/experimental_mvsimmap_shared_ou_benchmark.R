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

  install_lib <- Sys.getenv("MVSIMMAP_SHARED_BENCH_R_LIB", unset = "")
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
      repos = Sys.getenv("MVSIMMAP_SHARED_BENCH_CRAN_REPO", unset = "https://cloud.r-project.org"),
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

make_single_regime_simmap <- function(tree, regime = "A") {
  tree <- ape::reorder.phylo(tree, "cladewise")
  tree$maps <- lapply(tree$edge.length, function(x) structure(x, names = regime))
  tree$mapped.edge <- matrix(tree$edge.length, ncol = 1, dimnames = list(NULL, regime))
  class(tree) <- c("simmap", setdiff(class(tree), "simmap"))
  tree
}

drop_simmap_mapping <- function(tree) {
  base_tree <- ape::reorder.phylo(tree, order = "cladewise")
  class(base_tree) <- setdiff(class(base_tree), "simmap")
  base_tree$maps <- NULL
  base_tree$mapped.edge <- NULL
  base_tree
}

node_depths <- function(tree) {
  ape::node.depth.edgelength(tree)
}

node_tip_labels <- function(tree, node) {
  ape::extract.clade(tree, node)$tip.label
}

choose_regime_nodes <- function(tree,
                                target_B = 28L,
                                target_C = 18L,
                                min_B = 25L,
                                max_B = 30L,
                                min_C = 15L,
                                max_C = 20L) {
  internal_nodes <- seq_len(tree$Nnode) + ape::Ntip(tree)
  root_node <- ape::Ntip(tree) + 1L
  candidates <- setdiff(internal_nodes, root_node)
  depths <- node_depths(tree)

  info <- lapply(candidates, function(node) {
    tips <- node_tip_labels(tree, node)
    list(node = node, size = length(tips), depth = depths[node], tips = tips)
  })

  B_candidates <- Filter(function(x) x$size >= min_B && x$size <= max_B, info)
  if (!length(B_candidates)) {
    return(NULL)
  }

  score_B <- vapply(B_candidates, function(x) abs(x$size - target_B) + 0.25 * x$depth, numeric(1))
  B <- B_candidates[[which.min(score_B)]]

  C_candidates <- Filter(function(x) {
    x$size >= min_C &&
      x$size <= max_C &&
      length(intersect(x$tips, B$tips)) == 0L &&
      x$depth > B$depth
  }, info)
  if (!length(C_candidates)) {
    C_candidates <- Filter(function(x) {
      x$size >= min_C &&
        x$size <= max_C &&
        length(intersect(x$tips, B$tips)) == 0L
    }, info)
  }
  if (!length(C_candidates)) {
    return(NULL)
  }

  score_C <- vapply(C_candidates, function(x) abs(x$size - target_C) - 0.25 * x$depth, numeric(1))
  C <- C_candidates[[which.min(score_C)]]

  list(B = B, C = C)
}

make_three_regime_simmap <- function(seed = 20260330L,
                                     n_tips = 100L,
                                     target_B = 28L,
                                     target_C = 18L,
                                     min_B = 25L,
                                     max_B = 30L,
                                     min_C = 15L,
                                     max_C = 20L,
                                     max_seed_tries = 500L) {
  for (offset in 0:(max_seed_tries - 1L)) {
    current_seed <- seed + offset
    set.seed(current_seed)
    base_tree <- phytools::pbtree(n = n_tips, scale = 1)
    nodes <- choose_regime_nodes(
      base_tree,
      target_B = target_B,
      target_C = target_C,
      min_B = min_B,
      max_B = max_B,
      min_C = min_C,
      max_C = max_C
    )
    if (is.null(nodes)) next

    simmap <- make_single_regime_simmap(base_tree, regime = "A")
    simmap <- phytools::paintSubTree(simmap, node = nodes$B$node, state = "B", anc.state = "A", stem = TRUE)
    simmap <- phytools::paintSubTree(simmap, node = nodes$C$node, state = "C", anc.state = "A", stem = TRUE)
    return(list(
      tree = simmap,
      tree_seed = current_seed,
      node_B = nodes$B$node,
      node_C = nodes$C$node,
      size_B = nodes$B$size,
      size_C = nodes$C$size,
      depth_B = nodes$B$depth,
      depth_C = nodes$C$depth
    ))
  }
  stop("Could not find a suitable three-regime tree within the allowed seed search", call. = FALSE)
}

make_true_param <- function() {
  trait_names <- c("y1", "y2", "y3")
  theta <- rbind(
    root = c(0.0, 0.0, 0.0),
    "theta.B" = c(1.5, -1.0, 0.8),
    "theta.C" = c(-1.0, 1.2, 0.5)
  )
  colnames(theta) <- trait_names

  sigma_A <- matrix(
    c(0.050, 0.015, 0.010,
      0.015, 0.040, 0.012,
      0.010, 0.012, 0.030),
    3, 3, byrow = TRUE,
    dimnames = list(trait_names, trait_names)
  )
  sigma_shared <- matrix(
    c(0.035, 0.010, 0.006,
      0.010, 0.025, 0.008,
      0.006, 0.008, 0.020),
    3, 3, byrow = TRUE,
    dimnames = list(trait_names, trait_names)
  )
  alpha_shared <- diag(c(3.0, 2.0, 1.5), 3)
  dimnames(alpha_shared) <- list(trait_names, trait_names)

  list(
    theta = theta,
    sigma = list(A = sigma_A, ou_shared = sigma_shared),
    alpha = list(ou_shared = alpha_shared)
  )
}

fit_bm_baseline <- function(tree, data, model, method) {
  fit_primary <- mvBM(
    tree, data,
    model = model,
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
    model = model,
    method = method,
    optimization = "subplex",
    control = list(maxit = 20000),
    diagnostic = FALSE,
    echo = FALSE
  )
  fit_retry$benchmark.optimizer <- "subplex"
  if (fit_retry$AIC <= fit_primary$AIC) return(fit_retry)
  fit_primary$benchmark.optimizer <- "L-BFGS-B"
  fit_primary
}

fit_shared_model <- function(tree, data, method, decomp) {
  mvSIMMAP(
    tree, data,
    process = c(A = "BM", B = "OUM", C = "OUM"),
    process.groups = c(A = "A", B = "ou_shared", C = "ou_shared"),
    method = method,
    optimization = "L-BFGS-B",
    param = list(decomp = decomp, decompSigma = "cholesky"),
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
}

fit_separate_model <- function(tree, data, method, decomp) {
  mvSIMMAP(
    tree, data,
    process = c(A = "BM", B = "OUM", C = "OUM"),
    method = method,
    optimization = "L-BFGS-B",
    param = list(decomp = decomp, decompSigma = "cholesky"),
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
}

flatten_shared_fit <- function(fit, rep_id, seed) {
  out <- data.frame(rep = rep_id, seed = seed, stringsAsFactors = FALSE)
  for (row_name in rownames(fit$theta)) {
    for (trait_name in colnames(fit$theta)) {
      out[[paste0(gsub("\\.", "_", row_name), "_", trait_name)]] <- fit$theta[row_name, trait_name]
    }
  }

  for (matrix_name in c("A", "ou_shared")) {
    sigma <- fit$sigma[[matrix_name]]
    for (i in seq_len(nrow(sigma))) {
      for (j in i:ncol(sigma)) {
        out[[sprintf("sigma_%s_%d%d", matrix_name, i, j)]] <- sigma[i, j]
      }
    }
  }

  alpha <- fit$alpha$ou_shared
  for (i in seq_len(nrow(alpha))) {
    out[[sprintf("alpha_ou_shared_%d%d", i, i)]] <- alpha[i, i]
  }

  out$AIC_shared <- fit$AIC
  out$AICc_shared <- fit$AICc
  out$logLik_shared <- as.numeric(logLik(fit))
  out$nparam_shared <- fit$param$nparam
  out$shared_convergence <- fit$convergence
  out$shared_hessian <- fit$diagnostics$hessian$status
  out
}

truth_vector <- function(true_param) {
  c(
    root_y1 = true_param$theta["root", "y1"],
    root_y2 = true_param$theta["root", "y2"],
    root_y3 = true_param$theta["root", "y3"],
    theta_B_y1 = true_param$theta["theta.B", "y1"],
    theta_B_y2 = true_param$theta["theta.B", "y2"],
    theta_B_y3 = true_param$theta["theta.B", "y3"],
    theta_C_y1 = true_param$theta["theta.C", "y1"],
    theta_C_y2 = true_param$theta["theta.C", "y2"],
    theta_C_y3 = true_param$theta["theta.C", "y3"],
    sigma_A_11 = true_param$sigma$A[1, 1],
    sigma_A_12 = true_param$sigma$A[1, 2],
    sigma_A_13 = true_param$sigma$A[1, 3],
    sigma_A_22 = true_param$sigma$A[2, 2],
    sigma_A_23 = true_param$sigma$A[2, 3],
    sigma_A_33 = true_param$sigma$A[3, 3],
    sigma_ou_shared_11 = true_param$sigma$ou_shared[1, 1],
    sigma_ou_shared_12 = true_param$sigma$ou_shared[1, 2],
    sigma_ou_shared_13 = true_param$sigma$ou_shared[1, 3],
    sigma_ou_shared_22 = true_param$sigma$ou_shared[2, 2],
    sigma_ou_shared_23 = true_param$sigma$ou_shared[2, 3],
    sigma_ou_shared_33 = true_param$sigma$ou_shared[3, 3],
    alpha_ou_shared_11 = true_param$alpha$ou_shared[1, 1],
    alpha_ou_shared_22 = true_param$alpha$ou_shared[2, 2],
    alpha_ou_shared_33 = true_param$alpha$ou_shared[3, 3]
  )
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
      stringsAsFactors = FALSE
    )
  }))
}

summarize_model_comparison <- function(results) {
  model_cols <- c("AIC_shared", "AIC_separate", "AIC_BMM", "AIC_BM1")
  wins <- vapply(seq_len(nrow(results)), function(i) {
    model_cols[which.min(unlist(results[i, model_cols]))]
  }, character(1))
  weights <- t(apply(results[, model_cols], 1, function(x) {
    delta <- x - min(x)
    out <- exp(-0.5 * delta)
    out / sum(out)
  }))
  colnames(weights) <- c("weight_shared", "weight_separate", "weight_BMM", "weight_BM1")

  data.frame(
    metric = c(
      "mean_AIC_shared",
      "mean_AIC_separate",
      "mean_AIC_BMM",
      "mean_AIC_BM1",
      "mean_delta_AIC_separate_minus_shared",
      "mean_delta_AIC_BMM_minus_shared",
      "mean_delta_AIC_BM1_minus_shared",
      "mean_weight_shared",
      "mean_weight_separate",
      "mean_weight_BMM",
      "mean_weight_BM1",
      "shared_wins",
      "shared_wins_delta_gt_2_vs_separate",
      "shared_wins_delta_gt_10_vs_BM1"
    ),
    value = c(
      mean(results$AIC_shared),
      mean(results$AIC_separate),
      mean(results$AIC_BMM),
      mean(results$AIC_BM1),
      mean(results$AIC_separate - results$AIC_shared),
      mean(results$AIC_BMM - results$AIC_shared),
      mean(results$AIC_BM1 - results$AIC_shared),
      mean(weights[, "weight_shared"]),
      mean(weights[, "weight_separate"]),
      mean(weights[, "weight_BMM"]),
      mean(weights[, "weight_BM1"]),
      sum(wins == "AIC_shared"),
      sum((results$AIC_separate - results$AIC_shared) > 2),
      sum((results$AIC_BM1 - results$AIC_shared) > 10)
    )
  )
}

save_outputs <- function(output_dir, results, recovery, model_comparison, metadata) {
  if (!nzchar(output_dir)) return(invisible(NULL))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(results, file.path(output_dir, "replicate_results.csv"), row.names = FALSE)
  utils::write.csv(recovery, file.path(output_dir, "recovery_summary.csv"), row.names = FALSE)
  utils::write.csv(model_comparison, file.path(output_dir, "model_comparison_summary.csv"), row.names = FALSE)
  writeLines(metadata, file.path(output_dir, "run_metadata.txt"))
}

bootstrap_runtime(repo_root)
clean_staged_artifacts(repo_root)
pkgload::load_all(repo_root, quiet = TRUE, compile = TRUE, recompile = TRUE)

tree_seed <- env_int("MVSIMMAP_SHARED_BENCH_TREE_SEED", 20260330L)
n_tips <- env_int("MVSIMMAP_SHARED_BENCH_NTIPS", 100L)
nsim <- env_int("MVSIMMAP_SHARED_BENCH_REPS", 5L)
method <- env_chr("MVSIMMAP_SHARED_BENCH_METHOD", "rpf")
alpha_decomp <- env_chr("MVSIMMAP_SHARED_BENCH_ALPHA_DECOMP", "scalarPositive")
output_dir <- env_chr("MVSIMMAP_SHARED_BENCH_OUTPUT_DIR", "")
ncores <- env_int("MVSIMMAP_SHARED_BENCH_CORES", 4L)

if (!method %in% c("rpf", "inverse", "pseudoinverse")) {
  stop("MVSIMMAP_SHARED_BENCH_METHOD must be one of rpf, inverse, pseudoinverse", call. = FALSE)
}
if (!alpha_decomp %in% c("diagonalPositive", "scalarPositive")) {
  stop("MVSIMMAP_SHARED_BENCH_ALPHA_DECOMP must be one of diagonalPositive or scalarPositive", call. = FALSE)
}

benchmark <- make_three_regime_simmap(seed = tree_seed, n_tips = n_tips)
base_tree <- drop_simmap_mapping(benchmark$tree)
true_param <- make_true_param()
seeds <- seq.int(2001L, length.out = nsim)
ncores <- min(ncores, parallel::detectCores(logical = FALSE), nsim)

scaffold <- mvSIMMAP(
  benchmark$tree,
  data = NULL,
  process = c(A = "BM", B = "OUM", C = "OUM"),
  process.groups = c(A = "A", B = "ou_shared", C = "ou_shared"),
  optimization = "fixed",
  method = method,
  param = list(
    ntraits = 3,
    names_traits = c("y1", "y2", "y3"),
    decomp = alpha_decomp,
    decompSigma = "cholesky"
  ),
  echo = FALSE
)

Ylist <- lapply(seeds, function(seed) simulate(scaffold, nsim = 1, seed = seed, param = true_param))

run_one <- function(i) {
  Y <- Ylist[[i]]
  fit_shared <- fit_shared_model(benchmark$tree, Y, method = method, decomp = alpha_decomp)
  fit_separate <- fit_separate_model(benchmark$tree, Y, method = method, decomp = alpha_decomp)
  fit_bmm <- fit_bm_baseline(benchmark$tree, Y, model = "BMM", method = method)
  fit_bm1 <- fit_bm_baseline(base_tree, Y, model = "BM1", method = method)

  out <- flatten_shared_fit(fit_shared, rep_id = i, seed = seeds[i])
  out$AIC_separate <- fit_separate$AIC
  out$AICc_separate <- fit_separate$AICc
  out$logLik_separate <- as.numeric(logLik(fit_separate))
  out$nparam_separate <- fit_separate$param$nparam
  out$separate_convergence <- fit_separate$convergence
  out$separate_hessian <- fit_separate$diagnostics$hessian$status
  out$AIC_BMM <- fit_bmm$AIC
  out$AICc_BMM <- fit_bmm$AICc
  out$logLik_BMM <- as.numeric(logLik(fit_bmm))
  out$nparam_BMM <- fit_bmm$param$nparam
  out$BMM_convergence <- fit_bmm$convergence
  out$BMM_optimizer <- fit_bmm$benchmark.optimizer %||% "L-BFGS-B"
  out$AIC_BM1 <- fit_bm1$AIC
  out$AICc_BM1 <- fit_bm1$AICc
  out$logLik_BM1 <- as.numeric(logLik(fit_bm1))
  out$nparam_BM1 <- fit_bm1$param$nparam
  out$BM1_convergence <- fit_bm1$convergence
  out$BM1_optimizer <- fit_bm1$benchmark.optimizer %||% "L-BFGS-B"
  out$delta_AIC_separate_minus_shared <- fit_separate$AIC - fit_shared$AIC
  out$delta_AIC_BMM_minus_shared <- fit_bmm$AIC - fit_shared$AIC
  out$delta_AIC_BM1_minus_shared <- fit_bm1$AIC - fit_shared$AIC
  out
}

timing <- system.time({
  results <- do.call(
    rbind,
    parallel::mclapply(seq_along(seeds), run_one, mc.cores = ncores, mc.set.seed = FALSE)
  )
})

recovery <- summarize_recovery(results, truth_vector(true_param))
model_comparison <- summarize_model_comparison(results)

metadata <- c(
  paste("tree_seed:", benchmark$tree_seed),
  paste("n_tips:", n_tips),
  paste("nsim:", nsim),
  paste("method:", method),
  paste("alpha_decomp:", alpha_decomp),
  paste("node_B:", benchmark$node_B),
  paste("node_C:", benchmark$node_C),
  paste("size_B:", benchmark$size_B),
  paste("size_C:", benchmark$size_C),
  paste("ncores:", ncores),
  paste("elapsed_sec:", unname(timing[["elapsed"]]))
)

save_outputs(output_dir, results, recovery, model_comparison, metadata)

cat("TREE_SUMMARY\n")
cat("tree_seed", benchmark$tree_seed, "\n")
cat("total_tips", ape::Ntip(benchmark$tree), "\n")
cat("tree_height", max(node_depths(benchmark$tree)), "\n")
cat("node_B", benchmark$node_B, "\n")
cat("node_C", benchmark$node_C, "\n")
cat("size_B", benchmark$size_B, "\n")
cat("size_C", benchmark$size_C, "\n")
cat("depth_B", benchmark$depth_B, "\n")
cat("depth_C", benchmark$depth_C, "\n")
cat("method", method, "\n")
cat("alpha_decomp", alpha_decomp, "\n")
cat("ncores", ncores, "\n")
cat("elapsed_sec", unname(timing[["elapsed"]]), "\n")

cat("\nTRUE_THETA\n")
print(true_param$theta)

cat("\nTRUE_SIGMA_A\n")
print(true_param$sigma$A)

cat("\nTRUE_SIGMA_OU_SHARED\n")
print(true_param$sigma$ou_shared)

cat("\nTRUE_ALPHA_OU_SHARED\n")
print(true_param$alpha$ou_shared)

cat("\nMODEL_COMPARISON_BY_REPLICATE\n")
print(
  results[, c(
    "rep",
    "AIC_shared",
    "AIC_separate",
    "AIC_BMM",
    "AIC_BM1",
    "delta_AIC_separate_minus_shared",
    "delta_AIC_BMM_minus_shared",
    "delta_AIC_BM1_minus_shared",
    "shared_convergence",
    "separate_convergence",
    "shared_hessian",
    "separate_hessian"
  )],
  row.names = FALSE
)

cat("\nMODEL_COMPARISON_SUMMARY\n")
print(model_comparison, row.names = FALSE)

cat("\nRECOVERY_SUMMARY\n")
print(recovery, row.names = FALSE)

cat("\nSHARED_HESSIAN_TABLE\n")
print(table(results$shared_hessian))

cat("\nSEPARATE_HESSIAN_TABLE\n")
print(table(results$separate_hessian))

if (nzchar(output_dir)) {
  cat("\nOUTPUT_DIR\n")
  cat(normalizePath(output_dir), "\n")
}
