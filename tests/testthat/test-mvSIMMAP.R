make_single_regime_simmap <- function(tree, regime="A"){
    tree <- ape::reorder.phylo(tree, "cladewise")
    tree$maps <- lapply(tree$edge.length, function(x) structure(x, names=regime))
    tree$mapped.edge <- matrix(tree$edge.length, ncol=1, dimnames=list(NULL, regime))
    class(tree) <- c("simmap", setdiff(class(tree), "simmap"))
    tree
}

make_two_regime_simmap <- function(seed=1, n=8){
    set.seed(seed)
    repeat{
        tree <- ape::rtree(n)
        states <- setNames(rep(c("A", "B"), length.out=n), tree$tip.label)
        simmap <- suppressMessages(phytools::make.simmap(tree, states, model="ER", nsim=1))
        if(identical(sort(colnames(simmap$mapped.edge)), c("A", "B"))){
            return(simmap)
        }
    }
}

make_bm_ou_benchmark_simmap <- function(seed=20260330, n=50, min_derived_tips=15, target_derived_tips=28){
    set.seed(seed)
    tree <- phytools::pbtree(n=n, scale=1)
    internal_nodes <- seq_len(tree$Nnode) + ape::Ntip(tree)
    root_node <- ape::Ntip(tree) + 1L
    candidates <- setdiff(internal_nodes, root_node)
    sizes <- sapply(candidates, function(node) ape::Ntip(ape::extract.clade(tree, node)))
    eligible <- candidates[sizes >= min_derived_tips & sizes < n]
    if(!length(eligible)){
        stop("No eligible derived subtree found for the requested benchmark setup", call.=FALSE)
    }
    eligible_sizes <- sizes[sizes >= min_derived_tips & sizes < n]
    node <- eligible[which.min(abs(eligible_sizes - target_derived_tips))]
    simmap <- phytools::paintSubTree(tree, node=node, state="B", anc.state="A", stem=TRUE)
    list(
        tree=simmap,
        shift.node=node,
        derived.tips=sum(phytools::getStates(simmap, "tips") == "B")
    )
}

make_bm_ou_true_param <- function(ntraits=3){
    if(ntraits < 1L || ntraits > 3L){
        stop("This benchmark helper currently supports between 1 and 3 traits", call.=FALSE)
    }
    trait_names <- paste0("y", seq_len(ntraits))
    theta_B_base <- c(1.5, -1.0, 0.75)
    sigma_base <- matrix(
        c(0.040, 0.012, 0.008,
          0.012, 0.030, 0.009,
          0.008, 0.009, 0.025),
        3, 3, byrow=TRUE
    )
    alpha_base <- diag(c(3.0, 2.5, 2.0), 3)
    sigma <- sigma_base[seq_len(ntraits), seq_len(ntraits), drop=FALSE]
    alpha <- alpha_base[seq_len(ntraits), seq_len(ntraits), drop=FALSE]
    dimnames(sigma) <- list(trait_names, trait_names)
    dimnames(alpha) <- list(trait_names, trait_names)
    theta <- rbind(root=rep(0, ntraits), "theta.B"=theta_B_base[seq_len(ntraits)])
    colnames(theta) <- trait_names
    list(
        theta=theta,
        sigma=list(A=sigma, B=sigma),
        alpha=list(B=alpha)
    )
}

test_that("mvSIMMAP all-BM special case matches mvBM BMM likelihoods", {
    set.seed(1)
    tree <- ape::rtree(8)
    states <- setNames(sample(c("A", "B"), ape::Ntip(tree), replace=TRUE), tree$tip.label)
    simmap <- suppressMessages(phytools::make.simmap(tree, states, model="ER", nsim=1))
    X <- matrix(rnorm(ape::Ntip(simmap) * 2), ncol=2)
    rownames(X) <- simmap$tip.label

    sigma_A <- crossprod(matrix(c(0.4, 0.1, -0.2, 0.6), 2, 2))
    sigma_B <- crossprod(matrix(c(0.7, -0.1, 0.3, 0.5), 2, 2))
    theta <- c(0.25, -0.15)

    fit_mix <- mvSIMMAP(
        simmap, X,
        process=c(A="BM", B="BM"),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )
    fit_bmm <- mvBM(
        simmap, X,
        model="BMM",
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )

    par_mix <- c(mvMORPH:::sym.unpar(sigma_A), mvMORPH:::sym.unpar(sigma_B))
    par_theta <- c(par_mix, theta)

    expect_equal(
        fit_mix$llik(par_mix, root.mle=TRUE),
        fit_bmm$llik(par_mix, root.mle=TRUE),
        tolerance=1e-6
    )
    expect_equal(
        fit_mix$llik(par_theta, root.mle=FALSE),
        fit_bmm$llik(par_theta, root.mle=FALSE),
        tolerance=1e-6
    )
})

test_that("mvSIMMAP single-regime OU special case matches mvOU likelihoods", {
    set.seed(2)
    tree <- make_single_regime_simmap(ape::rtree(7), regime="A")
    X <- matrix(rnorm(ape::Ntip(tree) * 2), ncol=2)
    rownames(X) <- tree$tip.label

    alpha <- crossprod(matrix(c(0.8, 0.1, 0.0, 0.6), 2, 2))
    sigma <- crossprod(matrix(c(0.5, -0.2, 0.1, 0.4), 2, 2))
    theta <- c(0.1, -0.2, 0.5, 0.3)

    fit_mix <- mvSIMMAP(
        tree, X,
        process=c(A="OU"),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )
    fit_ou <- mvOU(
        tree, X,
        model="OUM",
        optimization="fixed",
        method="inverse",
        param=list(root=TRUE, vcv="fixedRoot"),
        echo=FALSE
    )

    par_process <- c(mvMORPH:::sym.unpar(alpha), mvMORPH:::sym.unpar(sigma))
    par_theta <- c(par_process, theta)

    expect_equal(
        fit_mix$llik(par_process, root.mle=TRUE),
        fit_ou$llik(par_process, root.mle=TRUE),
        tolerance=1e-6
    )
    expect_equal(
        fit_mix$llik(par_theta, root.mle=FALSE),
        fit_ou$llik(par_theta, root.mle=FALSE),
        tolerance=0.1
    )
})

test_that("mvSIMMAP uses separate OU parameter blocks for different painted regimes by default", {
    tree <- make_two_regime_simmap(seed=11, n=8)
    X <- matrix(rnorm(ape::Ntip(tree) * 2), ncol=2)
    rownames(X) <- tree$tip.label

    fit <- mvSIMMAP(
        tree, X,
        process=c(A="OU", B="OU"),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )

    expect_identical(names(fit$alpha), c("A", "B"))
    expect_identical(names(fit$sigma), c("A", "B"))
    expect_equal(fit$param$nprocess.groups, 2)
    expect_identical(fit$param$names_process_groups, c("A", "B"))
    expect_identical(rownames(fit$theta), c("root", "theta.A", "theta.B"))
})

test_that("mvSIMMAP shared OU groups collapse to one optimum for the likelihood", {
    tree <- make_two_regime_simmap(seed=12, n=8)
    X <- matrix(rnorm(ape::Ntip(tree) * 2), ncol=2)
    rownames(X) <- tree$tip.label

    fit <- mvSIMMAP(
        tree, X,
        process=c(A="OU", B="OU"),
        process.groups=c(A="ou_shared", B="ou_shared"),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )

    expect_identical(fit$process.groups, c(A="ou_shared", B="ou_shared"))
    expect_identical(names(fit$alpha), "ou_shared")
    expect_identical(names(fit$sigma), "ou_shared")
    expect_equal(fit$param$nprocess.groups, 1)
    expect_identical(fit$param$names_process_groups, "ou_shared")
    expect_identical(rownames(fit$theta), c("root", "theta.ou_shared"))
})

test_that("mvSIMMAP OUM can share OU dynamics while keeping painted-regime optima separate", {
    tree <- make_two_regime_simmap(seed=15, n=8)
    X <- matrix(rnorm(ape::Ntip(tree) * 2), ncol=2)
    rownames(X) <- tree$tip.label

    fit <- mvSIMMAP(
        tree, X,
        process=c(A="OUM", B="OUM"),
        process.groups=c(A="ou_shared", B="ou_shared"),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )

    expect_identical(fit$process.groups, c(A="ou_shared", B="ou_shared"))
    expect_identical(names(fit$alpha), "ou_shared")
    expect_identical(names(fit$sigma), "ou_shared")
    expect_equal(fit$param$nprocess.groups, 1)
    expect_identical(rownames(fit$theta), c("root", "theta.A", "theta.B"))
})

test_that("mvSIMMAP single-regime EB special case matches mvEB likelihoods", {
    set.seed(3)
    tree <- make_single_regime_simmap(ape::rtree(7), regime="A")
    X <- matrix(rnorm(ape::Ntip(tree) * 2), ncol=2)
    rownames(X) <- tree$tip.label

    sigma <- crossprod(matrix(c(0.6, 0.1, -0.2, 0.5), 2, 2))
    beta <- -0.2
    theta <- c(-0.4, 0.2)

    fit_mix <- mvSIMMAP(
        tree, X,
        process=c(A="EB"),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )
    fit_eb <- mvEB(
        tree, X,
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )

    par_mix <- c(mvMORPH:::sym.unpar(sigma), beta)
    par_mix_theta <- c(par_mix, theta)
    par_eb <- c(beta, mvMORPH:::sym.unpar(sigma))
    par_eb_theta <- c(par_eb, theta)

    expect_equal(
        fit_mix$llik(par_mix, root.mle=TRUE),
        fit_eb$llik(par_eb, root.mle=TRUE),
        tolerance=1e-6
    )
    expect_equal(
        fit_mix$llik(par_mix_theta, root.mle=FALSE),
        fit_eb$llik(par_eb_theta, root.mle=FALSE),
        tolerance=1e-6
    )
})

test_that("mvSIMMAP supports separate and shared EB process groups", {
    tree <- make_two_regime_simmap(seed=13, n=8)
    X <- matrix(rnorm(ape::Ntip(tree) * 2), ncol=2)
    rownames(X) <- tree$tip.label

    fit_sep <- mvSIMMAP(
        tree, X,
        process=c(A="EB", B="EB"),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )
    fit_shared <- mvSIMMAP(
        tree, X,
        process=c(A="EB", B="EB"),
        process.groups=c(A="eb_shared", B="eb_shared"),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )

    expect_identical(names(fit_sep$beta), c("A", "B"))
    expect_identical(names(fit_sep$sigma), c("A", "B"))
    expect_equal(fit_sep$param$nprocess.groups, 2)
    expect_identical(names(fit_shared$beta), "eb_shared")
    expect_identical(names(fit_shared$sigma), "eb_shared")
    expect_equal(fit_shared$param$nprocess.groups, 1)
})

test_that("mvSIMMAP rejects process groups that mix OU, BM, or EB families", {
    tree <- make_two_regime_simmap(seed=14, n=8)
    X <- matrix(rnorm(ape::Ntip(tree) * 2), ncol=2)
    rownames(X) <- tree$tip.label

    expect_error(
        mvSIMMAP(
            tree, X,
            process=c(A="OU", B="EB"),
            process.groups=c(A="shared", B="shared"),
            method="inverse",
            optimization="fixed",
            echo=FALSE
        ),
        "single process family"
    )

    expect_error(
        mvSIMMAP(
            tree, X,
            process=c(A="OU", B="OUM"),
            process.groups=c(A="shared", B="shared"),
            method="inverse",
            optimization="fixed",
            echo=FALSE
        ),
        "single process family"
    )
})

test_that("mvSIMMAP fixed objects use the mixed print method with a regime summary", {
    tree <- make_two_regime_simmap(seed=16, n=8)
    X <- matrix(rnorm(ape::Ntip(tree) * 2), ncol=2)
    rownames(X) <- tree$tip.label

    fit <- mvSIMMAP(
        tree, X,
        process=c(A="OU", B="OU"),
        process.groups=c(A="ou_shared", B="ou_shared"),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )

    expect_s3_class(fit, "mvmorph.mixed")
    out <- paste(capture.output(print(fit)), collapse="\n")
    expect_match(out, "Regime summary")
    expect_match(out, "process\\.group")
    expect_match(out, "theta\\.owner")
    expect_match(out, "ou_shared")
    expect_match(out, "No optimization performed")
})

test_that("mvSIMMAP print summary distinguishes OU and OUM theta ownership", {
    tree <- make_two_regime_simmap(seed=17, n=8)
    X <- matrix(rnorm(ape::Ntip(tree) * 2), ncol=2)
    rownames(X) <- tree$tip.label

    fit_ou <- mvSIMMAP(
        tree, X,
        process=c(A="OU", B="OU"),
        process.groups=c(A="ou_shared", B="ou_shared"),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )
    fit_oum <- mvSIMMAP(
        tree, X,
        process=c(A="OUM", B="OUM"),
        process.groups=c(A="ou_shared", B="ou_shared"),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )

    out_ou <- paste(capture.output(print(fit_ou)), collapse="\n")
    out_oum <- paste(capture.output(print(fit_oum)), collapse="\n")

    expect_match(out_ou, "theta\\.ou_shared")
    expect_match(out_oum, "theta\\.A")
    expect_match(out_oum, "theta\\.B")
})

test_that("mvSIMMAP summary method returns a compact structural summary", {
    tree <- make_two_regime_simmap(seed=18, n=8)
    X <- matrix(rnorm(ape::Ntip(tree) * 2), ncol=2)
    rownames(X) <- tree$tip.label

    fit <- mvSIMMAP(
        tree, X,
        process=c(A="OUM", B="OUM"),
        process.groups=c(A="ou_shared", B="ou_shared"),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )

    sum_fit <- summary(fit)
    expect_s3_class(sum_fit, "summary.mvmorph.mixed")
    expect_equal(sum_fit$dimensions$nprocess.groups, 1)
    expect_identical(sum_fit$blocks$theta, c("theta.A", "theta.B"))
    expect_identical(sum_fit$blocks$alpha, "ou_shared")
    expect_identical(sum_fit$blocks$sigma, "ou_shared")

    out <- paste(capture.output(print(sum_fit)), collapse="\n")
    expect_match(out, "mvSIMMAP summary")
    expect_match(out, "Parameter blocks")
    expect_match(out, "Regime summary")
    expect_match(out, "theta:\\ttheta\\.A, theta\\.B")
})

test_that("mvSIMMAP multivariate helpers keep fitted output interpretable", {
    tree <- make_two_regime_simmap(seed=181, n=8)
    X <- matrix(rnorm(ape::Ntip(tree) * 2), ncol=2)
    rownames(X) <- tree$tip.label

    fit <- mvSIMMAP(
        tree, X,
        process=c(A="BM", B="OU"),
        method="inverse",
        optimization="fixed",
        param=list(decomp="diagonalPositive", decompSigma="cholesky"),
        echo=FALSE
    )

    expect_identical(colnames(fit$theta), c("trait1", "trait2"))
    expect_identical(rownames(fit$sigma$A), c("trait1", "trait2"))
    expect_identical(colnames(fit$sigma$B), c("trait1", "trait2"))
    expect_identical(rownames(fit$alpha$B), c("trait1", "trait2"))

    expect_equal(coef(fit), fit$theta)

    pars <- parameters(fit)
    expect_s3_class(pars, "data.frame", exact=FALSE)
    expect_true(all(c("block", "group", "row", "col", "estimate", "parameter") %in% names(pars)))
    expect_true(any(pars$parameter == "theta[root,trait1]"))
    expect_true(any(pars$parameter == "sigma[A,trait1,trait2]"))
    expect_true(any(pars$parameter == "alpha[B,trait2,trait2]"))
    expect_equal(coef(fit, type="all"), pars)

    sum_fit <- summary(fit)
    out <- paste(capture.output(print(sum_fit)), collapse="\n")
    expect_match(out, "Multivariate overview")
    expect_match(out, "Sigma variances")
    expect_match(out, "Sigma correlations")
    expect_match(out, "OU pull summary")

    compact_out <- paste(capture.output(print(fit, compact=TRUE)), collapse="\n")
    expect_match(compact_out, "mvSIMMAP summary")
    expect_match(compact_out, "Multivariate overview")
})

test_that("mvSIMMAP supports a multivariate BM to OU simulation and recovery smoke path", {
    benchmark <- make_bm_ou_benchmark_simmap(seed=123, n=20, min_derived_tips=6, target_derived_tips=8)
    process <- c(A="BM", B="OU")
    true_param <- make_bm_ou_true_param(ntraits=2)

    scaffold <- mvSIMMAP(
        benchmark$tree, data=NULL,
        process=process,
        method="rpf",
        optimization="fixed",
        param=list(
            ntraits=2,
            names_traits=c("y1", "y2"),
            decomp="diagonalPositive",
            decompSigma="cholesky"
        ),
        echo=FALSE
    )

    X <- simulate(scaffold, nsim=1, seed=321, param=true_param)
    fit <- mvSIMMAP(
        benchmark$tree, X,
        process=process,
        method="rpf",
        optimization="L-BFGS-B",
        param=list(decomp="diagonalPositive", decompSigma="cholesky"),
        control=list(
            maxit=120,
            retry.unreliable=TRUE,
            retry.max=1,
            retry.jitter=0.05,
            retry.seed=1
        ),
        diagnostic=FALSE,
        echo=FALSE
    )

    expect_gte(benchmark$derived.tips, 6)
    expect_identical(fit$convergence, 0L)
    expect_identical(fit$diagnostics$hessian$status, "reliable")
    expect_identical(colnames(fit$theta), c("y1", "y2"))
    expect_equal(dim(fit$sigma$A), c(2L, 2L))
    expect_equal(dim(fit$alpha$B), c(2L, 2L))

    pars <- parameters(fit)
    expect_true(any(pars$parameter == "theta[theta.B,y1]"))
    expect_true(any(pars$parameter == "sigma[B,y1,y2]"))

    out <- paste(capture.output(print(fit, compact=TRUE)), collapse="\n")
    expect_match(out, "Multivariate overview")
    expect_match(out, "OU pull summary")
})

test_that("mvSIMMAP Hessian diagnostics distinguish reliable, flat, and saddle solutions", {
    reliable <- mvMORPH:::.mvSIMMAP_hessian_diagnostic(diag(c(2, 1)), tolerance=1e-6)
    flat <- mvMORPH:::.mvSIMMAP_hessian_diagnostic(diag(c(1, 1e-9)), tolerance=1e-6)
    saddle <- mvMORPH:::.mvSIMMAP_hessian_diagnostic(diag(c(1, -1e-2)), tolerance=1e-6)

    expect_identical(reliable$status, "reliable")
    expect_identical(reliable$code, 0L)
    expect_identical(flat$status, "flat")
    expect_identical(flat$code, 1L)
    expect_match(flat$label, "flat")
    expect_identical(saddle$status, "saddle")
    expect_identical(saddle$code, 1L)
    expect_match(saddle$label, "saddle")
})

test_that("mvSIMMAP can retry unreliable fits and reports the retained status", {
    set.seed(19)
    tree <- make_single_regime_simmap(ape::rtree(8), regime="A")
    X <- matrix(rnorm(ape::Ntip(tree) * 2), ncol=2)
    rownames(X) <- tree$tip.label

    fit <- mvSIMMAP(
        tree, X,
        process=c(A="OU"),
        method="inverse",
        optimization="L-BFGS-B",
        control=list(
            maxit=25,
            hessian.tol=1e9,
            retry.unreliable=TRUE,
            retry.max=1,
            retry.jitter=0.01,
            retry.seed=1
        ),
        diagnostic=FALSE,
        echo=FALSE
    )

    expect_true(is.list(fit$diagnostics))
    expect_identical(fit$diagnostics$retries$attempted, 1L)
    expect_identical(fit$diagnostics$hessian$status, "flat")
    expect_identical(fit$hess.values, 1L)

    sum_fit <- summary(fit)
    out <- paste(capture.output(print(sum_fit)), collapse="\n")
    expect_match(out, "Hessian:\\t\\s*nearly flat/boundary")
    expect_match(out, "Refits:\\t1 jittered restart")
})

test_that("mvSIMMAP print annotates EB rates at the BM boundary", {
    tree <- make_single_regime_simmap(ape::rtree(6), regime="A")
    X <- matrix(rnorm(ape::Ntip(tree)), ncol=1)
    rownames(X) <- tree$tip.label

    fit <- mvSIMMAP(
        tree, X,
        process=c(A="EB"),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )
    fit$beta[] <- 0

    out <- paste(capture.output(print(fit)), collapse="\n")
    expect_match(out, "BM boundary")
    expect_match(out, "process\\.group")
})

test_that("mvSIMMAP optimized fits expose a standards-compatible logLik object", {
    set.seed(20)
    tree <- make_single_regime_simmap(ape::rtree(8), regime="A")
    X <- matrix(rnorm(ape::Ntip(tree)), ncol=1)
    rownames(X) <- tree$tip.label

    fit <- mvSIMMAP(
        tree, X,
        process=c(A="OU"),
        method="inverse",
        optimization="L-BFGS-B",
        control=list(maxit=30),
        diagnostic=FALSE,
        echo=FALSE
    )

    expect_s3_class(fit$LogLik, "logLik")
    expect_identical(attr(fit$LogLik, "df"), as.integer(fit$param$nparam))
    expect_identical(attr(fit$LogLik, "nobs"), as.integer(fit$param$nobs))
    expect_equal(as.numeric(logLik(fit)), as.numeric(fit$LogLik))

    bic_expected <- -2 * as.numeric(fit$LogLik) +
        log(attr(fit$LogLik, "nobs")) * attr(fit$LogLik, "df")
    expect_equal(as.numeric(BIC(fit)), bic_expected)
})

test_that("simulate() can draw from a mixed mvSIMMAP scaffold with manual overrides", {
    tree <- make_two_regime_simmap(seed=21, n=8)

    scaffold <- mvSIMMAP(
        tree, data=NULL,
        process=c(A="OU", B="EB"),
        process.groups=c(A="ou_shared", B="eb_shared"),
        param=list(ntraits=2, names_traits=c("y1", "y2")),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )

    override <- list(
        theta=list(
            root=c(0.2, -0.1),
            "theta.ou_shared"=c(0.7, -0.4)
        ),
        sigma=list(
            ou_shared=diag(c(0.05, 0.08)),
            eb_shared=diag(c(0.09, 0.12))
        ),
        alpha=list(
            ou_shared=diag(c(1.2, 0.6))
        ),
        beta=c(eb_shared=-0.25)
    )

    spec <- mvMORPH:::.mvSIMMAP_extract_spec(scaffold)
    resolved <- mvMORPH:::.mvSIMMAP_simulation_values(scaffold, spec, param=override)
    built <- mvMORPH:::.mvSIMMAP_build_model(spec, resolved$values)
    expected_mean <- as.numeric(built$D %*% as.numeric(t(resolved$theta)))

    set.seed(101)
    expected <- matrix(
        mvMORPH:::rmvnorm_simul(n=1, mean=expected_mean, var=built$V, method="cholesky"),
        ncol=2
    )
    rownames(expected) <- spec$tree$tip.label
    colnames(expected) <- colnames(resolved$theta)

    actual <- simulate(scaffold, nsim=1, seed=101, param=override)

    expect_equal(actual, expected)
})

test_that("simulate() supports single-group scalar overrides and error inflation", {
    tree <- make_single_regime_simmap(ape::rtree(6), regime="A")

    scaffold <- mvSIMMAP(
        tree, data=NULL,
        process=c(A="OU"),
        param=list(ntraits=1, names_traits="y"),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )

    override <- list(
        theta=c(root=0.3, "theta.A"=-0.2),
        sigma=matrix(0.15, 1, 1),
        alpha=matrix(0.8, 1, 1)
    )
    error_mat <- matrix(seq_len(ape::Ntip(tree)) * 0.01, ncol=1,
                        dimnames=list(tree$tip.label, "y"))

    spec <- mvMORPH:::.mvSIMMAP_extract_spec(scaffold)
    resolved <- mvMORPH:::.mvSIMMAP_simulation_values(scaffold, spec, param=override)
    built <- mvMORPH:::.mvSIMMAP_build_model(spec, resolved$values)
    expected_mean <- as.numeric(built$D %*% as.numeric(t(resolved$theta)))
    V_error <- built$V
    diag(V_error) <- diag(V_error) + as.numeric(error_mat)

    set.seed(202)
    expected <- matrix(
        mvMORPH:::rmvnorm_simul(n=1, mean=expected_mean, var=V_error, method="cholesky"),
        ncol=1
    )
    rownames(expected) <- spec$tree$tip.label
    colnames(expected) <- colnames(resolved$theta)

    actual <- simulate(scaffold, nsim=1, seed=202, param=override, error=error_mat)

    expect_equal(actual, expected)
})

test_that("simulate() on mvSIMMAP scaffolds follows mvSIM return conventions", {
    tree <- make_two_regime_simmap(seed=22, n=8)

    scaffold_multi <- mvSIMMAP(
        tree, data=NULL,
        process=c(A="OUM", B="EB"),
        process.groups=c(A="ou_shared", B="eb_shared"),
        param=list(ntraits=2, names_traits=c("y1", "y2")),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )

    sims_multi <- simulate(scaffold_multi, nsim=3, seed=303)
    expect_true(is.list(sims_multi))
    expect_length(sims_multi, 3)
    expect_identical(dim(sims_multi[[1]]), c(ape::Ntip(tree), 2L))
    expect_identical(rownames(sims_multi[[1]]), tree$tip.label)
    expect_identical(colnames(sims_multi[[1]]), c("y1", "y2"))

    tree_uni <- make_single_regime_simmap(ape::rtree(6), regime="A")

    scaffold_uni <- mvSIMMAP(
        tree_uni, data=NULL,
        process=c(A="EB"),
        param=list(ntraits=1, names_traits="y"),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )

    sims_uni <- simulate(scaffold_uni, nsim=4, seed=404)
    expect_true(is.matrix(sims_uni))
    expect_identical(dim(sims_uni), c(ape::Ntip(tree_uni), 4L))
    expect_identical(rownames(sims_uni), tree_uni$tip.label)
})

test_that("mvSIMMAP fixed scaffolds can be created without data", {
    tree <- make_two_regime_simmap(seed=23, n=8)

    scaffold <- mvSIMMAP(
        tree, data=NULL,
        process=c(A="OU", B="EB"),
        process.groups=c(A="ou_shared", B="eb_shared"),
        param=list(ntraits=2, names_traits=c("trait1", "trait2")),
        method="inverse",
        optimization="fixed",
        echo=FALSE
    )

    expect_s3_class(scaffold, "mvmorph.mixed")
    expect_identical(dim(scaffold$theta), c(2L, 2L))
    expect_identical(rownames(scaffold$theta), c("root", "theta.ou_shared"))
    expect_identical(colnames(scaffold$theta), c("trait1", "trait2"))
    expect_identical(scaffold$param$ntraits, 2L)
    expect_identical(scaffold$param$nbspecies, ape::Ntip(tree))
})

test_that("mvSIMMAP still requires data outside fixed scaffold mode", {
    tree <- make_single_regime_simmap(ape::rtree(6), regime="A")

    expect_error(
        mvSIMMAP(
            tree, data=NULL,
            process=c(A="OU"),
            method="inverse",
            optimization="L-BFGS-B",
            echo=FALSE
        ),
        "unless optimization=\\\"fixed\\\""
    )
})
