make_single_regime_simmap <- function(tree, regime="A"){
    tree <- ape::reorder.phylo(tree, "cladewise")
    tree$maps <- lapply(tree$edge.length, function(x) structure(x, names=regime))
    tree$mapped.edge <- matrix(tree$edge.length, ncol=1, dimnames=list(NULL, regime))
    class(tree) <- c("simmap", setdiff(class(tree), "simmap"))
    tree
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
