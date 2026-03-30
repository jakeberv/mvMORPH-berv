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
