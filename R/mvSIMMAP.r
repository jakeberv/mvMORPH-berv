################################################################################
##                                                                            ##
##                           mvMORPH: mvSIMMAP                                ##
##                                                                            ##
## Experimental mixed-regime Gaussian models on SIMMAP trees                  ##
##                                                                            ##
################################################################################

.mvSIMMAP_reorder_tree <- function(tree, order="cladewise"){
    ind <- reorder.phylo(tree, order=order, index.only=TRUE)
    tree$edge <- tree$edge[ind, , drop=FALSE]
    tree$edge.length <- tree$edge.length[ind]
    if(!is.null(tree[["mapped.edge"]])){
        tree$mapped.edge <- tree$mapped.edge[ind, , drop=FALSE]
    }
    if(!is.null(tree[["maps"]])){
        tree$maps <- tree$maps[ind]
    }
    tree
}

.mvSIMMAP_scale_tree <- function(tree, scale.height=FALSE){
    if(!isTRUE(scale.height)) return(tree)
    maxHeight <- max(nodeHeights(tree))
    if(!is.finite(maxHeight) || maxHeight <= 0) return(tree)
    tree$edge.length <- tree$edge.length/maxHeight
    if(!is.null(tree[["mapped.edge"]])){
        tree$mapped.edge <- tree$mapped.edge/maxHeight
    }
    if(!is.null(tree[["maps"]])){
        tree$maps <- lapply(tree$maps, function(x){
            structure(as.numeric(x)/maxHeight, names=names(x))
        })
    }
    tree
}

.mvSIMMAP_safe_inverse <- function(x){
    inv <- try(solve(x), silent=TRUE)
    if(inherits(inv, "try-error")){
        pseudoinverse(x)
    }else{
        inv
    }
}

.mvSIMMAP_as_named_process <- function(process, regime_names){
    if(is.null(names(process))){
        if(length(process) != length(regime_names)){
            stop("The process vector must either be named by SIMMAP regime or have one entry per regime", call.=FALSE)
        }
        names(process) <- regime_names
    }
    process_names <- names(process)
    process <- toupper(as.character(process))
    names(process) <- process_names
    process <- process[regime_names]
    if(any(is.na(process))){
        stop("The process vector must include every regime present in tree$mapped.edge", call.=FALSE)
    }
    valid <- c("BM", "OU", "EB")
    if(any(!process %in% valid)){
        stop("Allowed process assignments are \"BM\", \"OU\", and \"EB\"", call.=FALSE)
    }
    structure(process, names=regime_names)
}

.mvSIMMAP_pick_list_value <- function(x, regime, regime_order){
    if(is.null(x)) return(NULL)
    if(!is.list(x)){
        return(x)
    }
    if(!is.null(names(x)) && regime %in% names(x)){
        return(x[[regime]])
    }
    idx <- match(regime, regime_order)
    if(is.na(idx) || idx > length(x)) return(NULL)
    x[[idx]]
}

.mvSIMMAP_pick_beta_value <- function(x, regime, regime_order){
    if(is.null(x)) return(NULL)
    if(is.list(x)){
        return(.mvSIMMAP_pick_list_value(x, regime, regime_order))
    }
    if(!is.null(names(x)) && regime %in% names(x)){
        return(x[[regime]])
    }
    if(length(x) == length(regime_order)){
        idx <- match(regime, regime_order)
        return(x[[idx]])
    }
    x[[1]]
}

.mvSIMMAP_sigma_start <- function(value, p, decompSigma, tree, data, fallback){
    if(is.null(value)) return(fallback)
    if(is.matrix(value)) return(startParamSigma(p, decompSigma, tree, data, guess=value))
    if(is.numeric(value)) return(as.numeric(value))
    stop("sigma starting values must be matrices, parameter vectors, or a list of such values", call.=FALSE)
}

.mvSIMMAP_alpha_start <- function(value, p, decomp, tree, fallback){
    if(is.null(value)) return(fallback)
    if(is.matrix(value)){
        switch(decomp,
            "cholesky" = sym.unpar(value),
            "diagonal" = diag(value),
            "diagonalPositive" = log(diag(value)),
            startParam(p, decomp, tree, hlife=Re(eigen(value)$values))
        )
    }else if(is.numeric(value)){
        as.numeric(value)
    }else{
        stop("alpha starting values must be matrices, parameter vectors, or a list of such values", call.=FALSE)
    }
}

.mvSIMMAP_symmetrize <- function(x){
    Re((x + t(x))/2)
}

.mvSIMMAP_eb_factor <- function(beta, duration, start_age=0){
    if(abs(beta) < sqrt(.Machine$double.eps)){
        duration
    }else{
        exp(beta * start_age) * (exp(beta * duration) - 1)/beta
    }
}

.mvSIMMAP_ou_transition <- function(alpha, sigma, duration){
    p <- nrow(alpha)
    if(duration <= 0){
        I_p <- diag(p)
        return(list(Phi=I_p, Theta=matrix(0, p, p), Q=matrix(0, p, p)))
    }
    eig <- eigen(alpha)
    invectors <- try(solve(eig$vectors), silent=TRUE)
    if(inherits(invectors, "try-error")){
        invectors <- pseudoinverse(eig$vectors)
    }
    lambda <- eig$values
    exp_diag <- diag(exp(-lambda * duration), p)
    Phi <- Re(eig$vectors %*% exp_diag %*% invectors)
    U <- invectors %*% sigma %*% t(invectors)
    G <- matrix(0+0i, p, p)
    for(i in seq_len(p)){
        for(j in seq_len(p)){
            denom <- lambda[i] + lambda[j]
            G[i, j] <- if(abs(denom) < sqrt(.Machine$double.eps)){
                duration * U[i, j]
            }else{
                ((1 - exp(-denom * duration))/denom) * U[i, j]
            }
        }
    }
    Q <- Re(eig$vectors %*% G %*% t(eig$vectors))
    I_p <- diag(p)
    list(
        Phi=Re(Phi),
        Theta=Re(I_p - Phi),
        Q=.mvSIMMAP_symmetrize(Q)
    )
}

.mvSIMMAP_transition <- function(type, sigma, duration, alpha=NULL, beta=NULL, eb_age=0){
    p <- nrow(sigma)
    I_p <- diag(p)
    if(duration <= 0){
        return(list(Phi=I_p, Theta=matrix(0, p, p), Q=matrix(0, p, p)))
    }
    switch(type,
        "BM" = {
            list(Phi=I_p, Theta=matrix(0, p, p), Q=.mvSIMMAP_symmetrize(sigma * duration))
        },
        "EB" = {
            list(
                Phi=I_p,
                Theta=matrix(0, p, p),
                Q=.mvSIMMAP_symmetrize(sigma * .mvSIMMAP_eb_factor(beta, duration, start_age=eb_age))
            )
        },
        "OU" = .mvSIMMAP_ou_transition(alpha=alpha, sigma=sigma, duration=duration),
        stop("Unknown process type in internal SIMMAP transition builder", call.=FALSE)
    )
}

.mvSIMMAP_edge_segments <- function(tree, edge_index){
    edge_map <- tree$maps[[edge_index]]
    if(is.null(edge_map) || length(edge_map) == 0L){
        return(list(regimes=character(0), times=numeric(0)))
    }
    times <- as.numeric(edge_map)
    edge_length <- tree$edge.length[edge_index]
    total <- sum(times)
    if(is.finite(total) && total > 0 && abs(total - edge_length) > 1e-10){
        times <- times * (edge_length/total)
    }
    list(regimes=names(edge_map), times=times)
}

.mvSIMMAP_build_design <- function(A_nodes, B_nodes, n, p){
    n_blocks <- 1 + length(B_nodes[[1]])
    D <- matrix(0, n * p, n_blocks * p)
    for(i in seq_len(n)){
        coeffs <- c(list(A_nodes[[i]]), B_nodes[[i]])
        for(block in seq_len(n_blocks)){
            mat <- Re(coeffs[[block]])
            cols <- ((block - 1) * p + 1):(block * p)
            for(trait in seq_len(p)){
                D[(trait - 1) * n + i, cols] <- mat[trait, ]
            }
        }
    }
    D
}

.mvSIMMAP_build_vcv <- function(A_nodes, V_nodes, invA_nodes, mrca_mat, n, p){
    I_p <- diag(p)
    V <- matrix(0, n * p, n * p)
    tip_rows <- lapply(seq_len(n), function(i) i + (0:(p - 1)) * n)
    for(i in seq_len(n)){
        rows_i <- tip_rows[[i]]
        for(j in seq_len(i)){
            rows_j <- tip_rows[[j]]
            anc <- if(i == j) i else mrca_mat[i, j]
            Pi <- if(i == anc) I_p else A_nodes[[i]] %*% invA_nodes[[anc]]
            Pj <- if(j == anc) I_p else A_nodes[[j]] %*% invA_nodes[[anc]]
            block <- Re(Pi %*% V_nodes[[anc]] %*% t(Pj))
            if(i == j) block <- .mvSIMMAP_symmetrize(block)
            V[rows_i, rows_j] <- block
            if(i != j){
                V[rows_j, rows_i] <- t(block)
            }
        }
    }
    V
}

.mvSIMMAP_unpack <- function(par, spec){
    idx <- 1L
    alpha <- list()
    sigma <- list()
    beta <- numeric(0)
    for(regime in spec$regime_names){
        type <- spec$process[[regime]]
        if(type == "OU"){
            alpha_par <- par[idx:(idx + spec$nalpha - 1)]
            idx <- idx + spec$nalpha
            alpha[[regime]] <- spec$alphafun(alpha_par)$A
        }
        sigma_par <- par[idx:(idx + spec$nsigma - 1)]
        idx <- idx + spec$nsigma
        sigma[[regime]] <- spec$sigmafun(sigma_par)
        if(type == "EB"){
            beta <- c(beta, par[[idx]])
            names(beta)[length(beta)] <- regime
            idx <- idx + 1L
        }
    }
    beta <- structure(as.numeric(beta), names=names(beta))
    list(alpha=alpha, sigma=sigma, beta=beta)
}

.mvSIMMAP_build_model <- function(spec, values){
    tree <- spec$tree
    n <- spec$n
    p <- spec$p
    n_total <- spec$n_total
    I_p <- diag(p)
    zero_p <- matrix(0, p, p)
    A_nodes <- vector("list", n_total)
    B_nodes <- vector("list", n_total)
    V_nodes <- vector("list", n_total)
    state_nodes <- vector("list", n_total)
    A_nodes[[spec$root_node]] <- I_p
    B_nodes[[spec$root_node]] <- if(length(spec$ou_regimes)){
        lapply(seq_along(spec$ou_regimes), function(x) zero_p)
    }else{
        list()
    }
    V_nodes[[spec$root_node]] <- zero_p
    state_nodes[[spec$root_node]] <- list(type=NULL, regime=NULL, eb_age=0)
    for(edge_index in seq_len(nrow(tree$edge))){
        parent <- tree$edge[edge_index, 1]
        child <- tree$edge[edge_index, 2]
        if(is.null(A_nodes[[parent]])){
            stop("Internal SIMMAP traversal error: parent node was not initialized", call.=FALSE)
        }
        A_curr <- A_nodes[[parent]]
        B_curr <- if(length(spec$ou_regimes)){
            lapply(B_nodes[[parent]], function(x) x)
        }else{
            list()
        }
        V_curr <- V_nodes[[parent]]
        state_curr <- state_nodes[[parent]]
        edge_segments <- spec$segments[[edge_index]]
        if(length(edge_segments$times)){
            for(seg_idx in seq_along(edge_segments$times)){
                duration <- edge_segments$times[[seg_idx]]
                if(duration <= 0) next
                regime <- edge_segments$regimes[[seg_idx]]
                type <- spec$process[[regime]]
                eb_age <- if(type == "EB" &&
                             identical(state_curr$type, "EB") &&
                             identical(state_curr$regime, regime)){
                    state_curr$eb_age
                }else{
                    0
                }
                trans <- .mvSIMMAP_transition(
                    type=type,
                    sigma=values$sigma[[regime]],
                    duration=duration,
                    alpha=if(type == "OU") values$alpha[[regime]] else NULL,
                    beta=if(type == "EB") values$beta[[regime]] else NULL,
                    eb_age=eb_age
                )
                A_curr <- trans$Phi %*% A_curr
                if(length(B_curr)){
                    B_curr <- lapply(B_curr, function(mat) trans$Phi %*% mat)
                }
                if(type == "OU"){
                    B_curr[[spec$ou_index[[regime]]]] <- B_curr[[spec$ou_index[[regime]]]] + trans$Theta
                }
                V_curr <- .mvSIMMAP_symmetrize(trans$Phi %*% V_curr %*% t(trans$Phi) + trans$Q)
                state_curr <- if(type == "EB"){
                    list(type=type, regime=regime, eb_age=eb_age + duration)
                }else{
                    list(type=type, regime=regime, eb_age=0)
                }
            }
        }
        A_nodes[[child]] <- A_curr
        B_nodes[[child]] <- B_curr
        V_nodes[[child]] <- V_curr
        state_nodes[[child]] <- state_curr
    }
    invA_nodes <- lapply(A_nodes, function(x){
        if(is.null(x)) NULL else .mvSIMMAP_safe_inverse(x)
    })
    list(
        D=.mvSIMMAP_build_design(A_nodes=A_nodes, B_nodes=B_nodes, n=n, p=p),
        V=.mvSIMMAP_build_vcv(
            A_nodes=A_nodes,
            V_nodes=V_nodes,
            invA_nodes=invA_nodes,
            mrca_mat=spec$mrca,
            n=n,
            p=p
        )
    )
}

.mvSIMMAP_make_theta_matrix <- function(theta, mean_names, trait_names){
    theta_mat <- matrix(theta, nrow=length(mean_names), byrow=TRUE)
    rownames(theta_mat) <- mean_names
    colnames(theta_mat) <- trait_names
    theta_mat
}

mvSIMMAP <- function(tree, data, process, error=NULL,
                     param=list(alpha=NULL, sigma=NULL, beta=NULL,
                                decomp=c("cholesky", "diagonal", "diagonalPositive"),
                                decompSigma=c("cholesky", "diagonal")),
                     method=c("rpf", "inverse", "pseudoinverse"),
                     scale.height=FALSE,
                     optimization=c("L-BFGS-B", "Nelder-Mead", "subplex"),
                     control=list(maxit=20000),
                     diagnostic=TRUE, echo=TRUE){

    if(missing(tree)) stop("The tree object is missing!", call.=FALSE)
    if(missing(data)) stop("You must provide a dataset along with your tree!", call.=FALSE)
    if(missing(process)) stop("You must provide a per-regime process assignment", call.=FALSE)
    if(!inherits(tree, "simmap")) stop("A tree of class \"simmap\" is required", call.=FALSE)
    if(is.null(tree[["mapped.edge"]]) || is.null(tree[["maps"]])){
        stop("The tree must include SIMMAP mapped regimes", call.=FALSE)
    }
    data <- as.matrix(data)
    if(!is.null(rownames(data))){
        if(any(tree$tip.label == rownames(data))){
            data <- data[tree$tip.label, , drop=FALSE]
        }else if(isTRUE(echo)){
            cat("row names of the data matrix must match tip names of your phylogeny!", "\n")
        }
    }else if(isTRUE(echo)){
        cat("species in the matrix are assumed to be in the same order as in the phylogeny, otherwise specify rownames of 'data'", "\n")
    }
    method <- method[1]
    if(!method %in% c("rpf", "inverse", "pseudoinverse")){
        stop("mvSIMMAP currently supports only the \"rpf\", \"inverse\", and \"pseudoinverse\" methods", call.=FALSE)
    }
    optimization <- optimization[1]
    p <- ncol(data)
    if(is.null(p)) p <- 1L
    n <- nrow(data)
    if(n != Ntip(tree)){
        stop("The number of rows in data must match the number of tips in the tree", call.=FALSE)
    }
    if(!is.null(error)){
        error <- as.vector(error)
        error[is.na(error)] <- 0
    }
    NA_val <- FALSE
    Indice_NA <- NULL
    if(any(is.na(data))){
        NA_val <- TRUE
        Indice_NA <- which(is.na(as.vector(data)))
    }
    tree <- .mvSIMMAP_scale_tree(tree, scale.height=scale.height)
    tree <- .mvSIMMAP_reorder_tree(tree, order="cladewise")
    regime_names <- colnames(tree$mapped.edge)
    process <- .mvSIMMAP_as_named_process(process, regime_names=regime_names)
    ou_regimes <- regime_names[process == "OU"]
    eb_regimes <- regime_names[process == "EB"]
    decomp <- if(is.null(param[["decomp"]])) "cholesky" else param$decomp[1]
    decompSigma <- if(is.null(param[["decompSigma"]])) "cholesky" else param$decompSigma[1]
    sigma_fallback <- startParamSigma(p, decompSigma, tree, data)
    alpha_fallback <- startParam(p, decomp, tree)
    maxHeight <- max(nodeHeights(tree))
    if(!is.finite(maxHeight) || maxHeight <= 0) maxHeight <- 1
    beta_default <- -1/maxHeight
    low_default <- log(1e-5)/maxHeight
    up_default <- 0
    beta_low <- if(is.null(param[["low"]])) {
        setNames(rep(low_default, length(eb_regimes)), eb_regimes)
    }else{
        vals <- param$low
        if(is.null(names(vals)) && length(vals) == 1L){
            setNames(rep(vals, length(eb_regimes)), eb_regimes)
        }else if(is.null(names(vals)) && length(vals) == length(eb_regimes)){
            setNames(vals, eb_regimes)
        }else{
            setNames(as.numeric(vals[eb_regimes]), eb_regimes)
        }
    }
    beta_up <- if(is.null(param[["up"]])) {
        setNames(rep(up_default, length(eb_regimes)), eb_regimes)
    }else{
        vals <- param$up
        if(is.null(names(vals)) && length(vals) == 1L){
            setNames(rep(vals, length(eb_regimes)), eb_regimes)
        }else if(is.null(names(vals)) && length(vals) == length(eb_regimes)){
            setNames(vals, eb_regimes)
        }else{
            setNames(as.numeric(vals[eb_regimes]), eb_regimes)
        }
    }
    if(anyNA(beta_low) || anyNA(beta_up)){
        stop("Named beta bounds must match the EB regimes in the SIMMAP tree", call.=FALSE)
    }
    if(length(eb_regimes) && any(beta_low > beta_up)){
        stop("Each EB lower bound must be less than or equal to its upper bound", call.=FALSE)
    }
    sigma_start <- vector("list", length(regime_names))
    names(sigma_start) <- regime_names
    for(regime in regime_names){
        raw <- .mvSIMMAP_pick_list_value(param[["sigma"]], regime, regime_names)
        sigma_start[[regime]] <- .mvSIMMAP_sigma_start(
            value=raw,
            p=p,
            decompSigma=decompSigma,
            tree=tree,
            data=data,
            fallback=sigma_fallback
        )
    }
    alpha_start <- vector("list", length(ou_regimes))
    names(alpha_start) <- ou_regimes
    for(regime in ou_regimes){
        raw <- .mvSIMMAP_pick_list_value(param[["alpha"]], regime, ou_regimes)
        alpha_start[[regime]] <- .mvSIMMAP_alpha_start(
            value=raw,
            p=p,
            decomp=decomp,
            tree=tree,
            fallback=alpha_fallback
        )
    }
    beta_start <- numeric(length(eb_regimes))
    names(beta_start) <- eb_regimes
    for(regime in eb_regimes){
        raw <- .mvSIMMAP_pick_beta_value(param[["beta"]], regime, eb_regimes)
        beta_start[[regime]] <- if(is.null(raw)) beta_default else as.numeric(raw)[1]
    }
    alphafun <- function(x) matrixParam(x, p, decomp)
    sigmafun <- function(x) symPar(x, decomp=decompSigma, p=p)
    nalpha <- length(alpha_fallback)
    nsigma <- length(sigma_fallback)
    start <- numeric(0)
    beta_idx <- integer(0)
    for(regime in regime_names){
        type <- process[[regime]]
        if(type == "OU"){
            start <- c(start, alpha_start[[regime]])
        }
        start <- c(start, sigma_start[[regime]])
        if(type == "EB"){
            beta_idx <- c(beta_idx, length(start) + 1L)
            start <- c(start, beta_start[[regime]])
        }
    }
    mean_names <- if(length(ou_regimes)){
        c("root", paste0("theta.", ou_regimes))
    }else{
        "root"
    }
    sizeD <- length(mean_names) * p
    spec <- list(
        tree=tree,
        process=process,
        regime_names=regime_names,
        ou_regimes=ou_regimes,
        ou_index=setNames(seq_along(ou_regimes), ou_regimes),
        segments=lapply(seq_len(nrow(tree$edge)), function(i) .mvSIMMAP_edge_segments(tree, i)),
        mrca=mrca(tree, full=FALSE),
        n=n,
        p=p,
        n_total=Ntip(tree) + tree$Nnode,
        root_node=Ntip(tree) + 1L,
        nalpha=nalpha,
        nsigma=nsigma,
        alphafun=alphafun,
        sigmafun=sigmafun
    )
    eval_loglik <- function(par, theta_mle=TRUE, theta=NULL){
        res <- try({
            values <- .mvSIMMAP_unpack(par, spec)
            built <- .mvSIMMAP_build_model(spec, values)
            loglik_mvmorph(
                data=data,
                V=built$V,
                D=built$D,
                n=n,
                k=p,
                error=error,
                method=method,
                sizeD=ncol(built$D),
                NA_val=NA_val,
                Indice_NA=Indice_NA,
                theta_mle=theta_mle,
                theta=theta
            )
        }, silent=TRUE)
        if(inherits(res, "try-error") || !is.finite(res$logl)){
            list(logl=-Inf, theta=rep(NA_real_, sizeD))
        }else{
            list(logl=as.numeric(res$logl), theta=as.numeric(res$anc))
        }
    }
    llfun <- function(par, ...){
        args <- list(...)
        if(is.null(args[["root.mle"]])) args$root.mle <- TRUE
        if(is.null(args[["theta"]])) args$theta <- FALSE
        if(isTRUE(args$root.mle)){
            result <- eval_loglik(par, theta_mle=TRUE)
            if(isTRUE(args$theta)){
                list(logl=result$logl, theta=result$theta)
            }else{
                result$logl
            }
        }else{
            theta_vec <- par[length(start) + seq_len(sizeD)]
            eval_loglik(par[seq_len(length(start))], theta_mle=FALSE, theta=theta_vec)$logl
        }
    }
    attr(llfun, "model") <- "SIMMAPmixed"
    attr(llfun, "alpha") <- length(ou_regimes) * nalpha
    attr(llfun, "sigma") <- length(regime_names) * nsigma
    attr(llfun, "beta") <- length(eb_regimes)
    attr(llfun, "theta") <- sizeD
    class(llfun) <- c("mvmorph.llik")
    theta_zero <- .mvSIMMAP_make_theta_matrix(
        theta=rep(0, sizeD),
        mean_names=mean_names,
        trait_names=if(is.null(colnames(data))) rep("", p) else colnames(data)
    )
    if(optimization == "fixed"){
        values <- .mvSIMMAP_unpack(start, spec)
        param$nbspecies <- n
        param$ntraits <- p
        param$nregimes <- length(regime_names)
        param$method <- method
        param$optimization <- optimization
        param$traits <- colnames(data)
        param$names_regimes <- regime_names
        param$process <- process
        param$model <- "SIMMAPmixed"
        param$decomp <- decomp
        param$decompSigma <- decompSigma
        param$mean.parameters <- mean_names
        param$alphafun <- alphafun
        param$sigmafun <- sigmafun
        results <- list(
            llik=llfun,
            theta=theta_zero,
            alpha=values$alpha,
            sigma=values$sigma,
            beta=values$beta,
            process=process,
            param=param
        )
        class(results) <- c("mvmorph")
        invisible(results)
        return(results)
    }
    objective <- function(par){
        logl <- eval_loglik(par, theta_mle=TRUE)$logl
        if(!is.finite(logl)) .Machine$double.xmax/100 else -logl
    }
    if(optimization == "subplex"){
        estim <- subplex(par=start, fn=objective, hessian=TRUE, control=control)
    }else if(optimization == "L-BFGS-B" && length(beta_idx)){
        lower <- rep(-Inf, length(start))
        upper <- rep(Inf, length(start))
        lower[beta_idx] <- beta_low[names(beta_start)]
        upper[beta_idx] <- beta_up[names(beta_start)]
        estim <- optim(
            par=start,
            fn=objective,
            method=optimization,
            lower=lower,
            upper=upper,
            hessian=TRUE,
            control=control
        )
    }else{
        estim <- optim(
            par=start,
            fn=objective,
            method=optimization,
            hessian=TRUE,
            control=control
        )
    }
    res.theta <- eval_loglik(estim$par, theta_mle=TRUE)$theta
    hessian_vals <- try(eigen(estim$hessian)$values, silent=TRUE)
    if(inherits(hessian_vals, "try-error")){
        hess.val <- 1
    }else{
        hess.val <- if(any(Re(hessian_vals) < 0)) 1 else 0
    }
    if(isTRUE(diagnostic)){
        if(estim$convergence == 0){
            cat("successful convergence of the optimizer", "\n")
        }else if(estim$convergence == 1){
            cat("maximum limit iteration has been reached, please consider increase maxit", "\n")
        }else{
            cat("convergence of the optimizer has not been reached, try simpler model", "\n")
        }
        if(hess.val == 0){
            cat("a reliable solution has been reached", "\n")
        }else{
            cat("unreliable solution has been reached, check hessian eigenvectors or try simpler model", "\n")
        }
    }
    LL <- -estim$value
    nparam <- length(estim$par) + sizeD
    nobs <- length(which(!is.na(data)))
    AIC <- -2 * LL + 2 * nparam
    AICc <- AIC + ((2 * nparam * (nparam + 1))/(nobs - nparam - 1))
    fitted_values <- .mvSIMMAP_unpack(estim$par, spec)
    trait_names <- colnames(data)
    if(is.null(trait_names)) trait_names <- rep("", p)
    theta.mat <- .mvSIMMAP_make_theta_matrix(
        theta=res.theta,
        mean_names=mean_names,
        trait_names=trait_names
    )
    sigma_named <- lapply(fitted_values$sigma, function(x){
        dimnames(x) <- list(trait_names, trait_names)
        x
    })
    alpha_named <- lapply(fitted_values$alpha, function(x){
        dimnames(x) <- list(trait_names, trait_names)
        x
    })
    if(isTRUE(echo)){
        cat("\n")
        cat("-- Summary results for the experimental SIMMAP mixed model --", "\n")
        cat("LogLikelihood:\t", LL, "\n")
        cat("AIC:\t", AIC, "\n")
        cat("AICc:\t", AICc, "\n")
        cat(nparam, "parameters", "\n")
        cat("Process assignment", "\n")
        cat("______________________", "\n")
        print(process)
        cat("\n")
        cat("Estimated mean parameters", "\n")
        cat("______________________", "\n")
        print(theta.mat)
        if(length(alpha_named)){
            cat("\n")
            cat("OU alpha matrices", "\n")
            cat("______________________", "\n")
            print(alpha_named)
        }
        cat("\n")
        cat("Regime sigma matrices", "\n")
        cat("______________________", "\n")
        print(sigma_named)
        if(length(fitted_values$beta)){
            cat("\n")
            cat("EB beta values", "\n")
            cat("______________________", "\n")
            print(fitted_values$beta)
        }
    }
    param$nparam <- nparam
    param$nbspecies <- n
    param$ntraits <- p
    param$nregimes <- length(regime_names)
    param$method <- method
    param$optimization <- optimization
    param$traits <- trait_names
    param$names_regimes <- regime_names
    param$process <- process
    param$model <- "SIMMAPmixed"
    param$decomp <- decomp
    param$decompSigma <- decompSigma
    param$mean.parameters <- mean_names
    param$alphafun <- alphafun
    param$sigmafun <- sigmafun
    param$opt <- estim
    results <- list(
        LogLik=LL,
        AIC=AIC,
        AICc=AICc,
        theta=theta.mat,
        alpha=alpha_named,
        sigma=sigma_named,
        beta=fitted_values$beta,
        process=process,
        convergence=estim$convergence,
        hess.values=hess.val,
        param=param,
        llik=llfun
    )
    class(results) <- c("mvmorph", "mvmorph.mixed")
    invisible(results)
}

print.mvmorph.mixed <- function(x, ...){
    cat("\n")
    cat("-- Summary results for the experimental SIMMAP mixed model --", "\n")
    cat("LogLikelihood:\t", x$LogLik, "\n")
    cat("AIC:\t", x$AIC, "\n")
    cat("AICc:\t", x$AICc, "\n")
    cat("Process assignment", "\n")
    cat("______________________", "\n")
    print(x$process)
    cat("\n")
    cat("Estimated mean parameters", "\n")
    cat("______________________", "\n")
    print(x$theta)
    if(length(x$alpha)){
        cat("\n")
        cat("OU alpha matrices", "\n")
        cat("______________________", "\n")
        print(x$alpha)
    }
    cat("\n")
    cat("Regime sigma matrices", "\n")
    cat("______________________", "\n")
    print(x$sigma)
    if(length(x$beta)){
        cat("\n")
        cat("EB beta values", "\n")
        cat("______________________", "\n")
        print(x$beta)
    }
}
