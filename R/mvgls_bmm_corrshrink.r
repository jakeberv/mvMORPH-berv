################################################################################
##                                                                            ##
##              mvMORPH: mvgls_bmm_corrshrink.r                               ##
##                                                                            ##
##   Experimental corr-shrink BMM fit for mvgls                               ##
##                                                                            ##
################################################################################

.mvgls_bmm_corrshrink_regimes <- function(tree){
    mapped_edge <- tree$mapped.edge
    mapped_sum <- rowSums(mapped_edge)
    edge_ratio <- rep(1, length(mapped_sum))
    nz <- mapped_sum > 0
    edge_ratio[nz] <- tree$edge.length[nz]/mapped_sum[nz]
    edge_ratio[!is.finite(edge_ratio)] <- 1
    mapped_edge <- mapped_edge * edge_ratio
    regime_names <- colnames(mapped_edge)
    regime_trees <- vector("list", length(regime_names))
    regime_cov <- vector("list", length(regime_names))
    names(regime_trees) <- regime_names
    names(regime_cov) <- regime_names
    for(i in seq_along(regime_names)){
        regime_trees[[i]] <- tree
        regime_trees[[i]]$edge.length <- mapped_edge[, i]
        regime_trees[[i]]$state <- regime_names[i]
        regime_cov[[i]] <- vcv.phylo(regime_trees[[i]])[tree$tip.label, tree$tip.label, drop=FALSE]
    }
    regime_cov
}

.mvgls_corrshrink_regime_cov <- function(object){
    if(!is.null(object$corrSt$regime_cov)){
        return(object$corrSt$regime_cov)
    }
    .mvgls_bmm_corrshrink_regimes(object$variables$tree)
}

.mvgls_corrshrink_omega <- function(object, sigma_object=object){
    regime_cov <- .mvgls_corrshrink_regime_cov(object)
    sigma_regime <- sigma_object$sigma$regime
    regime_names <- names(regime_cov)
    if(is.null(regime_names) || is.null(names(sigma_regime))){
        regime_cov <- unname(regime_cov)
        sigma_regime <- unname(sigma_regime)
    }else{
        sigma_regime <- sigma_regime[regime_names]
        if(any(vapply(sigma_regime, is.null, logical(1)))){
            stop("The corr-shrink sigma source is missing one or more regime covariance matrices")
        }
    }
    .Call(kroneckerSum,
        R=sigma_regime,
        C=regime_cov,
        Rrows=as.integer(object$dims$p),
        Crows=as.integer(object$dims$n),
        dimlist=as.integer(length(regime_cov))
    )
}

.mvgls_corrshrink_loglik <- function(Y, object, coefficients=object$coefficients, sigma_object=object){
    Y <- as.matrix(Y)
    if(!identical(dim(Y), dim(object$variables$Y))){
        stop("The response matrix does not match the fitted corr-shrink object dimensions")
    }
    coefficients <- as.matrix(coefficients)
    if(!identical(dim(coefficients), dim(object$coefficients))){
        stop("The coefficient matrix does not match the fitted corr-shrink object dimensions")
    }
    Omega <- .mvgls_corrshrink_omega(object, sigma_object=sigma_object)
    chol_omega <- try(chol(Omega), silent=TRUE)
    if(inherits(chol_omega, "try-error")) return(-Inf)
    residuals_vec <- as.numeric(Y - object$variables$X %*% coefficients)
    whitened <- try(backsolve(chol_omega, residuals_vec, transpose=TRUE), silent=TRUE)
    if(inherits(whitened, "try-error")) return(-Inf)
    logdet <- 2 * sum(log(diag(chol_omega)))
    -0.5 * (length(residuals_vec) * log(2 * pi) + logdet + sum(whitened^2))
}

.mvgls_corrshrink_refit <- function(object, Y, start=object$opt$par, echo=FALSE){
    if(is.null(object$model.frame)){
        stop("corr-shrink refits require the original model.frame on the fitted object")
    }
    refit_args <- list(
        formula=object$formula,
        data=object$model.frame,
        tree=object$variables$tree,
        model="BMM",
        method="LL",
        REML=FALSE,
        response=as.matrix(Y),
        bmm.structure="corrshrink",
        bmm.reference=object$reference_regime,
        grid.search=FALSE,
        start=start,
        echo=echo
    )
    if(!is.null(object$call$optimization)){
        refit_args$optimization <- if(is.character(object$call$optimization)) {
            object$call$optimization
        }else{
            eval(object$call$optimization, envir=environment(object$formula))
        }
    }
    do.call(mvgls, refit_args)
}

.mvgls_bmm_corrshrink_reference <- function(regime_names, bmm_reference=NULL){
    k <- length(regime_names)
    if(is.null(bmm_reference)){
        reference_index <- 1L
    }else if(is.numeric(bmm_reference)){
        reference_index <- as.integer(bmm_reference[1])
        if(length(reference_index) != 1L || is.na(reference_index) || reference_index < 1L || reference_index > k){
            stop("bmm.reference must be a valid regime index")
        }
    }else{
        reference_index <- match(as.character(bmm_reference[1]), regime_names)
        if(is.na(reference_index)){
            stop("bmm.reference must match one of the SIMMAP regime names")
        }
    }
    list(reference_index=reference_index, reference_regime=regime_names[reference_index])
}

.mvgls_bmm_corrshrink_start <- function(tree, Y, X, reference_index=1L){
    p <- ncol(Y)
    regime_names <- colnames(tree$mapped.edge)
    k <- length(regime_names)
    sigma_start <- startParamSigma(p, "cholesky", tree, Y)
    tip_values <- seq_len(Ntip(tree))
    index_tips <- tree$edge[, 2] %in% tip_values
    maps <- vapply(tree$maps[index_tips], function(x) names(x[length(x)]), character(1))
    if(length(unique(maps)) < k){
        scale_guess <- rep(mean(diag(rate_pic(tree, Y))), k)
    }else{
        scale_guess <- vapply(regime_names, function(map_name){
            dat_red <- which(maps == map_name)
            sp_to_remove <- tree$tip.label[!tree$tip.label %in% tree$tip.label[dat_red]]
            if(Ntip(tree) - length(sp_to_remove) <= 1 && length(sp_to_remove) > 0){
                sp_to_remove <- sp_to_remove[-sample(length(sp_to_remove), size=1)]
            }
            tree_red <- if(length(sp_to_remove) > 0) drop.tip(tree, sp_to_remove) else tree
            if(Ntip(tree_red) <= 1){
                sqrt(mean(apply(.rate_guess(tree, Y[tree$tip.label, , drop=FALSE], X[tree$tip.label, , drop=FALSE]), 2, var)))
            }else{
                sqrt(mean(apply(.rate_guess(tree_red, Y[tree_red$tip.label, , drop=FALSE], X[tree_red$tip.label, , drop=FALSE]), 2, var)))
            }
        }, numeric(1))
        scale_guess <- pmax(scale_guess^2, .Machine$double.eps)
    }
    ref_scale <- max(scale_guess[reference_index], .Machine$double.eps)
    rel_scales <- pmax(scale_guess[-reference_index]/ref_scale, .Machine$double.eps)
    c(sigma_start, sqrt(rel_scales), rep(stats::qlogis(0.95), k-1))
}

.mvgls_bmm_corrshrink_parameter_blocks <- function(p, regime_names, reference_index=1L){
    sigma_len <- p * (p + 1) / 2
    k <- length(regime_names)
    scale_idx <- if(k > 1) sigma_len + seq_len(k - 1) else integer(0)
    rho_idx <- if(k > 1) sigma_len + (k - 1) + seq_len(k - 1) else integer(0)
    free_regime_names <- if(k > 1) regime_names[-reference_index] else character(0)
    list(sigma_len=sigma_len, scale_idx=scale_idx, rho_idx=rho_idx, free_regime_names=free_regime_names)
}

.mvgls_bmm_corrshrink_start_pool <- function(base_start, p, regime_names, reference_index=1L){
    blocks <- .mvgls_bmm_corrshrink_parameter_blocks(p, regime_names, reference_index=reference_index)
    if(length(blocks$scale_idx) == 0L){
        start_table <- data.frame(
            start_id=1L,
            source="heuristic",
            scale_multiplier=NA_real_,
            rho_seed=NA_real_,
            stringsAsFactors=FALSE
        )
        return(list(starts=list(base_start), start_table=start_table))
    }
    scale_multipliers <- c(0.5, 1, 2)
    rho_seeds <- c(0.2, 0.5, 0.8, 0.95)
    starts <- vector("list", length(scale_multipliers) * length(rho_seeds))
    start_table <- data.frame(
        start_id=seq_along(starts),
        source="multistart",
        scale_multiplier=rep(scale_multipliers, each=length(rho_seeds)),
        rho_seed=rep(rho_seeds, times=length(scale_multipliers)),
        stringsAsFactors=FALSE
    )
    counter <- 1L
    for(scale_multiplier in scale_multipliers){
        for(rho_seed in rho_seeds){
            start_vec <- base_start
            start_vec[blocks$scale_idx] <- start_vec[blocks$scale_idx] * sqrt(scale_multiplier)
            start_vec[blocks$rho_idx] <- stats::qlogis(rho_seed)
            starts[[counter]] <- start_vec
            counter <- counter + 1L
        }
    }
    list(starts=starts, start_table=start_table)
}

.mvgls_bmm_corrshrink_unpack <- function(par, p, regime_names, reference_index=1L){
    blocks <- .mvgls_bmm_corrshrink_parameter_blocks(p, regime_names, reference_index=reference_index)
    sigma_len <- blocks$sigma_len
    k <- length(regime_names)
    theta_sigma <- par[seq_len(sigma_len)]
    theta_scale <- if(k > 1) par[blocks$scale_idx] else numeric(0)
    theta_rho <- if(k > 1) par[blocks$rho_idx] else numeric(0)
    scale <- rep(1, k)
    rho <- rep(1, k)
    if(k > 1){
        scale[-reference_index] <- theta_scale^2
        rho[-reference_index] <- stats::plogis(theta_rho)
    }
    names(scale) <- regime_names
    names(rho) <- regime_names
    list(theta_sigma=theta_sigma, scale=scale, rho=rho)
}

.mvgls_bmm_corrshrink_fit_components <- function(par, X, Y, regime_cov, reference_index=1L){
    n <- nrow(Y)
    p <- ncol(Y)
    trait_names <- colnames(Y)
    regime_names <- names(regime_cov)
    unpacked <- .mvgls_bmm_corrshrink_unpack(par, p, regime_names, reference_index=reference_index)
    R_base <- symPar(unpacked$theta_sigma, decomp="cholesky", p=p)
    V <- diag(diag(R_base))
    W <- R_base - V
    dimnames(R_base) <- list(trait_names, trait_names)
    dimnames(V) <- list(trait_names, trait_names)
    dimnames(W) <- list(trait_names, trait_names)
    C_var <- Reduce(`+`, Map(`*`, as.list(unpacked$scale), regime_cov))
    C_cov <- Reduce(`+`, Map(`*`, as.list(unpacked$scale * unpacked$rho), regime_cov))
    Omega <- .Call(kroneckerSum,
        R=list(V, W),
        C=list(C_var, C_cov),
        Rrows=as.integer(p),
        Crows=as.integer(n),
        dimlist=as.integer(2)
    )
    chol_omega <- try(chol(Omega), silent=TRUE)
    if(inherits(chol_omega, "try-error")) return(NULL)
    X_big <- kronecker(diag(p), X)
    y_vec <- as.numeric(Y)
    X_whitened <- try(backsolve(chol_omega, X_big, transpose=TRUE), silent=TRUE)
    y_whitened <- try(backsolve(chol_omega, y_vec, transpose=TRUE), silent=TRUE)
    if(inherits(X_whitened, "try-error") || inherits(y_whitened, "try-error")) return(NULL)
    XtX <- crossprod(X_whitened)
    XtY <- crossprod(X_whitened, y_whitened)
    beta_vec <- try(solve(XtX, XtY), silent=TRUE)
    if(inherits(beta_vec, "try-error")) return(NULL)
    residuals_whitened <- as.vector(y_whitened - X_whitened %*% beta_vec)
    logdet <- 2 * sum(log(diag(chol_omega)))
    nloglik <- 0.5 * (length(y_vec) * log(2 * pi) + logdet + sum(residuals_whitened^2))
    coefficients <- matrix(beta_vec, nrow=ncol(X), ncol=p, dimnames=list(colnames(X), colnames(Y)))
    fitted_values <- X %*% coefficients
    residuals_raw <- Y - fitted_values
    regime_vcv <- Map(function(scale, rho){
        scale * (V + rho * W)
    }, as.list(unpacked$scale), as.list(unpacked$rho))
    regime_vcv <- setNames(regime_vcv, regime_names)
    list(
        nloglik=nloglik,
        coefficients=coefficients,
        fitted=fitted_values,
        residuals=residuals_raw,
        omega=Omega,
        chol_omega=chol_omega,
        X_big=X_big,
        y_vec=y_vec,
        R_base=R_base,
        V=V,
        W=W,
        C_var=C_var,
        C_cov=C_cov,
        scale=unpacked$scale,
        rho=unpacked$rho,
        regime_vcv=regime_vcv,
        regime_cov=regime_cov
    )
}

.mvgls_bmm_corrshrink_fit_once <- function(start_par, objective, optimization, X, Y, regime_cov, reference_index=1L, hessian=FALSE){
    opt <- try(optim(
        par=start_par,
        fn=objective,
        method=optimization,
        hessian=isTRUE(hessian)
    ), silent=TRUE)
    fit <- if(inherits(opt, "try-error")) NULL else .mvgls_bmm_corrshrink_fit_components(opt$par, X=X, Y=Y, regime_cov=regime_cov, reference_index=reference_index)
    max_scale <- if(is.null(fit)) Inf else max(fit$scale)
    list(
        opt=if(inherits(opt, "try-error")) NULL else opt,
        fit=fit,
        nloglik=if(is.null(fit) || !is.finite(fit$nloglik)) Inf else fit$nloglik,
        max_scale=max_scale,
        par_norm=if(inherits(opt, "try-error")) Inf else sqrt(sum(opt$par^2))
    )
}

.mvgls_bmm_corrshrink_select_candidate <- function(candidates, tol=1e-6){
    if(length(candidates) == 0L) stop("No corr-shrink optimization candidates were provided")
    good_idx <- which(vapply(candidates, function(x){
        is.finite(x$nloglik) && !is.null(x$fit)
    }, logical(1)))
    if(length(good_idx) == 0L) stop("Optimization failed for all corr-shrink starting points")
    converged_idx <- good_idx[vapply(candidates[good_idx], function(x){
        !is.null(x$opt) && identical(x$opt$convergence, 0L)
    }, logical(1))]
    finite_idx <- if(length(converged_idx)) converged_idx else good_idx
    finite_scores <- vapply(candidates[finite_idx], function(x) x$nloglik, numeric(1))
    best_value <- min(finite_scores)
    tie_idx <- finite_idx[finite_scores <= best_value + tol]
    if(length(tie_idx) == 1L) return(tie_idx)
    tie_scale <- vapply(candidates[tie_idx], function(x) x$max_scale, numeric(1))
    tie_idx <- tie_idx[tie_scale == min(tie_scale)]
    if(length(tie_idx) == 1L) return(tie_idx)
    tie_norm <- vapply(candidates[tie_idx], function(x) x$par_norm, numeric(1))
    tie_idx[which.min(tie_norm)]
}

.mvgls_corrshrink_profile1d <- function(object, regime, param=c("rho", "scale"), grid, optimization=NULL, parameter=NULL){
    if(!is.null(parameter)) param <- parameter
    param <- match.arg(param)
    if(!isTRUE(!is.null(object$bmm.structure) && identical(object$bmm.structure, "corrshrink"))){
        stop("The profiling helper is only available for corr-shrink mvgls fits")
    }
    regime_names <- names(object$sigma$regime)
    regime_index <- match(regime, regime_names)
    if(is.na(regime_index)) stop("Unknown regime for corr-shrink profile")
    reference_index <- match(object$reference_regime, regime_names)
    if(is.na(reference_index)) stop("The corr-shrink fit is missing a valid reference_regime")
    if(regime_index == reference_index) stop("The reference regime is fixed and cannot be profiled")
    p <- object$dims$p
    blocks <- .mvgls_bmm_corrshrink_parameter_blocks(p, regime_names, reference_index=reference_index)
    free_regime_idx <- match(regime, blocks$free_regime_names)
    target_idx <- if(param == "scale") blocks$scale_idx[free_regime_idx] else blocks$rho_idx[free_regime_idx]
    fixed_transform <- if(param == "scale"){
        if(any(grid <= 0)) stop("Scale profile values must be strictly positive")
        sqrt(grid)
    }else{
        if(any(grid <= 0 | grid >= 1)) stop("Rho profile values must lie strictly between 0 and 1")
        stats::qlogis(grid)
    }
    objective_full <- function(par){
        fit <- .mvgls_bmm_corrshrink_fit_components(
            par,
            X=object$variables$X,
            Y=object$variables$Y,
            regime_cov=.mvgls_corrshrink_regime_cov(object),
            reference_index=reference_index
        )
        if(is.null(fit) || !is.finite(fit$nloglik)) return(1e25)
        fit$nloglik
    }
    if(is.null(optimization)){
        optimization <- if(is.null(object$call$optimization)) {
            "L-BFGS-B"
        }else if(is.character(object$call$optimization)) {
            object$call$optimization
        }else{
            eval(object$call$optimization, envir=environment(object$formula))
        }
    }
    profile_rows <- vector("list", length(grid))
    current_par <- object$opt$par
    free_idx <- setdiff(seq_along(current_par), target_idx)
    for(i in seq_along(grid)){
        current_par[target_idx] <- fixed_transform[i]
        opt <- optim(
            par=current_par[free_idx],
            fn=function(free_par){
                par_full <- current_par
                par_full[free_idx] <- free_par
                objective_full(par_full)
            },
            method=optimization
        )
        current_par[free_idx] <- opt$par
        fit <- .mvgls_bmm_corrshrink_fit_components(
            current_par,
            X=object$variables$X,
            Y=object$variables$Y,
            regime_cov=.mvgls_corrshrink_regime_cov(object),
            reference_index=reference_index
        )
        profile_rows[[i]] <- data.frame(
            regime=regime,
            param=param,
            fixed_value=grid[i],
            logLik=if(is.null(fit)) -Inf else -fit$nloglik,
            convergence=opt$convergence,
            max_scale=if(is.null(fit)) Inf else max(fit$scale),
            stringsAsFactors=FALSE
        )
        if(is.null(fit) || !is.finite(fit$nloglik)){
            current_par <- object$opt$par
        }
    }
    do.call(rbind, profile_rows)
}

.mvgls_corrshrink_profile <- function(object, parameter=c("rho", "scale"), regime, grid, optimization=NULL){
    .mvgls_corrshrink_profile1d(
        object=object,
        regime=regime,
        param=parameter,
        grid=grid,
        optimization=optimization
    )
}

.mvgls_bmm_corrshrink_profile <- function(object, parameter=c("rho", "scale"), regime, grid, optimization=NULL){
    .mvgls_corrshrink_profile(
        object=object,
        parameter=parameter,
        regime=regime,
        grid=grid,
        optimization=optimization
    )
}

.mvgls_bmm_corrshrink_object <- function(fit, formula, call_obj, model_frame, terms, xlevels, contrasts, variables, dims, start_values, method, optimization, tree, diagnostics, reference_regime){
    regime_names <- names(fit$scale)
    param <- as.vector(rbind(fit$scale, fit$rho))
    names(param) <- as.vector(rbind(paste0(regime_names, ".scale"), paste0(regime_names, ".rho")))
    results <- list(
        formula=formula,
        call=call_obj,
        model.frame=model_frame,
        coefficients=fit$coefficients,
        terms=terms,
        xlevels=xlevels,
        contrasts=contrasts,
        variables=variables,
        dims=dims,
        fitted=fit$fitted,
        logLik=-fit$opt$value,
        method=method,
        model="BMM",
        reference_regime=reference_regime,
        numIter=fit$opt$count[1],
        residuals=fit$residuals,
        sigma=list(
            Pinv=fit$R_base,
            P=solve(fit$R_base),
            S=NULL,
            regime=fit$regime_vcv
        ),
        tuning=NA,
        param=param,
        mserr=NA,
        start_values=start_values,
        corrSt=list(type="corrshrink", phy=tree, X=variables$X, Y=variables$Y, regime_cov=fit$regime_cov, reference_regime=reference_regime),
        penalty="LL",
        target="LL",
        REML=FALSE,
        FCI=NA,
        const_mtd=NA,
        opt=fit$opt,
        diagnostics=list(corrshrink=diagnostics),
        bmm.structure="corrshrink",
        df.free_cov=dims$p * (dims$p + 1) / 2,
        df.free_beta=length(fit$coefficients),
        df.free_model=2 * (length(regime_names) - 1)
    )
    results$df.free <- results$df.free_cov + results$df.free_beta + results$df.free_model
    results$regime.summary <- .mvgls_corrshrink_regime_summary(results)
    class(results) <- "mvgls"
    results
}

.fit_mvgls_bmm_corrshrink <- function(formula, call_obj, model_fr, X, Y, tree, terms, xlevels, contrasts, assign,
                                      qrx, method, REML, penalty, target, optimization, start,
                                      mserr, FCI, hessian, scale.height, bmm_reference=NULL){
    if(method != "LL") stop("The experimental corr-shrink BMM is available only with method=\"LL\"")
    if(isTRUE(REML)) stop("The experimental corr-shrink BMM does not support REML=TRUE")
    if(!is.null(mserr)) stop("The experimental corr-shrink BMM does not support measurement error")
    if(!identical(penalty, "RidgeArch")) stop("The experimental corr-shrink BMM does not support penalized settings")
    if(penalty == "EmpBayes") stop("The experimental corr-shrink BMM does not support EmpBayes")
    if(isTRUE(FCI)) stop("The experimental corr-shrink BMM does not support FCI")
    if(!inherits(tree, "simmap")) stop("A tree of class \"simmap\" is required for bmm.structure=\"corrshrink\"")
    if(scale.height) tree <- .scaleStruct(tree)
    regime_cov <- .mvgls_bmm_corrshrink_regimes(tree)
    reference <- .mvgls_bmm_corrshrink_reference(names(regime_cov), bmm_reference=bmm_reference)
    reference_index <- reference$reference_index
    base_start <- if(is.null(start)) .mvgls_bmm_corrshrink_start(tree, Y, X, reference_index=reference_index) else start
    objective <- function(par){
        fit <- .mvgls_bmm_corrshrink_fit_components(par, X=X, Y=Y, regime_cov=regime_cov, reference_index=reference_index)
        if(is.null(fit) || !is.finite(fit$nloglik)) return(1e25)
        fit$nloglik
    }
    regime_names <- names(regime_cov)
    start_pool <- if(is.null(start)) {
        .mvgls_bmm_corrshrink_start_pool(base_start, p=ncol(Y), regime_names=regime_names, reference_index=reference_index)
    }else{
        list(
            starts=list(base_start),
            start_table=data.frame(
                start_id=1L,
                source="user-provided",
                scale_multiplier=NA_real_,
                rho_seed=NA_real_,
                stringsAsFactors=FALSE
            )
        )
    }
    candidates <- lapply(start_pool$starts, function(start_par){
        .mvgls_bmm_corrshrink_fit_once(
            start_par=start_par,
            objective=objective,
            optimization=optimization,
            X=X,
            Y=Y,
            regime_cov=regime_cov,
            reference_index=reference_index,
            hessian=FALSE
        )
    })
    selected_idx <- .mvgls_bmm_corrshrink_select_candidate(candidates)
    selected <- candidates[[selected_idx]]
    opt <- selected$opt
    fit <- selected$fit
    if(isTRUE(as.logical(hessian))){
        opt <- optim(
            par=opt$par,
            fn=objective,
            method=optimization,
            hessian=TRUE
        )
        fit <- .mvgls_bmm_corrshrink_fit_components(opt$par, X=X, Y=Y, regime_cov=regime_cov, reference_index=reference_index)
    }
    if(!identical(opt$convergence, 0L)){
        warning("The experimental corr-shrink BMM optimizer reported convergence code ", opt$convergence)
    }
    if(is.null(fit) || !is.finite(fit$nloglik)) stop("Optimization failed for the experimental corr-shrink BMM fit")
    fit$opt <- opt
    start_table <- start_pool$start_table
    start_table$convergence <- vapply(candidates, function(x) if(is.null(x$opt)) NA_integer_ else x$opt$convergence, integer(1))
    start_table$nloglik <- vapply(candidates, function(x) x$nloglik, numeric(1))
    start_table$max_scale <- vapply(candidates, function(x) x$max_scale, numeric(1))
    start_table$selected <- start_table$start_id == selected_idx
    non_ref_rho <- fit$rho[-reference_index]
    diagnostics <- list(
        nstarts=nrow(start_table),
        selected_start_id=selected_idx,
        reference_regime=reference$reference_regime,
        start_table=start_table,
        boundary_rho=if(length(non_ref_rho) == 0L) FALSE else any(non_ref_rho < 0.02 | non_ref_rho > 0.98),
        pathological_scale=any(fit$scale > 20),
        max_scale=max(fit$scale)
    )
    dims <- list(
        n=nrow(Y),
        p=ncol(Y),
        m=qrx$rank,
        assign=assign,
        rank=qrx$rank,
        pivot=qrx$pivot,
        fullrank=(ncol(X) == qrx$rank)
    )
    variables <- list(Y=Y, X=X, tree=tree)
    .mvgls_bmm_corrshrink_object(
        fit=fit,
        formula=formula,
        call_obj=call_obj,
        model_frame=model_fr,
        terms=terms,
        xlevels=xlevels,
        contrasts=contrasts,
        variables=variables,
        dims=dims,
        start_values=start_pool$starts[[selected_idx]],
        method=method,
        optimization=optimization,
        tree=tree,
        diagnostics=diagnostics,
        reference_regime=reference$reference_regime
    )
}
