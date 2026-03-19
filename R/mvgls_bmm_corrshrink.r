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

.mvgls_bmm_corrshrink_start <- function(tree, Y, X){
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
    ref_scale <- max(scale_guess[1], .Machine$double.eps)
    rel_scales <- pmax(scale_guess[-1]/ref_scale, .Machine$double.eps)
    c(sigma_start, sqrt(rel_scales), rep(stats::qlogis(0.95), k-1))
}

.mvgls_bmm_corrshrink_unpack <- function(par, p, regime_names){
    sigma_len <- p * (p + 1) / 2
    k <- length(regime_names)
    theta_sigma <- par[seq_len(sigma_len)]
    theta_scale <- if(k > 1) par[sigma_len + seq_len(k - 1)] else numeric(0)
    theta_rho <- if(k > 1) par[sigma_len + (k - 1) + seq_len(k - 1)] else numeric(0)
    scale <- c(1, theta_scale^2)
    rho <- c(1, stats::plogis(theta_rho))
    names(scale) <- regime_names
    names(rho) <- regime_names
    list(theta_sigma=theta_sigma, scale=scale, rho=rho)
}

.mvgls_bmm_corrshrink_fit_components <- function(par, X, Y, regime_cov){
    n <- nrow(Y)
    p <- ncol(Y)
    trait_names <- colnames(Y)
    regime_names <- names(regime_cov)
    unpacked <- .mvgls_bmm_corrshrink_unpack(par, p, regime_names)
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

.mvgls_bmm_corrshrink_object <- function(fit, formula, call_obj, model_frame, terms, xlevels, contrasts, variables, dims, start_values, method, optimization, tree){
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
        corrSt=list(type="corrshrink", phy=tree, X=variables$X, Y=variables$Y, regime_cov=fit$regime_cov),
        penalty="LL",
        target="LL",
        REML=FALSE,
        FCI=NA,
        const_mtd=NA,
        opt=fit$opt,
        bmm.structure="corrshrink",
        df.free_cov=dims$p * (dims$p + 1) / 2,
        df.free_beta=length(fit$coefficients),
        df.free_model=2 * (length(regime_names) - 1)
    )
    results$df.free <- results$df.free_cov + results$df.free_beta + results$df.free_model
    class(results) <- "mvgls"
    results
}

.fit_mvgls_bmm_corrshrink <- function(formula, call_obj, model_fr, X, Y, tree, terms, xlevels, contrasts, assign,
                                      qrx, method, REML, penalty, target, optimization, start,
                                      mserr, FCI, hessian, scale.height){
    if(method != "LL") stop("The experimental corr-shrink BMM is available only with method=\"LL\"")
    if(isTRUE(REML)) stop("The experimental corr-shrink BMM does not support REML=TRUE")
    if(!is.null(mserr)) stop("The experimental corr-shrink BMM does not support measurement error")
    if(!identical(penalty, "RidgeArch")) stop("The experimental corr-shrink BMM does not support penalized settings")
    if(penalty == "EmpBayes") stop("The experimental corr-shrink BMM does not support EmpBayes")
    if(isTRUE(FCI)) stop("The experimental corr-shrink BMM does not support FCI")
    if(!inherits(tree, "simmap")) stop("A tree of class \"simmap\" is required for bmm.structure=\"corrshrink\"")
    if(scale.height) tree <- .scaleStruct(tree)
    regime_cov <- .mvgls_bmm_corrshrink_regimes(tree)
    opt_start <- if(is.null(start)) .mvgls_bmm_corrshrink_start(tree, Y, X) else start
    objective <- function(par){
        fit <- .mvgls_bmm_corrshrink_fit_components(par, X=X, Y=Y, regime_cov=regime_cov)
        if(is.null(fit) || !is.finite(fit$nloglik)) return(1e25)
        fit$nloglik
    }
    opt <- optim(
        par=opt_start,
        fn=objective,
        method=optimization,
        hessian=isTRUE(as.logical(hessian))
    )
    if(!identical(opt$convergence, 0L)){
        warning("The experimental corr-shrink BMM optimizer reported convergence code ", opt$convergence)
    }
    fit <- .mvgls_bmm_corrshrink_fit_components(opt$par, X=X, Y=Y, regime_cov=regime_cov)
    if(is.null(fit) || !is.finite(fit$nloglik)) stop("Optimization failed for the experimental corr-shrink BMM fit")
    fit$opt <- opt
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
        start_values=opt_start,
        method=method,
        optimization=optimization,
        tree=tree
    )
}
