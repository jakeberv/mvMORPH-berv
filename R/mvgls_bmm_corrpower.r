################################################################################
##                                                                            ##
##              mvMORPH: mvgls_bmm_corrpower.r                                ##
##                                                                            ##
##   Experimental corr-power BMM fit for mvgls                                ##
##                                                                            ##
################################################################################

.mvgls_corrpower_regime_cov <- function(object){
    if(!is.null(object$corrSt$regime_cov)){
        return(object$corrSt$regime_cov)
    }
    .mvgls_bmmcorr_regimes(object$variables$tree)
}

.mvgls_corrpower_omega <- function(object, sigma_object=object){
    regime_cov <- .mvgls_corrpower_regime_cov(object)
    sigma_regime <- sigma_object$sigma$regime
    regime_names <- names(regime_cov)
    if(is.null(regime_names) || is.null(names(sigma_regime))){
        regime_cov <- unname(regime_cov)
        sigma_regime <- unname(sigma_regime)
    }else{
        sigma_regime <- sigma_regime[regime_names]
        if(any(vapply(sigma_regime, is.null, logical(1)))){
            stop("The corr-power sigma source is missing one or more regime covariance matrices")
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

.mvgls_corrpower_loglik <- function(Y, object, coefficients=object$coefficients, sigma_object=object){
    Y <- as.matrix(Y)
    if(!identical(dim(Y), dim(object$variables$Y))){
        stop("The response matrix does not match the fitted corr-power object dimensions")
    }
    coefficients <- as.matrix(coefficients)
    if(!identical(dim(coefficients), dim(object$coefficients))){
        stop("The coefficient matrix does not match the fitted corr-power object dimensions")
    }
    Omega <- .mvgls_corrpower_omega(object, sigma_object=sigma_object)
    chol_omega <- try(chol(Omega), silent=TRUE)
    if(inherits(chol_omega, "try-error")) return(-Inf)
    residuals_vec <- as.numeric(Y - object$variables$X %*% coefficients)
    whitened <- try(backsolve(chol_omega, residuals_vec, transpose=TRUE), silent=TRUE)
    if(inherits(whitened, "try-error")) return(-Inf)
    logdet <- 2 * sum(log(diag(chol_omega)))
    -0.5 * (length(residuals_vec) * log(2 * pi) + logdet + sum(whitened^2))
}

.mvgls_corrpower_refit <- function(object, Y, start=object$opt$par, echo=FALSE){
    if(is.null(object$model.frame)){
        stop("corr-power refits require the original model.frame on the fitted object")
    }
    refit_args <- list(
        formula=object$formula,
        data=object$model.frame,
        tree=object$variables$tree,
        model="BMM",
        method="LL",
        REML=FALSE,
        response=as.matrix(Y),
        bmm.structure="corrpower",
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

.mvgls_bmm_corrpower_safe_correlation <- function(shared_cov, trait_names=NULL){
    shared_cov <- 0.5 * (shared_cov + t(shared_cov))
    jitter <- 0
    for(i in 0:6){
        cov_try <- if(jitter > 0) shared_cov + diag(jitter, nrow(shared_cov)) else shared_cov
        corr_try <- try(cov2cor(cov_try), silent=TRUE)
        if(!inherits(corr_try, "try-error")){
            eig_try <- eigen(corr_try, symmetric=TRUE)
            if(min(eig_try$values) > .Machine$double.eps){
                if(!is.null(trait_names)){
                    dimnames(cov_try) <- list(trait_names, trait_names)
                    dimnames(corr_try) <- list(trait_names, trait_names)
                }
                return(list(shared_cov=cov_try, baseline_correlation=corr_try, eig=eig_try))
            }
        }
        jitter <- if(jitter == 0) 1e-8 else jitter * 10
    }
    stop("failed to construct a positive-definite baseline correlation matrix for corr-power BMM")
}

.mvgls_bmm_corrpower_template <- function(theta_sigma, p, trait_names=NULL){
    shared_cov <- symPar(theta_sigma, decomp="cholesky", p=p)
    diag_values <- diag(shared_cov)
    D <- diag(sqrt(diag_values))
    corr_info <- .mvgls_bmm_corrpower_safe_correlation(shared_cov, trait_names=trait_names)
    if(!is.null(trait_names)){
        dimnames(D) <- list(trait_names, trait_names)
    }
    list(
        theta_sigma=theta_sigma,
        D=D,
        shared_cov=corr_info$shared_cov,
        baseline_correlation=corr_info$baseline_correlation,
        eig=corr_info$eig
    )
}

.mvgls_bmm_corrpower_matrix_power <- function(template, corr_power, trait_names=NULL){
    eig <- template$eig
    values <- pmax(eig$values, .Machine$double.eps)
    corr_raw <- eig$vectors %*% diag(values^corr_power, nrow=length(values)) %*% t(eig$vectors)
    corr_raw <- 0.5 * (corr_raw + t(corr_raw))
    corr_power_mat <- cov2cor(corr_raw)
    if(!is.null(trait_names)){
        dimnames(corr_power_mat) <- list(trait_names, trait_names)
    }
    corr_power_mat
}

.mvgls_bmm_corrpower_start <- function(tree, Y, X, reference_index=1L){
    p <- ncol(Y)
    regime_names <- colnames(tree$mapped.edge)
    k <- length(regime_names)
    sigma_start <- startParamSigma(p, "cholesky", tree, Y)
    safe_scale_guess <- function(tree_local, Y_local, X_local){
        residuals <- try(.rate_guess(
            tree_local,
            Y_local[tree_local$tip.label, , drop=FALSE],
            X_local[tree_local$tip.label, , drop=FALSE]
        ), silent=TRUE)
        if(inherits(residuals, "try-error") || any(!is.finite(residuals))){
            coef_guess <- try(pseudoinverse(X_local[tree_local$tip.label, , drop=FALSE]) %*% Y_local[tree_local$tip.label, , drop=FALSE], silent=TRUE)
            if(inherits(coef_guess, "try-error") || any(!is.finite(coef_guess))){
                residuals <- scale(Y_local[tree_local$tip.label, , drop=FALSE], center=TRUE, scale=FALSE)
            }else{
                residuals <- Y_local[tree_local$tip.label, , drop=FALSE] - X_local[tree_local$tip.label, , drop=FALSE] %*% coef_guess
            }
        }
        vars <- suppressWarnings(apply(residuals, 2, stats::var, na.rm=TRUE))
        vars <- vars[is.finite(vars) & vars > 0]
        if(!length(vars)){
            vars <- diag(stats::cov(Y_local[tree_local$tip.label, , drop=FALSE], use="pairwise.complete.obs"))
            vars <- vars[is.finite(vars) & vars > 0]
        }
        if(!length(vars)) return(1)
        sqrt(mean(vars))
    }
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
            if(Ntip(tree_red) <= 1) safe_scale_guess(tree, Y, X) else safe_scale_guess(tree_red, Y, X)
        }, numeric(1))
        scale_guess <- pmax(scale_guess^2, .Machine$double.eps)
    }
    ref_scale <- max(scale_guess[reference_index], .Machine$double.eps)
    rel_scales <- pmax(scale_guess[-reference_index] / ref_scale, .Machine$double.eps)
    c(
        sigma_start,
        log(rel_scales),
        rep(log(1), k - 1L)
    )
}

.mvgls_bmm_corrpower_parameter_blocks <- function(p, regime_names, reference_index=1L){
    sigma_len <- p * (p + 1) / 2
    k <- length(regime_names)
    scale_idx <- if(k > 1) sigma_len + seq_len(k - 1) else integer(0)
    corr_power_idx <- if(k > 1) sigma_len + length(scale_idx) + seq_len(k - 1) else integer(0)
    free_regime_names <- if(k > 1) regime_names[-reference_index] else character(0)
    list(
        sigma_len=sigma_len,
        scale_idx=scale_idx,
        corr_power_idx=corr_power_idx,
        free_regime_names=free_regime_names
    )
}

.mvgls_bmm_corrpower_regularization <- function(){
    lambda_scale <- getOption("mvMORPH.corrpower.lambda_scale", 0.05)
    lambda_corr_power <- getOption("mvMORPH.corrpower.lambda_corr_power", 0.05)
    lambda_scale <- suppressWarnings(as.numeric(lambda_scale[1]))
    lambda_corr_power <- suppressWarnings(as.numeric(lambda_corr_power[1]))
    if(!is.finite(lambda_scale) || lambda_scale < 0) lambda_scale <- 0.05
    if(!is.finite(lambda_corr_power) || lambda_corr_power < 0) lambda_corr_power <- 0.05
    list(
        lambda_scale=lambda_scale,
        lambda_corr_power=lambda_corr_power
    )
}

.mvgls_bmm_corrpower_regularization_from_object <- function(object){
    diag_info <- object$diagnostics$corrpower
    lambda_scale <- diag_info$regularization_scale
    lambda_corr_power <- diag_info$regularization_corr_power
    if(is.null(lambda_scale) || !is.finite(lambda_scale)) lambda_scale <- 0
    if(is.null(lambda_corr_power) || !is.finite(lambda_corr_power)) lambda_corr_power <- 0
    list(
        lambda_scale=lambda_scale,
        lambda_corr_power=lambda_corr_power
    )
}

.mvgls_bmm_corrpower_penalty <- function(par, p, regime_names, reference_index=1L,
                                         regularization=.mvgls_bmm_corrpower_regularization()){
    blocks <- .mvgls_bmm_corrpower_parameter_blocks(p, regime_names, reference_index=reference_index)
    penalty <- 0
    if(length(blocks$scale_idx) && regularization$lambda_scale > 0){
        penalty <- penalty + regularization$lambda_scale * sum(par[blocks$scale_idx]^2)
    }
    if(length(blocks$corr_power_idx) && regularization$lambda_corr_power > 0){
        penalty <- penalty + regularization$lambda_corr_power * sum(par[blocks$corr_power_idx]^2)
    }
    penalty
}

.mvgls_bmm_corrpower_start_pool <- function(base_start, p, regime_names, reference_index=1L){
    blocks <- .mvgls_bmm_corrpower_parameter_blocks(p, regime_names, reference_index=reference_index)
    if(length(blocks$scale_idx) == 0L){
        start_table <- data.frame(
            start_id=1L,
            source="heuristic",
            scale_multiplier=NA_real_,
            corr_power_seed=NA_real_,
            stringsAsFactors=FALSE
        )
        return(list(starts=list(base_start), start_table=start_table))
    }
    scale_multipliers <- c(0.5, 1, 2)
    corr_power_seeds <- c(0.25, 0.75, 1.0, 1.5)
    starts <- vector("list", length(scale_multipliers) * length(corr_power_seeds))
    start_table <- data.frame(
        start_id=seq_along(starts),
        source="multistart",
        scale_multiplier=rep(scale_multipliers, each=length(corr_power_seeds)),
        corr_power_seed=rep(corr_power_seeds, times=length(scale_multipliers)),
        stringsAsFactors=FALSE
    )
    counter <- 1L
    for(scale_multiplier in scale_multipliers){
        for(corr_power_seed in corr_power_seeds){
            start_vec <- base_start
            start_vec[blocks$scale_idx] <- start_vec[blocks$scale_idx] + log(scale_multiplier)
            start_vec[blocks$corr_power_idx] <- rep(log(corr_power_seed), length(blocks$corr_power_idx))
            starts[[counter]] <- start_vec
            counter <- counter + 1L
        }
    }
    list(starts=starts, start_table=start_table)
}

.mvgls_bmm_corrpower_unpack <- function(par, p, regime_names, reference_index=1L){
    blocks <- .mvgls_bmm_corrpower_parameter_blocks(p, regime_names, reference_index=reference_index)
    k <- length(regime_names)
    theta_sigma <- par[seq_len(blocks$sigma_len)]
    theta_scale <- if(length(blocks$scale_idx)) par[blocks$scale_idx] else numeric(0)
    theta_corr_power <- if(length(blocks$corr_power_idx)) par[blocks$corr_power_idx] else numeric(0)
    template <- .mvgls_bmm_corrpower_template(theta_sigma, p=p)
    scale <- rep(1, k)
    corr_power <- rep(1, k)
    if(length(blocks$scale_idx)){
        scale[-reference_index] <- exp(theta_scale)
        corr_power[-reference_index] <- exp(theta_corr_power)
    }
    names(scale) <- regime_names
    names(corr_power) <- regime_names
    list(
        theta_sigma=theta_sigma,
        scale=scale,
        corr_power=corr_power,
        template=template
    )
}

.mvgls_bmm_corrpower_fit_components <- function(par, X, Y, regime_cov, reference_index=1L){
    n <- nrow(Y)
    p <- ncol(Y)
    trait_names <- colnames(Y)
    regime_names <- names(regime_cov)
    unpacked <- .mvgls_bmm_corrpower_unpack(par, p, regime_names, reference_index=reference_index)
    template <- unpacked$template
    D <- template$D
    shared_cov <- template$shared_cov
    dimnames(shared_cov) <- list(trait_names, trait_names)
    dimnames(template$baseline_correlation) <- list(trait_names, trait_names)
    regime_vcv <- Map(function(scale, corr_power){
        corr_mat <- .mvgls_bmm_corrpower_matrix_power(template, corr_power, trait_names=trait_names)
        scale * (D %*% corr_mat %*% D)
    }, as.list(unpacked$scale), as.list(unpacked$corr_power))
    regime_vcv <- setNames(regime_vcv, regime_names)
    Omega <- .Call(kroneckerSum,
        R=regime_vcv,
        C=regime_cov,
        Rrows=as.integer(p),
        Crows=as.integer(n),
        dimlist=as.integer(length(regime_cov))
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
    list(
        nloglik=nloglik,
        coefficients=coefficients,
        fitted=fitted_values,
        residuals=residuals_raw,
        omega=Omega,
        chol_omega=chol_omega,
        X_big=X_big,
        y_vec=y_vec,
        shared_cov=shared_cov,
        baseline_correlation=template$baseline_correlation,
        scale=unpacked$scale,
        corr_power=unpacked$corr_power,
        regime_vcv=regime_vcv,
        regime_cov=regime_cov
    )
}

.mvgls_bmm_corrpower_fit_once <- function(start_par, objective, optimization, X, Y, regime_cov,
                                         reference_index=1L, hessian=FALSE, regularization=.mvgls_bmm_corrpower_regularization()){
    opt <- try(optim(
        par=start_par,
        fn=objective,
        method=optimization,
        hessian=isTRUE(hessian)
    ), silent=TRUE)
    regime_names <- names(regime_cov)
    p <- ncol(Y)
    fit <- if(inherits(opt, "try-error")) NULL else .mvgls_bmm_corrpower_fit_components(opt$par, X=X, Y=Y, regime_cov=regime_cov, reference_index=reference_index)
    penalty_value <- if(is.null(fit)) Inf else .mvgls_bmm_corrpower_penalty(
        if(inherits(opt, "try-error")) start_par else opt$par,
        p=p,
        regime_names=regime_names,
        reference_index=reference_index,
        regularization=regularization
    )
    max_scale <- if(is.null(fit)) Inf else max(fit$scale)
    max_corr_power <- if(is.null(fit)) Inf else max(fit$corr_power)
    list(
        opt=if(inherits(opt, "try-error")) NULL else opt,
        fit=fit,
        nloglik=if(is.null(fit) || !is.finite(fit$nloglik)) Inf else fit$nloglik,
        objective_value=if(is.null(fit) || !is.finite(fit$nloglik) || !is.finite(penalty_value)) Inf else fit$nloglik + penalty_value,
        penalty_value=penalty_value,
        max_scale=max_scale,
        max_corr_power=max_corr_power,
        par_norm=if(inherits(opt, "try-error")) Inf else sqrt(sum(opt$par^2))
    )
}

.mvgls_bmm_corrpower_select_candidate <- function(candidates, tol=1e-6){
    if(length(candidates) == 0L) stop("No corr-power optimization candidates were provided")
    good_idx <- which(vapply(candidates, function(x){
        is.finite(x$nloglik) && !is.null(x$fit)
    }, logical(1)))
    if(length(good_idx) == 0L) stop("Optimization failed for all corr-power starting points")
    converged_idx <- good_idx[vapply(candidates[good_idx], function(x){
        !is.null(x$opt) && identical(x$opt$convergence, 0L)
    }, logical(1))]
    finite_idx <- if(length(converged_idx)) converged_idx else good_idx
    finite_scores <- vapply(candidates[finite_idx], function(x) x$objective_value, numeric(1))
    best_value <- min(finite_scores)
    tie_idx <- finite_idx[finite_scores <= best_value + tol]
    if(length(tie_idx) == 1L) return(tie_idx)
    tie_scale <- vapply(candidates[tie_idx], function(x) x$max_scale, numeric(1))
    tie_idx <- tie_idx[tie_scale == min(tie_scale)]
    if(length(tie_idx) == 1L) return(tie_idx)
    tie_power <- vapply(candidates[tie_idx], function(x) x$max_corr_power, numeric(1))
    tie_idx <- tie_idx[tie_power == min(tie_power)]
    if(length(tie_idx) == 1L) return(tie_idx)
    tie_norm <- vapply(candidates[tie_idx], function(x) x$par_norm, numeric(1))
    tie_idx[which.min(tie_norm)]
}

.mvgls_corrpower_profile1d <- function(object, regime, param=c("corr_power", "scale"), grid, optimization=NULL, parameter=NULL){
    if(!is.null(parameter)) param <- parameter
    param <- match.arg(param, c("corr_power", "scale"))
    regime_names <- names(object$sigma$regime)
    regime_index <- match(regime, regime_names)
    if(is.na(regime_index)) stop("Unknown regime for corr-power profile")
    reference_index <- match(object$reference_regime, regime_names)
    if(is.na(reference_index)) stop("The corr-power fit is missing a valid reference_regime")
    if(regime_index == reference_index){
        stop("The reference regime parameters are fixed and cannot be profiled")
    }
    p <- object$dims$p
    blocks <- .mvgls_bmm_corrpower_parameter_blocks(p, regime_names, reference_index=reference_index)
    target_idx <- if(param == "scale"){
        blocks$scale_idx[match(regime, blocks$free_regime_names)]
    }else{
        blocks$corr_power_idx[match(regime, blocks$free_regime_names)]
    }
    if(any(grid <= 0)) stop("corr-power profile values must be strictly positive")
    objective_full <- function(par){
        fit <- .mvgls_bmm_corrpower_fit_components(
            par,
            X=object$variables$X,
            Y=object$variables$Y,
            regime_cov=.mvgls_corrpower_regime_cov(object),
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
    regularization <- .mvgls_bmm_corrpower_regularization_from_object(object)
    for(i in seq_along(grid)){
        opt <- optim(
            par=current_par[free_idx],
            fn=function(free_par){
                par_full <- current_par
                par_full[free_idx] <- free_par
                if(param == "scale"){
                    par_full[target_idx] <- log(grid[i])
                }else{
                    par_full[target_idx] <- log(grid[i])
                }
                objective_full(par_full) + .mvgls_bmm_corrpower_penalty(
                    par_full,
                    p=p,
                    regime_names=regime_names,
                    reference_index=reference_index,
                    regularization=regularization
                )
            },
            method=optimization
        )
        current_par[free_idx] <- opt$par
        if(param == "scale"){
            current_par[target_idx] <- log(grid[i])
        }else{
            current_par[target_idx] <- log(grid[i])
        }
        fit <- .mvgls_bmm_corrpower_fit_components(
            current_par,
            X=object$variables$X,
            Y=object$variables$Y,
            regime_cov=.mvgls_corrpower_regime_cov(object),
            reference_index=reference_index
        )
        profile_rows[[i]] <- data.frame(
            regime=regime,
            param=param,
            fixed_value=grid[i],
            logLik=if(is.null(fit)) -Inf else -fit$nloglik,
            convergence=opt$convergence,
            max_scale=if(is.null(fit)) Inf else max(fit$scale),
            max_corr_power=if(is.null(fit)) Inf else max(fit$corr_power),
            stringsAsFactors=FALSE
        )
        if(is.null(fit) || !is.finite(fit$nloglik)){
            current_par <- object$opt$par
        }
    }
    do.call(rbind, profile_rows)
}

.mvgls_corrpower_profile <- function(object, parameter=c("corr_power", "scale"), regime, grid, optimization=NULL){
    .mvgls_corrpower_profile1d(
        object=object,
        regime=regime,
        param=parameter,
        grid=grid,
        optimization=optimization
    )
}

.mvgls_bmm_corrpower_profile <- function(object, parameter=c("corr_power", "scale"), regime, grid, optimization=NULL){
    .mvgls_corrpower_profile(
        object=object,
        parameter=parameter,
        regime=regime,
        grid=grid,
        optimization=optimization
    )
}

.mvgls_bmm_corrpower_object <- function(fit, formula, call_obj, model_frame, terms, xlevels, contrasts, variables, dims, start_values, method, optimization, tree, diagnostics, reference_regime){
    regime_names <- names(fit$scale)
    param <- as.vector(rbind(fit$scale, fit$corr_power))
    names(param) <- as.vector(rbind(paste0(regime_names, ".scale"), paste0(regime_names, ".corr_power")))
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
        logLik=-fit$nloglik,
        method=method,
        model="BMM",
        reference_regime=reference_regime,
        numIter=if(is.null(fit$opt$count)) NA_integer_ else fit$opt$count[1],
        residuals=fit$residuals,
        sigma=list(
            Pinv=fit$shared_cov,
            P=.mvgls_bmmcorr_safe_inverse(fit$shared_cov),
            S=NULL,
            baseline_correlation=fit$baseline_correlation,
            regime=fit$regime_vcv
        ),
        tuning=NA,
        param=param,
        mserr=NA,
        start_values=start_values,
        corrSt=list(
            type="corrpower",
            phy=tree,
            X=variables$X,
            Y=variables$Y,
            regime_cov=fit$regime_cov,
            reference_regime=reference_regime
        ),
        penalty="LL",
        target="LL",
        REML=FALSE,
        FCI=NA,
        const_mtd=NA,
        opt=fit$opt,
        diagnostics=list(corrpower=diagnostics),
        bmm.structure="corrpower",
        df.free_cov=dims$p * (dims$p + 1) / 2,
        df.free_beta=length(fit$coefficients),
        df.free_model=2 * (length(regime_names) - 1),
        regime.summary=NULL
    )
    results$df.free <- results$df.free_cov + results$df.free_beta + results$df.free_model
    results$regime.summary <- .mvgls_bmmcorr_regime_summary(results)
    class(results) <- "mvgls"
    results
}

.fit_mvgls_bmm_corrpower <- function(formula, call_obj, model_fr, X, Y, tree, terms, xlevels, contrasts, assign,
                                     qrx, method, REML, penalty, target, optimization, start,
                                     mserr, FCI, hessian, scale.height, bmm_reference=NULL){
    if(method != "LL") stop("The experimental corr-power BMM is available only with method=\"LL\"")
    if(isTRUE(REML)) stop("The experimental corr-power BMM does not support REML=TRUE")
    if(!is.null(mserr)) stop("The experimental corr-power BMM does not support measurement error")
    if(!identical(penalty, "RidgeArch")) stop("The experimental corr-power BMM does not support penalized settings")
    if(penalty == "EmpBayes") stop("The experimental corr-power BMM does not support EmpBayes")
    if(isTRUE(FCI)) stop("The experimental corr-power BMM does not support FCI")
    if(!inherits(tree, "simmap")) stop("A tree of class \"simmap\" is required for bmm.structure=\"corrpower\"")
    if(scale.height) tree <- .scaleStruct(tree)
    regime_cov <- .mvgls_bmmcorr_regimes(tree)
    reference <- .mvgls_bmmcorr_reference(names(regime_cov), bmm_reference=bmm_reference)
    reference_index <- reference$reference_index
    regularization <- .mvgls_bmm_corrpower_regularization()
    base_start <- if(is.null(start)) .mvgls_bmm_corrpower_start(tree, Y, X, reference_index=reference_index) else start
    regime_names <- names(regime_cov)
    objective <- function(par){
        fit <- .mvgls_bmm_corrpower_fit_components(par, X=X, Y=Y, regime_cov=regime_cov, reference_index=reference_index)
        if(is.null(fit) || !is.finite(fit$nloglik)) return(1e25)
        fit$nloglik + .mvgls_bmm_corrpower_penalty(
            par,
            p=ncol(Y),
            regime_names=regime_names,
            reference_index=reference_index,
            regularization=regularization
        )
    }
    start_pool <- if(is.null(start)) {
        .mvgls_bmm_corrpower_start_pool(base_start, p=ncol(Y), regime_names=regime_names, reference_index=reference_index)
    }else{
        list(
            starts=list(base_start),
            start_table=data.frame(
                start_id=1L,
                source="user-provided",
                scale_multiplier=NA_real_,
                corr_power_seed=NA_real_,
                stringsAsFactors=FALSE
            )
        )
    }
    candidates <- lapply(start_pool$starts, function(start_par){
        .mvgls_bmm_corrpower_fit_once(
            start_par=start_par,
            objective=objective,
            optimization=optimization,
            X=X,
            Y=Y,
            regime_cov=regime_cov,
            reference_index=reference_index,
            regularization=regularization,
            hessian=FALSE
        )
    })
    selected_idx <- .mvgls_bmm_corrpower_select_candidate(candidates)
    selected <- candidates[[selected_idx]]
    opt <- selected$opt
    fit <- selected$fit
    opt_polish <- optim(
        par=opt$par,
        fn=objective,
        method=optimization,
        hessian=isTRUE(as.logical(hessian))
    )
    fit_polish <- .mvgls_bmm_corrpower_fit_components(opt_polish$par, X=X, Y=Y, regime_cov=regime_cov, reference_index=reference_index)
    if(!is.null(fit_polish) && is.finite(fit_polish$nloglik) &&
       (is.null(fit) || !is.finite(fit$nloglik) ||
        (fit_polish$nloglik + .mvgls_bmm_corrpower_penalty(
            opt_polish$par,
            p=ncol(Y),
            regime_names=regime_names,
            reference_index=reference_index,
            regularization=regularization
        )) <= (fit$nloglik + selected$penalty_value + 1e-8))){
        opt <- opt_polish
        fit <- fit_polish
    }
    if(!identical(opt$convergence, 0L)){
        warning("The experimental corr-power BMM optimizer reported convergence code ", opt$convergence)
    }
    if(is.null(fit) || !is.finite(fit$nloglik)) stop("Optimization failed for the experimental corr-power BMM fit")
    opt$penalized_value <- opt$value
    opt$gaussian_nloglik <- fit$nloglik
    fit$opt <- opt
    start_table <- start_pool$start_table
    start_table$convergence <- vapply(candidates, function(x) if(is.null(x$opt)) NA_integer_ else x$opt$convergence, integer(1))
    start_table$nloglik <- vapply(candidates, function(x) x$nloglik, numeric(1))
    start_table$objective_value <- vapply(candidates, function(x) x$objective_value, numeric(1))
    start_table$penalty_value <- vapply(candidates, function(x) x$penalty_value, numeric(1))
    start_table$max_scale <- vapply(candidates, function(x) x$max_scale, numeric(1))
    start_table$max_corr_power <- vapply(candidates, function(x) x$max_corr_power, numeric(1))
    start_table$selected <- start_table$start_id == selected_idx
    non_reference_power <- fit$corr_power[-reference_index]
    diagnostics <- list(
        nstarts=nrow(start_table),
        selected_start_id=selected_idx,
        reference_regime=reference$reference_regime,
        start_table=start_table,
        boundary_corr_power=any(non_reference_power < 0.05 | non_reference_power > 3),
        pathological_scale=any(fit$scale > 20),
        pathological_corr_power=any(non_reference_power > 4),
        regularization_scale=regularization$lambda_scale,
        regularization_corr_power=regularization$lambda_corr_power,
        max_scale=max(fit$scale),
        max_corr_power=max(fit$corr_power),
        corr_power_max=max(fit$corr_power)
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
    .mvgls_bmm_corrpower_object(
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
