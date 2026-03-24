################################################################################
##                                                                            ##
##                       mvMORPH: classes_methods.r                           ##
##                                                                            ##
##   Internal functions for the mvMORPH package                               ##
##                                                                            ##
##  Created by Julien Clavel - 31-07-2018                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@biologie.ens.fr)                 ##
##   require: phytools, ape, corpcor, subplex, spam, glassoFast, stats        ##
##                                                                            ##
################################################################################


# ------------ S3 Methods ------------------------------------------------- #
GIC <- function(object, ...) UseMethod("GIC")
EIC <- function(object, nboot=100L, nbcores=1L, ...) UseMethod("EIC")
corrstrength_diagnostics <- function(object, ...) UseMethod("corrstrength_diagnostics")
corrpower_diagnostics <- function(object, ...) UseMethod("corrpower_diagnostics")
`%||%` <- function(x, y) if(is.null(x)) y else x

.mvgls_bmmcorr_structure <- function(object){
    structure <- object$bmm.structure %||% NULL
    if(is.null(structure) || !structure %in% c("corrstrength", "corrpower")) return(NULL)
    structure
}

.mvgls_is_corrstrength <- function(object){
    identical(.mvgls_bmmcorr_structure(object), "corrstrength")
}

.mvgls_is_corrpower <- function(object){
    identical(.mvgls_bmmcorr_structure(object), "corrpower")
}

.mvgls_is_bmmcorr <- function(object){
    !is.null(.mvgls_bmmcorr_structure(object))
}

.mvgls_bmmcorr_label <- function(object=NULL, structure=NULL, capitalize=FALSE){
    structure <- structure %||% if(!is.null(object)) .mvgls_bmmcorr_structure(object) else NULL
    label <- switch(structure,
        corrstrength="corr-strength",
        corrpower="corr-power",
        "experimental correlation"
    )
    if(capitalize) {
        paste0(toupper(substring(label, 1L, 1L)), substring(label, 2L))
    } else {
        label
    }
}

.mvgls_bmmcorr_cor_param_name <- function(object=NULL, structure=NULL){
    structure <- structure %||% if(!is.null(object)) .mvgls_bmmcorr_structure(object) else "corrpower"
    switch(structure,
        corrpower="corr_power",
        corrstrength="corr_strength",
        "corr_strength"
    )
}

.mvgls_bmmcorr_regime_param <- function(object, regime_name, param_name){
    if(is.null(object$param)) return(NA_real_)
    key <- paste0(regime_name, ".", param_name)
    if(!key %in% names(object$param)) return(NA_real_)
    unname(object$param[[key]])
}

.mvgls_bmmcorr_primary_param_names <- function(object, include_reference=TRUE){
    if(is.null(object$sigma) || is.null(object$sigma$regime)) return(character(0))
    cor_param <- .mvgls_bmmcorr_cor_param_name(object)
    regime_names <- names(object$sigma$regime)
    if(is.null(regime_names)) regime_names <- character(0)
    if(!include_reference && !is.null(object$reference_regime)){
        regime_names <- setdiff(regime_names, object$reference_regime)
    }
    as.vector(rbind(paste0(regime_names, ".scale"), paste0(regime_names, ".", cor_param)))
}

.mvgls_bmmcorr_parse_parm <- function(object, parm){
    available <- .mvgls_bmmcorr_primary_param_names(object, include_reference=TRUE)
    if(is.null(parm)) return(available)
    parm <- as.character(parm)
    if(any(!parm %in% available)){
        stop("Unknown experimental BMM correlation parameter requested in 'parm'")
    }
    parm
}

.mvgls_bmmcorr_regime_covariances <- function(object, tree=object$variables$tree, internal=FALSE){
    regime_cov <- NULL
    same_tree <- inherits(tree, "phylo") &&
        inherits(object$variables$tree, "phylo") &&
        identical(tree$tip.label, object$variables$tree$tip.label)

    if(!internal && same_tree && !is.null(object$corrSt$regime_cov)){
        regime_cov <- object$corrSt$regime_cov
    }else{
        regime_cov <- vcvSplit(tree, internal=internal)
        names(regime_cov) <- colnames(tree$mapped.edge)
    }

    sigma_regime <- object$sigma$regime
    if(is.null(names(sigma_regime))){
        names(sigma_regime) <- names(regime_cov)
    }
    regime_names <- names(sigma_regime)
    regime_cov <- regime_cov[regime_names]
    if(any(vapply(regime_cov, is.null, logical(1)))){
        stop("Missing one or more regime covariance matrices for the ", .mvgls_bmmcorr_label(object), " fit")
    }
    regime_cov
}

.mvgls_bmmcorr_omega_from_tree <- function(object, tree=object$variables$tree, internal=FALSE){
    regime_cov <- .mvgls_bmmcorr_regime_covariances(object, tree=tree, internal=internal)
    sigma_regime <- object$sigma$regime[names(regime_cov)]
    .Call(kroneckerSum,
        R=sigma_regime,
        C=regime_cov,
        Rrows=as.integer(object$dims$p),
        Crows=as.integer(nrow(regime_cov[[1]])),
        dimlist=as.integer(length(regime_cov))
    )
}

.mvgls_bmmcorr_block_index <- function(base_index, block_size, p){
    as.integer(unlist(lapply(seq_len(p), function(i) base_index + (i - 1L) * block_size), use.names=FALSE))
}

.mvgls_bmmcorr_conditional_weights <- function(cov_train, rhs){
    chol_train <- try(chol(cov_train), silent=TRUE)
    if(!inherits(chol_train, "try-error")){
        tmp <- forwardsolve(t(chol_train), rhs)
        return(backsolve(chol_train, tmp))
    }
    pseudoinverse(cov_train) %*% rhs
}

.mvgls_bmmcorr_conditional_mean <- function(cov_cross, cov_train, residual_vec){
    weights <- .mvgls_bmmcorr_conditional_weights(cov_train, residual_vec)
    cov_cross %*% weights
}

.mvgls_bmmcorr_prediction_cov <- function(tree, object, target_names){
    if(is.null(target_names)) stop("You must provide species names to \"newdata\"")
    train_names <- rownames(object$variables$Y)
    if(is.null(train_names)) train_names <- object$variables$tree$tip.label
    if(any(!train_names %in% tree$tip.label)){
        stop("the provided tree must contain the training species used in the model fit")
    }
    if(any(!target_names %in% tree$tip.label)){
        stop("the \"newdata\" names does not matches names in the tree ")
    }
    if(any(target_names %in% train_names)){
        stop("tree-based ", .mvgls_bmmcorr_label(object), " prediction expects target tips not already used in the training data")
    }

    keep_tips <- unique(c(train_names, target_names))
    if(length(keep_tips) != length(tree$tip.label)){
        tree <- drop.tip(tree, setdiff(tree$tip.label, keep_tips))
    }
    train_names <- train_names[train_names %in% tree$tip.label]
    target_names <- target_names[target_names %in% tree$tip.label]

    omega <- .mvgls_bmmcorr_omega_from_tree(object, tree=tree, internal=FALSE)
    ntip <- Ntip(tree)
    target_idx <- match(target_names, tree$tip.label)
    train_idx <- match(train_names, tree$tip.label)
    target_block <- .mvgls_bmmcorr_block_index(target_idx, ntip, object$dims$p)
    train_block <- .mvgls_bmmcorr_block_index(train_idx, ntip, object$dims$p)

    list(
        omega_new_train = omega[target_block, train_block, drop=FALSE],
        omega_train = omega[train_block, train_block, drop=FALSE],
        train = train_names,
        target = target_names
    )
}

.mvgls_bmmcorr_normalized_residuals <- function(object){
    Omega <- .mvgls_bmmcorr_omega_from_tree(object, tree=object$variables$tree, internal=FALSE)
    chol_omega <- try(chol(Omega), silent=TRUE)
    if(inherits(chol_omega, "try-error")){
        stop("failed to compute normalized residuals because the ", .mvgls_bmmcorr_label(object), " covariance is not positive definite")
    }
    residual_vec <- as.numeric(object$residuals)
    whitened <- backsolve(chol_omega, residual_vec, transpose=TRUE)
    residuals <- matrix(whitened,
        nrow=object$dims$n,
        ncol=object$dims$p,
        dimnames=dimnames(object$residuals)
    )
    residuals
}

.mvgls_bmmcorr_predict_phylo <- function(object, X, sp_name, tree){
    rcov <- .mvgls_bmmcorr_prediction_cov(tree, object, sp_name)
    predicted_mean <- X %*% object$coefficients
    residual_vec <- as.numeric(object$residuals[rcov$train, , drop=FALSE])
    conditional <- .mvgls_bmmcorr_conditional_mean(rcov$omega_new_train, rcov$omega_train, residual_vec)
    conditional <- matrix(conditional,
        nrow=length(rcov$target),
        ncol=object$dims$p,
        dimnames=list(rcov$target, colnames(object$variables$Y))
    )
    predicted <- predicted_mean + conditional
    rownames(predicted) <- rcov$target
    predicted
}

.mvgls_bmmcorr_ancestral <- function(object, predicted_fit){
    tree <- object$variables$tree
    omega <- .mvgls_bmmcorr_omega_from_tree(object, tree=tree, internal=TRUE)
    ntot <- Ntip(tree) + Nnode(tree)
    tip_idx <- indiceTip(tree, object$dims$p)
    node_idx <- setdiff(seq_len(ntot * object$dims$p), tip_idx)
    residual_vec <- as.numeric(object$residuals[tree$tip.label, , drop=FALSE])
    conditional <- .mvgls_bmmcorr_conditional_mean(
        omega[node_idx, tip_idx, drop=FALSE],
        omega[tip_idx, tip_idx, drop=FALSE],
        residual_vec
    )
    recons_t <- matrix(conditional,
        nrow=Nnode(tree),
        ncol=object$dims$p,
        dimnames=list(NULL, colnames(object$variables$Y))
    ) + predicted_fit
    rownames(recons_t) <- paste("node_", Ntip(tree) + seq_len(Nnode(tree)), sep="")
    recons_t
}

.mvgls_bmm_nparam <- function(object){
    if(.mvgls_is_bmmcorr(object)){
        if(is.null(object$df.free)) stop(paste0(.mvgls_bmmcorr_label(object), " BMM fits must store object$df.free"))
        return(object$df.free)
    }
    
    p <- object$dims$p
    if(object$model=="BM") (length(object$start_values)-1) + length(object$coefficients) + p*(p + 1)/2 else length(object$start_values) + length(object$coefficients) + p*(p + 1)/2
}

.mvgls_bmmcorr_loglik_guard <- function(object){
    if(.mvgls_is_bmmcorr(object) && object$method!="LL"){
        stop(.mvgls_bmmcorr_label(object, capitalize=TRUE), " BMM is only implemented for fits with method=\"LL\"")
    }
    if(.mvgls_is_bmmcorr(object) && isTRUE(object$REML)){
        stop(.mvgls_bmmcorr_label(object, capitalize=TRUE), " BMM is only available for ML fits")
    }
}

.mvgls_bmmcorr_diag_key <- function(object){
    .mvgls_bmmcorr_structure(object)
}

.mvgls_bmmcorr_loglik <- function(Y, object, coefficients=object$coefficients, sigma_object=object){
    switch(.mvgls_bmmcorr_structure(object),
        corrstrength=.mvgls_corrstrength_loglik(Y, object=object, coefficients=coefficients, sigma_object=sigma_object),
        corrpower=.mvgls_corrpower_loglik(Y, object=object, coefficients=coefficients, sigma_object=sigma_object),
        stop("Experimental BMM correlation log-likelihood is unavailable for this fit")
    )
}

.mvgls_bmmcorr_refit <- function(object, Y, start=object$opt$par, echo=FALSE){
    switch(.mvgls_bmmcorr_structure(object),
        corrstrength=.mvgls_corrstrength_refit(object, Y, start=start, echo=echo),
        corrpower=.mvgls_corrpower_refit(object, Y, start=start, echo=echo),
        stop("Experimental BMM correlation refit helper is unavailable for this fit")
    )
}

.mvgls_bmmcorr_profile <- function(object, parameter, regime, grid, optimization=NULL){
    switch(.mvgls_bmmcorr_structure(object),
        corrstrength=.mvgls_corrstrength_profile(object=object, parameter=parameter, regime=regime, grid=grid, optimization=optimization),
        corrpower=.mvgls_corrpower_profile(object=object, parameter=parameter, regime=regime, grid=grid, optimization=optimization),
        stop("Experimental BMM correlation profile helper is unavailable for this fit")
    )
}

.mvgls_bmmcorr_boundary_indicator <- function(refit, cor_est){
    cor_param <- .mvgls_bmmcorr_cor_param_name(refit)
    if(!is.finite(cor_est)) return(NA_real_)
    if(.mvgls_is_corrpower(refit)){
        return(as.numeric(cor_est < 0.05 || cor_est > 3))
    }
    corr_strength_max <- refit$diagnostics$corrstrength$corr_strength_max %||% NA_real_
    as.numeric(cor_est < 0.05 || (!is.na(corr_strength_max) && cor_est > 0.95 * corr_strength_max))
}

.mvgls_bmmcorr_diagnostics <- function(object){
    diagnostics <- NULL
    diag_key <- .mvgls_bmmcorr_diag_key(object)
    if(!is.null(object$diagnostics) && is.list(object$diagnostics)){
        diagnostics <- object$diagnostics[[diag_key]]
        if(is.null(diagnostics)) diagnostics <- object$diagnostics
    }
    if(is.null(diagnostics) || !is.list(diagnostics)) return(NULL)

    start_table <- diagnostics$start_table
    selected_row <- NULL
    if(is.data.frame(start_table) && nrow(start_table) > 0L){
        if("selected" %in% names(start_table) && any(start_table$selected %in% TRUE, na.rm = TRUE)){
            selected_row <- start_table[start_table$selected %in% TRUE, , drop=FALSE][1, , drop=FALSE]
        }else if("start_id" %in% names(start_table) && !is.null(diagnostics$selected_start_id)){
            selected_idx <- which(start_table$start_id == diagnostics$selected_start_id)
            if(length(selected_idx) > 0L) selected_row <- start_table[selected_idx[1], , drop=FALSE]
        }
    }

    selected_start_id <- diagnostics$selected_start_id
    if(is.null(selected_start_id) && !is.null(selected_row) && "start_id" %in% names(selected_row)){
        selected_start_id <- selected_row$start_id[[1]]
    }

    nstarts <- diagnostics$nstarts
    if(is.null(nstarts) && is.data.frame(start_table)) nstarts <- nrow(start_table)

    max_scale <- diagnostics$max_scale
    if(is.null(max_scale) && !is.null(selected_row) && "max_scale" %in% names(selected_row)){
        max_scale <- selected_row$max_scale[[1]]
    }
    if((is.null(max_scale) || !is.numeric(max_scale) || !is.finite(max_scale)) && !is.null(object$param)){
        scale_names <- grep("\\.scale$", names(object$param), value=TRUE)
        if(length(scale_names) > 0L){
            max_scale <- max(as.numeric(object$param[scale_names]), na.rm=TRUE)
            if(!is.finite(max_scale)) max_scale <- NULL
        }
    }

    cor_label <- .mvgls_bmmcorr_cor_param_name(object)
    boundary_field <- paste0("boundary_", cor_label)
    pathological_cor_field <- paste0("pathological_", cor_label)
    boundary_cor <- diagnostics[[boundary_field]]
    pathological_scale <- diagnostics$pathological_scale
    pathological_cor <- diagnostics[[pathological_cor_field]]
    reference_regime <- diagnostics$reference_regime
    if(is.null(reference_regime) && !is.null(object$reference_regime)){
        reference_regime <- object$reference_regime
    }

    list(
        nstarts = nstarts,
        selected_start_id = selected_start_id,
        reference_regime = reference_regime,
        max_scale = max_scale,
        boundary_corr_strength = boundary_cor,
        boundary_corr_power = diagnostics$boundary_corr_power %||% if(identical(cor_label, "corr_power")) boundary_cor else NULL,
        cor_label = cor_label,
        pathological_scale = pathological_scale,
        pathological_cor = pathological_cor,
        pathological_corr_power = diagnostics$pathological_corr_power %||% if(identical(cor_label, "corr_power")) pathological_cor else NULL,
        regularization_scale = diagnostics$regularization_scale %||% NULL,
        regularization_corr_power = diagnostics$regularization_corr_power %||% NULL,
        diagnostics_label = .mvgls_bmmcorr_label(object, capitalize=TRUE),
        start_table = start_table
    )
}

.mvgls_print_bmmcorr_diagnostics <- function(object, digits = max(3L, getOption("digits") - 3L)){
    info <- .mvgls_bmmcorr_diagnostics(object)
    if(is.null(info)) return(invisible(NULL))

    cat(info$diagnostics_label, " diagnostics:\n", sep="")
    if(!is.null(info$nstarts)) cat("nstarts:", info$nstarts, "\n")
    if(!is.null(info$selected_start_id)) cat("selected start:", info$selected_start_id, "\n")
    if(!is.null(info$reference_regime)) cat("reference regime:", info$reference_regime, "\n")
    if(!is.null(info$max_scale)) cat("max scale:", round(info$max_scale, digits = digits), "\n")
    if(!is.null(info$boundary_corr_strength)) cat("boundary", info$cor_label, "flag:", info$boundary_corr_strength, "\n")
    if(!is.null(info$pathological_scale)) cat("pathological scale flag:", info$pathological_scale, "\n")
    if(!is.null(info$pathological_cor)) cat("pathological", info$cor_label, "flag:", info$pathological_cor, "\n")
    cat("\n")
    invisible(info)
}

.mvgls_bmmcorr_regime_summary <- function(object){
    if(!.mvgls_is_bmmcorr(object)) return(NULL)
    if(is.null(object$sigma) || is.null(object$sigma$regime)) return(NULL)

    regime_names <- names(object$sigma$regime)
    if(is.null(regime_names)) regime_names <- seq_along(object$sigma$regime)
    p <- object$dims$p
    cor_param <- .mvgls_bmmcorr_cor_param_name(object)

    rows <- lapply(regime_names, function(regime_name){
        Sigma <- object$sigma$regime[[regime_name]]
        if(is.null(Sigma)) return(NULL)
        diag_vals <- diag(Sigma)
        cov_vals <- Sigma[upper.tri(Sigma)]
        cor_vals <- try({
            corr <- stats::cov2cor(Sigma)
            corr[upper.tri(corr)]
        }, silent=TRUE)
        if(inherits(cor_vals, "try-error")) cor_vals <- rep(NA_real_, length(cov_vals))

        data.frame(
            reference = identical(regime_name, object$reference_regime),
            scale = .mvgls_bmmcorr_regime_param(object, regime_name, "scale"),
            corr_strength = .mvgls_bmmcorr_regime_param(object, regime_name, cor_param),
            mean_rate = sum(diag_vals)/p,
            mean_variance = mean(diag_vals),
            mean_covariance = if(length(cov_vals)) mean(cov_vals) else 0,
            mean_correlation = if(length(cor_vals)) mean(cor_vals) else 0,
            mean_abs_correlation = if(length(cor_vals)) mean(abs(cor_vals)) else 0,
            row.names = regime_name,
            stringsAsFactors = FALSE
        )
    })

    out <- do.call(rbind, rows)
    rownames(out) <- regime_names
    names(out)[names(out) == "corr_strength"] <- cor_param
    out
}

.mvgls_bmmcorr_metric_summary <- function(values, level=0.95){
    values <- as.numeric(values)
    values <- values[is.finite(values)]
    if(!length(values)){
        return(list(low=NA_real_, high=NA_real_, mean=NA_real_, n=0L))
    }
    alpha <- (1 - level)/2
    list(
        low=unname(stats::quantile(values, probs=alpha, names=FALSE, type=7)),
        high=unname(stats::quantile(values, probs=1 - alpha, names=FALSE, type=7)),
        mean=mean(values),
        n=length(values)
    )
}

.mvgls_bmmcorr_default_profile_grid <- function(estimate, parameter, points=9L){
    points <- as.integer(points[1])
    if(is.na(points) || points < 2L) points <- 9L
    if(identical(parameter, "scale")){
        center <- max(as.numeric(estimate), 1e-4)
        lo <- max(center/4, 1e-4)
        hi <- max(center * 4, lo * 1.1)
        exp(seq(log(lo), log(hi), length.out=points))
    }else if(identical(parameter, "corr_power")){
        center <- max(as.numeric(estimate), 0.05)
        lo <- max(0.01, center/2)
        hi <- max(center * 2, lo * 1.5)
        unique(round(exp(seq(log(lo), log(hi), length.out=points)), 6))
    }else{
        center <- max(as.numeric(estimate), 0.02)
        lo <- max(0.01, center - 0.6)
        hi <- max(center + 0.6, lo * 1.5)
        unique(round(seq(lo, hi, length.out=points), 6))
    }
}

.mvgls_bmmcorr_profile_interval <- function(object, regime, parameter, level=0.95, profile_points=9L){
    cor_param <- .mvgls_bmmcorr_cor_param_name(object)
    parameter <- match.arg(parameter, c("scale", cor_param))
    estimate <- .mvgls_bmmcorr_regime_param(object, regime, parameter)
    if(identical(regime, object$reference_regime)){
        return(list(
            low=estimate,
            high=estimate,
            estimate=estimate,
            table=data.frame(
                regime=regime,
                param=parameter,
                fixed_value=estimate,
                logLik=as.numeric(object$logLik),
                convergence=0L,
                max_scale=max(object$regime.summary$scale, na.rm=TRUE),
                stringsAsFactors=FALSE
            ),
            error=NULL
        ))
    }
    grid <- .mvgls_bmmcorr_default_profile_grid(estimate, parameter, points=profile_points)
    prof <- tryCatch(
        .mvgls_bmmcorr_profile(
            object=object,
            parameter=parameter,
            regime=regime,
            grid=grid
        ),
        error=function(e) e
    )
    if(inherits(prof, "error")){
        return(list(low=NA_real_, high=NA_real_, estimate=estimate, table=NULL, error=conditionMessage(prof)))
    }
    finite <- prof[is.finite(prof$logLik), , drop=FALSE]
    if(!nrow(finite)){
        return(list(low=NA_real_, high=NA_real_, estimate=estimate, table=prof, error=NULL))
    }
    cutoff <- max(finite$logLik) - stats::qchisq(level, df=1)/2
    inside <- finite$fixed_value[finite$logLik >= cutoff]
    list(
        low=if(length(inside)) min(inside) else NA_real_,
        high=if(length(inside)) max(inside) else NA_real_,
        estimate=estimate,
        table=prof,
        error=NULL
    )
}

.mvgls_bmmcorr_refit_with_anchor <- function(object, anchor, echo=FALSE){
    if(is.null(object$model.frame)){
        stop(.mvgls_bmmcorr_label(object, capitalize=TRUE), " diagnostics require the original model.frame on the fitted object")
    }
    refit_args <- list(
        formula=object$formula,
        data=object$model.frame,
        tree=object$variables$tree,
        model="BMM",
        method="LL",
        REML=FALSE,
        response=as.matrix(object$variables$Y),
        bmm.structure=.mvgls_bmmcorr_structure(object),
        bmm.reference=anchor,
        grid.search=FALSE,
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

.mvgls_bmmcorr_bootstrap_draws <- function(object, nboot=100L, nbcores=1L){
    cor_param <- .mvgls_bmmcorr_cor_param_name(object)
    regime_names <- names(object$sigma$regime)
    if(is.null(regime_names)) regime_names <- character(0)
    draws <- .parallel_mapply(function(i){
        tryCatch({
            Yp <- simulate(object, nsim=1)
            refit <- .mvgls_bmmcorr_refit(object, Yp, start=object$opt$par)
            diag <- .mvgls_bmmcorr_diagnostics(refit)
            param_rows <- lapply(regime_names, function(regime){
                scale_est <- .mvgls_bmmcorr_regime_param(refit, regime, "scale")
                cor_est <- .mvgls_bmmcorr_regime_param(refit, regime, cor_param)
                data.frame(
                    draw=i,
                    regime=regime,
                    reference=identical(regime, refit$reference_regime),
                    scale=scale_est,
                    corr_strength=cor_est,
                    boundary=.mvgls_bmmcorr_boundary_indicator(refit, cor_est),
                    runaway_scale=as.numeric(is.finite(scale_est) && scale_est > 20),
                    stringsAsFactors=FALSE
                )
            })
            param_df <- do.call(rbind, param_rows)
            names(param_df)[names(param_df) == "corr_strength"] <- cor_param

            regime_summary <- refit$regime.summary
            if(!is.null(regime_summary)){
                regime_summary$draw <- i
                regime_summary$regime <- rownames(regime_summary)
                rownames(regime_summary) <- NULL
            }
            list(parameter_draws=param_df, regime_summary_draws=regime_summary)
        }, error=function(e){
            NULL
        })
    }, seq_len(as.integer(nboot)), mc.cores=getOption("mc.cores", nbcores))
    draws <- Filter(Negate(is.null), draws)
    if(!length(draws)){
        empty_param <- data.frame(
            draw=integer(0),
            regime=character(0),
            reference=logical(0),
            scale=numeric(0),
            stringsAsFactors=FALSE
        )
        empty_param[[cor_param]] <- numeric(0)
        empty_param$boundary <- numeric(0)
        empty_param$runaway_scale <- numeric(0)
        empty_summary <- data.frame()
        return(list(parameter_draws=empty_param, regime_summary_draws=empty_summary, n_success=0L))
    }
    param_draws <- do.call(rbind, lapply(draws, `[[`, "parameter_draws"))
    regime_summary_list <- Filter(Negate(is.null), lapply(draws, `[[`, "regime_summary_draws"))
    regime_summary_draws <- if(length(regime_summary_list)) do.call(rbind, regime_summary_list) else data.frame()
    list(parameter_draws=param_draws, regime_summary_draws=regime_summary_draws, n_success=length(draws))
}

.mvgls_bmmcorr_parameter_summary <- function(object, level=0.95, profile_points=9L, bootstrap=NULL){
    cor_param <- .mvgls_bmmcorr_cor_param_name(object)
    regime_names <- names(object$sigma$regime)
    if(is.null(regime_names)) regime_names <- character(0)
    rows <- lapply(regime_names, function(regime){
        lapply(c("scale", cor_param), function(parameter){
            estimate <- .mvgls_bmmcorr_regime_param(object, regime, parameter)
            prof <- if(profile_points > 0L) {
                .mvgls_bmmcorr_profile_interval(object, regime, parameter, level=level, profile_points=profile_points)
            }else{
                list(low=NA_real_, high=NA_real_, estimate=estimate, table=NULL, error=NULL)
            }
            boot_low <- boot_high <- boot_mean <- boundary_rate <- runaway_rate <- NA_real_
            boot_success <- 0L
            if(!is.null(bootstrap) && is.data.frame(bootstrap$parameter_draws) && nrow(bootstrap$parameter_draws) > 0L){
                sub <- bootstrap$parameter_draws[bootstrap$parameter_draws$regime == regime, , drop=FALSE]
                metric <- .mvgls_bmmcorr_metric_summary(sub[[parameter]], level=level)
                boot_low <- metric$low
                boot_high <- metric$high
                boot_mean <- metric$mean
                boot_success <- metric$n
                if(identical(parameter, cor_param)) boundary_rate <- mean(sub$boundary, na.rm=TRUE)
                if(identical(parameter, "scale")) runaway_rate <- mean(sub$runaway_scale, na.rm=TRUE)
            }
            data.frame(
                label=paste0(regime, ".", parameter),
                regime=regime,
                parameter=parameter,
                reference=identical(regime, object$reference_regime),
                estimate=estimate,
                profile_low=prof$low,
                profile_high=prof$high,
                bootstrap_low=boot_low,
                bootstrap_high=boot_high,
                bootstrap_mean=boot_mean,
                boundary_rate=boundary_rate,
                runaway_scale_rate=runaway_rate,
                boot_success=boot_success,
                stringsAsFactors=FALSE
            )
        })
    })
    do.call(rbind, unlist(rows, recursive=FALSE))
}

.mvgls_bmmcorr_regime_summary_diagnostics <- function(object, level=0.95, bootstrap=NULL){
    info <- if(!is.null(object$regime.summary)) object$regime.summary else .mvgls_bmmcorr_regime_summary(object)
    if(is.null(info)) return(NULL)
    metrics <- setdiff(colnames(info), c("reference", "scale", .mvgls_bmmcorr_cor_param_name(object)))
    rows <- lapply(rownames(info), function(regime){
        lapply(metrics, function(metric){
            estimate <- info[regime, metric]
            boot_low <- boot_high <- boot_mean <- NA_real_
            boot_success <- 0L
            if(!is.null(bootstrap) && is.data.frame(bootstrap$regime_summary_draws) && nrow(bootstrap$regime_summary_draws) > 0L){
                sub <- bootstrap$regime_summary_draws[bootstrap$regime_summary_draws$regime == regime, , drop=FALSE]
                metric_summary <- .mvgls_bmmcorr_metric_summary(sub[[metric]], level=level)
                boot_low <- metric_summary$low
                boot_high <- metric_summary$high
                boot_mean <- metric_summary$mean
                boot_success <- metric_summary$n
            }
            data.frame(
                label=paste0(regime, ".", metric),
                regime=regime,
                metric=metric,
                reference=info[regime, "reference"],
                estimate=estimate,
                bootstrap_low=boot_low,
                bootstrap_high=boot_high,
                bootstrap_mean=boot_mean,
                boot_success=boot_success,
                stringsAsFactors=FALSE
            )
        })
    })
    do.call(rbind, unlist(rows, recursive=FALSE))
}

.mvgls_bmmcorr_anchor_sensitivity <- function(object, anchors=NULL){
    if(identical(anchors, FALSE)) return(NULL)
    regime_names <- names(object$sigma$regime)
    if(is.null(regime_names)) regime_names <- character(0)
    if(is.null(anchors)) anchors <- regime_names
    anchors <- intersect(as.character(anchors), regime_names)
    if(!length(anchors)) return(NULL)

    fits <- lapply(anchors, function(anchor){
        tryCatch(.mvgls_bmmcorr_refit_with_anchor(object, anchor=anchor, echo=FALSE), error=function(e) e)
    })
    names(fits) <- anchors

    anchor_rows <- lapply(names(fits), function(anchor){
        fit <- fits[[anchor]]
        if(inherits(fit, "error")){
            return(data.frame(
                anchor=anchor,
                success=FALSE,
                logLik=NA_real_,
                AIC=NA_real_,
                GIC=NA_real_,
                selected_start=NA_integer_,
                max_scale=NA_real_,
                boundary_corr_strength=NA,
                pathological_scale=NA,
                error=conditionMessage(fit),
                stringsAsFactors=FALSE
            ))
        }
        diag <- fit$diagnostics[[.mvgls_bmmcorr_diag_key(fit)]]
        data.frame(
            anchor=fit$reference_regime,
            success=TRUE,
            logLik=as.numeric(fit$logLik),
            AIC=AIC(fit)$AIC,
            GIC=GIC(fit)$GIC,
            selected_start=diag$selected_start_id %||% NA_integer_,
            max_scale=diag$max_scale %||% NA_real_,
            boundary_corr_strength=diag$boundary_corr_strength %||% NA,
            pathological_scale=diag$pathological_scale %||% NA,
            error=NA_character_,
            stringsAsFactors=FALSE
        )
    })
    anchor_summary <- do.call(rbind, anchor_rows)
    rownames(anchor_summary) <- NULL

    ok_names <- names(fits)[vapply(fits, inherits, logical(1), "mvgls")]
    pairwise_cov <- NULL
    if(length(ok_names) >= 2L){
        pair_rows <- list()
        for(pair in utils::combn(ok_names, 2L, simplify=FALSE)){
            left <- fits[[pair[1]]]
            right <- fits[[pair[2]]]
            common_regimes <- intersect(names(left$sigma$regime), names(right$sigma$regime))
            for(regime in common_regimes){
                delta <- left$sigma$regime[[regime]] - right$sigma$regime[[regime]]
                pair_rows[[length(pair_rows) + 1L]] <- data.frame(
                    anchor_left=pair[1],
                    anchor_right=pair[2],
                    regime=regime,
                    max_abs_diff=max(abs(delta)),
                    frobenius_diff=sqrt(sum(delta^2)),
                    stringsAsFactors=FALSE
                )
            }
        }
        pairwise_cov <- do.call(rbind, pair_rows)
        rownames(pairwise_cov) <- NULL
    }

    success_summary <- anchor_summary[anchor_summary$success, , drop=FALSE]
    loglik_spread <- if(nrow(success_summary) > 1L) diff(range(success_summary$logLik)) else 0
    acceptance <- list(
        logLik_spread=loglik_spread,
        logLik_spread_ok=is.finite(loglik_spread) && loglik_spread < 5,
        no_pathological_scale=!any(as.logical(success_summary$pathological_scale), na.rm=TRUE)
    )

    list(anchor_summary=anchor_summary, pairwise_cov=pairwise_cov, acceptance=acceptance, fits=fits)
}

.mvgls_bmmcorr_format_value <- function(x, digits=3){
    if(is.logical(x)) return(ifelse(is.na(x), "NA", ifelse(x, "TRUE", "FALSE")))
    if(is.character(x)) return(ifelse(is.na(x), "NA", x))
    x <- as.numeric(x)
    ifelse(is.na(x), "NA", formatC(x, digits=digits, format="f"))
}

corrstrength_diagnostics.mvgls <- function(object, level=0.95, nboot=100L, nbcores=1L, profile_points=9L, anchors=NULL, include_bootstrap_draws=FALSE, include_anchor_fits=FALSE, ...){
    if(!.mvgls_is_corrstrength(object)){
        stop("corrstrength_diagnostics() is only implemented for mvgls fits with bmm.structure=\"corrstrength\"")
    }
    .mvgls_bmmcorr_loglik_guard(object)
    nboot <- as.integer(nboot[1])
    if(is.na(nboot) || nboot < 0L) stop("'nboot' must be a non-negative integer")
    profile_points <- as.integer(profile_points[1])
    if(is.na(profile_points) || profile_points < 0L) stop("'profile_points' must be a non-negative integer")

    bootstrap <- if(nboot > 0L) .mvgls_bmmcorr_bootstrap_draws(object, nboot=nboot, nbcores=nbcores) else NULL
    parameter_summary <- .mvgls_bmmcorr_parameter_summary(object, level=level, profile_points=profile_points, bootstrap=bootstrap)
    regime_summary <- .mvgls_bmmcorr_regime_summary_diagnostics(object, level=level, bootstrap=bootstrap)
    anchor_sensitivity <- .mvgls_bmmcorr_anchor_sensitivity(object, anchors=anchors)

    results <- list(
        parameter_summary=parameter_summary,
        regime_summary=regime_summary,
        anchor_summary=if(!is.null(anchor_sensitivity)) anchor_sensitivity$anchor_summary else NULL,
        anchor_pairwise_cov=if(!is.null(anchor_sensitivity)) anchor_sensitivity$pairwise_cov else NULL,
        acceptance=if(!is.null(anchor_sensitivity)) anchor_sensitivity$acceptance else NULL,
        settings=list(level=level, nboot=nboot, nbcores=nbcores, profile_points=profile_points, anchors=anchors)
    )
    if(include_bootstrap_draws && !is.null(bootstrap)){
        results$bootstrap_draws <- bootstrap
    }
    if(include_anchor_fits && !is.null(anchor_sensitivity)){
        results$anchor_fits <- anchor_sensitivity$fits
    }
    class(results) <- "corrstrength_diagnostics"
    results
}

corrpower_diagnostics.mvgls <- function(object, level=0.95, nboot=100L, nbcores=1L, profile_points=9L, anchors=NULL, include_bootstrap_draws=FALSE, include_anchor_fits=FALSE, ...){
    if(!.mvgls_is_corrpower(object)){
        stop("corrpower_diagnostics() is only implemented for mvgls fits with bmm.structure=\"corrpower\"")
    }
    .mvgls_bmmcorr_loglik_guard(object)
    nboot <- as.integer(nboot[1])
    if(is.na(nboot) || nboot < 0L) stop("'nboot' must be a non-negative integer")
    profile_points <- as.integer(profile_points[1])
    if(is.na(profile_points) || profile_points < 0L) stop("'profile_points' must be a non-negative integer")

    bootstrap <- if(nboot > 0L) .mvgls_bmmcorr_bootstrap_draws(object, nboot=nboot, nbcores=nbcores) else NULL
    parameter_summary <- .mvgls_bmmcorr_parameter_summary(object, level=level, profile_points=profile_points, bootstrap=bootstrap)
    regime_summary <- .mvgls_bmmcorr_regime_summary_diagnostics(object, level=level, bootstrap=bootstrap)
    anchor_sensitivity <- .mvgls_bmmcorr_anchor_sensitivity(object, anchors=anchors)
    diag_info <- .mvgls_bmmcorr_diagnostics(object)
    acceptance <- if(!is.null(anchor_sensitivity)) anchor_sensitivity$acceptance else NULL
    if(is.null(acceptance)) acceptance <- list()
    acceptance$selected <- diag_info$selected_start_id %||% NA_integer_
    acceptance$boundary_corr_power <- diag_info$boundary_corr_power %||% NA
    acceptance$pathological_scale <- diag_info$pathological_scale %||% NA
    acceptance$regularization_scale <- diag_info$regularization_scale %||% NA
    acceptance$regularization_corr_power <- diag_info$regularization_corr_power %||% NA

    results <- list(
        parameter_summary=parameter_summary,
        regime_summary=regime_summary,
        anchor_summary=if(!is.null(anchor_sensitivity)) anchor_sensitivity$anchor_summary else NULL,
        anchor_pairwise_cov=if(!is.null(anchor_sensitivity)) anchor_sensitivity$pairwise_cov else NULL,
        acceptance=acceptance,
        settings=list(level=level, nboot=nboot, nbcores=nbcores, profile_points=profile_points, anchors=anchors)
    )
    if(include_bootstrap_draws && !is.null(bootstrap)){
        results$bootstrap_draws <- bootstrap
    }
    if(include_anchor_fits && !is.null(anchor_sensitivity)){
        results$anchor_fits <- anchor_sensitivity$fits
    }
    class(results) <- "corrpower_diagnostics"
    results
}

print.corrstrength_diagnostics <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
    cat("Corr-strength diagnostics\n")
    if(!is.null(x$settings)){
        cat("level:", round(x$settings$level, digits = digits), "\n")
        cat("bootstrap replicates:", x$settings$nboot, "\n")
        cat("profile points:", x$settings$profile_points, "\n\n")
    }
    if(is.data.frame(x$parameter_summary)){
        display <- x$parameter_summary
        numeric_cols <- vapply(display, is.numeric, logical(1))
        display[numeric_cols] <- lapply(display[numeric_cols], round, digits = digits)
        cat("Parameter summary:\n")
        print(display, row.names=FALSE)
        cat("\n")
    }
    if(is.data.frame(x$regime_summary)){
        display <- x$regime_summary
        numeric_cols <- vapply(display, is.numeric, logical(1))
        display[numeric_cols] <- lapply(display[numeric_cols], round, digits = digits)
        cat("Derived regime summary:\n")
        print(display, row.names=FALSE)
        cat("\n")
    }
    if(is.data.frame(x$anchor_summary)){
        display <- x$anchor_summary
        numeric_cols <- vapply(display, is.numeric, logical(1))
        display[numeric_cols] <- lapply(display[numeric_cols], round, digits = digits)
        cat("Anchor sensitivity:\n")
        print(display, row.names=FALSE)
        cat("\n")
    }
    if(!is.null(x$acceptance)){
        cat("Acceptance checks:\n")
        cat("logLik spread < 5:", x$acceptance$logLik_spread_ok,
            "(spread =", .mvgls_bmmcorr_format_value(x$acceptance$logLik_spread, digits = digits), ")\n")
        cat("no pathological scale flags:", x$acceptance$no_pathological_scale, "\n")
    }
    invisible(x)
}

print.corrpower_diagnostics <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
    cat("Corr-power diagnostics\n")
    if(!is.null(x$settings)){
        cat("level:", round(x$settings$level, digits = digits), "\n")
        cat("bootstrap replicates:", x$settings$nboot, "\n")
        cat("profile points:", x$settings$profile_points, "\n\n")
    }
    if(is.data.frame(x$parameter_summary)){
        display <- x$parameter_summary
        numeric_cols <- vapply(display, is.numeric, logical(1))
        display[numeric_cols] <- lapply(display[numeric_cols], round, digits = digits)
        cat("parameter_summary:\n")
        print(display, row.names=FALSE)
        cat("\n")
    }
    if(is.data.frame(x$regime_summary)){
        display <- x$regime_summary
        numeric_cols <- vapply(display, is.numeric, logical(1))
        display[numeric_cols] <- lapply(display[numeric_cols], round, digits = digits)
        cat("regime_summary:\n")
        print(display, row.names=FALSE)
        cat("\n")
    }
    if(is.data.frame(x$anchor_summary)){
        display <- x$anchor_summary
        numeric_cols <- vapply(display, is.numeric, logical(1))
        display[numeric_cols] <- lapply(display[numeric_cols], round, digits = digits)
        cat("anchor_summary:\n")
        print(display, row.names=FALSE)
        cat("\n")
    }
    if(!is.null(x$acceptance)){
        cat("Acceptance checks:\n")
        cat("logLik spread < 5:", x$acceptance$logLik_spread_ok,
            "(spread =", .mvgls_bmmcorr_format_value(x$acceptance$logLik_spread, digits = digits), ")\n")
        cat("no pathological scale flags:", x$acceptance$no_pathological_scale, "\n")
    }
    invisible(x)
}

.mvgls_print_bmmcorr_regime_summary <- function(object, digits = max(3L, getOption("digits") - 3L)){
    info <- if(!is.null(object$regime.summary)) object$regime.summary else .mvgls_bmmcorr_regime_summary(object)
    if(is.null(info)) return(invisible(NULL))

    display <- info
    numeric_cols <- vapply(display, is.numeric, logical(1))
    display[numeric_cols] <- lapply(display[numeric_cols], round, digits = digits)

    cat("Regime summary:\n")
    print(display)
    cat("\n")
    invisible(info)
}

# Legacy internal aliases retained for archived corr-strength benchmarks.
.mvgls_corrstrength_regime_param <- .mvgls_bmmcorr_regime_param
.mvgls_corrstrength_primary_param_names <- .mvgls_bmmcorr_primary_param_names
.mvgls_corrstrength_parse_parm <- .mvgls_bmmcorr_parse_parm
.mvgls_corrstrength_regime_covariances <- .mvgls_bmmcorr_regime_covariances
.mvgls_corrstrength_omega_from_tree <- .mvgls_bmmcorr_omega_from_tree
.mvgls_corrstrength_prediction_cov <- .mvgls_bmmcorr_prediction_cov
.mvgls_corrstrength_normalized_residuals <- .mvgls_bmmcorr_normalized_residuals
.mvgls_corrstrength_predict_phylo <- .mvgls_bmmcorr_predict_phylo
.mvgls_corrstrength_ancestral <- .mvgls_bmmcorr_ancestral
.mvgls_corrstrength_loglik_guard <- .mvgls_bmmcorr_loglik_guard
.mvgls_corrstrength_diagnostics <- .mvgls_bmmcorr_diagnostics
.mvgls_print_corrstrength_diagnostics <- .mvgls_print_bmmcorr_diagnostics
.mvgls_corrstrength_regime_summary <- .mvgls_bmmcorr_regime_summary
.mvgls_corrstrength_profile_interval <- .mvgls_bmmcorr_profile_interval
.mvgls_corrstrength_bootstrap_draws <- .mvgls_bmmcorr_bootstrap_draws
.mvgls_corrstrength_parameter_summary <- .mvgls_bmmcorr_parameter_summary
.mvgls_corrstrength_regime_summary_diagnostics <- .mvgls_bmmcorr_regime_summary_diagnostics
.mvgls_corrstrength_anchor_sensitivity <- .mvgls_bmmcorr_anchor_sensitivity
.mvgls_corrstrength_format_value <- .mvgls_bmmcorr_format_value
.mvgls_print_corrstrength_regime_summary <- .mvgls_print_bmmcorr_regime_summary

# ------------------------------------------------------------------------- #
# BIC.mvgls                                                                 #
# options: object, ...                                                      #
# S3 method - Bayesian Information Criterion - generic from stats           #
# ------------------------------------------------------------------------- #

BIC.mvgls <- function(object, ...){
    
    # retrieve arguments
    args <- list(...)
    if(is.null(args[["REML"]])) args$forceREML <- TRUE else args$forceREML <- args$REML
    .mvgls_bmmcorr_loglik_guard(object)
    if(object$REML & args$forceREML==FALSE) LL <- .reml_to_ml(object) else LL <- object$logLik

    if(.mvgls_is_bmmcorr(object)){
        n <- object$dims$n
        nparam <- .mvgls_bmm_nparam(object)
        BIC <- -2*LL + log(n)*nparam

        results <- list(LogLikelihood=LL, BIC=BIC, nparam=nparam, k=log(n))
        class(results) <- c("bic.mvgls","bic")
        return(results)
    }
    
    # TODO generalize to mvXX functions
    if(object$method=="LL"){
        p <- object$dims$p
        n <- object$dims$n # should take n or n*p?

        nparam = if(object$model=="BM") (length(object$start_values)-1) + length(object$coefficients) + p*(p + 1)/2 else length(object$start_values) + length(object$coefficients) + p*(p + 1)/2
        # BIC
        BIC = -2*LL+ log(n)*nparam
    }else if(object$method=="EmpBayes"){
        p <- object$dims$p
        n <- object$dims$n # should take n or n*p?

        nparam = if(object$model=="BM") (length(object$start_values)-1) + length(object$coefficients) else length(object$start_values) + length(object$coefficients)
        # BIC
        BIC = -2*LL+ log(n)*nparam
    }else{
        stop("BIC works only for models fit by Maximum Likelihood (method=\"LL\")")
    }
    
    # return the results
    results <- list(LogLikelihood=LL, BIC=BIC, nparam=nparam, k=log(n))
    class(results) <- c("bic.mvgls","bic")
    return(results)
}

# ------------------------------------------------------------------------- #
# AIC.mvgls                                                                 #
# options: object, ..., k = 2                                               #
# S3 method - Akaike Information Criterion - generic from stats             #
# ------------------------------------------------------------------------- #

AIC.mvgls <- function(object, ..., k = 2){
    
    # retrieve arguments
    args <- list(...)
    if(is.null(args[["REML"]])) args$forceREML <- TRUE else args$forceREML <- args$REML
    .mvgls_bmmcorr_loglik_guard(object)
    if(object$REML & args$forceREML==FALSE) LL <- .reml_to_ml(object) else LL <- object$logLik

    if(.mvgls_is_bmmcorr(object)){
        nparam <- .mvgls_bmm_nparam(object)
        AIC <- -2*LL + k*nparam

        results <- list(LogLikelihood=LL, AIC=AIC, nparam=nparam, k=k)
        class(results) <- c("aic.mvgls","aic")
        return(results)
    }
  
    if(object$method=="LL"){
        p <- object$dims$p
        nparam = if(object$model=="BM") (length(object$start_values)-1) + length(object$coefficients) + p*(p + 1)/2 else length(object$start_values) + length(object$coefficients) + p*(p + 1)/2
        # AIC
        AIC = -2*LL+k*nparam
    }else if(object$method=="EmpBayes"){
        nparam = if(object$model=="BM") (length(object$start_values)-1) + length(object$coefficients)  else length(object$start_values) + length(object$coefficients)
        # AIC
        AIC = -2*LL+k*nparam
    }else{
        stop("AIC works only for models fit by Maximum Likelihood (method=\"LL\")")
    }
    
    # return the results
    results <- list(LogLikelihood=LL, AIC=AIC, nparam=nparam, k=k)
    class(results) <- c("aic.mvgls","aic")
    return(results)
}

# ------------------------------------------------------------------------- #
# GIC.mvgls                                                                 #
# options: model,...                                                        #
# S3 method from "RPANDA"  package                                          #
# ------------------------------------------------------------------------- #
GIC.mvgls <- function(object, ...){
    
    # retrieve arguments
    args <- list(...)
    if(is.null(args[["eigSqm"]])) eigSqm <- TRUE else eigSqm <- args$eigSqm
    if(is.null(args[["REML"]])) args$forceREML <- TRUE else args$forceREML <- args$REML # force to true to follow what has been done in the first paper?

    if(.mvgls_is_bmmcorr(object)){
        .mvgls_bmmcorr_loglik_guard(object)
        if(!is.null(args[["eigSqm"]]) && !isTRUE(args$eigSqm)) warning("The 'eigSqm' argument is ignored for ", .mvgls_bmmcorr_label(object), " BMM fits")
        
        LL <- object$logLik
        n <- object$dims$n
        p <- object$dims$p
        bias_cov <- if(is.null(object$df.free_cov)) NA_real_ else object$df.free_cov
        beta_df <- if(is.null(object$df.free_beta)) 0 else object$df.free_beta
        mod.par <- if(is.null(object$df.free_model)) 0 else object$df.free_model
        bias <- if(is.null(object$df.free)) bias_cov + beta_df + mod.par else object$df.free
        GIC <- -2*LL + 2*bias
        
        results <- list(LogLikelihood=LL, GIC=GIC, p=p, n=n, bias=bias, bias_cov=bias_cov, args=args)
        class(results) <- c("gic.mvgls","gic")
        return(results)
    }
     
    method <- object$method
    penalty <- object$penalty
    target <- object$target
    n <- object$dims$n
    p <- object$dims$p
    m <- object$dims$m
    tuning <- object$tuning
    P <- object$sigma$P # The precision matrix
    Pi <- object$sigma$Pinv # the covariance matrix
    S <- object$sigma$S # the sample estimate
    if(penalty!="EmpBayes") Target <- .targetM(S=S, targM=target, penalty=penalty)
    beta <- object$coefficients
    
    if(eigSqm){ # to follow the scheme in RPANDA
        sqM1 <- .sqM1(object$corrSt$phy)
        if(!is.null(object$corrSt$diagWeight)){
            w <- 1/object$corrSt$diagWeight
            Y <- crossprod(sqM1, matrix(w*object$variables$Y, nrow=n))
            X <- crossprod(sqM1, matrix(w*object$variables$X, nrow=n))
        }else{
            X <- crossprod(sqM1, object$variables$X)
            Y <- crossprod(sqM1, object$variables$Y)
        }
        residuals <- Y - X%*%beta
    }else{
        residuals <- residuals(object, type="normalized")
        X <- object$corrSt$X
        Y <- object$corrSt$Y
    }
    
    if(object$model=="BM"){
       mod.par=0
    }else if(object$model=="BMM"){
       mod.par=(ncol(object$corrSt$phy$mapped.edge)-1) # should we consider k parameters or k-1 (i.e. relative scaling to the first group)
    }else{
       mod.par=1
    }
    if(is.numeric(object$mserr)) mod.par = mod.par + 1 # already included in the covariance matrix structure?
    if(object$REML & args$forceREML==TRUE) ndimCov = n - m else ndimCov = n
    # Nominal loocv
    XtX <- pseudoinverse(crossprod(X))
    # hat matrix
    h <- diag(X%*%pseudoinverse(X))
    
    # check for hat score of 1 (e.g. MANOVA design)
    nloo <- 1:n
    nloo <- nloo[!h+1e-8>=1]
    nC = length(nloo)
    
    if(penalty=="RidgeArch"){
        
        # First and second derivative of the functional (we can use patterned matrix to target some matrix elements)
        # We use the Kronecker-vec identity to speed up the computations
        T1 <- sapply(nloo, function(i){
            Sk <- tcrossprod(residuals[i,]) ;
            VSV <- 0.5*(Pi - (1-tuning)*Sk - tuning*Target);
            VSV2 <- 0.5*(Pi - Sk);
            sum(VSV * 2*(P%*%VSV2%*%P))
        })
        
        df = sum(T1)/nC
        sigma_df <- df
        
    }else if(penalty=="LASSO" | penalty=="LL"){
        
        # LASSO or ML
        Tf2 <- function(S, P) {
            I <- ifelse(P==0,0,1) ;
            t(.vec(S*I))%*%.vec(P%*%(S*I)%*%P)
        }
        
        sigma_df <- (1/(2*nC))*sum(sapply(nloo, function(i){ Tf2(tcrossprod(residuals[i,]) , P)})) - (1/2)*Tf2(S,P)
        
    }else if(penalty=="RidgeAlt"){
        # Alternative Ridge
        eig <- eigen(Pi)
        V <- eig$vectors
        d <- eig$values
        H <- (1/(0.5*(kronecker(d,d)+tuning)))
        
        # 2) First derivative of the functional
        T1 <- sapply(nloo, function(i){
            Sk <- tcrossprod(residuals[i,]) ;
            VSV <- .vec(crossprod(V, (0.5*(Pi - (Sk - tuning*Target) - tuning*P))%*%V));
            VSV2 <- .vec(crossprod(V, (0.5*(Pi - Sk))%*%V));
            sum(VSV * (H*VSV2))
        })
        
        df = sum(T1)/nC
        sigma_df <- df
    }else{
        sigma_df <- 0
    }
    
    # Number of parameters for the root state:
    # The Information matrix from the Hessian and gradients scores
    T2 <- sapply(nloo, function(i){
        gradient <- (X[i,])%*%t(P%*%t(Y[i,]-X[i,]%*%beta))
        sum(gradient * (XtX%*%gradient%*%Pi))
    })
    beta_df <- sum(T2)
    
    if( min(m, sum(object$dims$assign!=0))>1 & args$forceREML==FALSE) warning("GIC criterion with multiple predictors has not been fully tested. Please use it with care and consider EIC or simulations instead")
    
    # LogLikelihood (minus)
    if(penalty=="EmpBayes"){
         if(object$REML & args$forceREML==FALSE) llik <- -.reml_to_ml(object) else llik <- -object$logLik
         
         # compute the bias term with effective df (+1 for the hyper-parameter/regularization term)
         bias <- (beta_df+mod.par+1)
         
    }else{
        DP <- as.numeric(determinant(Pi)$modulus)
        if(object$REML==TRUE & args$forceREML==FALSE) Ccov <- as.numeric(object$corrSt$det -     determinant(crossprod(object$corrSt$X))$modulus + object$corrSt$const) else Ccov <- as.numeric(object$corrSt$det)
        if(object$REML==TRUE & args$forceREML==FALSE) S <- (S*(n-m))/ndimCov # want quadratic product with REML estimate for P
        llik <- 0.5 * (ndimCov*p*log(2*pi) + p*Ccov + ndimCov*DP + ndimCov*sum(S*P))
        
        # compute the bias term with effective df
        bias <- (sigma_df+beta_df+mod.par)
   }
    
    GIC <- 2*llik + 2*bias
    
    # return the results
    results <- list(LogLikelihood=-llik, GIC=GIC, p=p, n=n, bias=bias, bias_cov=sigma_df, args=args)
    class(results) <- c("gic.mvgls","gic")
    return(results)
}

# ------------------------------------------------------------------------- #
# EIC.mvgls                                                                 #
# options: object, nboot, nbcores, ...                                      #
# S3 method - Extended/Efron Information Criterion                          #
# ------------------------------------------------------------------------- #

EIC.mvgls <- function(object, nboot=100L, nbcores=1L, ...){
    
    # retrieve arguments
    args <- list(...)
    if(.mvgls_is_bmmcorr(object)){
        .mvgls_bmmcorr_loglik_guard(object)
        if(!is.null(args[["eigSqm"]])) warning("The 'eigSqm' argument is ignored for ", .mvgls_bmmcorr_label(object), " BMM fits")
        if(!is.null(args[["REML"]])) warning("The 'REML' argument is ignored for ", .mvgls_bmmcorr_label(object), " BMM fits")
        if(is.null(args[["restricted"]])) restricted <- FALSE else restricted <- args$restricted
        
        llik <- .mvgls_bmmcorr_loglik(object$variables$Y, object=object)
        if(!is.finite(llik)) stop("The ", .mvgls_bmmcorr_label(object), " BMM likelihood could not be evaluated for the empirical fit")
        bias <- .parallel_mapply(function(i){
            tryCatch({
                Yp <- simulate(object, nsim=1)
                objectBoot <- .mvgls_bmmcorr_refit(object, Yp, start=object$opt$par)
                
                ll1 <- .mvgls_bmmcorr_loglik(objectBoot$variables$Y, object=objectBoot)
                if(!is.finite(ll1)) stop("The bootstrap experimental BMM fit produced a non-finite log-likelihood")
                
                if(restricted){
                    ll2 <- .mvgls_bmmcorr_loglik(Yp, object=object, coefficients=objectBoot$coefficients)
                    ll4 <- .mvgls_bmmcorr_loglik(object$variables$Y, object=object, coefficients=object$coefficients,
                        sigma_object=objectBoot)
                }else{
                    ll2 <- .mvgls_bmmcorr_loglik(Yp, object=object, coefficients=object$coefficients)
                    ll4 <- .mvgls_bmmcorr_loglik(object$variables$Y, object=object, coefficients=objectBoot$coefficients,
                        sigma_object=objectBoot)
                }
                if(!is.finite(ll2) || !is.finite(ll4)) stop("Cross-evaluated experimental BMM log-likelihood was non-finite")
                
                (ll1 - ll2) + (llik - ll4)
            }, error=function(e){
                NA_real_
            })
        }, 1:nboot, mc.cores = getOption("mc.cores", nbcores))
        bias <- as.numeric(bias)
        if(sum(is.na(bias)) > 0L) warning("There were multiple issues/aborted estimations in the bootstrapped samples")
        bias <- bias[!is.na(bias)]
        nboot_eff <- length(bias)
        if(nboot_eff == 0L) stop("All experimental BMM EIC bootstrap replicates failed")
        
        pboot <- mean(bias)
        EIC <- -2*llik + 2*pboot
        se <- if(nboot_eff > 1L) sd(bias)/sqrt(nboot_eff) else 0
        
        results <- list(EIC=EIC, bias=bias, LogLikelihood=llik, se=se, p=object$dims$p, n=object$dims$n)
        class(results) <- c("eic.mvgls","eic")
        return(results)
    }
    if(is.null(args[["eigSqm"]])) eigSqm <- TRUE else eigSqm <- args$eigSqm
    if(is.null(args[["restricted"]])) restricted <- FALSE else restricted <- args$restricted
    if(is.null(args[["REML"]])) args$forceREML <- FALSE else args$forceREML <- args$REML
    if(is.null(args[["method"]])) args$method <- object$method
    
    # retrieve data to simulate bootstrap samples
    beta <- object$coefficients
    if(eigSqm){ # to follow the scheme in RPANDA
      sqM1 <- .sqM1(object$corrSt$phy)
      if(!is.null(object$corrSt$diagWeight)){
        w <- 1/object$corrSt$diagWeight
        Y <- crossprod(sqM1, matrix(w*object$variables$Y, nrow=object$dims$n))
        X <- crossprod(sqM1, matrix(w*object$variables$X, nrow=object$dims$n))
      }else{
        X <- crossprod(sqM1, object$variables$X)
        Y <- crossprod(sqM1, object$variables$Y)
      }
      residuals <- Y - X%*%beta
    }else{
      residuals <- residuals(object, type="normalized")
      X <- object$corrSt$X
      Y <- object$corrSt$Y
    }
    
    N = nrow(Y)
    p = object$dims$p
    if(object$REML) m = object$dims$m else m = 0
    if(object$REML & args$forceREML==TRUE) ndimCov = object$dims$n - m else ndimCov = object$dims$n
    tuning <- object$tuning
    target <- object$target
    penalty <- object$penalty
    
    # Mute unecessary options from EmpBayes
    if(object$method=="EmpBayes"){
        object$MMSE <- quote(FALSE)
        object$FCI <- quote(FALSE)
        }
    
    # Weight matrix (OU, etc)
    if(is.null(object$corrSt$diagWeight)){
      diagWeight <- 1; is_weight = FALSE
    }else{
      diagWeight <- object$corrSt$diagWeight; is_weight = TRUE
      diagWeightInv <- 1/diagWeight
    }
    Dsqrt <- .pruning_general(object$corrSt$phy, trans=FALSE, inv=FALSE)$sqrtM # return warning message if n-ultrametric tree is used with OU?
    # TODO (change to allow n-ultrametric and OU) > just need to standardize the data by the weights
    # if(object$model=="OU" & !is.ultrametric(object$variables$tree)) stop("The EIC method does not handle yet non-ultrametric trees with OU processes")
    
    DsqrtInv <- .pruning_general(object$corrSt$phy, trans=FALSE, inv=TRUE)$sqrtM
    modelPerm <- object$call
    modelPerm$formula <- object$formula
    modelPerm$tree <- quote(object$variables$tree)
    if(!is.null(object$model.frame)) modelPerm$data <- quote(object$model.frame)
    modelPerm$grid.search <- quote(FALSE)
    modelPerm$start <- quote(object$opt$par)
    
    # Mean and residuals for the model
    MeanNull <- object$variables$X%*%beta
    
    
    # Estimate the bias term
    D1 <- function(objectBoot, objectFit, ndimCov, m, p, sqM, Ccov2){ # LL(Y*|param*) - LL(Y*| param)
      
      # Y*|param*
      residualsBoot <- residuals(objectBoot, type="normalized")
      
      # For boot "i" LL1(Y*|param*)
      if(objectFit$REML==TRUE & args$forceREML==FALSE) Ccov1 <- as.numeric(objectBoot$corrSt$det - determinant(crossprod(objectBoot$corrSt$X))$modulus + objectBoot$corrSt$const) else Ccov1 <- as.numeric(objectBoot$corrSt$det)
      
      # compute LL1(Y*|param*)
      llik1 <- .llik_fn(object_boot=objectBoot, object_emp=objectFit, residualsBoot=residualsBoot,
                       method=args$method, Ccov=Ccov1, n=ndimCov, m=m,
                       p=p, v=p+1, bias_type="D1_1")
      # Y*|param
      #if(!restricted) residualsBoot <- objectBoot$corrSt$Y - objectBoot$corrSt$X%*%objectFit$coefficients # does not account for the phylo model of the original fit
      if(!restricted){
        if(is_weight){
          residualsBoot <- crossprod(sqM, (objectBoot$variables$Y - objectBoot$variables$X%*%objectFit$coefficients)*diagWeightInv)
        }else{
          residualsBoot <- crossprod(sqM, objectBoot$variables$Y - objectBoot$variables$X%*%objectFit$coefficients)
          
        }
      }
      
      # For boot "i" LL2(Y*|param)
      # if(objectFit$REML==TRUE & args$forceREML==FALSE) Ccov2 <- as.numeric(objectFit$corrSt$det - determinant(crossprod(objectFit$corrSt$X))$modulus + objectFit$corrSt$const) else Ccov2 <- as.numeric(objectFit$corrSt$det)
     
      llik2 <- .llik_fn(object_boot=objectBoot, object_emp=objectFit, residualsBoot=residualsBoot,
                       method=args$method, Ccov=Ccov2, n=ndimCov, m=m,
                       p=p, v=p+1, bias_type="D1_2")
      
      # Return the difference in LL for D1
      return(llik1 - llik2)
    }
    
    D3 <- function(objectBoot, objectFit, loglik, ndimCov, m, p, method){ # LL(Y|param) - LL(Y| param*)
      
      # Y|param*
      if(!restricted) {
        sqM_temp <- .pruning_general(objectBoot$corrSt$phy, trans=FALSE, inv=TRUE)$sqrtM
        if(is_weight){
          residualsBoot <- try(crossprod(sqM_temp, (objectFit$variables$Y - objectFit$variables$X%*%objectBoot$coefficients)/objectBoot$corrSt$diagWeight), silent=TRUE)
        } else {
          residualsBoot <- try(crossprod(sqM_temp, objectFit$variables$Y - objectFit$variables$X%*%objectBoot$coefficients), silent=TRUE)
          
        }
      }else{ residualsBoot <- objectFit$corrSt$Y - objectFit$corrSt$X%*%objectFit$coefficients}
      
        
      # For boot "i" LL2(Y|param*)
      if(objectFit$REML==TRUE & args$forceREML==FALSE) Ccov1 <- as.numeric(objectBoot$corrSt$det - determinant(crossprod(objectBoot$corrSt$X))$modulus + objectBoot$corrSt$const) else Ccov1 <- as.numeric(objectBoot$corrSt$det)
      
      
      # compute LL2(Y|param*)
      llik2 <- .llik_fn(object_boot=objectBoot, object_emp=objectFit, residualsBoot=residualsBoot,
                       method=args$method, Ccov=Ccov1, n=ndimCov, m=m,
                       p=p, v=p+1, bias_type="D1_1")
      
      # Return the difference in LL for D3
      return(loglik - llik2)
    }
    
    # Estimate EIC: LL+bias
    
    # Maximum Likelihood
    if(object$REML==TRUE & args$forceREML==FALSE) Ccov <- as.numeric(object$corrSt$det - determinant(crossprod(object$corrSt$X))$modulus + object$corrSt$const) else Ccov <- as.numeric(object$corrSt$det)
    
    if(args$method=="EmpBayes"){
        
        # df for the Matrix T distribution
        v = p+1
        
        if(object$target=="Variance"){
          target = colSums(residuals^2)*(1/(ndimCov-m))*object$tuning
          SigS2 <- svd(t(residuals)*sqrt(1/(target*(v-p))), nu=0, nv=0)$d^2
          detSig <- sum(log(target*(v-p)))
        }else{
          tuning = mean(colSums(residuals^2)*(1/(ndimCov-m)))*object$tuning
          SigS2 <- svd(residuals*sqrt(1/((v-p)*tuning)), nu=0, nv=0)$d^2
          detSig <- p*log((v-p)*tuning) # note, with default df, v-p=1
        }
        
        Kdet <- 0.5*(v+ndimCov+p-1)*sum(log(1+SigS2))
        
        llik <-  lmvgamma((v+ndimCov+p-1)/2, p) - lmvgamma((v+p-1)/2, p) - 0.5*(ndimCov*p*log(pi)) - 0.5*p*Ccov - 0.5*ndimCov*detSig - Kdet
        
    }else{
        Gi <- try(chol(object$sigma$Pinv), silent=TRUE)
        if(inherits(Gi, 'try-error')) return("error")
        quadprod <- sum(backsolve(Gi, t(residuals), transpose = TRUE)^2)
        detValue <- sum(2*log(diag(Gi)))
        llik <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov + ndimCov*detValue + quadprod)
    }
    
    
    # Estimate parameters on bootstrap samples
    bias <- .parallel_mapply(function(i){
      
      # generate bootstrap sample: TODO check degenerate case when all the resampled values are identical?
      Yp <- MeanNull + Dsqrt%*%(residuals[sample(N, replace=TRUE),])*diagWeight # sampling with replacement for bootstrap
      rownames(Yp) <- rownames(object$variables$Y)
      
      modelPerm$response <- quote(Yp);
      estimModelNull <- eval(modelPerm);
      d1res <- D1(objectBoot=estimModelNull, objectFit=object, ndimCov=ndimCov, m=m, p=p, sqM=DsqrtInv, Ccov2=Ccov)
      d3res <- D3(objectBoot=estimModelNull, objectFit=object, loglik=llik, ndimCov=ndimCov, m=m, p=p)
      d1res+d3res
      
    }, 1:nboot, mc.cores = getOption("mc.cores", nbcores))
    
    ## Sometimes, warnings are added to the object (string) and then any mathematical operation is prevented
    if(class(bias)=='list' & length(bias)>1) bias <- bias[[1]]
    
    # check for errors first?
    bias <- .check_samples(bias)
    nboot_eff <- length(bias)
    
    # compute the EIC
    pboot <- mean(bias)
    EIC <- -2*llik + 2*pboot
    
    # standard-error
    se <- sd(bias)/sqrt(nboot_eff)
    
    # concatenate the results
    results <- list(EIC=EIC, bias=bias, LogLikelihood=llik, se=se, p=p, n=N)
    class(results) <- c("eic.mvgls","eic")
    
    return(results)
  }

# ------------------------------------------------------------------------- #
# .reml_to_ml                                                               #
# options: object, ...                                                      #
# Compute the log-likelihood using reml estimates                           #
# ------------------------------------------------------------------------- #
.reml_to_ml <- function(object, ...){
    
    # parameters
    ndimCov = object$dims$n
    p = object$dims$p
    
    # Maximum Likelihood determinant from REML fit
    Ccov <- as.numeric(object$corrSt$det - determinant(crossprod(object$corrSt$X))$modulus + object$corrSt$const)
    
    # residuals
    residuals <- object$corrSt$Y - object$corrSt$X%*%object$coefficients
    
    # switch between methods
    if(object$method=="EmpBayes"){
        
        # df for the Matrix T distribution
        v = p+1
        
        if(object$target=="Variance"){
          target = colSums(residuals^2)*(1/(ndimCov-object$dims$m))*object$tuning
          SigS2 <- svd(t(residuals)*sqrt(1/(target*(v-p))), nu=0, nv=0)$d^2
          detSig <- sum(log(target*(v-p)))
        }else{
          tuning = mean(colSums(residuals^2)*(1/(ndimCov-object$dims$m)))*object$tuning
          SigS2 <- svd(residuals*sqrt(1/((v-p)*tuning)), nu=0, nv=0)$d^2
          detSig <- p*log((v-p)*tuning) # note, with default df, v-p=1
        }
        
        Kdet <- 0.5*(v+ndimCov+p-1)*sum(log(1+SigS2))
        
        llik <- lmvgamma((v+ndimCov+p-1)/2, p) - lmvgamma((v+p-1)/2, p) - 0.5*(ndimCov*p*log(pi)) - 0.5*p*Ccov - 0.5*ndimCov*detSig - Kdet
        
        
    }else{
        # using the Matrix normal distribution
        Gi <- try(chol(object$sigma$Pinv), silent=TRUE)
        if(inherits(Gi, 'try-error')) return("error")
        quadprod <- sum(backsolve(Gi, t(residuals), transpose = TRUE)^2)
        detValue <- sum(2*log(diag(Gi)))
        llik <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov + ndimCov*detValue + quadprod)
    }
    
    return(llik)
}

# ------------------------------------------------------------------------- #
# .llik_fn_eic                                                              #
# options: object_boot, object_emp, residualsBoot, method, Ccov, n, p, v,   #
# bias_type, ...                                                            #
# Compute the log-likelihood with bootstrap or empirical model fit in EIC   #
# ------------------------------------------------------------------------- #
.llik_fn <- function(object_boot, object_emp, residualsBoot, method, Ccov, n, m, p, v, bias_type, ...){
  
  if(method=="EmpBayes"){
    
    
    # Switch between various evaluation of the bootstrapped samples
    switch(bias_type,
           "D1_1"={ # Y*|param* or for Y|param*
             
             if(object_emp$target=="Variance"){
               target = colSums(residualsBoot^2)*(1/(n-m))*object_boot$tuning
               SigS2 <- svd(t(residualsBoot)*sqrt(1/(target*(v-p))), nu=0, nv=0)$d^2
               detSig <- sum(log(target*(v-p)))
             }else{
                 tuning = mean(colSums(residualsBoot^2)*(1/(n-m)))*object_boot$tuning
               SigS2 <- svd(residualsBoot*sqrt(1/((v-p)*tuning)), nu=0, nv=0)$d^2
               detSig <- p*log((v-p)*tuning) # note, with default df, v-p=1
             }
             
           },
           "D1_2"={ # Y*|param
             
             if(object_emp$target=="Variance"){
               target = colSums(residualsBoot^2)*(1/(n-m))*object_emp$tuning
               SigS2 <- svd(t(residualsBoot)*sqrt(1/(target*(v-p))), nu=0, nv=0)$d^2
               detSig <- sum(log(target*(v-p)))
             }else{
                 tuning = mean(colSums(residualsBoot^2)*(1/(n-m)))*object_emp$tuning
               SigS2 <- svd(residualsBoot*sqrt(1/((v-p)*tuning)), nu=0, nv=0)$d^2
               detSig <- p*log((v-p)*tuning) # note, with default df, v-p=1
             }
             
           })
    
    # compute the log-likelihood
    Kdet <- 0.5*(v+n+p-1)*sum(log(1+SigS2))
    llik <- lmvgamma((v+n+p-1)/2, p) - lmvgamma((v+p-1)/2, p) - 0.5*(n*p*log(pi)) - 0.5*p*Ccov - 0.5*n*detSig - Kdet
    
  }else{
    
    # Gaussian distribution
    # Switch between various evaluation of the bootstrapped samples
    switch(bias_type,
           "D1_1"={ # Y*|param* or for Y|param*
             
             Gi1 <- try(chol(object_boot$sigma$Pinv), silent=TRUE)
             if(inherits(Gi1, 'try-error')) return("error")
             quadprod <- sum(backsolve(Gi1, t(residualsBoot), transpose = TRUE)^2)
             detValue <- sum(2*log(diag(Gi1)))
             
             
           },
           "D1_2"={ # Y*|param
             
             Gi2 <- try(chol(object_emp$sigma$Pinv), silent=TRUE)
             if(inherits(Gi2, 'try-error')) return("error")
             quadprod <- sum(backsolve(Gi2, t(residualsBoot), transpose = TRUE)^2)
             detValue <- sum(2*log(diag(Gi2)))
             
           })
    
    # compute the log-likelihood
    llik <- -0.5 * (n*p*log(2*pi) + p*Ccov + n*detValue + quadprod)
  }
  return(llik)
}


# ------------------------------------------------------------------------- #
# fitted.values.mvgls  / fitted.mvgls                                       #
# options: object, ...                                                      #
# S3 method from "stats" package                                            #
# ------------------------------------------------------------------------- #
fitted.mvgls <- function(object, ...){
    return(object$fitted)
}

# ------------------------------------------------------------------------- #
# residuals.mvgls                                                           #
# options: type = c("response","normalized")                                #
# S3 method "mvgls" class                                                   #
# ------------------------------------------------------------------------- #
residuals.mvgls <- function(object, type=c("response","normalized"), ...){
    type <- match.arg(type)[1]
    if(.mvgls_is_bmmcorr(object) && type=="normalized"){
        .mvgls_bmmcorr_loglik_guard(object)
        residuals <- .mvgls_bmmcorr_normalized_residuals(object)
    }else if(type=="response"){
        residuals <- object$residuals
    }else{
        residuals <- .mvGLS(object$corrSt)$residuals
    }
    return(residuals)
}

# ------------------------------------------------------------------------- #
# vcov.mvgls                                                                #
# options: object, ...                                                      #
# S3 method from "stats" package                                            #
# ------------------------------------------------------------------------- #
vcov.mvgls <- function(object, ...){
    args <- list(...)
    if(is.null(args[["type"]])) type <- "coef" else type <- args$type
    
    if(.mvgls_is_bmmcorr(object)){
        switch(type,
        "covariance"={return(object$sigma$Pinv)},
        "precision"={return(object$sigma$P)},
        "coef"={
            stop("coefficient covariance is not implemented for ", .mvgls_bmmcorr_label(object), " BMM fits")
        })
    }
    
    switch(type,
    "covariance"={return(object$sigma$Pinv)}, # inverse of the precision matrix
    "precision"={return(object$sigma$P)},     # precision matrix
    "coef"={
        XtX <- solve(crossprod(object$corrSt$X))
        covBeta <- kronecker(object$sigma$Pinv, XtX)
        rownames(covBeta) <- colnames(covBeta) <- rep(attr(object$variables$X,"dimnames")[[2]], object$dims$p)
       
        return(covBeta)})
}

# ------------------------------------------------------------------------- #
# coef.mvgls     / coefficients.mvgls                                       #
# options: object, ...                                                      #
# S3 method from "stats" package                                            #
# ------------------------------------------------------------------------- #
coef.mvgls <- function(object, ...){
    
    coeffs <- object$coefficients
    rownames(coeffs) <- attr(object$variables$X,"dimnames")[[2]]
    
    return(coeffs)
}

.mvgls_oum_coefficient_blocks <- function(object){
    coeffs <- coef(object)
    if(!identical(object$model, "OUM")){
        return(list(combined=coeffs, regime_optima=NULL, regression_effects=NULL))
    }

    regime_names <- colnames(object$variables$regimes)
    if(is.null(regime_names) || !length(regime_names)){
        return(list(combined=coeffs, regime_optima=coeffs, regression_effects=NULL))
    }

    regime_idx <- match(regime_names, rownames(coeffs), nomatch=0L)
    regime_idx <- regime_idx[regime_idx > 0L]
    if(!length(regime_idx)){
        return(list(combined=coeffs, regime_optima=coeffs, regression_effects=NULL))
    }

    regression_idx <- setdiff(seq_len(nrow(coeffs)), regime_idx)

    list(
        combined=coeffs,
        regime_optima=coeffs[regime_idx, , drop=FALSE],
        regression_effects=if(length(regression_idx)) coeffs[regression_idx, , drop=FALSE] else NULL
    )
}

.mvgls_print_coefficient_matrix <- function(coeffs, digits, heading="Coefficients:"){
    if(is.null(coeffs) || nrow(coeffs) == 0L) return(invisible(NULL))

    if(ncol(coeffs) < 10L){
        cat(heading, "\n", sep="")
        print.default(format(coeffs, digits = digits), print.gap = 2L, quote = FALSE)
    } else {
        cat(sub(":$", "", heading), " (truncated):\n", sep="")
        coefHead <- coeffs[, 1:10, drop=FALSE]
        print(coefHead, digits = digits, quote = FALSE)
        cat("Use \"coef\" to display all the coefficients\n")
    }
    cat("\n")
    invisible(NULL)
}

.mvgls_print_coefficients <- function(object, digits){
    coef_blocks <- .mvgls_oum_coefficient_blocks(object)
    if(identical(object$model, "OUM") && !is.null(coef_blocks$regime_optima)){
        .mvgls_print_coefficient_matrix(coef_blocks$regime_optima, digits, heading="Regime optima (theta):")
        if(!is.null(coef_blocks$regression_effects)){
            .mvgls_print_coefficient_matrix(coef_blocks$regression_effects, digits, heading="Regression effects (beta):")
        }
    }else{
        .mvgls_print_coefficient_matrix(coef_blocks$combined, digits, heading="Coefficients:")
    }
    invisible(NULL)
}

# ------------------------------------------------------------------------- #
# confint.mvgls                                                             #
# options: object, parm, level, method, ...                                 #
# S3 method from "stats" package                                            #
# ------------------------------------------------------------------------- #
confint.mvgls <- function(object, parm=NULL, level=0.95, method=c("profile","bootstrap","both"), nboot=100L, nbcores=1L, profile_points=9L, ...){
    method <- match.arg(method)
    if(!.mvgls_is_bmmcorr(object)){
        stop("confint() is only implemented for mvgls fits with bmm.structure=\"corrstrength\" or \"corrpower\"")
    }
    .mvgls_bmmcorr_loglik_guard(object)

    diag_obj <- if(.mvgls_is_corrpower(object)){
        corrpower_diagnostics(
            object,
            level=level,
            nboot=if(identical(method, "profile")) 0L else nboot,
            nbcores=nbcores,
            profile_points=if(identical(method, "bootstrap")) 0L else profile_points,
            anchors=FALSE
        )
    }else{
        corrstrength_diagnostics(
            object,
            level=level,
            nboot=if(identical(method, "profile")) 0L else nboot,
            nbcores=nbcores,
            profile_points=if(identical(method, "bootstrap")) 0L else profile_points,
            anchors=FALSE
        )
    }

    summary_df <- diag_obj$parameter_summary
    requested <- .mvgls_bmmcorr_parse_parm(object, parm)
    summary_df <- summary_df[match(requested, summary_df$label), , drop=FALSE]
    rownames(summary_df) <- summary_df$label

    if(identical(method, "profile")){
        out <- cbind(lower=summary_df$profile_low, upper=summary_df$profile_high)
    }else if(identical(method, "bootstrap")){
        out <- cbind(lower=summary_df$bootstrap_low, upper=summary_df$bootstrap_high)
    }else{
        out <- summary_df[, c("estimate", "profile_low", "profile_high", "bootstrap_low", "bootstrap_high"), drop=FALSE]
    }
    rownames(out) <- summary_df$label
    out
}

# ------------------------------------------------------------------------- #
# logLik.mvgls                                                              #
# options: object, ...                                                      #
# S3 method from "stats" package                                            #
# ------------------------------------------------------------------------- #
logLik.mvgls<-function(object,...){
    
    if(.mvgls_is_bmmcorr(object)){
        .mvgls_bmmcorr_loglik_guard(object)
        return(object$logLik)
    }
    
    if(object$method=="LL" | object$method=="EmpBayes"){
        LL = object$logLik # it's already the LL returned
    }else{
        # param
        n <- object$dims$n
        p <- object$dims$p
        m <- object$dims$rank
        if(object$REML) ndimCov = n - m else ndimCov = n
        DP <- as.numeric(determinant(object$sigma$Pi)$modulus)
        Ccov <- object$corrSt$det
        LL <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov + ndimCov*DP + ndimCov*sum(object$sigma$S*object$sigma$P))
    }
    return(LL)
}

# ------------ S3 Printing Methods ----------------------------------------- #

# Generic S3 print for linear models in R stats library (R core team).
print.mvgls <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
    # model call
    cat("\nCall:\n",
    paste(deparse(x$call), sep = "", collapse = "\n"), "\n\n", sep = "")
    
    # loocv or LL
    meth <- ifelse(x$REML, "REML", "ML")
    if(x$method=="LL" | x$method=="EmpBayes"){
        if(inherits(x, "mvols")) cat("\nOrdinary least squares fit by",meth,"\n") else cat("\nGeneralized least squares fit by",meth,"\n")
        if(x$REML) cat("Log-restricted-likelihood:",round(x$logLik, digits=digits), "\n\n") else cat("Log-likelihood:",round(x$logLik, digits=digits), "\n\n")
    }else{
        if(inherits(x, "mvols")) cat("\nOrdinary least squares fit by penalized",meth,"\n") else cat("\nGeneralized least squares fit by penalized",meth,"\n")
        if(x$REML){
            cat("LOOCV of the log-restricted-likelihood:",round(x$logLik, digits=digits), "\n\n")
        }else{
            cat("LOOCV of the log-likelihood:",round(x$logLik, digits=digits), "\n\n")
        }
    }
    
    
    # Model parameters
    cat("\nParameter estimate(s):\n")
    if(!any(is.na(x$param))){
        switch(x$model,
        "OU"={ cat("alpha:",round(x$param, digits=digits),"\n\n")},
        "OUM"={ cat("alpha:",round(x$param, digits=digits),"\n\n")},
        "EB"={ cat("r:",round(x$param, digits=digits),"\n\n")},
        "lambda"={cat("lambda:",round(x$param, digits=digits),"\n\n")},
        "BMM"={print(round(x$param, digits=digits)); cat("\n")},
        cat("parameter(s):",round(x$param, digits=digits),"\n\n")
        )
    }

    if(.mvgls_is_bmmcorr(x)){
        .mvgls_print_bmmcorr_regime_summary(x, digits = digits)
        .mvgls_print_bmmcorr_diagnostics(x, digits = digits)
    }
    
    # Regularization parameter
    if(!is.na(x$tuning)){
        cat("Regularization parameter (gamma):", round(x$tuning, digits=digits), "\n\n")
    }
    
    # size of the evolutionary covariance matrix
    cat("\nCovariance matrix of size:",x$dims$p,"by",x$dims$p,"\n")
    cat("for",x$dims$n,"observations","\n\n")
    
    # coefficients of the linear model
    .mvgls_print_coefficients(x, digits)
    invisible(x)
}


# Generic S3 print for linear models in R stats library (R core team).
print.summary.mvgls <- function(x, digits = max(3, getOption("digits") - 3), ...){
   
    # model call
    cat("\nCall:\n",
    paste(deparse(x$call), sep = "", collapse = "\n"), "\n\n", sep = "")
    
    # loocv or LL
    meth <- ifelse(x$REML, "REML", "ML")
    
    if(x$method=="LL" | x$method=="EmpBayes"){
        if(x$GLS) cat("\nGeneralized least squares fit by",meth,"\n") else cat("\nOrdinary least squares fit by",meth,"\n")
        print(x$results.fit,  quote = FALSE )
    }else{
        
        if(x$GLS) cat("\nGeneralized least squares fit by penalized",meth,"\n") else cat("\nOrdinary least squares fit by penalized",meth,"\n")
        print(x$results.fit,  quote = FALSE )
    }
    
    
    # Model parameters
    cat("\nParameter estimate(s):\n")
    if(!any(is.na(x$param))){
        switch(x$model,
        "OU"={ cat("alpha:",round(x$param, digits=digits),"\n\n")},
        "OUM"={ cat("alpha:",round(x$param, digits=digits),"\n\n")},
        "EB"={ cat("r:",round(x$param, digits=digits),"\n\n")},
        "lambda"={cat("lambda:",round(x$param, digits=digits),"\n\n")},
        "BMM"={print(round(x$param, digits=digits)); cat("\n")},
        cat("parameter(s):",round(x$param, digits=digits),"\n\n")
        )
    }

    if(.mvgls_is_bmmcorr(x)){
        .mvgls_print_bmmcorr_regime_summary(x, digits = digits)
        .mvgls_print_bmmcorr_diagnostics(x, digits = digits)
    }
    
    # Regularization parameter
    if(!is.na(x$tuning)){
        cat("Regularization parameter (gamma):", round(x$tuning, digits=digits), "\n\n")
    }
    
    # Error parameter
    if(is.numeric(x$mserr)){
        cat("Nuisance parameter (error variance):", round(x$mserr, digits=digits), "\n\n")
    }
    
    # size of the evolutionary covariance matrix
    cat("\nCovariance matrix of size:",x$dims$p,"by",x$dims$p,"\n")
    cat("for",x$dims$n,"observations","\n\n")
    
    # coefficients of the linear model
    .mvgls_print_coefficients(x, digits)
    invisible(x)
}

summary.mvgls <- function(object, ...){
    
    # param
    n <- object$dims$n
    p <- object$dims$p
    m <- object$dims$rank
    if(object$REML) ndimCov = n - m else ndimCov = n
    
    # loocv or LL
    meth <- ifelse(object$REML, "REML", "ML")

    if(.mvgls_is_bmmcorr(object)){
        .mvgls_bmmcorr_loglik_guard(object)
        LL <- object$logLik
        nparam <- .mvgls_bmm_nparam(object)
        AIC <- -2*LL + 2*nparam
        GIC <- -2*LL + 2*nparam
        results.fit <- data.frame("AIC"=AIC, "GIC"=GIC, "logLik"=LL, row.names = " ")
        
        object$results.fit <- results.fit
        object$regime.summary <- .mvgls_bmmcorr_regime_summary(object)
        object$GLS <- if(inherits(object,"mvols")) FALSE else TRUE
        class(object) <- c("summary.mvgls","mvgls")
        return(object)
    }
    
    if(object$method=="LL"){
        LL = object$logLik
        nparam = if(object$model=="BM") (length(object$start_values)-1) + length(object$coefficients) + p*(p + 1)/2 else length(object$start_values) + length(object$coefficients) + p*(p + 1)/2 
        # AIC
        AIC = -2*LL+2*nparam
        # GIC
        GIC = GIC(object)$GIC
        
        results.fit <- data.frame("AIC"=AIC, "GIC"=GIC, "logLik"=LL, row.names = " ")
        
    }else if(object$method=="EmpBayes"){
        
        LL = object$logLik
        nparam = if(object$model=="BM") (length(object$start_values)-1) + length(object$coefficients)  else length(object$start_values) + length(object$coefficients) # the covariance matrices were marginalized
        
        # AIC
        AIC = -2*LL+2*nparam
        # GIC
        GIC = GIC(object)$GIC
        
        results.fit <- data.frame("AIC"=AIC, "GIC"=GIC, "logLik"=LL, row.names = " ")
        
    }else{
        # LogLikelihood (minus)
        DP <- as.numeric(determinant(object$sigma$Pi)$modulus)
        Ccov <- object$corrSt$det
        LL <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov + ndimCov*DP + ndimCov*sum(object$sigma$S*object$sigma$P))
        # GIC
        GIC = GIC(object)$GIC
        results.fit <- data.frame("GIC"=GIC, "logLik"=LL, row.names = " ")
    }
    
    
    object$results.fit <- results.fit
    object$GLS <- if(inherits(object,"mvols")) FALSE else TRUE
    class(object) <- c("summary.mvgls","mvgls")
    object
}

# AIC printing options
print.aic.mvgls<-function(x,...){
    cat("\n")
    message("-- Akaike Information Criterion --","\n")
    cat("AIC:",x$AIC,"| Log-likelihood",x$LogLikelihood,"\n")
    cat("\n")
}

# BIC printing options
print.bic.mvgls<-function(x,...){
    cat("\n")
    message("-- Bayesian Information Criterion --","\n")
    cat("BIC:",x$BIC,"| Log-likelihood",x$LogLikelihood,"\n")
    cat("\n")
}

# GIC printing options
print.gic.mvgls<-function(x,...){
    cat("\n")
    message("-- Generalized Information Criterion --","\n")
    cat("GIC:",x$GIC,"| Log-likelihood",x$LogLikelihood,"\n")
    cat("\n")
}

# EIC printing options
print.eic.mvgls<-function(x,...){
    cat("\n")
    message("-- Extended Information Criterion --","\n")
    cat("EIC:",x$EIC,"| +/-",3.92*x$se,"| Log-likelihood",x$LogLikelihood,"\n")
    cat("\n")
}

# ------------------------------------------------------------------------- #
# .parallel_mapply wrapper switch options for parallel calculation          #
# options: ...                                                              #
#                                                                           #
# ------------------------------------------------------------------------- #

.parallel_mapply <- function(FUN,..., MoreArgs = NULL, mc.style = "ETA", mc.substyle = NA,
                       mc.cores = getOption("mc.cores", 2L),
                       ignore.interactive = getOption("ignore.interactive", F),
                       mc.preschedule = TRUE, mc.set.seed = TRUE, mc.cleanup = TRUE, verbose=TRUE){
    
    # Check if MacOS with new silicon chip is used to turn to safer socket parallel computing rather than forking
    is_macos <- Sys.info()[["sysname"]] == "Darwin"
    is_apple_silicon <- is_macos && grepl("arm64|aarch64", R.version$arch)
    
    # switch depending on the options
    if(is_apple_silicon){
        
        if(verbose){
            # Should switch to pbapply in the long term ? TODO
            cl <- makeCluster(mc.cores)
            result <- pblapply(..., FUN, cl=cl)
            stopCluster(cl)
            return(simplify2array(result))
        }else{
            # Should switch to pbapply in the long term ? TODO
            cl <- makeCluster(mc.cores)
            result <- parLapply(cl = cl, ..., fun=FUN)
            stopCluster(cl)
            return(simplify2array(result))
        }
        
    }else{
    
        if(verbose){
            return(pbmcmapply(FUN, ..., MoreArgs = MoreArgs, mc.style = mc.style, mc.substyle = mc.substyle,    mc.cores = mc.cores,
                              ignore.interactive = ignore.interactive, mc.preschedule = mc.preschedule,     mc.set.seed = mc.set.seed, mc.cleanup = mc.cleanup))
        }else{
            return(mcmapply(FUN, ..., MoreArgs = MoreArgs, mc.cores = mc.cores,
                              mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed, mc.cleanup =  mc.cleanup))
        }
    }
}

# ------------------------------------------------------------------------- #
# .pruning_general wrapper switch options for OLS vs GLS and various models #
# options: tree, inv, scaled, trans, check                                  #
#                                                                           #
# ------------------------------------------------------------------------- #

.pruning_general <- function(tree, inv=TRUE, scaled=TRUE, trans=TRUE, check=TRUE){
    if(inherits(tree, "phylOLS")){
        n <- Ntip(tree)
        if((sum(tree$edge.length) - n)<=.Machine$double.eps){
            # Return the determinant
            det <- 0
            sqrtMat <- diag(n)
        }else{
            descendent <- tree$edge[,2]
            extern <- (descendent <= n)
            if(inv) sqrt_phy <- 1/sqrt(tree$edge.length[extern]) else sqrt_phy <- sqrt(tree$edge.length[extern])
            sqrtMat <- diag(sqrt_phy)
            # Return the determinant => variance terms of the 'star' tree
            det <- sum(2*log(sqrt_phy))
        }
        return(list(sqrtMat=sqrtMat, det=det))
    }else{
        return(pruning(tree, inv=inv, scaled=scaled, trans=trans, check=check))
    }
}
                        
# ------------------------------------------------------------------------- #
# print option for MANOVA tests  (output borrowed from "car" package)       #
# options: x, digits, ...                                                   #
#                                                                           #
# ------------------------------------------------------------------------- #


print.manova.mvgls <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  
  # select the appropriate output
  if(x$param){
    
    if(x$type=="I") cat("Sequential MANOVA Tests:",x$test,"test statistic","\n")
    if(x$type=="II") cat("Type II MANOVA Tests:",x$test,"test statistic","\n")
    if(x$type=="III") cat("Type III MANOVA Tests:",x$test,"test statistic","\n")
    if(x$type=="glh") cat("General Linear Hypothesis Test:",x$test,"test statistic","\n")
    if(x$type=="glhrm") cat("General Linear Hypothesis Test (repeated measures design):",x$test,"test statistic","\n")
    
    signif <- sapply(x$pvalue, function(i) if(i<0.001){"***"}else if(i<0.01){
      "**"}else if(i<0.05){"*"}else if(i<0.1){"."}else{""})
    
    table_results <- data.frame(Df=x$Df, stat=x$stat, approxF=x$approxF, numDf=x$NumDf, denDf=x$DenDf, pval=x$pvalue, signif=signif)
    if(x$type!="glh" & x$type!="glhrm"){
        if(x$type=="III") rownames(table_results) <- x$terms[unique(x$dims$assign)+1] else rownames(table_results) <- x$terms[unique(x$dims$assign)]
    }else{
        rownames(table_results) <- "Contrasts L"
    }
    colnames(table_results) <- c("Df", "test stat", "approx F", "num Df", "den Df", "Pr(>F)", "")
    print(table_results, digits = digits, ...)
    cat("---","\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1","\n")
    
    
    
  }else{ # permutation methods
    
    if(x$type=="I") cat("Sequential MANOVA Tests with",x$nperm,"permutations:",x$test,"test statistic","\n")
    if(x$type=="II") cat("Type II MANOVA Tests with",x$nperm,"permutations:",x$test,"test statistic","\n")
    if(x$type=="III") cat("Type III MANOVA Tests with",x$nperm,"permutations:",x$test,"test statistic","\n")
    if(x$type=="glh")  cat("General Linear Hypothesis Test with",x$nperm,"permutations:",x$test,"test statistic","\n")
    if(x$type=="glhrm")  cat("General Linear Hypothesis Test (repeated measures design) with",x$nperm,"permutations:",x$test,"test statistic","\n")
    signif <- sapply(x$pvalue, function(i) if(i<0.001){"***"}else if(i<0.01){
      "**"}else if(i<0.05){"*"}else if(i<0.1){"."}else{""})
    
    table_results <- data.frame(stat=x$stat, pval=x$pvalue, signif=signif)
    if(x$type!="glh" & x$type!="glhrm"){
        if(x$type=="III") rownames(table_results) <- x$terms[unique(x$dims$assign)+1] else rownames(table_results) <- x$terms[unique(x$dims$assign)]
    }else{
        rownames(table_results) <- "Contrasts L"
    }
    colnames(table_results) <- c("Test stat", "Pr(>Stat)", "")
    print(table_results, digits = digits, ...)
    cat("---","\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1","\n")
    
  }
  
}
                     
                    
# ------------------------------------------------------------------------- #
# plot option for MANOVA tests  (distribution of the test statistic)        #
# options: x, ...                                                           #
#                                                                           #
# ------------------------------------------------------------------------- #

plot.manova.mvgls <- function(x,...){
    
    args <- list(...)
    if(is.null(args[["density"]])) density = FALSE else density = args$density
    if(is.null(args[["breaks"]])) breaks = 50 else breaks = args$breaks
    
    nterms <- length(x$terms)
    
    if(x$param==TRUE){
        
        for(i in 1:nterms){
            df_mod <- x
            d1=df_mod$NumDf[i]
            d2=df_mod$DenDf[i]
            
        curve(df(x, df1=d1, df2=d2), 0, qf(0.9999, d1, d2), las=1, 
           main=paste("F test:",x$terms[i]), 
           xlab=paste("F value","(",round(x$approxF[i],3),")","p-value :", round(x$pvalue[i],3)), ylab="density" );
            abline(v=x$approxF[i], col="red")
        }
        
    }else{
    
    
    # plot histogram with permuted statistics
    for(i in 1:nterms){
      
      if(density){
          plot(density(x$nullstat[,i]), main=paste("Statistic distribution:",x$terms[i]),xlab=paste(x$test,"(",round(x$stat[i],3),")","p-value :",
                round(x$pvalue[i],3)), las=1, xlim=range(c(x$nullstat[,i],x$stat[i])))
      }else{
          hist(x$nullstat[,i], freq=FALSE, main=paste("Statistic distribution:",x$terms[i]),
        xlab=paste(x$test,"(",round(x$stat[i],3),")","p-value :",
                round(x$pvalue[i],3)), las=1, breaks=breaks, border=NA, col="lightgrey", xlim=range(c(x$nullstat[,i],x$stat[i]))); 
      }
        abline(v=x$stat[i], col="red", lwd=2)
        }
    }

}

# ------------------------------------------------------------------------- #
# plot option for Pairwise tests  (distribution of the test statistic)      #
# options: x, ...                                                           #
#                                                                           #
# ------------------------------------------------------------------------- #

plot.pairs.mvgls <- function(x,...){
  
  args <- list(...)
  if(is.null(args[["density"]])) density = FALSE else density = args$density
  if(is.null(args[["breaks"]])) breaks = 50 else breaks = args$breaks
  
  nterms <- nrow(x$L)
  if(is.null(rownames(x$L))) namesContrasts <- "contrasts L" else namesContrasts <- rownames(x$L)
  
  if(x$param==TRUE){
    
    for(i in 1:nterms){
      df_mod <- x
      d1=df_mod$NumDf[i]
      d2=df_mod$DenDf[i]
      
      curve(df(x, df1=d1, df2=d2), 0, qf(0.9999, d1, d2), las=1,
            main=paste("F test:", namesContrasts[i]),
            xlab=paste("F value","(",round(x$approxF[i],3),")","p-value :", round(x$pvalue[i],3)), ylab="density" );
      abline(v=x$approxF[i], col="red")
    }
    
  }else{
    
    
    # plot histogram with permuted statistics
    for(i in 1:nterms){
      
      if(density){
        plot(density(x$nullstat[,i]), main=paste("Statistic distribution:",x$terms[i]),xlab=paste(x$test,"(",round(x$stat[i],3),")","p-value :",
                                                                                                  round(x$pvalue[i],3)), las=1, xlim=range(c(x$nullstat[,i],x$stat[i])))
      }else{
        hist(x$nullstat[,i], main=paste("Statistic distribution:",namesContrasts[i]),
             xlab=paste(x$test,"(",round(x$stat[i],3),")","p-value :",
                        round(x$pvalue[i],3)), las=1, breaks=breaks, border=NA, col="lightgrey", xlim=range(c(x$nullstat[,i],x$stat[i])));
      }
      abline(v=x$stat[i], col="red", lwd=2)
    }
  }
  
}



# ------------------------------------------------------------------------- #
# plot.mvgls                                                                #
# options: x, term, ..., fitted=TRUE                                        #
#                                                                           #
# ------------------------------------------------------------------------- #
plot.mvgls <- function(x, term, ..., fitted=FALSE, residuals=FALSE){
    
    if(missing(term)){
        term <- which(attr(x$variables$X,"dimnames")[[2]]!="(Intercept)")[1]
        term <- attr(x$variables$X,"dimnames")[[2]][term]
    }
    
    if(!is.numeric(term) & !term%in%attr(x$variables$X,"dimnames")[[2]]) stop("Unknown predictor name.","\n")
    
    # based on Drake & Klingenberg 2008 shape score
    betas <- coefficients(x)[term,,drop=TRUE]
    standBeta <- betas %*% sqrt(solve(crossprod(betas)))
    if(residuals) scoreVar <- (x$residuals)%*% standBeta else scoreVar <- (x$variables$Y)%*% standBeta
    
    # plot
    plot(scoreVar ~ x$variables$X[,term], xlab=term, ylab="mvScore", ...)
    
    # plot predictions on the same space?
    if(fitted){
        scoreVar2 <- (x$fitted ) %*% standBeta
        points(scoreVar2 ~ x$variables$X[,term], col="red", pch=16)
    }
    
    # loess on the residuals?
    if(residuals){
        abline(h=0, lty=2)
        scores_residuals <- data.frame(score=scoreVar, xvar=x$variables$X[,term])
        loess_fit <- loess(score ~ xvar, data=scores_residuals)
        xseq <- seq(from=min(scores_residuals$xvar), to=max(scores_residuals$xvar), length=80)
        pred <- predict(loess_fit, newdata=data.frame(xvar=xseq))
        lines(pred~xseq, col="red", xpd=FALSE)
    }
    
    
    results <- list(scores = scoreVar, standBeta=standBeta, betas=betas, term=term)
    invisible(results)
}

# ------------------------------------------------------------------------- #
# predict.mvgls                                                             #
# options: object, newdata, ...                                             #
#                                                                           #
# ------------------------------------------------------------------------- #
.mvgls_oum_predictor_matrix <- function(object, X_formula=NULL, tree=object$variables$tree, rows=NULL, include_internal=FALSE){
    root <- object$root %||% "stationary"
    root_std <- if(identical(root, "stationary")) 1L else 0L
    nterm <- if(include_internal) Ntip(tree) + Nnode(tree) else Ntip(tree)
    W <- .mvgls_oum_weight_matrix(tree, object$param, root=root, std=root_std, nterm=nterm)

    if(!is.null(rows)){
        if(any(!rows %in% rownames(W))) stop("Could not match OUM regimes to the requested rows")
        W <- W[rows, , drop=FALSE]
    }else if(!is.null(X_formula) && nrow(X_formula) != nrow(W)){
        stop("OUM prediction requires row names matching the tree tips or nodes")
    }

    assign_formula <- object$dims$formula.assign %||% integer(0)
    if(is.null(X_formula) || length(assign_formula) == 0L) return(W)

    cov_idx <- which(assign_formula != 0)
    if(length(cov_idx) == 0L) return(W)
    if(nrow(X_formula) != nrow(W)) stop("The OUM covariate design does not match the regime design")

    cbind(W, X_formula[, cov_idx, drop=FALSE])
}

predict.mvgls <- function(object, newdata, ...){
    
    args <- list(...)
    # if "tree" is provided
    if(!is.null(args[["tree"]])){
        if(!inherits(args$tree, "phylo")) stop("the provided tree is not of class \"phylo\" ") else tree <- args$tree
        if(!is.data.frame(newdata)) stop("the \"newdata\" should be a data.frame object with column names matching predictors names, and row names matching names in the tree ")
    } else tree <- NULL
    if(is.null(args[["na.action"]])) na.action <- na.pass else na.action <- args$na.action
    
    # check if newdata is provided
    if(missing(newdata) || is.null(newdata)) {
        X <- object$variables$X # simply return fitted values when newdata is empty
    }else{
        
        Terms <- delete.response(object$terms)
        # as in "stats v3.3.0"
        m <- model.frame(Terms, newdata, xlev = object$xlevels, na.action = na.action)
        
        # check the arguments
        if(!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
        X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
        
        # FIXME allow lists
        predictors_names <- rownames(newdata)
        if(object$model=="OUM"){
            tree_use <- if(is.null(tree)) object$variables$tree else tree
            if(is.null(predictors_names)){
                if(nrow(X) == Ntip(tree_use)){
                    predictors_names <- tree_use$tip.label
                }else{
                    stop("OUM prediction requires row names on \"newdata\" matching tree tips")
                }
            }
            X <- .mvgls_oum_predictor_matrix(object, X_formula=X, tree=tree_use, rows=predictors_names)
        }
    }
    
    
    # GLS/OLS prediction
    if(is.null(tree)){
        predicted <- X%*%object$coefficients # simply return fitted values when newdata is empty
    }else if(.mvgls_is_bmmcorr(object)){
        .mvgls_bmmcorr_loglik_guard(object)
        if(missing(newdata) || is.null(newdata)) stop(.mvgls_bmmcorr_label(object, capitalize=TRUE), " phylogenetic prediction requires a \"newdata\" data.frame with row names matching target tips")
        if(is.null(rownames(newdata))) stop(.mvgls_bmmcorr_label(object, capitalize=TRUE), " phylogenetic prediction requires row names on \"newdata\" matching target tips")
        predicted <- .mvgls_bmmcorr_predict_phylo(object, X, predictors_names, tree)
    }else{
        rcov <- .resid_cov_phylo(tree, object, predictors_names)
        predicted <- X%*%object$coefficients + rcov$w%*%solve(rcov$Vt)%*%object$residuals[rcov$train,,drop=FALSE] # FIXME account for (multivariate) variance scaling? Rao & Toutenberg => no just the correlation structure
    }
    
    return(predicted)
}

# ------------------------------------------------------------------------- #
# .resid_cov_phylo                                                          #
# options: tree, object, sp_name, ...                                       #
#                                                                           #
# ------------------------------------------------------------------------- #
.resid_cov_phylo <- function(tree, object, sp_name, ...){
    
    if(is.null(sp_name)) stop("You must provide species names to \"newdata\"")
    if(any(!sp_name%in%tree$tip.label)) stop("the \"newdata\" names does not matches names in the tree ")
    train_sample <- tree$tip.label[!tree$tip.label%in%sp_name]
    
    # check first that species in the training sample are the same as in the model fit object
    if(any(!train_sample%in%object$corrSt$phy$tip.label)) train_sample <- object$corrSt$phy$tip.label
    
    # helper to obtain the covariances between data used in a model and newdata
    switch(object$model,
    "BM"={ V <- vcv.phylo(tree)},
    "OU"={
        V <- .Call("mvmorph_covar_ou_fixed", A=vcv.phylo(tree), alpha=as.double(object$param), sigma=1, PACKAGE="mvMORPH")
        rownames(V) <- colnames(V) <- tree$tip.label
    },
    "OUM"={
        V <- .Call("mvmorph_covar_ou_fixed", A=vcv.phylo(tree), alpha=as.double(object$param), sigma=1, PACKAGE="mvMORPH")
        rownames(V) <- colnames(V) <- tree$tip.label
    },
    "EB"={ V <- vcv.phylo(.transformPhylo(tree, model="EB", param=object$param)) },
    "lambda"={ V <- vcv.phylo(.transformPhylo(tree, model="lambda", param=object$param)) },
    #FIXME -- add BMM
    "BMM"={stop("BMM model is not handled yet. Please contact the author for further assistance.")},
    )
    
    # If error=TRUE, we add it to the covariance matrix here
    #if(!is.na(object$mserr)) diag(V) = diag(V) + object$mserr
    
    # Build the covariance matrices
    w <- V[sp_name, train_sample, drop=FALSE]
    Vt <- V[train_sample, train_sample, drop=FALSE]
    
    # return the covariances
    results <- list(w=w, Vt=Vt, train=train_sample)
    return(results)
}

# ------------------------------------------------------------------------- #
# .transformPhylo                                                           #
# options: phy, model, param, ...                                           #
#                                                                           #
# ------------------------------------------------------------------------- #
.transformPhylo <- function(phy, model, param, ...){
    
    # precomputations
    n <- Ntip(phy)
    parent <- phy$edge[,1]
    descendent <- phy$edge[,2]
    extern <- (descendent <= n)
    
    switch(model,
    "EB"={
        if (param!=0){
            distFromRoot <- node.depth.edgelength(phy)
            phy$edge.length = (exp(param*distFromRoot[descendent])-exp(param*distFromRoot[parent]))/param
        }
    },
    "lambda"={
        # Pagel's lambda tree transformation
        if(param!=1) {
            root2tipDist <- node.depth.edgelength(phy)[1:n] # for non-ultrametric trees. The 'up' limit should be exactly 1 to avoid singularity issues
            phy$edge.length <- phy$edge.length * param
            phy$edge.length[extern] <- phy$edge.length[extern] + (root2tipDist * (1-param))
        }
    },)
    
    return(phy)
}

# ------------------------------------------------------------------------- #
# print option for effect/association of multivariate tests                 #
# options: x, digits, ...                                                   #
#                                                                           #
# ------------------------------------------------------------------------- #

print.effects.mvgls <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
    
    if(x$adjusted & x$parametric==TRUE){
        tatsuoka <- attr(x$adjusted, "tatsuoka")
        cat("## Multivariate measure(s) of association ##","\n")
        if(tatsuoka) cat("## Tatsuoka' bias adjustment              ##","\n") else cat("## Serlin' bias adjustment                ##","\n")
        print(x$effect, digits=digits)
    }else{
        cat("## Multivariate measure(s) of association ##","\n")
        print(x$effect, digits=digits)
        if(x$adjusted) cat("## Note: bias is empirically adjusted     ##","\n")
    }
    if(any(x$effect<0)) message("## Values < 0 represent no association    ##","\n")
}

# ------------------------------------------------------------------------- #
# ancestral.mvgls                                                           #
# options: object, ...                                                      #
#                                                                           #
# ------------------------------------------------------------------------- #


## S3 Method for ancestral states estimation
ancestral <- function(object, ...) UseMethod("ancestral")

# core function
ancestral.mvgls <- function(object, ...){
    
    # arguments
    args <- list(...)
    
    # extract objects
    if(!inherits(object,"mvgls")){
        
        # wrapper to "estim"
        if(is.null(args[["data"]]) | is.null(args[["tree"]])) stop("Need a \"tree\" object and a new \"data\" matrix to predict ancestral states. See ?estim")
        estim(tree = args$tree, data = args$data, object = object, asr=TRUE)
        
    }else if(inherits(object,"mvgls")){
        # If regression, must provide a regressor for the ancestral (nodes) states
        # check if newdata is provided
        if(any(object$dims$assign>0)){
            
            if(is.null(args[["na.action"]])) na.action <- na.pass else na.action <- args$na.action
            if(is.null(args[["newdata"]])) stop("Regression model. You must provide a new \"dataset\" of predictors for each nodes. See also ?predict")
            
            Terms <- delete.response(object$terms)
            # as in "stats v3.3.0"
            m <- model.frame(Terms, args$newdata, xlev = object$xlevels, na.action = na.action)
            
            # check the arguments
            if(!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
            X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
            if(object$model=="OUM"){
                if(nrow(X) != Nnode(object$variables$tree)){
                    stop("OUM ancestral reconstruction requires one row of predictors per internal node")
                }
                node_rows <- paste("node_", Ntip(object$variables$tree) + seq_len(Nnode(object$variables$tree)), sep="")
                X <- .mvgls_oum_predictor_matrix(object, X_formula=X, tree=object$variables$tree, rows=node_rows, include_internal=TRUE)
            }
            predicted_fit <- X %*% object$coefficients
            if(.mvgls_is_bmmcorr(object) && nrow(predicted_fit) != Nnode(object$variables$tree)){
                stop(.mvgls_bmmcorr_label(object, capitalize=TRUE), " ancestral state estimation requires one row of predictors per internal node")
            }
            
        }else{
            # Just use the grand mean - i.e. the ancestral states at the root
            if(.mvgls_is_bmmcorr(object)){
                root_mean <- object$variables$X[1, , drop=FALSE] %*% object$coefficients
                predicted_fit <- matrix(root_mean,
                    nrow=Nnode(object$variables$tree),
                    ncol=object$dims$p,
                    byrow=TRUE
                )
                colnames(predicted_fit) <- colnames(object$variables$Y)
            }else if(object$model=="OUM"){
                node_rows <- paste("node_", Ntip(object$variables$tree) + seq_len(Nnode(object$variables$tree)), sep="")
                predicted_fit <- .mvgls_oum_predictor_matrix(object, tree=object$variables$tree, rows=node_rows, include_internal=TRUE) %*% object$coefficients
            }else{
                predicted_fit <- object$variables$X %*% object$coefficients
                predicted_fit <- predicted_fit[-1,,drop=TRUE]
            }
        }

        if(.mvgls_is_bmmcorr(object)){
            .mvgls_bmmcorr_loglik_guard(object)
            return(.mvgls_bmmcorr_ancestral(object, predicted_fit))
        }
        
        # start estimating ancestral states using GLS
        n <- object$dims$n
        p <- object$dims$p
        
        # covariance for the nodes
        if(!is.null(object$corrSt$diagWeight)){
            V<-.Call("mvmorph_covar_ou_fixed", A=.vcvPhyloInternal(object$variables$tree), alpha=as.double(object$param), sigma=1, PACKAGE="mvMORPH")
        }else{
            V <- .vcvPhyloInternal(object$corrSt$phy)
        }
        indice <- (1:n)
        AY <- V[-indice,indice]
        vY <- V[indice,indice]
        
        # states at the nodes
        residuals_fit <- object$residuals
        recons_t <- (AY%*%pseudoinverse(vY)%*%residuals_fit)+predicted_fit
        colnames(recons_t) = colnames(object$variables$Y)
        rownames(recons_t) = paste("node_",n+1:Nnode(object$variables$tree), sep="")
        #class(recons_t) = "anc.mvgls"
        return(recons_t)
        
    }else{
        stop("only works with \"mvgls\" class objects. See ?mvgls, or use instead \"estim\" function")
    }
    
}
