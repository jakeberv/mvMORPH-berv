\name{mvSIMMAP}
\alias{mvSIMMAP}
\title{
Experimental mixed-regime Gaussian models on SIMMAP trees
}
\description{
Fits an experimental mixed-regime model on a painted SIMMAP tree where each mapped regime is assigned one of four process families: Brownian motion (BM), single-optimum Ornstein-Uhlenbeck (OU), multiple-optimum Ornstein-Uhlenbeck (OUM), or Early-Burst (EB).
}
\usage{
mvSIMMAP(tree, data, process, process.groups = NULL, error = NULL,
         param = list(alpha = NULL, sigma = NULL, beta = NULL,
                      decomp = c("cholesky", "diagonal", "diagonalPositive"),
                      decompSigma = c("cholesky", "diagonal")),
         method = c("rpf", "inverse", "pseudoinverse"),
         scale.height = FALSE,
         optimization = c("L-BFGS-B", "Nelder-Mead", "subplex"),
         control = list(maxit = 20000),
         diagnostic = TRUE, echo = TRUE)
}
\arguments{
  \item{tree}{
Phylogenetic tree of class \code{"simmap"} with mapped regimes in \code{tree$maps} and \code{tree$mapped.edge}.
}
  \item{data}{
Matrix or data frame with species in rows and continuous traits in columns.
}
  \item{process}{
Named character vector assigning each SIMMAP regime to one of \code{"BM"}, \code{"OU"}, \code{"OUM"}, or \code{"EB"}. If unnamed, it must have one entry per column of \code{tree$mapped.edge}. \code{"OU"} uses one optimum per process group, whereas \code{"OUM"} keeps a separate optimum for each painted regime.
}
  \item{process.groups}{
Optional named character vector controlling which painted regimes share dynamic process parameters. By default, each regime forms its own group. Regimes assigned the same group label must share the same process family and will share the fitted \code{sigma}, \code{alpha}, and/or \code{beta} parameters for that family. For \code{"OU"}, sharing a group also shares the optimum, so the grouped regimes are treated as a single OU regime for the likelihood. For \code{"OUM"}, grouped regimes still keep separate painted-regime optima while sharing the OU dynamics.
}
  \item{error}{
Optional measurement-error variances supplied in the same format used by the other \code{mvMORPH} fitters.
}
  \item{param}{
List of starting values and tuning options. \code{sigma} can be a matrix/vector reused across process groups or a list keyed by process group or regime. \code{alpha} is used only for OU or OUM groups. \code{beta} is used only for EB groups. \code{low} and \code{up} optionally specify lower and upper bounds for EB \code{beta} values. The root state and OU optima are estimated analytically by GLS rather than supplied as starting values.
}
  \item{method}{
Likelihood evaluation method. The current implementation supports only dense methods: \code{"rpf"}, \code{"inverse"}, and \code{"pseudoinverse"}.
}
  \item{scale.height}{
Whether the tree should be scaled to unit height before fitting.
}
  \item{optimization}{
Optimization routine. The special value \code{"fixed"} returns only the log-likelihood function and starting-value parameterization.
}
  \item{control}{
Control list passed to the optimizer. In addition to the optimizer's own entries, \code{mvSIMMAP} recognizes \code{hessian.tol} for the post-fit Hessian tolerance, \code{retry.unreliable} to enable automatic jittered restarts when the first fit converges to an unreliable or nearly flat solution, \code{retry.max} for the number of such restarts, \code{retry.jitter} for the relative amount of start-value perturbation, and \code{retry.seed} for reproducible restart draws. These special entries are stripped before the remaining control list is sent to the optimizer.
}
  \item{diagnostic}{
Whether convergence diagnostics are printed.
}
  \item{echo}{
Whether a fitted-model summary is printed.
}
}
\details{
This function implements a dense likelihood-based mixed-regime extension for painted SIMMAP trees. Each mapped branch segment is treated as a local Gaussian transition and the transitions are composed along each lineage. BM, OU, OUM, and EB regimes can therefore coexist on the same painted tree.

At present, this is the only public mixed-process SIMMAP entry point in the fork; the planned \code{mvgls} integration remains a separate design problem.

The current implementation is experimental and is designed around \code{mvMORPH}'s existing dense likelihood machinery. The underlying modeling idea was taken from the branchwise mixed-Gaussian formulation developed in \code{PCMBase}/\code{PCMFit}: branch-local Gaussian transitions are composed along the tree to handle changes in process type and parameters among mapped segments. This implementation is not a code port of \code{PCMBase}/\code{PCMFit}, and it does not yet implement their linear-time pruning algorithm. Consequently, \code{mvSIMMAP} is best suited to moderate tree sizes.

By default, each painted regime gets its own dynamic parameter block. Supplying \code{process.groups} allows multiple painted regimes to share the same dynamic OU, OUM, BM, or EB parameters. The fitted mean parameters always include a root-state vector. For \code{"OU"}, each process group contributes a single optimum vector, so grouped painted regimes are collapsed into one OU regime for the likelihood. For \code{"OUM"}, each painted regime contributes its own optimum vector even when several regimes share the same OU dynamics. BM and EB regimes contribute only covariance/process parameters.

Printed model objects show the estimated parameter matrices in detail. Calling \code{summary()} on a fitted \code{mvSIMMAP} object gives a more compact structural view listing the painted regimes, process groups, theta ownership, and the names of the estimated parameter blocks. EB beta values printed at or extremely near zero are annotated as the BM boundary for readability.

The post-fit Hessian check is slightly more informative than the older binary \code{mvMORPH} diagnostic used in some other fitters. Tiny negative eigenvalues inside \code{hessian.tol} are treated as numerical noise, nearly zero eigenvalues are reported as a flat or boundary solution, and clearly negative eigenvalues are reported as a possible saddle point. When \code{retry.unreliable = TRUE}, the function can automatically try a small number of jittered restarts and keep the best retained solution.
}
\value{
\item{LogLik }{The log-likelihood of the fitted model, stored as a standard \code{"logLik"} object with \code{df} and \code{nobs} attributes for compatibility with generic comparison tools such as \code{BIC()}.}
\item{AIC }{Akaike Information Criterion.}
\item{AICc }{Sample size-corrected AIC.}
\item{theta }{Estimated mean parameters, including the root state and any OU/OUM optimum vectors.}
\item{alpha }{Named list of OU/OUM pull matrices by process group.}
\item{sigma }{Named list of process-group-specific rate matrices.}
\item{beta }{Named vector of EB beta values by process group.}
\item{process }{The normalized regime-to-process assignment vector used in the fit.}
\item{process.groups }{The normalized regime-to-process-group assignment vector used in the fit.}
\item{convergence }{Convergence status of the optimizing function.}
\item{hess.values }{Compatibility-friendly binary Hessian diagnostic: \code{0} for reliable and \code{1} otherwise.}
\item{diagnostics }{Detailed post-fit diagnostics, including Hessian status, eigenvalue summary, and any automatic restart information.}
\item{param }{List of fit metadata and internal parameterization functions.}
\item{llik }{The corresponding log-likelihood function.}
}
\references{
Clavel J., Escarguel G., Merceron G. 2015. mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data. Methods in Ecology and Evolution, 6(11):1311-1319.

Mitov V., Bartoszek K., Asimomitis G., Stadler T. 2019. Fast likelihood calculation for multivariate Gaussian phylogenetic models with shifts. Theoretical Population Biology. 131:66-78.

Mitov V., Bartoszek K., Stadler T. PCMBase and PCMFit project documentation. The mixed-Gaussian branchwise formulation in these projects was the main conceptual inspiration for the present \code{mvSIMMAP} implementation.
}
\author{
Julien Clavel; experimental SIMMAP mixed-model extension in this fork.
}
\seealso{
\code{\link{mvBM}}
\code{\link{mvOU}}
\code{\link{mvEB}}
\code{\link{mvSHIFT}}
\code{\link{mvSIM}}
\code{\link{mvgls}}
}
\examples{
\donttest{
set.seed(1)
tree <- phytools::make.simmap(ape::rtree(20),
                              setNames(sample(c("A", "B"), 20, replace = TRUE),
                                       paste0("t", seq_len(20))),
                              model = "ER", nsim = 1)
X <- mvSIM(tree, model = "BMM",
           param = list(sigma = list(diag(2) * 0.05, diag(2) * 0.2),
                        theta = c(0, 0),
                        ntraits = 2))

fit <- mvSIMMAP(
  tree, X,
  process = c(A = "OU", B = "BM"),
  param = list(
    sigma = list(A = diag(2) * 0.05, B = diag(2) * 0.2),
    alpha = list(A = diag(2))
  )
)

fit_shared_ou <- mvSIMMAP(
  tree, X,
  process = c(A = "OU", B = "OU"),
  process.groups = c(A = "ou_shared", B = "ou_shared"),
  param = list(
    sigma = list(ou_shared = diag(2) * 0.1),
    alpha = list(ou_shared = diag(2))
  )
)

fit_shared_dynamics <- mvSIMMAP(
  tree, X,
  process = c(A = "OUM", B = "OUM"),
  process.groups = c(A = "ou_shared", B = "ou_shared"),
  param = list(
    sigma = list(ou_shared = diag(2) * 0.1),
    alpha = list(ou_shared = diag(2))
  )
)
}
}
\keyword{ Ornstein Uhlenbeck }
\keyword{ Brownian Motion }
\keyword{ Early-Burst }
\keyword{ SIMMAP }
\keyword{ Regimes }
