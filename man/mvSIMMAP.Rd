\name{mvSIMMAP}
\alias{mvSIMMAP}
\title{
Experimental mixed-regime Gaussian models on SIMMAP trees
}
\description{
Fits an experimental mixed-regime model on a painted SIMMAP tree where each mapped regime is assigned one of three process families: Brownian motion (BM), Ornstein-Uhlenbeck (OU), or Early-Burst (EB).
}
\usage{
mvSIMMAP(tree, data, process, error = NULL,
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
Named character vector assigning each SIMMAP regime to one of \code{"BM"}, \code{"OU"}, or \code{"EB"}. If unnamed, it must have one entry per column of \code{tree$mapped.edge}.
}
  \item{error}{
Optional measurement-error variances supplied in the same format used by the other \code{mvMORPH} fitters.
}
  \item{param}{
List of starting values and tuning options. \code{sigma} can be a matrix/vector reused across regimes or a list keyed by regime. \code{alpha} is used only for OU regimes. \code{beta} is used only for EB regimes. \code{low} and \code{up} optionally specify lower and upper bounds for EB \code{beta} values. The root state and OU optima are estimated analytically by GLS rather than supplied as starting values.
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
Control list passed to the optimizer.
}
  \item{diagnostic}{
Whether convergence diagnostics are printed.
}
  \item{echo}{
Whether a fitted-model summary is printed.
}
}
\details{
This function implements a dense likelihood-based mixed-regime extension for painted SIMMAP trees. Each mapped branch segment is treated as a local Gaussian transition and the transitions are composed along each lineage. BM, OU, and EB regimes can therefore coexist on the same painted tree.

The current implementation is experimental and is designed around \code{mvMORPH}'s existing dense likelihood machinery. It is conceptually close to the branchwise mixed-Gaussian strategy used by \code{PCMBase}/\code{PCMFit}, but it does not yet implement their linear-time pruning algorithm. Consequently, \code{mvSIMMAP} is best suited to moderate tree sizes.

The fitted mean parameters always include a root-state vector. Each OU regime contributes an additional regime-specific optimum vector. BM and EB regimes contribute only covariance/process parameters.
}
\value{
\item{LogLik }{The log-likelihood of the fitted model.}
\item{AIC }{Akaike Information Criterion.}
\item{AICc }{Sample size-corrected AIC.}
\item{theta }{Estimated mean parameters, including the root state and any OU regime optima.}
\item{alpha }{Named list of OU pull matrices by regime.}
\item{sigma }{Named list of regime-specific rate matrices.}
\item{beta }{Named vector of EB beta values.}
\item{process }{The normalized regime-to-process assignment vector used in the fit.}
\item{convergence }{Convergence status of the optimizing function.}
\item{hess.values }{Basic Hessian reliability diagnostic.}
\item{param }{List of fit metadata and internal parameterization functions.}
\item{llik }{The corresponding log-likelihood function.}
}
\references{
Clavel J., Escarguel G., Merceron G. 2015. mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data. Methods in Ecology and Evolution, 6(11):1311-1319.

Mitov V., Bartoszek K., Asimomitis G., Stadler T. 2019. Fast likelihood calculation for multivariate Gaussian phylogenetic models with shifts. Theoretical Population Biology. 131:66-78.
}
\author{
Julien Clavel; experimental SIMMAP mixed-model extension in this branch.
}
\seealso{
\code{\link{mvBM}}
\code{\link{mvOU}}
\code{\link{mvEB}}
\code{\link{mvSHIFT}}
\code{\link{mvSIM}}
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
}
}
\keyword{ Ornstein Uhlenbeck }
\keyword{ Brownian Motion }
\keyword{ Early-Burst }
\keyword{ SIMMAP }
\keyword{ Regimes }
