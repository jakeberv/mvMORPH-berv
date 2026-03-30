\name{mvSIMMAP}
\alias{mvSIMMAP}
\title{
Experimental mixed-regime Gaussian models on SIMMAP trees
}
\description{
Fits a mixed-regime Gaussian model on a painted SIMMAP tree where each mapped regime is assigned one of four process families: Brownian motion (BM), single-optimum Ornstein-Uhlenbeck (OU), multiple-optimum Ornstein-Uhlenbeck (OUM), or Early-Burst (EB).
}
\usage{
mvSIMMAP(tree, data, process, process.groups = NULL, error = NULL,
         param = list(alpha = NULL, sigma = NULL, beta = NULL,
                      decomp = c("scalarPositive", "diagonalPositive", "diagonal", "cholesky"),
                      decompSigma = c("cholesky", "diagonal")),
         method = c("rpf", "inverse", "pseudoinverse"),
         scale.height = FALSE,
         optimization = c("L-BFGS-B", "Nelder-Mead", "subplex", "fixed"),
         control = list(maxit = 20000),
         diagnostic = TRUE, echo = TRUE)
}
\arguments{
  \item{tree}{
Phylogenetic tree of class \code{"simmap"} with mapped regimes in \code{tree$maps} and \code{tree$mapped.edge}.
}
  \item{data}{
Matrix or data frame with species in rows and continuous traits in columns. When row names are present they must match the tree tip labels and will be reordered to the phylogeny. NA values are supported through the dense methods currently exposed by \code{mvSIMMAP}. When \code{optimization="fixed"}, \code{data} can also be \code{NULL} to build a simulation scaffold without observed traits; in that case provide \code{param$ntraits} or enough matrix-valued parameter information for \code{mvSIMMAP} to infer the trait dimension.
}
  \item{process}{
Named character vector assigning each painted SIMMAP regime to one of \code{"BM"}, \code{"OU"}, \code{"OUM"}, or \code{"EB"}. If unnamed, it must have one entry per column of \code{tree$mapped.edge}. This argument chooses the process family for each painted regime; it does not by itself control whether regimes share parameter values. \code{"OU"} uses one optimum per process group, whereas \code{"OUM"} keeps a separate optimum for each painted regime.
}
  \item{process.groups}{
Optional named character vector controlling which painted regimes share dynamic process parameters. By default, each regime forms its own group. Regimes assigned the same group label must share the same process family and will share the fitted \code{sigma}, \code{alpha}, and/or \code{beta} parameters for that family. This is the main user-facing way to impose equality among regimes in \code{mvSIMMAP}. For \code{"OU"}, sharing a group also shares the optimum, so the grouped regimes are treated as a single OU regime for the likelihood. For \code{"OUM"}, grouped regimes still keep separate painted-regime optima while sharing the OU dynamics.
}
  \item{error}{
Optional measurement-error variances. For fitting, these should correspond to the same tip and trait layout used by \code{data}; matrix-like inputs can be aligned by matching row names to tip labels and column names to trait names. Missing error entries are treated as zero.
}
  \item{param}{
List of starting values and tuning options. \code{sigma} can be a single matrix or parameter vector reused across process groups, or a list keyed by process group or painted regime. \code{alpha} is used only for OU or OUM groups and follows the same conventions. Unnamed lists are interpreted positionally using the current process-group or regime order, so named inputs are safer when several groups are present. \code{beta} is used only for EB groups and can be supplied as a scalar, vector, or list under the same naming rules. \code{low} and \code{up} optionally specify lower and upper bounds for EB \code{beta} values; they may be scalars, positional vectors matching the EB process-group order, or named vectors keyed by EB process group. The root state and OU optima are estimated analytically by GLS rather than supplied as starting values. By default, each process group gets its own full dynamic block; use \code{process.groups} when several regimes should share a full \code{sigma}, \code{alpha}, or \code{beta} block. The default OU/OUM alpha parameterization is \code{"scalarPositive"}, which fits one positive pull parameter per OU/OUM process group for both \code{"OU"} and \code{"OUM"} models. \code{mvSIMMAP} does not currently expose a public \code{root} or \code{vcv} switch; the present root treatment corresponds to a fixed-root-style mean-parameter formulation. For fixed scaffolds built with \code{data=NULL}, \code{param$ntraits} and optional \code{param$names_traits} can be used to declare the trait layout explicitly.
}
  \item{method}{
Likelihood evaluation method. The current implementation supports only dense methods: \code{"rpf"}, \code{"inverse"}, and \code{"pseudoinverse"}. These methods evaluate the same likelihood for the same model and differ only in the linear algebra used to compute it. \code{"rpf"} is the recommended default, \code{"inverse"} is mainly a slower explicit-matrix cross-check, and \code{"pseudoinverse"} can be safer for numerically awkward cases at higher computational cost.
}
  \item{scale.height}{
Whether the tree should be scaled to unit height before fitting.
}
  \item{optimization}{
Optimization routine. The special value \code{"fixed"} skips numerical optimization and returns a scaffold containing the log-likelihood function together with the resolved starting-value parameterization. This mode is intended mainly for simulation-first workflows, and it is the only mode that allows \code{data=NULL}.
}
  \item{control}{
Control list passed to the optimizer. In addition to the optimizer's own entries, \code{mvSIMMAP} recognizes \code{hessian.tol} for the post-fit Hessian tolerance, \code{retry.unreliable} to enable automatic jittered restarts when the first fit converges to an unreliable or nearly flat solution, \code{retry.max} for the number of such restarts, \code{retry.jitter} for the relative amount of start-value perturbation, and \code{retry.seed} for reproducible restart draws. These special entries are stripped before the remaining control list is sent to the optimizer, so they should be thought of as \code{mvSIMMAP}-specific controls rather than generic optimizer arguments.
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

The current root treatment is fixed-root-style: \code{mvSIMMAP} always includes one global root-state vector in the mean design and estimates that root analytically by GLS together with the OU/OUM optima. There is not yet a public analogue of the broader \code{mvOU} \code{root=} or \code{vcv=} options such as \code{"randomRoot"} or \code{"stationary"}. In the pure OU/OUM overlap cases covered by the tests, the present semantics correspond to \code{mvOU(..., root=TRUE, vcv="fixedRoot")}.

In practice, the default \code{mvSIMMAP()} fit is per-regime: each painted regime forms its own process group and therefore gets its own dynamic parameter block unless \code{process.groups} is used to tie regimes together. The public interface does not yet expose an \code{mvBM}-style partial constraint layer for sharing only correlations, only variances, shared eigenvectors, or proportional rate matrices across groups.

The \code{method} argument changes only how that same dense likelihood is evaluated. \code{"rpf"} reuses the package's factorization-based dense machinery and is the best default choice when the fit is numerically well behaved. \code{"inverse"} uses an explicit inverse/determinant calculation and is mainly useful as a slower numerical cross-check. \code{"pseudoinverse"} uses a generalized inverse and can be more forgiving near singularity, but it is typically the slowest of the three. None of these method choices change the model interpretation, parameterization, or set of fitted parameters.

The documented decomposition choices for \code{mvSIMMAP} are the values shown in the usage section: \code{decomp = c("scalarPositive", "diagonalPositive", "diagonal", "cholesky")} for OU/OUM \code{alpha} blocks and \code{decompSigma = c("cholesky", "diagonal")} for \code{sigma} blocks. The default \code{"scalarPositive"} option constrains each OU/OUM process group to a single positive pull parameter so that \code{alpha = a I}; this works for both \code{"OU"} and \code{"OUM"} groups and is useful when full traitwise OU pull strengths appear too weakly identified. \code{"diagonalPositive"} is the next richer option when per-trait pull strengths are needed, and \code{"cholesky"} should be treated as an opt-in richer structure rather than the starting point for routine multivariate fits. Other helper-level decompositions exist elsewhere in \code{mvMORPH}, but they should be treated as internal unless and until they are explicitly documented for \code{mvSIMMAP}.

A practical starting recipe for new multivariate fits is therefore \code{method = "rpf"} with the default \code{decomp = "scalarPositive"}, only moving to richer \code{alpha} structures if the scientific question or a clear AIC gain justifies it.

Compared with the older standalone fitters, \code{mvSIMMAP} intentionally exposes a smaller surface for now. It does not currently provide public \code{root}, \code{vcv}, \code{precalc}, \code{smean}, \code{trend}, \code{sparse}, \code{pic}, or \code{univarpf} options, and it does not yet expose the user-defined matrix-index constraint layer available in some \code{mvBM}/\code{mvOU} workflows.

For pure special-case overlap, the standalone mixed fitter recovers familiar models already present elsewhere in the package: pure BM/BMM assignments compare against \code{mvBM}, and pure OU/OUM assignments compare against \code{mvOU}. The main difference is that \code{mvSIMMAP} keeps the common fixed-root-style mean structure described above rather than exposing the broader root/VCV menu of \code{mvOU}.

Printed model objects show the estimated parameter matrices in detail. Calling \code{summary()} on a fitted \code{mvSIMMAP} object gives a more compact structural view listing the painted regimes, process groups, theta ownership, and the names of the estimated parameter blocks. EB beta values printed at or extremely near zero are annotated as the BM boundary for readability.

For multivariate fits, \code{summary()} also adds a compact overview of process-group variances, cross-trait correlations, and traitwise OU pull summaries. Both the full print and compact summary now report the root treatment and mean-structure layout explicitly. \code{print(fit, compact=TRUE)} is a shorthand for printing that summary view directly.

The post-fit Hessian check is slightly more informative than the older binary \code{mvMORPH} diagnostic used in some other fitters. Tiny negative eigenvalues inside \code{hessian.tol} are treated as numerical noise, nearly zero eigenvalues are reported as a flat or boundary solution, and clearly negative eigenvalues are reported as a possible saddle point. When \code{retry.unreliable = TRUE}, the function can automatically try a small number of jittered restarts and keep the best retained solution.

Fixed or optimized \code{mvSIMMAP} objects can also be used as scaffolds for \code{simulate()}. The simulator reuses the fitted tree, \code{process}, and \code{process.groups} structure from the object and, when needed, accepts a \code{param=list(...)} override with replacement values for \code{theta}, \code{sigma}, \code{alpha}, and \code{beta}. This mixed-model simulation path is currently tied to the scaffold's original tree rather than supporting a new \code{tree=} override.

For simulation-first workflows, \code{mvSIMMAP(tree, data=NULL, optimization="fixed", ...)} can be used as a lightweight scaffold constructor. Internally this creates a zero-valued placeholder dataset only to establish dimensions and default parameter blocks; the resulting object is intended mainly for \code{simulate()} rather than for likelihood-based inference. In that mode, \code{param$ntraits} or compatible matrix-valued parameter information is required so the trait dimension can be inferred unambiguously.

For EB groups, \code{beta} values at or near zero correspond to the BM boundary. The \code{low} and \code{up} entries in \code{param} only affect those EB \code{beta} parameters and are ignored for BM, OU, and OUM groups.

\code{coef(fit)} returns the estimated mean-parameter matrix, whereas \code{parameters(fit)} returns a long-form table of all fitted parameter elements with block, group, row, and column labels that are easier to use in simulation studies and post-processing scripts.
}
\value{
\item{LogLik }{The log-likelihood of the fitted model, stored as a standard \code{"logLik"} object with \code{df} and \code{nobs} attributes for compatibility with generic comparison tools such as \code{BIC()}.}
\item{AIC }{Akaike Information Criterion.}
\item{AICc }{Sample size-corrected AIC.}
\item{theta }{Estimated mean parameters. This matrix always includes a \code{root} row, plus one row per shared \code{"OU"} process group and/or one row per painted \code{"OUM"} regime, depending on \code{process} and \code{process.groups}.}
\item{alpha }{Named list of OU/OUM pull matrices by process group.}
\item{sigma }{Named list of process-group-specific rate matrices.}
\item{beta }{Named vector of EB beta values by process group.}
\item{process }{The normalized regime-to-process assignment vector used in the fit.}
\item{process.groups }{The normalized regime-to-process-group assignment vector used in the fit.}
\item{convergence }{Convergence status of the optimizing function.}
\item{hess.values }{Compatibility-friendly binary Hessian diagnostic: \code{0} for reliable and \code{1} otherwise.}
\item{diagnostics }{Detailed post-fit diagnostics, including Hessian status, eigenvalue summary, and any automatic restart information.}
\item{param }{List of fit metadata and internal parameterization functions, including the resolved root-treatment and mean-structure labels used for printing and summaries.}
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
\code{\link{parameters}}
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
    alpha = list(A = 1)
  )
)

fit_shared_ou <- mvSIMMAP(
  tree, X,
  process = c(A = "OU", B = "OU"),
  process.groups = c(A = "ou_shared", B = "ou_shared"),
  param = list(
    sigma = list(ou_shared = diag(2) * 0.1),
    alpha = list(ou_shared = 1)
  )
)

fit_shared_dynamics <- mvSIMMAP(
  tree, X,
  process = c(A = "OUM", B = "OUM"),
  process.groups = c(A = "ou_shared", B = "ou_shared"),
  param = list(
    sigma = list(ou_shared = diag(2) * 0.1),
    alpha = list(ou_shared = 1)
  )
)

scaffold <- mvSIMMAP(
  tree, data = NULL,
  process = c(A = "BM", B = "OU"),
  optimization = "fixed",
  param = list(ntraits = 2, names_traits = c("y1", "y2"))
)

Ysim <- simulate(
  scaffold, nsim = 2,
  param = list(
    theta = rbind(root = c(0, 0), theta.B = c(1, -1)),
    sigma = list(A = diag(2) * 0.05, B = diag(2) * 0.05),
    alpha = list(B = 1.5)
  )
)
}
}
\keyword{ Ornstein Uhlenbeck }
\keyword{ Brownian Motion }
\keyword{ Early-Burst }
\keyword{ SIMMAP }
\keyword{ Regimes }
