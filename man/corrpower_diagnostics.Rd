\name{corrpower_diagnostics}
\alias{corrpower_diagnostics}
\alias{corrpower_diagnostics.mvgls}
\alias{confint.mvgls}
\title{
Diagnostics For Experimental Corr-Power mvgls Fits
}
\description{
User-facing diagnostics for experimental \code{mvgls} fits with \code{model="BMM"} and \code{bmm.structure="corrpower"}. This is the correlation-power sibling of the earlier experimental corr-strength interface.
}
\usage{
corrpower_diagnostics(object, ...)

\method{corrpower_diagnostics}{mvgls}(object, level = 0.95, nboot = 100L,
  nbcores = 1L, profile_points = 9L, anchors = NULL,
  include_bootstrap_draws = FALSE, include_anchor_fits = FALSE, ...)

\method{confint}{mvgls}(object, parm = NULL, level = 0.95,
  method = c("profile", "bootstrap", "both"),
  nboot = 100L, nbcores = 1L, profile_points = 9L, ...)
}
\arguments{
\item{object}{
A fitted \code{mvgls} object with \code{model="BMM"} and \code{bmm.structure="corrpower"}.
}

\item{parm}{
Optional character vector naming the regime-specific raw parameters to report in \code{confint()}. Valid names are of the form \code{"<regime>.scale"} and \code{"<regime>.corr_power"}.
}

\item{level}{
Confidence level used for the profile and bootstrap intervals.
}

\item{method}{
For \code{confint()}, whether to return profile intervals, bootstrap intervals, or both in a single data.frame.
}

\item{nboot}{
Number of parametric-bootstrap replicates used for the bootstrap summaries.
}

\item{nbcores}{
Number of cores used for the parametric-bootstrap refits.
}

\item{profile_points}{
Number of grid points used for the one-dimensional profile calculations.
}

\item{anchors}{
Optional character vector of regimes to use for anchor-sensitivity refits. The default (\code{NULL}) profiles all fitted regimes as anchors. Use \code{FALSE} to skip anchor-sensitivity refits.
}

\item{include_bootstrap_draws}{
Logical; whether to include the raw bootstrap draws in the returned diagnostics object.
}

\item{include_anchor_fits}{
Logical; whether to include the fitted anchor-specific \code{mvgls} objects in the returned diagnostics object.
}

\item{...}{
Unused, for compatibility.
}
}
\details{
\code{corrpower_diagnostics()} summarizes the current experimental corr-power fit along three axes:

\enumerate{
\item raw regime parameters (\code{scale} and \code{corr_power}) with profile and/or bootstrap summaries,
\item derived regime summaries such as mean rate, variance, and correlation based on \code{sigma$regime},
\item anchor-sensitivity summaries obtained by refitting the same response matrix under alternative \code{bmm.reference} choices.
}

The returned object has class \code{"corrpower_diagnostics"} and contains data.frames named \code{parameter_summary}, \code{regime_summary}, \code{anchor_summary}, and \code{anchor_pairwise_cov}, plus a small \code{acceptance} list summarizing the seeded anchor-sensitivity checks.

\code{confint()} is a lighter-weight front end returning only the regime-level raw parameter intervals.
}
\value{
\code{corrpower_diagnostics()} returns an object of class \code{"corrpower_diagnostics"}.

\code{confint()} returns either a two-column matrix with \code{lower}/\code{upper} bounds, or a data.frame with \code{estimate}, profile bounds, and bootstrap bounds when \code{method="both"}.
}
\seealso{
\code{\link{mvgls}}
\code{\link{EIC}}
\code{\link{GIC}}
\code{\link{corrstrength_diagnostics}}
}
\examples{
\donttest{
tree <- phytools::pbtree(n = 20)
states <- setNames(sample(c("A", "B"), size = length(tree$tip.label), replace = TRUE), tree$tip.label)
simmap <- suppressMessages(phytools::make.simmap(tree, states, model = "ER", nsim = 1))
Y <- mvSIM(simmap, model = "BMM", nsim = 1,
           param = list(sigma = list(A = diag(2), B = 1.5 * diag(2)),
                        theta = c(0, 0), ntraits = 2))
if (is.list(Y)) Y <- Y[[1]]
fit_corrpower <- mvgls(Y ~ 1, data = list(Y = Y), tree = simmap, model = "BMM",
                       method = "LL", REML = FALSE,
                       bmm.structure = "corrpower")

\dontrun{
diag_fit <- corrpower_diagnostics(fit_corrpower, nboot = 20, profile_points = 7)
print(diag_fit)
confint(fit_corrpower, method = "both", nboot = 20, profile_points = 7)
}
}
}
