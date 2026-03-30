\name{parameters}
\alias{parameters}
\alias{parameters.mvmorph.mixed}
\alias{coef.mvmorph.mixed}
\title{
Extract coefficients and parameter tables from experimental mvSIMMAP fits
}
\usage{
parameters(object, ...)

\method{parameters}{mvmorph.mixed}(object, ...)

\method{coef}{mvmorph.mixed}(object, type = c("theta", "all"), ...)
}
\arguments{
  \item{object}{
A fitted or fixed \code{"mvmorph.mixed"} object returned by \code{\link{mvSIMMAP}}.
}
  \item{type}{
For \code{coef.mvmorph.mixed}, either \code{"theta"} to return the mean-parameter matrix or \code{"all"} to return the same long-form parameter table as \code{parameters()}.
}
  \item{...}{
Reserved for future extensions.
}
}
\value{
\code{coef(object)} returns the labeled mean-parameter matrix stored in \code{object$theta}. \code{parameters(object)} returns a data frame with one row per scalar parameter element and columns \code{block}, \code{group}, \code{row}, \code{col}, \code{estimate}, and \code{parameter}.
}
\description{
Convenience extractors for experimental mixed-regime \code{mvSIMMAP} fits.
}
\details{
The \code{parameters()} table is designed to make multivariate fits easier to inspect and summarize programmatically. Matrix-valued blocks such as \code{sigma} and \code{alpha} are flattened into one row per element while preserving process-group and trait labels.
}
\seealso{
\code{\link{mvSIMMAP}}
}
