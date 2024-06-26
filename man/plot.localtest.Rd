% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.localtest.R
\name{plot.localtest}
\alias{plot.localtest}
\title{Plot of the result of the local permutation test}
\usage{
\method{plot}{localtest}(x, col = "grey", cols = "lightgrey", col2 = "red", ...)
}
\arguments{
\item{x}{An object of class  \code{localtest}}

\item{col}{Color of the significant points}

\item{cols}{Color of the linear network. If applicable.}

\item{col2}{Color of the non-significant points}

\item{...}{additional unused argument}
}
\description{
This function plots the result of the local permutation test performed with
\link{localtest} on either a \code{stp} or \code{stlp} object. It highlights the points of the background pattern \code{X},
 which exhibit local differences in the second-order
 structure with respect to \code{Z}, according to the previously performed test.
 The remaining points of \code{X} are also represented.

 It also shows the underlying linear network, if the local test has been applied
 to point patterns occurring on the same linear network, that is, if \link{localtest}
 has been applied to a \code{stlp} object.
}
\examples{

\donttest{

set.seed(2)
X <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)},
            par = c(.005, 5))
Z <- rstpp(lambda = 30)

test <- localtest(X, Z, method = "K", k = 3) 

plot(test)

}


}
\references{
D’Angelo, N., Adelfio, G., and Mateu, J. (2021). Assessing local differences between the spatio-temporal second-order structure of two point patterns occurring on the same linear network. Spatial Statistics, 45, 100534.

Siino, M., Rodríguez‐Cortés, F. J., Mateu, J. ,and Adelfio, G. (2018). Testing for local structure in spatiotemporal point pattern data. Environmetrics, 29(5-6), e2463.
}
\seealso{
\link{localtest},  \link{print.localtest},
\link{summary.localtest}
}
\author{
Nicoletta D'Angelo
}
