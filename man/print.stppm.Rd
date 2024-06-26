% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.stppm.R
\name{print.stppm}
\alias{print.stppm}
\title{Print of a fitted spatio-temporal Poisson process model}
\usage{
\method{print}{stppm}(x, ...)
}
\arguments{
\item{x}{An object of class \code{stppm}}

\item{...}{additional unused argument}
}
\description{
The function prints the main information of the fitted model.
}
\examples{

set.seed(2)
pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6))
inh1 <- stppm(pin, formula = ~ x)

inh1



}
\references{
D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
}
\seealso{
\link{stppm}, \link{print.stppm},
\link{plot.stppm}
}
\author{
Nicoletta D'Angelo
}
