\name{coef.stppm}
\alias{coef.stppm}
\title{Extract the fitted coefficients of a spatio-temporal Poisson process model}
\usage{
\method{coef}{stppm}(object, ...)
}
\arguments{
\item{object}{An object of class \code{stppm}}

\item{...}{additional unused argument}
}
\value{
A list containing the fitted coefficients
}
\description{
Extract the fitted coefficients of a spatio-temporal Poisson process model
}
\examples{

\donttest{
## Homogeneous
set.seed(2)
ph <- rstpp(lambda = 200, nsim = 1, seed = 2, verbose = TRUE)
hom1 <- stppm(ph, formula = ~ 1)

coef(hom1)

## Inhomogeneous
pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6),
nsim = 1, seed = 2, verbose = TRUE)
inh1 <- stppm(pin, formula = ~ x)

coef(inh1)

}

}
\references{
D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
}
\seealso{
\link{stppm}, \link{print.stppm}, \link{summary.stppm}, \link{plot.stppm}
}
\author{
Nicoletta D'Angelo
}
