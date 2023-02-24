\name{coef.locstppm}
\alias{coef.locstppm}
\title{Extract the fitted coefficients of a local spatio-temporal Poisson process model}
\usage{
\method{coef}{locstppm}(object, ...)
}
\arguments{
\item{object}{An object of class \code{locstppm}}

\item{...}{additional unused argument}
}
\value{
A list containing the global and local fitted coefficients
}
\description{
Extract the fitted coefficients of a local spatio-temporal Poisson process model
}
\examples{

\donttest{
set.seed(2)
pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6),
             nsim = 1, seed = 2, verbose = TRUE)
inh00_local <- locstppm(pin, formula = ~ 1)
inh01_local <- locstppm(pin, formula = ~ x)

coef(inh00_local)
coef(inh01_local)

}

}
\references{
D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
}
\seealso{
\link{locstppm}, \link{print.locstppm}, \link{summary.locstppm}, \link{plot.locstppm}
}
\author{
Nicoletta D'Angelo
}
