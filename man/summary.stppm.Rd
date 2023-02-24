\name{summary.stppm}
\alias{summary.stppm}
\title{Summary of a fitted spatio-temporal Poisson process model}
\usage{
\method{summary}{stppm}(object, ...)
}
\arguments{
\item{object}{An object of class \code{stppm}}

\item{...}{additional unused argument}
}
\description{
The function summarises the main information of the fitted model.
}
\value{No return value, called for side effects}
\examples{
\donttest{
## Homogeneous
set.seed(2)
ph <- rstpp(lambda = 200, nsim = 1, seed = 2, verbose = TRUE)
hom1 <- stppm(ph, formula = ~ 1)

summary(hom1)

# Homogeneous Poisson process
# with Intensity: 202.093
#
# Estimated coefficients:
#   (Intercept)
# 5.309

## Inhomogeneous
pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6),
nsim = 1, seed = 2, verbose = TRUE)
inh1 <- stppm(pin, formula = ~ x)

summary(inh1)

# Inhomogeneous Poisson process
# with Trend: ~x
#
# Estimated coefficients:
#   (Intercept)           x
# 2.180       5.783
#
}


}
\references{
D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
}
\seealso{
\link{stppm}, \link{print.stppm},
\link{plot.stppm}, \link{coef.stppm}
}
\author{
Nicoletta D'Angelo
}
