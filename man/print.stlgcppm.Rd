\name{print.stlgcppm}
\alias{print.stlgcppm}
\title{Print of a fitted LGCP model}
\usage{
\method{print}{stlgcppm}(x, ...)
}
\arguments{
\item{x}{An object of class \code{stlgcppm}}

\item{...}{additional unused argument}
}
\description{
The function prints the main information on the fitted model.
  In this case of local parameters (both first- and second-order),
  the summary function contains information on
 their distributions.
}
\value{No return value, called for side effects}
\examples{
\donttest{
# Example with complex seismic point pattern
data("greececatalog")

lgcp1 <- stlgcppm(greececatalog, formula = ~ 1, first = "global", second = "global")
lgcp2 <- stlgcppm(greececatalog, formula = ~ x, first = "local", second = "global")

lgcp1


# Joint minimum contrast fit
# for a log-Gaussian Cox process with
# global first-order intensity and
# global second-order intensity
# --------------------------------------------------
#   Homogeneous Poisson process
# with Intensity: 0.0064
#
# Estimated coefficients of the first-order intensity:
#   (Intercept)
# -5.052
# --------------------------------------------------
#   Covariance function: separable
#
# Estimated coefficients of the second-order intensity:
#   sigma  alpha   beta
# 10.984  0.224 47.076
# --------------------------------------------------
#   Model fitted in 0.873 minutes


lgcp2


# Joint minimum contrast fit
# for a log-Gaussian Cox process with
# local first-order intensity and
# global second-order intensity
# --------------------------------------------------
#   Inhomogeneous Poisson process
# with Trend: ~x
#
# Summary of estimated coefficients of the first-order intensity
# (Intercept)           x
# Min.   :-6.282   Min.   :-0.96831
# 1st Qu.:-2.387   1st Qu.:-0.36685
# Median : 2.122   Median :-0.25871
# Mean   : 2.052   Mean   :-0.26309
# 3rd Qu.: 4.569   3rd Qu.:-0.07325
# Max.   :17.638   Max.   : 0.10269
# --------------------------------------------------
#   Covariance function: separable
#
# Estimated coefficients of the second-order intensity:
#   sigma  alpha   beta
# 2.612  0.001 36.415
# --------------------------------------------------
#   Model fitted in 3.503 minutes

}


}
\references{
D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.

Siino, M., Adelfio, G., and Mateu, J. (2018). Joint second-order parameter estimation for spatio-temporal log-Gaussian Cox processes. Stochastic environmental research and risk assessment, 32(12), 3525-3539.
}
\seealso{
\link{stlgcppm}, \link{print.stlgcppm}, \link{localsummary.stlgcppm},
\link{plot.stlgcppm}, \link{localplot.stlgcppm}
}
\author{
Nicoletta D'Angelo and Giada Adelfio
}
