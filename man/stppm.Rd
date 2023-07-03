\name{stppm}
\alias{stppm}
\title{Fit a Poisson process model to a spatio-temporal point pattern}
\usage{
stppm(X, formula, verbose = TRUE, mult = 4)
}
\arguments{
\item{X}{A \code{stp} object}

\item{formula}{An object of class \code{"formula"}:
a symbolic description of the model to be fitted.
The current version only supports formulas depending on the spatial and temporal coordinates:
\code{x}, \code{y}, \code{t}.}

\item{verbose}{Default to \code{TRUE}}

\item{mult}{The multiplicand of the number of data points,
 for setting the number of dummy
points to generate for the quadrature scheme}
}
\value{
An object of class \code{stppm}. A list of
\describe{
\item{\code{IntCoefs}}{The fitted coefficients}
\item{\code{X}}{The \code{stp} object provided as input}
\item{\code{nX}}{The number of points in \code{X}}
\item{\code{I}}{Vector indicating which points are dummy or data}
\item{\code{y_resp}}{The response variable of the model fitted to the quadrature scheme}
\item{\code{formula}}{The formula provided as input}
\item{\code{l}}{Fitted intensity}
\item{\code{mod_global}}{The \code{glm} object of the model fitted to the quadrature scheme}
\item{\code{newdata}}{The data used to fit the model, without the dummy points}
\item{\code{time}}{Time elapsed to fit the model, in minutes}
}
}
\description{
This function fits a Poisson process model to an observed spatio-temporal
point pattern stored in a \code{stp} object.

Estimation is performed through the fitting of a \code{glm} using a spatio-temporal
version of the quadrature scheme by Berman and Turner (1992).

}

\examples{

\donttest{
## Homogeneous
set.seed(2)
ph <- rstpp(lambda = 200, nsim = 1, seed = 2, verbose = TRUE)
hom1 <- stppm(ph, formula = ~ 1)

## Inhomogeneous
pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6),
nsim = 1, seed = 2, verbose = TRUE)
inh1 <- stppm(pin, formula = ~ x)
}

}
\references{
Berman, M. and Turner, T. R. (1992). Approximating point process likelihoods with glim. Journal of the Royal Statistical Society: Series C (Applied Statistics), 41(1):31–38

D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.

}
\seealso{
\link{locstppm}, \link{AIC.stppm}, \link{BIC.stppm}
}
\author{
Nicoletta D'Angelo
}
