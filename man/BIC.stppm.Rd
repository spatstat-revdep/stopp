\name{BIC.stppm}
\alias{BIC.stppm}
\title{Bayesian Information Criterion for spatio-temporal Poisson Point Process Model}
\usage{
\method{BIC}{stppm}(object, ...)
}
\arguments{
\item{object}{An abject of class \code{stppm}}

\item{...}{additional unused argument}
}
\value{
BIC value
}
\description{
This function returns the \eqn{BIC = k\log{n} - 2 \log(\hat{L})} of a point process
model fitted through the
function \link{stppm} applied to an observed
spatio-temporal point pattern of class \code{stp}.

As the model returned by \link{stppm} is fitted through a quadrature scheme,
the log-likelihood is computed trough the quantity:
\eqn{- \log{L(\hat{\theta}; \boldsymbol{x})} = \frac{D}{2} + \sum_{j = 1}^{n}I_j\log{w_j}+n(\boldsymbol{x})}
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


BIC(hom1)
BIC(inh1)


}


}
\references{
Baddeley, A., Rubak, E., and Turner, R. (2015). Spatial point patterns: methodology and applications with R. CRC press.
}
\seealso{
\link{stppm}, \link{AIC.stppm}
}
\author{
Nicoletta D'Angelo
}
