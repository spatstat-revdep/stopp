\name{AIC.stppm}
\alias{AIC.stppm}
\title{Akaike's Information Criterion for spatio-temporal Poisson Point Process Model}
\usage{
\method{AIC}{stppm}(object, ..., k = 2)
}
\arguments{
\item{object}{An abject of class \code{stppm}}

\item{...}{additional unused argument}

\item{k}{numeric, the penalty per parameter to be used; the default is k = 2}
}
\value{
AIC value
}
\description{
This function returns the \eqn{AIC = 2k - 2 \log(\hat{L})} of a point process
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

AIC(hom1)
AIC(inh1)

}


}
\references{
Baddeley, A., Rubak, E., and Turner, R. (2015). Spatial point patterns: methodology and applications with R. CRC press.
}
\seealso{
\link{stppm}, \link{BIC.stppm}
}
\author{
Nicoletta D'Angelo
}
