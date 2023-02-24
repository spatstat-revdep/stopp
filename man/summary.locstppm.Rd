\name{summary.locstppm}
\alias{summary.locstppm}
\title{Summary of a fitted local spatio-temporal Poisson process model}
\usage{
\method{summary}{locstppm}(object, ...)
}
\arguments{
\item{object}{An object of class \code{locstppm}}

\item{...}{additional unused argument}
}
\description{
The function summarises the main information on the distribution of the parameters
 of a fitted local spatio-temporal Poisson process model.
}
\value{No return value, called for side effects}
\examples{
\donttest{
set.seed(2)
pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6),
             nsim = 1, seed = 2, verbose = TRUE)
inh00_local <- locstppm(pin, formula = ~ 1)
inh01_local <- locstppm(pin, formula = ~ x)

summary(inh00_local)
summary(inh01_local)

}


}
\references{
D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
}
\seealso{
\link{locstppm}, \link{print.locstppm},
\link{plot.locstppm}
}
\author{
Nicoletta D'Angelo
}
