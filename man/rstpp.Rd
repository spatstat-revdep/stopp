\name{rstpp}
\alias{rstpp}
\title{Simulate homogeneous and inhomogeneous spatio-temporal Poisson point patterns}
\usage{
rstpp(
  lambda = 500,
  nsim = 1,
  seed = 2,
  verbose = TRUE,
  par = NULL,
  minX = 0,
  maxX = 1,
  minY = 0,
  maxY = 1,
  minT = 0,
  maxT = 1
)
}
\arguments{
\item{lambda}{Expected number of points to simulate}

\item{nsim}{Number of patterns to simulate. Default to 1.}

\item{seed}{Seed to set, if ones wished to reproduce the analyses}

\item{verbose}{Default to \code{TRUE}}

\item{par}{Parameters of the reference intensity}

\item{minX}{Minimum of x coordinate range}

\item{maxX}{Maximum of x coordinate range}

\item{minY}{Minimum of y coordinate range}

\item{maxY}{Maximum of y coordinate range}

\item{minT}{Minimum of t coordinate range}

\item{maxT}{Maximum of t coordinate range}
}
\value{
A \code{stp} object
}
\description{
This function creates a \code{stp} object, simulating a  spatio-temporal point pattern
 following either an
homogeneous or inhomogeneous intensity
}
\examples{

# homogeneous Poisson processes
h1 <- rstpp(lambda = 500, nsim = 1, seed = 2, verbose = TRUE)

h2 <- rstpp(lambda = 500, nsim = 1, seed = 2, verbose = TRUE, minX = 0,
             maxX = 2, minY = 3, maxY = 5, minT = 1, maxT = 9)

h3 <- rstpp(lambda = 900, nsim = 3, seed = 2, verbose = TRUE)

# inhomogeneous Poisson process
inh <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6),
             nsim = 1, seed = 2, verbose = TRUE)

}
\seealso{
\link{stppm}, \link{AIC.stppm}, \link{BIC.stppm}
}
\author{
Nicoletta D'Angelo
}
