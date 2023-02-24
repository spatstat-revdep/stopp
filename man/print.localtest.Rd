\name{print.localtest}
\alias{print.localtest}
\title{Print of the result of the permutation local test}
\usage{
\method{print}{localtest}(x, ...)
}
\arguments{
\item{x}{An object of class  \code{localtest}}

\item{...}{additional unused argument}
}
\description{
It prints the main information on the result of the local permutation test
 performed with \link{localtest}  on either a \code{stp} or \code{stlp} object:
 whether the local test was run on point patterns lying on a linear network or not;
 the number of points in the background \code{X} and alternative \code{Z} patterns;
 the number of points in \code{X} which exhibit local differences in the second-order
 structure with respect to \code{Z}, according to the performed test.
}
\value{No return value, called for side effects}
\examples{

\dontrun{

# background pattern
set.seed(12345)
X <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(.05, 4),
           nsim = 1, seed = 2, verbose = TRUE)

# alternative pattern
set.seed(12345)
Z <- rstpp(lambda = 25, nsim = 1, seed = 2, verbose = TRUE)

# run the local test
test <- localtest(X, Z, method = "K", k = 9)

test

# Test for local differences between two
# spatio-temporal point patterns
# --------------------------------------
#   Backgound pattern X: 17
# Alternative pattern Z: 20
#
# 1 significant points at alpha = 0.05


# Linear networks

# background pattern
set.seed(12345)
X <- rETASlp(cat = NULL, params = c(0.078915 / 1.95, 0.003696,  0.013362,  1.2,
                                        0.424466,  1.164793),
                 betacov = 0.5, m0 = 2.5, b = 1.0789, tmin = 0, t.lag = 200,
                 xmin = 600, xmax = 2200, ymin = 4000, ymax = 5300,
                 iprint = TRUE, covdiag = FALSE, covsim = FALSE, L = chicagonet)

# alternative pattern, on the same linear network
l <- 20 / (spatstat.geom::volume(chicagonet) * (200 - 25))
set.seed(12345)
stlppPOIS <- stlnpp::rpoistlpp(lambda = l, a = 25, b = 200, L = chicagonet)
Z <- as.stlp(stlppPOIS)

# run the local test
test <- localtest(X, Z, method = "K", k = 9)

test

# Test for local differences between two
# spatio-temporal point patterns on a linear network
# --------------------------------------------------
#   Backgound pattern X: 31
# Alternative pattern Z: 22
#
# 19 significant points at alpha = 0.05

}


}
\references{
D’Angelo, N., Adelfio, G., and Mateu, J. (2021). Assessing local differences between the spatio-temporal second-order structure of two point patterns occurring on the same linear network. Spatial Statistics, 45, 100534.

Siino, M., Rodríguez‐Cortés, F. J., Mateu, J. ,and Adelfio, G. (2018). Testing for local structure in spatiotemporal point pattern data. Environmetrics, 29(5-6), e2463.
}
\seealso{
\link{localtest}, \link{summary.localtest},
\link{plot.localtest}
}
\author{
Nicoletta D'Angelo
}
