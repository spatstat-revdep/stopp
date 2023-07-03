\name{localtest}
\alias{localtest}
\title{Test of local structure for spatio-temporal point processes}
\usage{
localtest(X, Z, method = "K", k, alpha = 0.05, verbose = TRUE)
}
\arguments{
\item{X}{Background spatio-temporal point pattern. Usually, the most clustered
between \code{X} and \code{Z}. Must be either a \code{stp} or \code{stlp} object.}

\item{Z}{Other spatio-temporal point pattern. Must also be of the same class as \code{X}.}

\item{method}{Character string indicating which version of LISTA function to use:
either
 \code{"K"} or \code{"g"}.
If \code{"K"}, the local spatio-temporal K-function is used to run the test.
If \code{"g"}, the local spatio-temporal pair correlation function is used.}

\item{k}{Number of permutations}

\item{alpha}{Significance level}

\item{verbose}{If TRUE (default) the progress of the test is printed}
}
\value{
A list of class \code{localtest}, containing
\describe{
\item{\code{p}}{A vector of p-values, one for each of the points in \code{X}}
\item{\code{X}}{The background spatio-temporal point pattern given in input}
\item{\code{Z}}{The alternative spatio-temporal point pattern given in input}
\item{\code{alpha}}{The threshold given in input}
\item{\code{Xsig}}{A \code{stp} object storing the resulting significant points}
\item{\code{Xnosig}}{A \code{stp} object storing the resulting non-significant points}
\item{\code{id}}{The ids of the resulting significan points}
}
}
\description{
This function performs the permutation test of the local structure for
spatio-temporal point pattern data, proposed in Siino et al. (2018), as well as for
spatio-temporal point pattern data occurring on the same linear network,
 following
D'Angelo et al. (2021).
}
\details{
The test detects local differences between two observed spatio-temporal point patterns occurring on the same space-time region.

 If the test is  performed for spatio-temporal point patterns as in
   Siino et al. (2018), that is, on an object of class \code{stp}, the LISTA
   functions employed are the local functions of
   Adelfio et al. (2020), documented in
 \link{KLISTAhat} and \link{LISTAhat} of the \code{stpp} package (Gabriel et al, 2013).

 If the function is applied to a \code{stlp} object, that is, on two spatio-temporal
 point patterns observed on the same linear network \code{L}, the LISTA function
used are the ones proposed in D'Angelo et al. (2021), documented
in  \link{localSTLKinhom} and \link{localSTLginhom}.

Details on the test and its performance are found in Siino et al. (2018) and
D'Angelo et al. (2021), for Euclidean and network spaces, respectively.
}
\examples{



\donttest{
# background pattern
set.seed(12345)
X <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(.05, 4),
           nsim = 1, seed = 2, verbose = TRUE)

# alternative pattern
set.seed(12345)
Z <- rstpp(lambda = 25, nsim = 1, seed = 2, verbose = TRUE)


# run the local test
test <- localtest(X, Z, method = "K", k = 9)

## Example on a linear network

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



}


}
\references{
Adelfio, G., Siino, M., Mateu, J., and Rodríguez-Cortés, F. J. (2020). Some properties of local weighted second-order statistics for spatio-temporal point processes. Stochastic Environmental Research and Risk Assessment, 34(1), 149-168.

D’Angelo, N., Adelfio, G., and Mateu, J. (2021). Assessing local differences between the spatio-temporal second-order structure of two point patterns occurring on the same linear network. Spatial Statistics, 45, 100534.

Gabriel, E., Rowlingson, B. S., and Diggle, P. J. (2013). stpp: An R Package for Plotting, Simulating and Analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53(2), 1–29. https://doi.org/10.18637/jss.v053.i02

Siino, M., Rodríguez‐Cortés, F. J., Mateu, J. ,and Adelfio, G. (2018). Testing for local structure in spatiotemporal point pattern data. Environmetrics, 29(5-6), e2463.
}
\seealso{
\link{print.localtest}, \link{summary.localtest}, \link{plot.localtest}
}
\author{
Nicoletta D'Angelo and Marianna Siino
}
