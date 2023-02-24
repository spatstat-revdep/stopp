\name{localdiag}
\alias{localdiag}
\title{Local diagnostics of spatio-temporal point process models}
\usage{
localdiag(X, intensity, p = 0.95)
}
\arguments{
\item{X}{Either a \code{stp} or a \code{stlp} object}

\item{intensity}{A vector of intensity values, of the same length as the number
of point in \code{X}}

\item{p}{The percentile to consider as threshold for the outlying points.
Default to 0.95.}
}
\value{
A list object of class \code{localdiag}, containing
\describe{
\item{\code{X}}{The \code{stp} object provided as input}
\item{\code{listas}}{The LISTA functions, in a list object}
\item{\code{ids}}{The ids of the points identified as outlying}
\item{\code{x2}}{A vector with the individual contributions to the Chi-squared statistics,
 normalized}
\item{\code{p}}{The percentile considered}
}
}
\description{
This function performs local diagnostics of a model fitted for the
first-order intensity of a spatio-temporal point pattern, by means of  the
local spatio-temporal inhomogeneous K-function (Adelfio et al, 2020)
 documented by the function
 \link[stpp]{KLISTAhat}  of the \code{stpp} package (Gabriel et al, 2013).

The function can also perform local diagnostics of a model fitted for the
first-order intensity of an  spatio-temporal point pattern on a linear network,
 by means of  the
local spatio-temporal inhomogeneous K-function on linear networks
 (D'Angelo et al, 2021)
 documented by the function
 \link{localSTLKinhom}.

In both cases, it returns the   points identified as outlying following the diagnostics
procedure on individual points of an observed
point pattern, as introduced in
Adelfio et al. (2020), and applied in D'Angelo et al. (2022) for the linear
network case

See the section 'Details'.

The points resulting from the local diagnostic procedure provided by this
function can be inspected via the \link{plot}, \link{print}, \link{summary},
and \link{infl}
functions.
}
\details{
Adelfio et al. (2020) derived the expectation of the local
inhomogeneous spatio-temporal K-function, under the Poisson case:
\eqn{\mathbb{E}[\hat{K}^i(r,h) ]= \pi r^ 2 h.}

Moreover, they found that when the local estimator is weighted by the true
 intensity function,
 its
expectation, \eqn{\mathbb{E}[\hat{K}_{I}^i(r,h)]}, is the same as the expectation of
\eqn{\hat{K}^i(r,h)}.

These results motivate the usage of such local estimator
\eqn{\hat{K}_{I}^i(r,h)} as a diagnostic tool for general spatio-temporal
point processes for assessing the goodness-of-fit of spatio-temporal
point processes of any generic first-order
intensity function \eqn{\lambda}.

Indeed, if the estimated intensity function
used for weighting in our proposed LISTA functions is
the true one, then the LISTA functions should behave as the
corresponding ones of a homogeneous Poisson process,
resulting in small discrepancies between the two.

Therefore, this function computes such discrepancies
by means of the \eqn{\chi_i^2} values, obtained following the expression
\deqn{
  \chi_i^2=\int_L \int_T \Bigg(
    \frac{\big(\hat{K}^i_{I}(r,h)- \mathbb{E}[\hat{K}^i(r,h) ]
    \big)^2}{\mathbb{E}[\hat{K}^i(r,h) ]}
    \Bigg) \text{d}h \text{d}r ,
  }
    one for each point in the point pattern.

Basically, departures of the LISTA functions \eqn{\hat{K}_{I}^i(r,h)} from
the Poisson expected value \eqn{rh} directly suggest the unsuitability of
the intensity function \eqn{\lambda(\cdot)} used in the weighting of the
LISTA functions for that specific point. This can be referred to as  an *outlying point*.


Given that D'Angelo et al. (2022) proved the same results for the network case,
that is,
\eqn{\mathbb{E}[\hat{K}_{L}^i(r,h) ]= rh} and
 \eqn{\mathbb{E}[\hat{K}_{L,I}^i(r,h) ]=\mathbb{E}[\hat{K}_{L}^i(r,h) ] }
when \eqn{\hat{K}_{L,I}^i(r,h) } is weighted by the true intensity function,
we implemented the same above-mentioned diagnostics procedure to work on
intensity functions fitted on spatio-temporal point patterns occurring on
linear networks.

Note that the Euclidean procedure is implemented by the
local K-functions of
   Adelfio et al. (2020), documented in
 \link{KLISTAhat} of the \code{stpp} package (Gabriel et al, 2013).
The network case uses  the local K-functions on networks (D'Angelo et al., 2021),
 documented
in  \link{localSTLKinhom}.
}
\examples{

\donttest{

#load data
set.seed(12345)
id <- sample(1:nrow(etasFLP::catalog.withcov), 200)
cat <- etasFLP::catalog.withcov[id, ]
stp1 <- stp(cat[, 5:3])

#fit two competitor models
# and extract the fitted spatio-temporal intensity

lETAS <- etasFLP::etasclass(cat.orig = cat, magn.threshold = 2.5, magn.threshold.back = 3.9,
  mu = 0.3, k0 = 0.02, c = 0.015, p = 1.01, gamma = 0, d = 1,
  q = 1.5, params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE,
  TRUE), formula1 = "time ~  magnitude- 1",
  declustering = TRUE,
  thinning = FALSE, flp = TRUE, ndeclust = 15, onlytime = FALSE,
  is.backconstant = FALSE, sectoday = FALSE, usenlm = TRUE,
  compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36)$l

lPOIS <- etasFLP::etasclass(cat.orig = cat, magn.threshold = 2.5, magn.threshold.back = 3.9,
  mu = 0.3, k0 = 0.02, c = 0.015, p = 1.01, gamma = 0, d = 1,
  q = 1.5, params.ind = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
  FALSE), formula1 = "time ~  magnitude- 1",
  declustering = TRUE,
  thinning = FALSE, flp = TRUE, ndeclust = 15, onlytime = FALSE,
  is.backconstant = FALSE, sectoday = FALSE, usenlm = TRUE,
  compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36)$l

# let's identify the outlying points at a .9 percentile

resETAS <- localdiag(stp1, lETAS, p = .9)
resPOIS <- localdiag(stp1, lPOIS, p = .9)


#### Network case

set.seed(12345)
stlp1 <- rETASlp(cat = NULL, params = c(0.078915 / 2, 0.003696,  0.013362,  1.2,
                                             0.424466,  1.164793),
                      betacov = 0.5, m0 = 2.5, b = 1.0789, tmin = 0, t.lag = 200,
                      xmin = 600, xmax = 2200, ymin = 4000, ymax = 5300,
                      iprint = TRUE, covdiag = FALSE, covsim = FALSE, L = chicagonet)

res <- localdiag(stlp1, intensity = density(as.stlpp(stlp1),
                                                at = "points"), p = .65)


}


}
\references{
Adelfio, G., Siino, M., Mateu, J., and Rodríguez-Cortés, F. J. (2020). Some properties of local weighted second-order statistics for spatio-temporal point processes. Stochastic Environmental Research and Risk Assessment, 34(1), 149-168.

D’Angelo, N., Adelfio, G.  and Mateu, J. (2022) Local inhomogeneous second-order characteristics for spatio-temporal point processes on linear networks. Stat Papers. https://doi.org/10.1007/s00362-022-01338-4

Gabriel, E., Rowlingson, B. S., and Diggle, P. J. (2013). stpp: An R Package for Plotting, Simulating and Analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53(2), 1–29. https://doi.org/10.18637/jss.v053.i02
}
\seealso{
\link{infl}, \link{plot.localdiag}, \link{print.localdiag},
\link{summary.localdiag},
\link{globaldiag}
}
\author{
Nicoletta D'Angelo and Giada Adelfio
}
