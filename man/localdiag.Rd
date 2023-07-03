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
The Euclidean procedure is implemented by the
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
