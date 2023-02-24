\name{globaldiag}
\alias{globaldiag}
\title{Global diagnostics of a spatio-temporal point process first-order intensity}
\usage{
globaldiag(X, intensity, samescale = TRUE)
}
\arguments{
\item{X}{A \code{stp} object}

\item{intensity}{A vector of intensity values, of the same length as the number
of point in \code{X}}

\item{samescale}{Logical value. It indicates whether to plot the observed
and the theoretical K-function in the same or
different scale. Default to \code{TRUE}.}
}
\value{Character returning the sum of squared differences.
The function plots three panels: the observed K-function, as returned by \link{STLKinhom};
the theoretical one; their difference. The function also prints the sum of
 squared differences between the observed and theoretical
K-function on the console.
}
\description{
This function performs global diagnostics of a model fitted for the
first-order intensity of a  spatio-temporal point pattern, by returning
the plots of the inhomogeneous K-function weighted by the
provided intensity to diagnose, its theoretical value,
and their difference.

If applied to a \code{stp} object, it resorts to  the
spatio-temporal inhomogeneous K-function (Gabriel and Diggle, 2009)
 documented by the function
 \link{STIKhat}  of the \code{stpp} package (Gabriel et al, 2013).

 If applied to a \code{stlp} object, it uses the
spatio-temporal inhomogeneous K-function on a linear network (Moradi and Mateu, 2020)
 documented by the function
 \link{STLKinhom}  of the \code{stlnpp} package (Moradi et al., 2020).
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

globaldiag(stp1, lETAS)


globaldiag(stp1, lPOIS)


globaldiag(stp1, lPOIS, samescale = FALSE)

## Network case

set.seed(12345)
stlp1 <- rETASlp(cat = NULL, params = c(0.078915 / 2, 0.003696,  0.013362,  1.2,
                                             0.424466,  1.164793),
                      betacov = 0.5, m0 = 2.5, b = 1.0789, tmin = 0, t.lag = 200,
                      xmin = 600, xmax = 2200, ymin = 4000, ymax = 5300,
                      iprint = TRUE, covdiag = FALSE, covsim = FALSE, L = chicagonet)

globaldiag(stlp1, intensity = density(as.stlpp(stlp1),
at = "points"))

globaldiag(stlp1, intensity = density(as.stlpp(stlp1),
                                          at = "points"), samescale = FALSE)




}



}
\references{
Adelfio, G., Siino, M., Mateu, J., and Rodríguez-Cortés, F. J. (2020). Some properties of local weighted second-order statistics for spatio-temporal point processes. Stochastic Environmental Research and Risk Assessment, 34(1), 149-168.

D’Angelo, N., Adelfio, G.  and Mateu, J. (2022) Local inhomogeneous second-order characteristics for spatio-temporal point processes on linear networks. Stat Papers. https://doi.org/10.1007/s00362-022-01338-4

Gabriel, E., and Diggle, P. J. (2009). Second‐order analysis of inhomogeneous spatio‐temporal point process data. Statistica Neerlandica, 63(1), 43-51.

Gabriel, E., Rowlingson, B. S., & Diggle, P. J. (2013). stpp: An R Package for Plotting, Simulating and Analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53(2), 1–29. https://doi.org/10.18637/jss.v053.i02

Moradi M, Cronie O, and Mateu J (2020). stlnpp: Spatio-temporal analysis of point patterns on linear networks.

Moradi, M. M., and Mateu, J. (2020). First-and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics, 29(3), 432-443.
}
\seealso{
\link{infl}, \link{plot.localdiag}, \link{print.localdiag},
\link{summary.localdiag}, \link{localdiag},
}
\author{
Nicoletta D'Angelo and Giada Adelfio
}
