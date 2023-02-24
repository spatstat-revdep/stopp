\name{plot.localdiag}
\alias{plot.localdiag}
\title{Plot of the local diagnostics' result on a spatio-temporal point process model}
\usage{
\method{plot}{localdiag}(x, marg = TRUE, col = "grey", col2 = "red", cols = "lightgrey", ...)
}
\arguments{
\item{x}{A \code{localdiag} object}

\item{marg}{Default to \code{TRUE}. If \code{marg = F}, only the
spatio-temporal point pattern is plotted}

\item{col}{Color of the outlying points}

\item{col2}{Color of the network (if applicable)}

\item{cols}{Color of the non-outlying points}

\item{...}{additional unused argument}
}
\description{
This function plots the result of the local diagnostics performed with
\link{localdiag} on either a \code{stp} or \code{stlp} object.
 It highlights the  points of the analysed spatio-temporal point pattern \code{X}
 which are identified as outlying by the
 previously performed local diagnostics; the remaining points of \code{X}
 are also represented.

 It also shows the underlying linear network, if the local diagnostics has been applied
 to point patterns occurring on the same linear network, that is, if \link{localdiag}
 has been applied to a \code{stlp} object.
}
\value{No return value, called for side effects}
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

plot(resETAS)
plot(resPOIS)

plot(resETAS, marg = FALSE)
plot(resPOIS, marg = FALSE)

}


}
\references{
Adelfio, G., Siino, M., Mateu, J., and Rodríguez-Cortés, F. J. (2020). Some properties of local weighted second-order statistics for spatio-temporal point processes. Stochastic Environmental Research and Risk Assessment, 34(1), 149-168.

D’Angelo, N., Adelfio, G.  and Mateu, J. (2022) Local inhomogeneous second-order characteristics for spatio-temporal point processes on linear networks. Stat Papers. https://doi.org/10.1007/s00362-022-01338-4
}
\seealso{
\link{infl}, \link{print.localdiag}, \link{summary.localdiag}
}
\author{
Nicoletta D'Angelo and Giada Adelfio
}
