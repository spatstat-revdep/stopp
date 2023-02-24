\name{infl}
\alias{infl}
\title{Display outlying LISTA functions}
\usage{
infl(x, id = NULL)
}
\arguments{
\item{x}{An object of class \code{localdiag}}

\item{id}{The id of the LISTA to display.
Default is set to the ids identified and stored in the \code{localdiag} object}
}
\description{
This function works on the objects of class \code{localdiag}, as returned by
\link{localdiag}, plotting the identified 'outlying'
LISTA functions. These correspond to the influential points in the fitting
of the model provided by \link{localdiag}
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

# let's look at the outlying LISTA functions selected by localdiag() ...

infl(resETAS)
infl(resPOIS)

#... and at a some specific LISTA

infl(resETAS, id = c(75))
infl(resETAS, id = c(2, 4, 6))
infl(resPOIS, id = 1:6)

}


}
\references{
Adelfio, G., Siino, M., Mateu, J., and Rodríguez-Cortés, F. J. (2020). Some properties of local weighted second-order statistics for spatio-temporal point processes. Stochastic Environmental Research and Risk Assessment, 34(1), 149-168.

D’Angelo, N., Adelfio, G.  and Mateu, J. (2022) Local inhomogeneous second-order characteristics for spatio-temporal point processes on linear networks. Stat Papers. https://doi.org/10.1007/s00362-022-01338-4
}
\seealso{
\link{localdiag}, \link{plot.localdiag},
 \link{print.localdiag}, \link{summary.localdiag}
}
\author{
Nicoletta D'Angelo and Giada Adelfio
}
