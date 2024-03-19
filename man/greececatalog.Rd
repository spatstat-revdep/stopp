\docType{data}
\name{greececatalog}
\alias{greececatalog}
\title{Catalog of Greek earthquakes}
\format{
A \code{stp} object for a spatio-temporal point pattern with 1111 points
}
\usage{
data(greececatalog)
}
\description{
A dataset in \code{stp} format containing the catalog of Greek earthquakes
of magnitude at least 4.0 from year 2005 to year 2014. 
Data come from the Hellenic Unified Seismic Network (H.U.S.N.).
}
\details{
The variables are as follows:
\itemize{
  \item x. longitude, ranging from 20.02 to 27.98
  \item y. latitude, ranging from 33.75 to 40.45
  \item t. time, ranging from 38354, 42000
}
}
\examples{

data(greececatalog)

plot(greececatalog)

}
\references{
D’Angelo, N., Siino, M., D’Alessandro, A., and Adelfio, G. (2022). Local spatial log-Gaussian Cox processes for seismic data. AStA Advances in Statistical Analysis, 1-39.

D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.

Gabriel, E., Rodriguez-Cortes, F., Coville, J., Mateu, J., and Chadoeuf, J. (2022). Mapping the intensity function of a non-stationary point process in unobserved areas. Stochastic Environmental Research and Risk Assessment, 1-17.

Siino, M., Adelfio, G., Mateu, J., Chiodi, M., and D’alessandro, A. (2017). Spatial pattern analysis using hybrid models: an application to the Hellenic seismicity. Stochastic Environmental Research and Risk Assessment, 31(7), 1633-1648.
}
\author{
Nicoletta D'Angelo
}
\keyword{datasets}
