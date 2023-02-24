\name{stlgcppm}
\alias{stlgcppm}
\title{Fit a log-Gaussian Cox process model to a spatio-temporal point pattern}
\usage{
stlgcppm(
  X,
  formula,
  verbose = TRUE,
  cov = "separable",
  first = "local",
  second = "local",
  mult = 4,
  hs = "global",
  npx0 = 10,
  npt0 = 10,
  itnmax = 100
)
}
\arguments{
\item{X}{A stp object}

\item{formula}{An object of class \code{formula}: a symbolic description of the first-order intensity to be fitted.
The current version only supports formulas depending on the spatial and temporal coordinates:
\code{x}, \code{y}, \code{t}.}

\item{verbose}{Default to TRUE}

\item{cov}{Covariance function to be fitted for the second-order intensity function.
Default to \code{separable}. Other options are \code{gneiting} and \code{iaco-cesare}".}

\item{first}{Character string indicating whether to fit a first-order intensity function
with global or local parameters:
either \code{local} (default) or \code{global}.}

\item{second}{Character string indicating whether to fit a second-order intensity function
with global or local parameters:
either \code{local} (default) or \code{global}.}

\item{mult}{The multiplicand of the number of data points,
 for setting the number of dummy
points to generate for the quadrature scheme}

\item{hs}{Character string indicating whether to select fixed or variable bandwidths
for the kernel weights to be used in the log-likelihood.
In any of those cases, the well-supported rule-of-thumb for choosing the
bandwidth of a Gaussian kernel density estimator is employed.
If \code{hs = "global"} (default), a fixed bandwidth is selected.
If \code{hs = "local"}, an individual badwidth is selected for each point in the
pattern  \code{X}.}

\item{npx0}{A positive integer representing the spatial distance to np-th closest event.
Used in the computation of the local bandwidth.
Suitable values are in the range from 10 (default) to 100.}

\item{npt0}{A positive integer representing the temporal distance to np-th closest event.
Used in the computation of the local bandwidth.
Suitable values are in the range from 10 (default) to 100.}

\item{itnmax}{Maximum number of iterations to run in the optimization procedure
for the estiamtion of the second-order intensity parameters.}
}
\value{
A list of the class \code{stlgcppm}, containing
\describe{
\item{\code{IntCoefs}}{The fitted coefficients of the first-order intensity function}
\item{\code{CovCoefs}}{The fitted coefficients of the second-order intensity function}
\item{\code{X}}{The stp objectt provided as input}
\item{\code{formula}}{The formula provided as input}
\item{\code{cov}}{A string with the chosen covariance type}
\item{\code{l}}{Fitted first-order intensity}
\item{\code{mu}}{Mean function of the random intensity}
\item{\code{mod_global}}{The glm object of the model fitted to the quadrature scheme
for the first-order intensity parameters estimation}
\item{\code{newdata}}{The data used to fit the model, without the dummy points}
\item{\code{time}}{Time elapsed to fit the model, in minutes}
}
}
\description{
This function estimates a local log-Gaussian Cox process (LGCP), following the
*locally weighted minimum contrast* procedure introduced in
D'Angelo et al. (2023).

Three covariances are available: separable exponential,
Gneiting, and De Iaco-Cesare.

If both \code{first} and \code{second} arguments are set to \code{global}, a log-Gaussian
Cox process is fitted by means of the *joint minimum contrast* procedure proposed in
Siino et al. (2018).
}
\details{
Following the inhomogeneous specification in Diggle et al. (2013),
we consider LGCPs with intensity
\deqn{
\Lambda(\textbf{u},t)=\lambda(\textbf{u},t)\exp(S(\textbf{u},t)).
}

Following Siino et al. (2018), the second-order
parameters \eqn{\boldsymbol{\psi}} are found by minimising
\deqn{  M_J\{ \boldsymbol{\psi}\}=\int_{h_0}^{h_{max}} \int_{r_0}^{r_{max}} \phi(r,h) \{\nu[\hat{J}(r,h)]-\nu[J(r,h;\boldsymbol{\psi})]\}^2 \text{d}r \text{d}h,
} where \eqn{\phi(r, h)} is a weight that depends on the space-time
distance and \eqn{\nu} is a transformation function.

They suggest \eqn{\phi(r,h)=1} and \eqn{\nu} as
the identity function, while \eqn{r_{max}} and \eqn{h_{max}} are selected as 1/4
of the maximum observable spatial and temporal distances.

Following D'Angelo et al. (2023), we can fit a localised version of the LGCP,
that is,  obtain a
vector of parameters \eqn{\boldsymbol{\psi}_i} for each point \eqn{i}, by
minimising \deqn{
  M_{J,i}\{ \boldsymbol{\psi}_i \}=\int_{h_0}^{h_{max}}\int_{r_0}^{r_{max}}
   \phi(r,h) \{ \nu[\bar{J}_i(r,h)]-\nu[J(r,h;\boldsymbol{\psi})]\}^2 \text{d}r \text{d}h,
} where \eqn{ \bar{J}_i(r,h)= \frac{\sum_{i=1}^{n}\hat{J}_i(r,h)w_i}{\sum_{i=1}^{n}w_i}} is the average of the local functions
\eqn{\hat{J}_i(r,h)}, weighted by some point-wise kernel estimates.

  In particular, we consider \eqn{\hat{J}_i(\cdot)} as the local
spatio-temporal pair correlation function (Gabriel et al, 2013) documented in \link{LISTAhat}.
}
\examples{

\donttest{
# Example with complex seismic point pattern
data("greececatalog")

##

# If both first and second arguments are set to "global", a log-Gaussian
# Cox process is fitted by means of the joint minimum contrast.

lgcp1 <- stlgcppm(greececatalog, formula = ~ 1, first = "global", second = "global")

 ##

 # If first = "local", local parameters for the first-order intensity are
 # provided. In this case, the summary function contains information on
 # their distributions.

 lgcp2 <- stlgcppm(greececatalog, formula = ~ x, first = "local", second = "global")

 ##

# Finally, if second = "local" (default option), the model with local
 #parameters for the covariance is fitted.

 lgcp3 <- stlgcppm(greececatalog, formula = ~ 1, first = "global", second = "local")

 ##

#  It is also possible to fit local parameters for both the first and
#  second-order intensity, but we abstain from addressing this case here.


}


}
\references{
Baddeley, A. (2017). Local composite likelihood for spatial point processes. Spatial Statistics, 22, 261-295.

D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.

Diggle, P. J., Moraga, P., Rowlingson, B., and Taylor, B. M. (2013). Spatial and spatio-temporal log-gaussian cox processes: extending the geostatistical paradigm. Statistical Science, 28(4):542–563.

Gabriel, E., Rowlingson, B. S., and Diggle, P. J. (2013). stpp: An R Package for Plotting, Simulating and Analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53(2), 1–29. https://doi.org/10.18637/jss.v053.i02

Siino, M., Adelfio, G., and Mateu, J. (2018). Joint second-order parameter estimation for spatio-temporal log-Gaussian Cox processes. Stochastic environmental research and risk assessment, 32(12), 3525-3539.
}
\seealso{
\link{print.stlgcppm}, \link{summary.stlgcppm}, \link{localsummary.stlgcppm},
 \link{plot.stlgcppm}, \link{localplot.stlgcppm}
}
\author{
Nicoletta D'Angelo, Giada Adelfio, and Marianna Siino
}
