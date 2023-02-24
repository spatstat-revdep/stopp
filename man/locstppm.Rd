\name{locstppm}
\alias{locstppm}
\title{Fit a local Poisson process model to a spatio-temporal point pattern}
\usage{
locstppm(
  X,
  formula,
  verbose = TRUE,
  mult = 4,
  hs = "global",
  npx0 = 10,
  npt0 = 10
)
}
\arguments{
\item{X}{A \code{stp} object}

\item{formula}{An object of class \code{"formula"}: a symbolic description of
 the model to be fitted.
The current version only supports formulas depending on the spatial and temporal coordinates:
\code{x}, \code{y}, \code{t}.}

\item{verbose}{Default to \code{TRUE}}

\item{mult}{The multiplicand of the number of data points,
 for setting the number of dummy
points to generate for the quadrature scheme}

\item{hs}{Character string indicating whether to select fixed or variable bandwidths
for the kernel weights to be used in the log-likelihood.
In any of those cases, the well-supported rule-of-thumb for choosing the
bandwidth of a Gaussian kernel density estimator is employed.
If \code{hs = "global"} (default), a fixed bandwidth is selected.
If \code{hs = "local"}, an individual bandwidth is selected for each point in the
pattern  \code{X}.}

\item{npx0}{Number of lags for the space grid period for variable bandwidths kernel}

\item{npt0}{Number of lags for the time period for variable bandwidths kernel}
}
\value{
An object of class \code{locstppm}. A list of
\describe{
\item{\code{IntCoefs}}{The fitted global coefficients}
\item{\code{IntCoefs_local}}{The fitted local coefficients}
\item{\code{X}}{The \code{stp} objectt provided as input}
\item{\code{nX}}{The number of points in \code{X}}
\item{\code{I}}{Vector indicating which points are dummy or data}
\item{\code{y_resp}}{The response variable of the model fitted to the quadrature scheme}
\item{\code{formula}}{The formula provided as input}
\item{\code{l}}{Fitted intensity through the global parameters}
\item{\code{l_local}}{Fitted intensity through the local parameters}
\item{\code{mod_global}}{The \code{glm} object of the model fitted to the quadrature scheme}
\item{\code{newdata}}{The data used to fit the model, without the dummy points}
\item{\code{time}}{Time elapsed to fit the model, in minutes}
}
}
\description{
This function fits a Poisson process model to an observed spatio-temporal
point pattern stored in a \code{stp} object, that is, a Poisson model with
a set of parameters \eqn{\theta_i} for each point \eqn{i}.

We assume that the template model is a Poisson process, with a parametric
intensity or rate function \eqn{\lambda(\textbf{u}, t; \theta_i)}  with space
and time locations \eqn{\textbf{u} \in
W,  t \in T} and parameters \eqn{\theta_i \in \Theta.}

Estimation is performed through the fitting of a \code{glm} using a localized
version of the quadrature scheme by Berman and Turner (1992), firstly introduced
in the purely spatial context by Baddeley (2017), and in the spatio-temporal
framework by D'Angelo et al. (2023).

See the 'Details' section.
}
\details{
The local log-likelihood associated with the  spatio-temporal location
\eqn{(\textbf{v},s)} is given by
\deqn{
\log L((\textbf{v},s);\theta) = \sum_i w_{\sigma_s}(\textbf{u}_i - \textbf{v}) w_{\sigma_t}(t_i - s)
\lambda(\textbf{u}_i, t_i; \theta)  - \int_W \int_T
\lambda(\textbf{u}, t; \theta) w_{\sigma_s}(\textbf{u}_i - \textbf{v}) w_{\sigma_t}(t_i - s) \text{d}t \text{d}u
}
where \eqn{w_{\sigma_s}} and \eqn{w_{\sigma_t}} are weight functions, and
\eqn{\sigma_s, \sigma_t > 0} are the smoothing bandwidths. It is not
necessary to assume that \eqn{w_{\sigma_s}} and \eqn{w_{\sigma_t}}
 are probability densities. For simplicity, we shall consider only kernels of fixed
bandwidth, even though spatially adaptive kernels could also be used.

Note that if the template model is the homogeneous Poisson process with intensity
\eqn{\lambda}, then the local
likelihood estimate \eqn{\hat{\lambda}(\textbf{v}, s)}
reduces to the kernel estimator of the point process intensity with
kernel proportional to \eqn{w_{\sigma_s}w_{\sigma_t}}.

We now use an  approximation similar to
\eqn{
\log L(\theta) \approx
\sum_j
a_j
(y_j \log \lambda(\textbf{u}_j, t_j; \theta) - \lambda(\textbf{u}_j, t_j; \theta))
+
  \sum_j
a_j,}
but for the local log-likelihood associated
with each desired location \eqn{(\textbf{v},s) \in W \times T}, that is:
  \deqn{
\log L((\textbf{v},s); \theta) \approx
\sum_j
w_j(\textbf{v},s)a_j
(y_j \log \lambda(\textbf{u}_j,t_j; \theta) - \lambda(\textbf{u}_j,t_j; \theta))
+
  \sum_j
w_j(\textbf{v},s)a_j   ,}
where \eqn{w_j(\textbf{v},s) = w_{\sigma_s}(\textbf{v} - \textbf{u}_j)
w_{\sigma_t}(s - t_j)}.

Basically, for each
desired location \eqn{(\textbf{v},s)},
 we replace the vector of quadrature weights \eqn{a_j} by
 \eqn{a_j(\textbf{v},s)= w_j(\textbf{v},s)a_j} where
\eqn{w_j (\textbf{v},s) = w_{\sigma_s}(\textbf{v} - \textbf{u}_j)w_{\sigma_t}(s - t_j)},
 and use the GLM software to fit the Poisson regression.

The local likelihood is defined at any location \eqn{(\textbf{v},s)} in continuous space.
In practice it is sufficient to
consider a grid of points \eqn{(\textbf{v},s)}.

We refer to D'Angelo et al. (2023) for further discussion on bandwidth selection
and on computational costs.
}
\examples{

\donttest{
set.seed(2)
pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6),
nsim = 1, seed = 2, verbose = TRUE)
inh00_local <- locstppm(pin, formula = ~ 1)
inh01_local <- locstppm(pin, formula = ~ x)
}


}
\references{
Baddeley, A. (2017). Local composite likelihood for spatial point processes. Spatial Statistics, 22, 261-295.

D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
}
\seealso{
\link{stppm}
}
\author{
Nicoletta D'Angelo
}
