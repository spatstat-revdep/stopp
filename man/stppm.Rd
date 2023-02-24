\name{stppm}
\alias{stppm}
\title{Fit a Poisson process model to a spatio-temporal point pattern}
\usage{
stppm(X, formula, verbose = TRUE, mult = 4)
}
\arguments{
\item{X}{A \code{stp} object}

\item{formula}{An object of class \code{"formula"}:
a symbolic description of the model to be fitted.
The current version only supports formulas depending on the spatial and temporal coordinates:
\code{x}, \code{y}, \code{t}.}

\item{verbose}{Default to \code{TRUE}}

\item{mult}{The multiplicand of the number of data points,
 for setting the number of dummy
points to generate for the quadrature scheme}
}
\value{
An object of class \code{stppm}. A list of
\describe{
\item{\code{IntCoefs}}{The fitted coefficients}
\item{\code{X}}{The \code{stp} object provided as input}
\item{\code{nX}}{The number of points in \code{X}}
\item{\code{I}}{Vector indicating which points are dummy or data}
\item{\code{y_resp}}{The response variable of the model fitted to the quadrature scheme}
\item{\code{formula}}{The formula provided as input}
\item{\code{l}}{Fitted intensity}
\item{\code{mod_global}}{The \code{glm} object of the model fitted to the quadrature scheme}
\item{\code{newdata}}{The data used to fit the model, without the dummy points}
\item{\code{time}}{Time elapsed to fit the model, in minutes}
}
}
\description{
This function fits a Poisson process model to an observed spatio-temporal
point pattern stored in a \code{stp} object.

We assume that the template model is a Poisson process, with a parametric
intensity or rate function \eqn{\lambda(\textbf{u}, t; \theta)}  with space
and time locations \eqn{\textbf{u} \in
W,  t \in T} and parameters \eqn{\theta \in \Theta.}


Estimation is performed through the fitting of a \code{glm} using a spatio-temporal
version of the quadrature scheme by Berman and Turner (1992).

See the 'Details' section.
}
\details{
The log-likelihood of the template model is
\deqn{\log L(\theta) = \sum_i
\lambda(\textbf{u}_i, t_i; \theta) - \int_W\int_T
\lambda(\textbf{u}, t; \theta) \text{d}t\text{d}u}
up to an additive constant, where the sum is over all points \eqn{\textbf{u}_i}
 in the spatio-temporal point process \eqn{X}.

Following Berman and Turner (1992), we use a finite quadrature approximation
to the log-likelihood. Renaming the data points as \eqn{\textbf{x}_1,\dots ,
\textbf{x}_n} with \eqn{(\textbf{u}_i,t_i) = \textbf{x}_i} for \eqn{i = 1, \dots , n},
then generate \eqn{m}  additional 'dummy points' \eqn{(\textbf{u}_{n+1},t_{n+1})
\dots , (\textbf{u}_{m+n},t_{m+n})} to
form a set of \eqn{n + m} quadrature points (where \eqn{m > n}).

Then we determine quadrature weights \eqn{a_1, \dots , a_m}
so that integrals in the log-likelihood can be approximated by a Riemann sum
\deqn{ \int_W \int_T \lambda(\textbf{u},t;\theta)\text{d}t\text{d}u \approx \sum_{k = 1}^{n + m}a_k\lambda(\textbf{u}_{k},t_{k};\theta)}
where \eqn{a_k} are the quadrature weights such that
\eqn{\sum_{k = 1}^{n + m}a_k = l(W \times T)} where \eqn{l} is the Lebesgue measure.

Then the log-likelihood  of the template model can be approximated by
\deqn{ \log L(\theta)   \approx \sum_i \log \lambda(\textbf{x}_i; \theta) +\sum_j(1 - \lambda(\textbf{u}_j,t_j; \theta))a_j=\sum_je_j \log \lambda(\textbf{u}_j, t_j; \theta) + (1 - \lambda(\textbf{u}_j, t_j; \theta))a_j}
where \eqn{e_j = 1\{j \leq n\}} is the indicator that equals \eqn{1} if
\eqn{u_j} is a data point. Writing \eqn{y_j = e_j/a_j} this becomes
\deqn{
\log L(\theta) \approx
\sum_j
a_j
(y_j \log \lambda(\textbf{u}_j, t_j; \theta) - \lambda(\textbf{u}_j, t_j; \theta))
+
  \sum_j
a_j.}

Apart from the constant \eqn{\sum_j a_j}, this expression is formally equivalent
 to the weighted log-likelihood of
a Poisson regression model with responses \eqn{y_j} and means
\eqn{\lambda(\textbf{u}_j,t_j; \theta) = \exp(\theta Z(\textbf{u}_j,t_j) +
 B(\textbf{u}_j,t_j))}.

  This is
maximised by this function by using standard GLM software.


In detail, we define the spatio-temporal quadrature scheme by considering a
spatio-temporal
partition of \eqn{W \times T} into cubes \eqn{C_k} of equal volume \eqn{\nu},
 assigning the weight \eqn{a_k=\nu/n_k}
 to each quadrature point (dummy or data) where \eqn{n_k} is the number of
 points that lie in the same cube as the point \eqn{u_k} (Raeisi et al, 2021).

The number of dummy points should be sufficient for an accurate estimate of the
likelihood. Following Baddeley et al. (2000) and Raeisi et al. (2021),
we start with a number of dummy points \eqn{m \approx 4 n}, increasing it until
\eqn{\sum_k a_k = l(W \times T)}.
}
\examples{

\donttest{
## Homogeneous
set.seed(2)
ph <- rstpp(lambda = 200, nsim = 1, seed = 2, verbose = TRUE)
hom1 <- stppm(ph, formula = ~ 1)

## Inhomogeneous
pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6),
nsim = 1, seed = 2, verbose = TRUE)
inh1 <- stppm(pin, formula = ~ x)
}

}
\references{
Baddeley, A. J., Møller, J., and Waagepetersen, R. (2000). Non-and semi-parametric estimation of interaction in inhomogeneous point patterns. Statistica Neerlandica, 54(3):329–350

Berman, M. and Turner, T. R. (1992). Approximating point process likelihoods with glim. Journal of the Royal Statistical Society: Series C (Applied Statistics), 41(1):31–38

D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.

Raeisi, M., Bonneu, F., and Gabriel, E. (2021). A spatio-temporal multi-scale model for geyer saturation point process: application to forest fire occurrences. Spatial Statistics, 41:100492.
}
\seealso{
\link{locstppm}, \link{AIC.stppm}, \link{BIC.stppm}
}
\author{
Nicoletta D'Angelo
}
