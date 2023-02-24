\name{localSTLKinhom}
\alias{localSTLKinhom}
\title{Local inhomogeneous Spatio-temporal K-functions on a linear network}
\usage{
localSTLKinhom(
  x,
  lambda = lambda,
  normalize = FALSE,
  r = NULL,
  t = NULL,
  nxy = 10
)
}
\arguments{
\item{x}{A realisation of a spatio-temporal point processes on a linear network in \code{stlp} format}

\item{lambda}{values of estimated intensity.}

\item{normalize}{normalization factor to be considered.}

\item{r}{values of argument r where K-function will be evaluated. optional.}

\item{t}{values of argument t where K-function will be evaluated. optional.}

\item{nxy}{pixel array dimensions. optional.}
}
\value{
A list of objects of class \code{sumstlpp} (Moradi and Mateu, 2020).
}
\description{
The functions \code{localSTLKinhom} and \code{localSTLginhom} implement the
inhomogeneous LISTA functions proposed in D'Angelo et al. (2022).

The homogeneous K-function and pair correlation functions, in
D'Angelo et al. (2021), can be obtained easily with \code{localSTLKinhom} and
\code{localSTLginhom}, by imputing a lambda vector of constant intensity
values, the same for each point.

These local functions are the local counterparts of \link{STLKinhom} and
\link{STLginhom} by
 Moradi and Mateu (2020), available in the \code{stlnpp} package (Moradi et al., 2020).
Basically, we avoid summing up all
the points as in the global statistics counterparts, and we denote the
individual contribution to the global statistics with the index \eqn{i}.
}
\details{
The *local spatio-temporal inhomogeneous*
  K-function for the i-th event \eqn{(\boldsymbol{u}_i,t_i)} on a linear network
is \deqn{\hat{K}^i_{L,I}(r,h)=\frac{1}{ \vert L  \vert  \vert T  \vert}\sum_{(\boldsymbol{u}_i,t_i)\ne (\boldsymbol{v},s)} \frac{I\{ d_L(\boldsymbol{u}_i,\boldsymbol{v})<r,\vert t_i-s\vert <h\} }{\hat{\lambda}(\boldsymbol{u}_i,t_i)\hat{\lambda}(\boldsymbol{v},s)M((\boldsymbol{u}_i,t_i),d_L(\boldsymbol{u}_i,\boldsymbol{v}),\vert t_i-s\vert )},}
 and the corresponding *local pair correlation function* (pcf)
  \deqn{\hat{g}^i_{L,I}(r,h)=\frac{1}{ \vert L  \vert  \vert T  \vert}\sum_{(\boldsymbol{u}_i,t_i)\ne (\boldsymbol{v},s)} \frac{\kappa( d_L(\boldsymbol{u}_i,\boldsymbol{v})-r)\kappa(\vert t_i-s\vert -h) }{\hat{\lambda}(\boldsymbol{u}_i,t_i)\hat{\lambda}(\boldsymbol{v},s)M((\boldsymbol{u}_i,t_i),d_L(\boldsymbol{u}_i,\boldsymbol{v}),\vert t_i-s\vert )},}
with
  \deqn{D(X) = \frac{n-1}{ \vert L  \vert  \vert T  \vert}\sum_{i=1}^n\sum_{i \ne j}}
  normalization factor. This leads to the unbieased estimators \eqn{\frac{1}{D(X)}\hat{K}^i_{L,I}(r,h)} and
\eqn{\frac{1}{D(X)}\hat{g}^i_{L,I}(r,h)} (D'Angelo et al., 2022) .

The homogeneous versions (D'Angelo et al., 2021) can be obtained by weighting the second-order
summary statistics (either K or pcf) by a constant intensity
\eqn{\hat{\lambda}=n/( \vert L \vert \vert T \vert)}, giving
  \deqn{\hat{K}_L^i(r,h)=\frac{1}{\hat{\lambda}^{2} \vert L  \vert  \vert T  \vert}\sum_{(\boldsymbol{u}_i,t_i)\ne (\boldsymbol{v},s)} \frac{I\{ d_L(\boldsymbol{u}_i,\boldsymbol{v})<r,\vert t_i-s\vert <h\} }{M((\boldsymbol{u}_i,t_i),d_L(\boldsymbol{u}_i,\boldsymbol{v}),\vert t_i-s\vert )},}
and
 \deqn{\hat{g}_L^i(r,h)=\frac{1}{\hat{\lambda}^{2} \vert L  \vert  \vert T  \vert}\sum_{(\boldsymbol{u}_i,t_i)\ne (\boldsymbol{v},s)} \frac{\kappa( d_L(\boldsymbol{u}_i,\boldsymbol{v})-r)\kappa(\vert t_i-s\vert -h) }{M((\boldsymbol{u}_i,t_i),d_L(\boldsymbol{u}_i,\boldsymbol{v}),\vert t_i-s\vert )}.}
}
\examples{

\donttest{
set.seed(10)
X <- stlnpp::rpoistlpp(.2, a = 0, b = 5, L = stlnpp::easynet)
X
lambda <- density(X, at = "points")
x <- as.stlp(X)
k <- localSTLKinhom(x, lambda = lambda, normalize = TRUE)
#select an individual point
j = 1
k[[j]]
#plot the lista function and compare it with its theoretical value
inhom <- list(x = k[[j]]$r, y = k[[j]]$t, z = k[[j]]$Kinhom)
theo <- list(x = k[[j]]$r, y = k[[j]]$t, z = k[[j]]$Ktheo)
diff <- list(x = k[[j]]$r, y = k[[j]]$t, z = k[[j]]$Kinhom - k[[j]]$Ktheo)
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(1, 3))
fields::image.plot(inhom, main= "Kinhom", col = hcl.colors(12, "YlOrRd", rev = FALSE),
xlab = "Spatial distance", ylab = "Temporal distance")
fields::image.plot(theo, main = "Ktheo", col = hcl.colors(12, "YlOrRd", rev = FALSE),
xlab = "Spatial distance", ylab = "Temporal distance")
fields::image.plot(diff, main = "Kinhom - Ktheo", col = hcl.colors(12, "YlOrRd", rev = FALSE),
xlab = "Spatial distance", ylab = "Temporal distance")
par(oldpar)
}


}
\references{
D’Angelo, N., Adelfio, G., and Mateu, J. (2021). Assessing local differences between the spatio-temporal second-order structure of two point patterns occurring on the same linear network. Spatial Statistics, 45, 100534.

D’Angelo, N., Adelfio, G.,  and Mateu, J. (2022). Local inhomogeneous second-order characteristics for spatio-temporal point processes on linear networks. Stat Papers. https://doi.org/10.1007/s00362-022-01338-4

Moradi M, Cronie O, and Mateu J (2020). stlnpp: Spatio-temporal analysis of point patterns on linear networks.

Moradi, M. M., and Mateu, J. (2020). First-and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics, 29(3), 432-443.
}
\seealso{
\link{localSTLginhom}, \link{STLKinhom}, \link{STLginhom}
}
\author{
Nicoletta D'Angelo
}
