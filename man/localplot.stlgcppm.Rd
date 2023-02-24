\name{localplot.stlgcppm}
\alias{localplot.stlgcppm}
\title{Plot the coefficients of a fitted local LGCP model}
\usage{
localplot.stlgcppm(x, par = TRUE)
}
\arguments{
\item{x}{An object of class \code{stlgcppm}}

\item{par}{Default to \code{TRUE}.}
}
\description{
The function plots the local estimates.
In the case of local covariance parameters, the function
displays the local estimates of the chosen covariance function.
}
\value{No return value, called for side effects}
\examples{
\donttest{
# Example with complex seismic point pattern
data("greececatalog")

lgcp2 <- stlgcppm(greececatalog, formula = ~ x, first = "local", second = "global")
localplot.stlgcppm(lgcp2)

}




}
\references{
D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
}
\seealso{
\link{stlgcppm}, \link{print.stlgcppm}, \link{summary.stlgcppm},
\link{localsummary.stlgcppm}, \link{plot.stlgcppm}
}
\author{
Nicoletta D'Angelo and Giada Adelfio
}
