\name{as.stp}
\alias{as.stp}
\title{Convert a \code{stpp} object to a \code{stp} object}
\usage{
as.stp(x)
}
\arguments{
\item{x}{An object of class \code{stpp}}
}
\value{
An object of class \code{stp}
}
\description{
This function converts an object of class \code{stpp} to an object of class \code{stp}
}
\examples{

set.seed(12345)
rpp1 <- stpp::rpp(lambda = 200, replace = FALSE)

class(as.stp(rpp1$xyt))

}
\references{
Gabriel, E., Rowlingson, B. S., & Diggle, P. J. (2013). stpp: An R Package for Plotting, Simulating and Analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53(2), 1–29. https://doi.org/10.18637/jss.v053.i02
}
\seealso{
\link{stp}, \link{print.stp}, \link{summary.stp},
\link{plot.stp}, \link{as.stpp}
}
\author{
Nicoletta D'Angelo
}
