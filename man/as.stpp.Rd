\name{as.stpp}
\alias{as.stpp}
\title{Convert a \code{stp} object to a \code{stpp} object}
\usage{
as.stpp(x)
}
\arguments{
\item{x}{An object of class \code{stp}}
}
\value{
An object of class \code{stpp}
}
\description{
This function converts an object of class \code{stp} to an object of class \code{stpp}
}
\examples{

set.seed(12345)
rpp1 <- stpp::rpp(lambda = 200, replace = FALSE)
df0 <- cbind(rpp1$xyt[, 1], rpp1$xyt[, 2], rpp1$xyt[, 3])
stp1 <- stp(df0)

class(as.stpp(stp1))

}
\references{
Gabriel, E., Rowlingson, B. S., & Diggle, P. J. (2013). stpp: An R Package for Plotting, Simulating and Analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53(2), 1–29. https://doi.org/10.18637/jss.v053.i02
}
\seealso{
\link{stp}, \link{print.stp}, \link{summary.stp},
\link{plot.stp}, \link{as.stp}
}
\author{
Nicoletta D'Angelo
}
