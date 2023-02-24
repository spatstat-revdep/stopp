\name{plot.stp}
\alias{plot.stp}
\title{Plot a \code{stp} object}
\usage{
\method{plot}{stp}(x, tcum = FALSE, marg = TRUE, col = 1, ...)
}
\arguments{
\item{x}{An object of class \code{stp}}

\item{tcum}{Default to \code{FALSE}. If \code{TRUE},
the temporal point pattern is displayed cumulatively.
A barplot is automatically plotted if there are repeated
counts (typically with discrete times).}

\item{marg}{Default to \code{TRUE}. If \code{FALSE}, only the
spatio-temporal point pattern is plotted.}

\item{col}{The color of the points. Default to \code{"black"}}

\item{...}{additional unused argument}
}
\description{
This function plots the point pattern stored in the stp object given in input,
in a three panel plot representing the 3Dplot of the coordinates, and the
marginal spatial and temporal coordinates.
}
\value{No return value, called for side effects}
\examples{

set.seed(12345)
rpp1 <- stpp::rpp(lambda = 200, replace = FALSE)
df0 <- cbind(rpp1$xyt[, 1], rpp1$xyt[, 2], rpp1$xyt[, 3])
stp1 <- stp(df0)

#plot
plot(stp1)

#cumulative time occurrances
plot(stp1, tcum = TRUE)

#change color of points
plot(stp1, col = "blue")

#display only in space-time
plot(stp1, marg = FALSE)

#discrete times
stp2 <- stp(cbind(rpp1$xyt[, 1], rpp1$xyt[, 2], round(rpp1$xyt[, 3] * 10)))
plot(stp2)

}
\seealso{
\link{stp}, \link{print.stp}, \link{summary.stp},
\link{as.stp}, \link{as.stpp}
}
\author{
Nicoletta D'Angelo
}
