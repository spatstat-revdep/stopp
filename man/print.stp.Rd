\name{print.stp}
\alias{print.stp}
\title{Print a \code{stp} object}
\usage{
\method{print}{stp}(x, ...)
}
\arguments{
\item{x}{An object of class \code{stp}}

\item{...}{additional unused argument}
}
\value{No return value, called for side effects}
\description{
It prints the main information on the spatio-temporal point pattern
stored in the \code{stp} object: the number of points; the enclosing spatial window;
the temporal time period.
}
\examples{

set.seed(12345)
rpp1 <- stpp::rpp(lambda = 200, replace = FALSE)
df0 <- cbind(rpp1$xyt[, 1], rpp1$xyt[, 2], rpp1$xyt[, 3])
stp1 <- stp(df0)

stp1

# Spatio-temporal point pattern
# 208 points
# Enclosing window: rectangle = [0.0011366, 0.9933775] x [0.0155277,
#                                                          0.9960438] units
# Time period: [0.004, 0.997]

}
\seealso{
\link{stp}, \link{summary.stp}, \link{plot.stp},
\link{as.stp}, \link{as.stpp}
}
\author{
Nicoletta D'Angelo
}
