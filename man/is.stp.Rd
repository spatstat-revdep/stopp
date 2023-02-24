\name{is.stp}
\alias{is.stp}
\title{Test whether an object is a spatio-temporal point pattern}
\usage{
is.stp(x)
}
\arguments{
\item{x}{Any object}
}
\value{
\code{TRUE} is \code{x} is a spatio-temporal point pattern, otherwise \code{FALSE}
}
\description{
Test whether an object is a spatio-temporal point pattern
}
\details{
This function tests whether the object \code{x} is a spatio-temporal
point pattern object of class \code{stp}. See \link{stp} for details of this class
}
\examples{
set.seed(12345)
rpp1 <- stpp::rpp(lambda = 200, replace = FALSE)
is.stp(rpp1)

df0 <- cbind(rpp1$xyt[, 1], rpp1$xyt[, 2], rpp1$xyt[, 3])
stp1 <- stp(df0)

is.stp(stp1)

}
\seealso{
\link{stp}, \link{summary.stp}, \link{print.stp}, \link{plot.stp}
}
\author{
Nicoletta D'Angelo
}
