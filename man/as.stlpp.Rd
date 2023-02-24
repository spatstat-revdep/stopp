\name{as.stlpp}
\alias{as.stlpp}
\title{Convert a \code{stlp} object to a \code{stlpp} object}
\usage{
as.stlpp(x)
}
\arguments{
\item{x}{An object of class \code{stlp}}
}
\value{
An object of class \code{stlpp}
}
\description{
This function converts an object of class \code{stlp} to an object of class \code{stlpp}
}
\examples{

\donttest{

set.seed(12345)
rpp1 <- stpp::rpp(lambda = 200, replace = FALSE)
df0 <- cbind(rpp1$xyt[, 1], rpp1$xyt[, 2], rpp1$xyt[, 3])
stlp1 <- stp(df0, chicagonet)

class(as.stlpp(stlp1))

}

}
\references{
Moradi M, Cronie O, Mateu J (2020). stlnpp: Spatio-temporal analysis of point patterns on linear networks.
}
\seealso{
\link{stp}, \link{summary.stlp}, \link{print.stlp},
\link{as.stlp}, \link{plot.stlp}
}
\author{
Nicoletta D'Angelo
}
