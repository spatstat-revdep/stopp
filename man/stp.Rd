\name{stp}
\alias{stp}
\title{Create \code{stp} and \code{stlp} objects for point patterns storage}
\usage{
stp(df, L)
}
\arguments{
\item{df}{A matrix with three columns, containing to two space and the
temporal coordinates}

\item{L}{Optional. The linear network of class \code{linnet}.
If provided, the function returns a \code{stlp} object.}
}
\value{
An \code{stp} or \code{stlpp} object, depending on whether or not an object of class
linnet is provided for the \code{L} argument.
}
\description{
This function creates a \code{stp} object as a dataframe with three columns:
\code{x}, \code{y}, and \code{t}.
If also the linear network \code{L}, of class linnet, is provided, a \code{stlp}
object is created instead.
}
\examples{

# stp

set.seed(12345)
rpp1 <- stpp::rpp(lambda = 200, replace = FALSE)
df0 <- cbind(rpp1$xyt[, 1], rpp1$xyt[, 2], rpp1$xyt[, 3])
stp1 <- stp(df0)



# stlp

set.seed(12345)
stlpp1 <- stlnpp::rpoistlpp(.2, a = 0, b = 5, L = stlnpp::easynet)
df0 <- cbind(stlpp1$data$x, stlpp1$data$y, stlpp1$data$t)
L0 <- stlpp1$domain
stlp1 <- stp(df0, L0)

}
\seealso{
\link{summary.stp}, \link{print.stp}, \link{plot.stp}

\link{stppm}, \link{AIC.stppm}, \link{BIC.stppm},
\link{print.stp}, \link{summary.stp}, \link{plot.stp},
\link{as.stp}, \link{as.stpp}, \link{print.stlp},
\link{summary.stlp}, \link{plot.stlp}, \link{as.stlp},  \link{as.stlpp}
}
\author{
Nicoletta D'Angelo
}
