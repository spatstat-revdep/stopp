\name{summary.stlp}
\alias{summary.stlp}
\title{Summary of a \code{stlp} object}
\usage{
\method{summary}{stlp}(object, ...)
}
\arguments{
\item{object}{An object of class \code{stlp}}

\item{...}{additional unused argument}
}
\description{
It prints the main information on the spatio-temporal point pattern on a linear
network stored in the \code{stlp} object: the number of points;
vertices and lines of the linear network; the enclosing spatial window;
the temporal time period.
}
\value{No return value, called for side effects}
\examples{

set.seed(12345)
stlpp1 <- stlnpp::rpoistlpp(.2, a = 0, b = 5, L = stlnpp::easynet)
df0 <- cbind(stlpp1$data$x, stlpp1$data$y, stlpp1$data$t)
L0 <- stlpp1$domain
stlp1 <- stp(df0, L0)

summary(stlp1)

}
\seealso{
\link{stp}, \link{plot.stlp}, \link{print.stlp},
\link{as.stlpp}, \link{as.stlp}
}
\author{
Nicoletta D'Angelo
}
