\name{print.stlp}
\alias{print.stlp}
\title{Print a \code{stlp} object}
\usage{
\method{print}{stlp}(x, ...)
}
\arguments{
\item{x}{An object of class \code{stlp}}

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

\donttest{

set.seed(12345)
stlpp1 <- stlnpp::rpoistlpp(.2, a = 0, b = 5, L = stlnpp::easynet)
df0 <- cbind(stlpp1$data$x, stlpp1$data$y, stlpp1$data$t)
L0 <- stlpp1$domain
stlp1 <- stp(df0, L0)

stlp1

# Spatio-temporal point pattern on a linear network
# 43 points
# Linear network with 19 vertices and 26 lines
# Enclosing window: rectangle = [-0.01, 5.1] x [-0.01, 5.1] units
# Time period: [0.043, 4.93]

}

}
\seealso{
\link{stp}, \link{plot.stlp}, \link{summary.stlp},
\link{as.stlpp}, \link{as.stlp}
}
\author{
Nicoletta D'Angelo
}
