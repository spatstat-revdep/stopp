\name{plot.stlp}
\alias{plot.stlp}
\title{Plot a \code{stlp} object}
\usage{
\method{plot}{stlp}(x, tcum = FALSE, marg = TRUE, col = 1, cols = "grey", ...)
}
\arguments{
\item{x}{An object of class \code{stp}}

\item{tcum}{Default to \code{FALSE}.
If \code{TRUE}, the temporal point pattern is displayed cumulatively.
A barplot is automatically plotted if there are repeated
counts (typically with discrete times).}

\item{marg}{Default to \code{TRUE}. If \code{FALSE}, only the
spatio-temporal point pattern is plotted.}

\item{col}{The color of the points. Default to \code{"black"}}

\item{cols}{The color of the linear network. Default to \code{"grey"}}

\item{...}{additional unused argument}
}
\description{
This function plots the point pattern on a linear network
 stored in the \code{stlp} object given in input,
in a three panel plot representing the \link{plot3D} of the coordinates, and the
marginal spatial and temporal coordinates.
}
\value{No return value, called for side effects}
\examples{

set.seed(12345)
stlpp1 <- stlnpp::rpoistlpp(.2, a = 0, b = 5, L = stlnpp::easynet)
df0 <- cbind(stlpp1$data$x, stlpp1$data$y, stlpp1$data$t)
L0 <- stlpp1$domain
stlp1 <- stp(df0, L0)

#plot
plot(stlp1)

#cumulative time occurrances
plot(stlp1, tcum = TRUE)

#change color of points
plot(stlp1, col = "blue")

#change color of network
plot(stlp1, cols = "magenta")

#display only in space-time
plot(stlp1, marg = FALSE)

#discrete times
stp2 <- stp(cbind(stlpp1$data$x, stlpp1$data$y, round(stlpp1$data$t)))
plot(stp2)

}
\seealso{
\link{stp}, \link{summary.stlp}, \link{print.stlp},
\link{as.stlpp}, \link{as.stlp}
}
\author{
Nicoletta D'Angelo
}
