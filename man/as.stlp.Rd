\name{as.stlp}
\alias{as.stlp}
\title{Convert a \code{stlpp} object to a \code{stlp} object}
\usage{
as.stlp(x)
}
\arguments{
\item{x}{An object of class \code{stlpp}}
}
\value{
An object of class \code{stlp}
}
\description{
This function converts an object of class \code{stlpp} to an object of class \code{stlp}
}
\examples{

set.seed(12345)
stlpp1 <- stlnpp::rpoistlpp(.2, a = 0, b = 5, L = stlnpp::easynet)

class(as.stlp(stlpp1))

}
\references{
Moradi M, Cronie O, Mateu J (2020). stlnpp: Spatio-temporal analysis of point patterns on linear networks.
}
\seealso{
\link{stp}, \link{summary.stlp}, \link{print.stlp},
\link{as.stlpp}, \link{plot.stlp}
}
\author{
Nicoletta D'Angelo
}
