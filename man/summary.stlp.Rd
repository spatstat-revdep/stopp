% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.stlp.R
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
\examples{


set.seed(2)                       
df_net <- data.frame(cbind(runif(100, 0, 0.85), runif(100, 0, 0.85), runif(100)))

stlp1 <- stp(df_net, L = chicagonet)
summary(stlp1)


}
\seealso{
\link{stp}, \link{plot.stlp}, \link{print.stlp}
}
\author{
Nicoletta D'Angelo
}
