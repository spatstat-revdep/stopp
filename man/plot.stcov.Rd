% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.stcov.R
\name{plot.stcov}
\alias{plot.stcov}
\title{Plot a \code{stcov} object}
\usage{
\method{plot}{stcov}(x, ...)
}
\arguments{
\item{x}{An object of class \code{stcov}}

\item{...}{additional unused argument}
}
\description{
This function plots the covariate stored in the \code{stcov} object given in input,
in a three panel plot representing the 3Dplot of the coordinates, and the
covariate values.
}
\examples{


set.seed(2)
df <- data.frame(runif(100), runif(100), runif(100), rpois(100, 15))

cov <- stcov(df, interp = FALSE)

plot(cov)

}
\seealso{
\link{stcov}
}
\author{
Nicoletta D'Angelo
}
