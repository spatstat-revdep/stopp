% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.stlpm.R
\name{summary.stlpm}
\alias{summary.stlpm}
\title{Summary of a \code{stlpm} object}
\usage{
\method{summary}{stlpm}(object, ...)
}
\arguments{
\item{object}{An object of class \code{stlpm}}

\item{...}{additional unused argument}
}
\description{
It prints the summary statistics of the spatio-temporal coordinates and the marks
of  the spatio-temporal point pattern
stored in the \code{stlpm} object.
}
\examples{

set.seed(2)
df <- data.frame(x = runif(100, 0, 0.8), y = runif(100, 0, 0.8),
 t = runif(100), m = rpois(100, 15))

stlpm1 <- stpm(df, L = chicagonet)

summary(stlpm1)

}
\author{
Nicoletta D'Angelo
}
