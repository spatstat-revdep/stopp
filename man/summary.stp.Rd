\name{summary.stp}
\alias{summary.stp}
\title{Summary of a \code{stp} object}
\usage{
\method{summary}{stp}(object, ...)
}
\arguments{
\item{object}{An object of class \code{stp}}

\item{...}{additional unused argument}
}
\description{
It prints the summary statistics of the spatial and temporal coordinates
of  the spatio-temporal point pattern
stored in the \code{stp} object.
}
\value{No return value, called for side effects}
\examples{

set.seed(12345)
rpp1 <- stpp::rpp(lambda = 200, replace = FALSE)
df0 <- cbind(rpp1$xyt[, 1], rpp1$xyt[, 2], rpp1$xyt[, 3])
stp1 <- stp(df0)

summary(stp1)

#       x                  y                 t
# Min.   :0.001137   Min.   :0.01553   Min.   :0.004254
# 1st Qu.:0.290510   1st Qu.:0.27316   1st Qu.:0.266398
# Median :0.559250   Median :0.53679   Median :0.524548
# Mean   :0.536354   Mean   :0.52031   Mean   :0.513976
# 3rd Qu.:0.793498   3rd Qu.:0.78442   3rd Qu.:0.758390
# Max.   :0.993378   Max.   :0.99604   Max.   :0.996996

}
\seealso{
\link{stp}, \link{print.stp}, \link{plot.stp},
\link{as.stp}, \link{as.stpp}
}
\author{
Nicoletta D'Angelo
}
