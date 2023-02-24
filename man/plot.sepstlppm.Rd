\name{plot.sepstlppm}
\alias{plot.sepstlppm}
\title{Plot of the fitted intensity of a separable spatio-temporal Poisson model}
\usage{
\method{plot}{sepstlppm}(x,
do.points = TRUE,
par = TRUE,
...)
}
\arguments{
\item{x}{An object of class \code{sepstlppm}}

\item{do.points}{Add points to plot}

\item{par}{Default to \code{TRUE}. If \code{par=FALSE},
the user is asked for input, before a new figure is drawn.}

\item{...}{additional unused argument}
}
\description{
The function plots the fitted intensity, displayed both in space and in space and time.
}
\value{No return value, called for side effects}
\examples{
\donttest{

mod1 <- sepstlppm(valenciacrimes[1:2500, ], spaceformula = ~x,
timeformula = ~ crime_hour + week_day, L = valencianet)

plot.sepstlppm(mod1)

}
}
\author{
Nicoletta D'Angelo
}
