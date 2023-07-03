\name{sepstlppm}
\alias{sepstlppm}
\title{Fit a separable spatio-temporal Poisson process model on a linear network}
\usage{
sepstlppm(x, spaceformula, timeformula, L)
}
\arguments{
\item{x}{A dataframe}

\item{spaceformula}{A formula for the spatial component. See \link{lppm} for details}

\item{timeformula}{A formula for the temporal component. It fits a log-linear model with the \link{glm} function}

\item{L}{A linear network of class \code{linnet}}
}
\value{
An object of class \code{sepstlppm}
}
\description{
Fit a separable spatio-temporal Poisson process model on a linear network
}
\examples{
\donttest{
mod1 <- sepstlppm(valenciacrimes[1:2500, ], spaceformula = ~x,
timeformula = ~ crime_hour + week_day, L = valencianet)

}
}
