\name{sepstppm}
\alias{sepstppm}
\title{Fit a separable spatio-temporal Poisson process model}
\usage{
sepstppm(x, spaceformula, timeformula)
}
\arguments{
\item{x}{A dataframe}

\item{spaceformula}{A formula for the spatial component. See \link{ppm} for details}

\item{timeformula}{A formula for the temporal component. It fits a log-linear model with the \link{glm} function}
}
\value{
An object of class \code{sepstppm}
}
\description{
Fit a separable spatio-temporal Poisson process model
}
\examples{
\donttest{
df1 <- valenciacrimes[valenciacrimes$x < 210000 & valenciacrimes$x > 206000
& valenciacrimes$y < 4377000 & valenciacrimes$y > 4373000, ]

mod1 <- sepstppm(df1, spaceformula = ~x * y,
                timeformula = ~ crime_hour + week_day)
}



}
