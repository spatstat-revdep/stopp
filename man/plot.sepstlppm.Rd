% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.sepstlppm.R
\name{plot.sepstlppm}
\alias{plot.sepstlppm}
\title{Plot of the fitted intensity of a separable spatio-temporal Poisson model}
\usage{
\method{plot}{sepstlppm}(x, do.points = TRUE, par = TRUE, ...)
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
\examples{

crimesub <- stpm(valenciacrimes$df[101:200, ],
           names = colnames(valenciacrimes$df)[-c(1:3)],
L = valencianet)

mod1 <- sepstlppm(crimesub, spaceformula = ~x ,
                  timeformula = ~ day)
                  
plot(mod1)

}
\author{
Nicoletta D'Angelo
}
