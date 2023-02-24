\name{plot.sepstppm}
\alias{plot.sepstppm}
\title{Plot of the fitted intensity of a separable spatio-temporal Poisson model}
\usage{
\method{plot}{sepstppm}(
  x,
  scaler = c("silverman", "IQR", "sd", "var"),
  do.points = TRUE,
  print.bw = FALSE,
  zap = 1e-05,
  par = TRUE,
  sig = NULL,
  ...
)
}
\arguments{
\item{x}{An object of class \code{sepstppm}}

\item{scaler}{Optional. Controls the value for a scalar representation of the
spatial scale of the data.
Either a character string, \code{"silverman"} (default), \code{"IQR"},
 \code{"sd"}, or \code{"var"};
or positive numeric value(s). See \link{OS}.}

\item{do.points}{Add points to plot}

\item{print.bw}{It prints the estimated oversmoothing (\link{OS}) bandwidth selector}

\item{zap}{Noise threshold factor (default to 0.00001). A numerical value greater than or equal to 1.
If the range of pixel values is less than \code{zap * .Machine$double.eps},
 the image will be treated as constant. This avoids displaying images which
 should be constant but contain small numerical errors.}

\item{par}{Default to \code{TRUE}. If \code{par=FALSE},
the user is asked for input, before a new figure is drawn.}

\item{sig}{Smoothing bandwidth for spatial representation}

\item{...}{additional unused argument}
}
\description{
The function plots the fitted intensity, displayed both in space and in space and time.
}
\value{No return value, called for side effects}
\examples{
\donttest{

df1 <- valenciacrimes[valenciacrimes$x < 210000 & valenciacrimes$x > 206000
& valenciacrimes$y < 4377000 & valenciacrimes$y > 4373000, ]

mod1 <- sepstppm(df1, spaceformula = ~x * y,
                 timeformula = ~ crime_hour + week_day)

plot.sepstppm(mod1)

}
}
\author{
Nicoletta D'Angelo
}
