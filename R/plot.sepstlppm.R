#' Plot of the fitted intensity of a separable spatio-temporal Poisson model
#'
#' The function plots the fitted intensity, displayed both in space and in space and time.
#'
#' @param x An object of class \code{sepstlppm}
#' @param do.points Add points to plot
#' @param par Default to \code{TRUE}. If \code{par=FALSE},
#' the user is asked for input, before a new figure is drawn.
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#'
#' @examples
#' \dontrun{
#'
#' mod1 <- sepstlppm(valenciacrimes[1:2500, ], spaceformula = ~x,
#' timeformula = ~ crime_hour + week_day, L = valencianet)
#'
#' plot.sepstlppm(mod1)
#'
#' }
plot.sepstlppm <- function(x,
                          do.points = TRUE,
                          par = TRUE,
                          ...){

  if (!inherits(x, "sepstlppm")) stop("x should an object of class sepstlppm")

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  mark_int <- x$l

  if(par == TRUE){
    par(mfrow = c(1, 2))
    par(mar = c(5, 4, 4, 2) + 0.1 - c(4, 1 , 1, 1))
    plot(x$spaceint_plot, col = grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
         main = c("Intensity in space \n Density Kernel Smoothing"))
    if(do.points == TRUE){points(x$x$x, x$x$y, pch = ".")}

    par(mar = c(5, 4, 4, 2) + 0.1 - c(4, 1 , 1, -0.5))

    plot3D::scatter3D(x$x$x, x$x$y, x$x$t,
                      theta = - 45, phi = 20,
                      col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                          range = range(mark_int)),
                                 "stuff")$outputs,
                      ticktype = "detailed", pch = 19,
                      colvar = mark_int,
                      xlab="x",ylab="y",zlab="t",
                      main = c("Intensity in space-time \n Pointwise computation"))
    par(mar = c(5, 4, 4, 2) + 0.1)
  } else {
    par(mfrow = c(1, 1))
    par(ask = TRUE)
    plot(x$spaceint_plot, col = grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
         main = c("Intensity in space \n Density Kernel Smoothing"))
    if(do.points == TRUE){points(x$x$x, x$x$y, pch = ".")}
    par(ask = TRUE)
    par(mar = c(5, 4, 4, 2) + 0.1 - c(4, 1 , 1, -0.5))
    plot3D::scatter3D(x$x$x, x$x$y, x$x$t,
                      theta = - 45, phi = 20,
                      col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                          range = range(mark_int)),
                                 "stuff")$outputs,
                      ticktype = "detailed", pch = 19,
                      colvar = mark_int,
                      xlab="x",ylab="y",zlab="t",
                      main = c("Intensity in space-time \n Pointwise computation"))
    par(mar = c(5, 4, 4, 2) + 0.1)
    par(ask = FALSE)
  }
}
