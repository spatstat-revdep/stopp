#' Plot a \code{stp} object
#'
#' This function plots the point pattern stored in the stp object given in input,
#' in a three panel plot representing the 3Dplot of the coordinates, and the
#' marginal spatial and temporal coordinates.
#'
#' @param x An object of class \code{stp}
#' @param tcum If \code{TRUE} (default option),
#' the temporal point pattern is displayed cumulatively.
#' A barplot is automatically plotted if there are repeated
#' counts (typically with discrete times).
#' @param marg Default to \code{TRUE}. If \code{FALSE}, only the
#' spatio-temporal point pattern is plotted.
#' @param col The color of the points. Default to \code{"black"}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stp}, \link{print.stp}, \link{summary.stp}
#'
#' @examples
#'
#' \dontrun{
#'
#' set.seed(2)
#' df <- data.frame(cbind(runif(100), runif(100), runif(100)))
#' 
#' stp1 <- stp(df)
#' #plot
#' plot(stp1)
#' 
#' #cumulative time occurrances
#' plot(stp1, tcum = FALSE)
#' 
#' #change color of points
#' plot(stp1, col = "blue")
#' 
#' #display only in space-time
#' plot(stp1, marg = FALSE)
#' 
#' #discrete times
#' set.seed(2)
#' stp2 <- stp(data.frame(cbind(runif(100), runif(100), round(runif(100) * 100))))
#' plot(stp2)
#' 
#' }
#'
plot.stp <- function(x, tcum = TRUE, marg = TRUE, col = 1, ...)
{
  if(!inherits(x,"stp")) stop("class(x) must be stp")

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if(marg == TRUE) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(1, 3))
  }

  plot3D::scatter3D(x$df$x, x$df$y, x$df$t,
                    theta = - 45, phi = 20, pch = 20,
                    ticktype = "detailed",
                    col = col,
                    xlab = "x", ylab = "y", zlab = "t",
                    main = c("Spatio-temporal point pattern"))
  if(marg == TRUE) {
    plot(x$df$x ,x$df$y, main = "Spatial point pattern", xlab = "x", ylab = "y",
         pch = 20, col = col)
    if(any(table(x$df$t) > 1)){
      barplot(table(x$df$t), main = "Temporal point pattern", xlab = "t")
    } else {
      if(tcum == TRUE){
        plot(sort(x$df$t), sort(cumsum(x$df$t)), main = "Temporal point pattern",
             xlab = "t", ylab = "", yaxt = "n", type = "l", pch = 20, col = col)
      } else {
        plot(x$df$t, rep(1, length(x$df$t)), main = "Temporal point pattern",
             xlab = "t", ylab = "", yaxt = "n", pch = 20, col = col)
      }
    }
  }
}
