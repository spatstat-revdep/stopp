#' Plot a \code{stlp} object
#'
#' This function plots the point pattern on a linear network
#'  stored in the \code{stlp} object given in input,
#' in a three panel plot representing the \link{plot3D} of the coordinates, and the
#' marginal spatial and temporal coordinates.
#'
#' @param x An object of class \code{stp}
#' @param tcum Default to \code{FALSE}.
#' If \code{TRUE}, the temporal point pattern is displayed cumulatively.
#' A barplot is automatically plotted if there are repeated
#' counts (typically with discrete times).
#' @param marg Default to \code{TRUE}. If \code{FALSE}, only the
#' spatio-temporal point pattern is plotted.
#' @param col The color of the points. Default to \code{"black"}
#' @param cols The color of the linear network. Default to \code{"grey"}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stp}, \link{summary.stlp}, \link{print.stlp},
#' \link{as.stlpp}, \link{as.stlp}
#'
#' @examples
#'
#' \dontrun{
#'
#' set.seed(12345)
#' stlpp1 <- rpoistlpp(.2, a = 0, b = 5, L = easynet)
#' df0 <- cbind(stlpp1$data$x, stlpp1$data$y, stlpp1$data$t)
#' L0 <- stlpp1$domain
#' stlp1 <- stp(df0, L0)
#'
#' #plot
#' plot(stlp1)
#'
#' #cumulative time occurrances
#' plot(stlp1, tcum = TRUE)
#'
#' #change color of points
#' plot(stlp1, col = "blue")
#'
#' #change color of network
#' plot(stlp1, cols = "magenta")
#'
#' #display only in space-time
#' plot(stlp1, marg = FALSE)
#'
#' #discrete times
#' stp2 <- stp(cbind(stlpp1$data$x, stlpp1$data$y, round(stlpp1$data$t)))
#' plot(stp2)
#'
#' }
plot.stlp <- function(x, tcum = FALSE, marg = TRUE, col = 1, cols = "grey", ...)
{
  if(!any(class(x) == "stlp")) stop("class(x) must be stlp")

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if(marg == TRUE) par(mfrow = c(1, 3)) else par(mfrow = c(1, 1))

  plot3D::scatter3D(x$df$x, x$df$y, x$df$t,
                    theta = - 45, phi = 20,  pch = 16,
                    ticktype = "detailed",
                    col = col,
                    xlab = "x", ylab = "y", zlab = "t",
                    main = c("Spatio-temporal point pattern"))

  if(marg == TRUE) {
    plot(x$L, col = cols, main = "Spatial point pattern")
    title(xlab = "x", ylab = "y")
    axis(1)
    axis(2)
    box()
    points(x$df$x, x$df$y, pch = 16, col = col)
    if(any(table(x$df$t) > 1)){
      barplot(table(x$df$t), main = "Temporal point pattern", xlab = "t")
    } else {
      if(tcum == TRUE){
        plot(sort(x$df$t), sort(cumsum(x$df$t)), main = "Temporal point pattern",
             xlab = "t", ylab = "", yaxt = "n", type = "l", pch = 16, col = col)
      }else{
        plot(x$df$t, rep(1, length(x$df$t)), main = "Temporal point pattern",
             xlab = "t", ylab = "", yaxt = "n", pch = 16, col = col)
      }

    }
  }
}
