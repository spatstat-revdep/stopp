#' Plot a \code{stcov} object
#'
#' This function plots the covariate stored in the \code{stcov} object given in input,
#' in a three panel plot representing the 3Dplot of the coordinates, and the
#' covariate values.
#'
#' @param x An object of class \code{stcov}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stcov}
#'
#' @examples
#'
#'
#' set.seed(2)
#' df <- data.frame(runif(100), runif(100), runif(100), rpois(100, 15))
#'
#' cov <- stcov(df, interp = FALSE)
#' 
#' plot(cov)
#'
plot.stcov <- function(x, ...){

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  names <- colnames(x$df)
      plot3D::scatter3D(x$df[, 1], x$df[, 2], x$df[, 3],
                        theta = - 45, phi = 20, pch = 20,
                        ticktype = "detailed",
                        colvar = x$df[, 4],
                        xlab = "x", ylab = "y", zlab = "t",
                        main = paste(names[4]),
                        col = grDevices::hcl.colors(12, 
                                                    "ag_Sunset", rev = TRUE))
    
  
}
