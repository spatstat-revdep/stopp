#' Plot a \code{stlpm} object
#'
#' This function plots the covariate stored in the stcov object given in input,
#' in a three panel plot representing the 3Dplot of the coordinates, and the
#' mark values.
#'
#' @param x An object of class \code{stpm}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stppm}
#'
#' @examples
#'
#' set.seed(2)
#' df <- data.frame(x = runif(100, 0, 0.8), y = runif(100, 0, 0.8), t = runif(100), m = rpois(100, 15))
#'
#' stlpm1 <- stpm(df, L = chicagonet)
#'
#' plot(stlpm1)
#' 
#'
plot.stlpm <- function(x, ...){

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  id <- !sapply(x$df, is.numeric)
  
  for(k in which(id)){
    x$df[, k] <- as.factor(x$df[, k])
  }
  
  m <- ncol(x$df) - 3
  if(m == 1) {
    oldpar <- par(mfrow = c(1, 1))
    on.exit(par(oldpar))
  } else {
    a <- as.numeric(max(spatstat.utils::primefactors(m)))
    b <- ceiling(m / a)
    oldpar <- par(mfrow = c(b, a))
    on.exit(par(oldpar))
  }
  names <- colnames(x$df)
  for(i in 4:(m + 3)){
    if(!id[i]){
      plot3D::scatter3D(x$df[, 1], x$df[, 2], x$df[, 3],
                        theta = - 45, phi = 20, pch = 20,
                        ticktype = "detailed",
                        colvar = x$df[, i],
                        xlab = "x", ylab = "y", zlab = "t",
                        main = paste(names[i]))
    } else {
      plot3D::scatter3D(x$df[, 1], x$df[, 2], x$df[, 3],
                        theta = - 45, phi = 20, pch = 20,
                        ticktype = "detailed",
                        colvar = as.integer(x$df[, i]),
                        col = 1:length(levels(x$df[, i])) + 4,
                        xlab = "x", ylab = "y", zlab = "t",
                        colkey = list(at = unique(as.integer(x$df[, i])), 
                                      labels = c(levels(x$df[, i]))),
                        main = paste(names[i]))
    }
    
  }
  
}
