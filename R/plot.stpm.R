#' Plot a \code{stpm} object
#'
#' This function plots the marked point pattern stored in the \code{stpm} object
#'  given in input,
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
#' \dontrun{
#'
#' df <- data.frame(cbind(runif(100), runif(100), runif(100), rpois(100, 15),
#' rpois(100, 30)))
#'
#' stpm1 <- stpm(df)
#'
#' plot(stpm1)
#' 
#' ## Categorical marks
#' 
#' dfA <- data.frame(x = runif(100), y = runif(100), t = runif(100), 
#'                   m1 = rnorm(100), m2 = rep(c("C"), times = 100))
#' dfB <- data.frame(x = runif(50), y = runif(50), t = runif(50), 
#'                   m1 = rnorm(25), m2 = rep(c("D"), times = 50))
#' 
#' stpm2 <- stpm(rbind(dfA, dfB), names = c("continuous", "dichotomous"))
#' 
#' plot(stpm2)
#' }
#'
plot.stpm <- function(x, ...){
  if(!inherits(x,"stpm")) stop("class(x) must be stpm")
  
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
