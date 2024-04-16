#' Display LISTA functions
#'
#' @description 
#' This function works on the objects of class \code{lista}, 
#' as returned by
#' \code{localSTLKinhom} or \code{localSTLginhom}, plotting the specified
#' LISTA functions. 
#' 
#' @param x An object of class \code{lista} 
#' @param id The id of the LISTA to display
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo 
#'
#' @seealso
#' \link{localSTLKinhom}, \link{localSTLginhom}
#'
#'
#' @examples
#'
#' 
#' set.seed(2)
#' df_net <- data.frame(x = runif(25, 0, 0.85), y = runif(25, 0, 0.85), t = runif(25))
#' stlp1 <- stp(df_net, L = chicagonet)
#' lambda <- rep(diff(range(stlp1$df$x)) * diff(range(stlp1$df$y))
#'  * diff(range(stlp1$df$t)) / spatstat.geom::volume(stlp1$L),
#' nrow(stlp1$df))
#'
#' k <- localSTLKinhom(stlp1, lambda = lambda, normalize = TRUE)
#' 
#' plot(k, id = 1:9)
#'
#'
#'
#' @references
#' D’Angelo, N., Adelfio, G., and Mateu, J. (2021). Assessing local differences between the spatio-temporal second-order structure of two point patterns occurring on the same linear network. Spatial Statistics, 45, 100534.
#'
#' D’Angelo, N., Adelfio, G.,  and Mateu, J. (2022). Local inhomogeneous second-order characteristics for spatio-temporal point processes on linear networks. Stat Papers. https://doi.org/10.1007/s00362-022-01338-4
#'
plot.lista <- function(x, id, ...){
  
  if(missing(id)) stop("Specify which lista function to plot")
  
  n <- length(id)
  
  if(n == 1) {
    oldpar <- par(mfrow = c(1, 1))
    on.exit(par(oldpar))
  } else {
    a <- as.numeric(max(spatstat.utils::primefactors(n)))
    b <- n / a
    oldpar <- par(mfrow = c(b, a))
    on.exit(par(oldpar))
  }
  nd <- length(x[[1]]$r)
  nt <- length(x[[1]]$t)
    for(j in id){
      fields::image.plot(x[[j]][[1]], main = paste("id", j), xlab = "r",
                      ylab = "h", 
                         col = hcl.colors(12, "YlOrRd", rev = TRUE),
                      axes = FALSE)
      axis(1, at = seq(0, 1, l = nd), labels = round(x[[j]]$r, 3))
      axis(2, at = seq(0, 1, l = nt), labels = round(x[[j]]$t, 3))
      box()
    }
  
}
