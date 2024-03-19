#' Display LISTA functions
#'
#' @description 
#' This function works on the objects of class \code{lista}, 
#' as returned by
#' \code{localSTLKinhom} or \code{localSTLginhom}, plotting the specified
#' LISTA functions. 
#' 
#' @param x An object of class \code{lista} 
#' @param id The id of the LISTA to display.
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
#' \donttest{
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
#' display(k, id = 1:8)
#'
#' }
#'
#'
#' @references
#' D’Angelo, N., Adelfio, G., and Mateu, J. (2021). Assessing local differences between the spatio-temporal second-order structure of two point patterns occurring on the same linear network. Spatial Statistics, 45, 100534.
#'
#' D’Angelo, N., Adelfio, G.,  and Mateu, J. (2022). Local inhomogeneous second-order characteristics for spatio-temporal point processes on linear networks. Stat Papers. https://doi.org/10.1007/s00362-022-01338-4
#'
display <- function(x, id){
  
  if(!inherits(x,"lista")) stop("class(x) must be lista")
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
  
    for(j in id){
      inhom <- list(x = x[[j]]$r, y = x[[j]]$t, z = x[[j]][[1]])
      fields::image.plot(inhom, main = paste("id", j), col = hcl.colors(12, "YlOrRd", rev = TRUE))
    }
  
}
