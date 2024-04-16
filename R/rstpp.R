#' Simulate homogeneous and inhomogeneous spatio-temporal Poisson point patterns
#'
#' This function creates a \code{stp} object, simulating a  spatio-temporal point pattern
#'  following either an
#' homogeneous or inhomogeneous intensity
#'
#' @param lambda Expected number of points to simulate
#' @param par Parameters of the reference intensity
#' @param nsim Number of patterns to simulate. Default to 1.
#' @param verbose Default to \code{FALSE}
#' @param minX Minimum of x coordinate range
#' @param maxX Maximum of x coordinate range
#' @param minY Minimum of y coordinate range
#' @param maxY Maximum of y coordinate range
#' @param minT Minimum of t coordinate range
#' @param maxT Maximum of t coordinate range
#'
#' @return A \code{stp} object
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stppm}
#'
#'
#' @examples
#'
#' # homogeneous Poisson processes
#' set.seed(2)
#' h1 <- rstpp(lambda = 500)
#'
#' set.seed(2)
#' h2 <- rstpp(lambda = 500, minX = 0,
#'              maxX = 2, minY = 3, maxY = 5, minT = 1, maxT = 9)
#'
#' set.seed(2)
#' h3 <- rstpp(lambda = 900, nsim = 3, verbose = TRUE)
#'
#' # inhomogeneous Poisson process
#' set.seed(2)
#' inh <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6))
#'
rstpp <- function(lambda = 500, nsim = 1, verbose = FALSE,  par = NULL,
                  minX = 0, maxX = 1, minY = 0, maxY = 1, minT = 0, maxT = 1){
  
  
  if(is.numeric(lambda)){
    par <- log(lambda)
    lambda <- function(x,y,t, a) {exp(a[1])}
  }
  
  if(nsim != 1){pp0 <- list(l = nsim)}
  
  for(i in 1:nsim){
    if(verbose == T) progressreport(i, nsim)
    lam = lambda(x = 1, y = 1, t = 1, a = par)
    candn = rpois(1, lam)
    candx = runif(candn, min = minX, max = maxX)
    candy = runif(candn, min = minY, max = maxY)
    candt = runif(candn, min = minT, max = maxT)
    d = runif(candn)
    lam2 = lambda(x = candx, y = candy, t = candt, a = par)
    lmax = max(lam2)
    keep = (d < lam2 / lmax)
    n = sum(keep)
    t1 = candt[keep]
    t = t1[order(t1)]
    lon = candx[keep][order(t1)]
    lat = candy[keep][order(t1)]
    if(nsim != 1){
      pp0[[i]] <- stp(cbind(lon, lat, t))
    } else {
      pp0 <- stp(cbind(lon, lat, t))
    }
    
  }
  
  return(pp0)
}
