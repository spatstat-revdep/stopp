#' Simulate homogeneous and inhomogeneous spatio-temporal Poisson point patterns on linear networks
#'
#' This function creates a \code{stlp} object, simulating a  spatio-temporal point pattern on
#' a linear network
#'  following either an
#' homogeneous or inhomogeneous intensity
#'
#' @param lambda Expected number of points to simulate
#' @param par Parameters of the reference intensity
#' @param nsim Number of patterns to simulate. Default to 1.
#' @param seed Seed to set, if ones wished to reproduce the analyses
#' @param verbose Default to \code{FALSE}
#' @param minX Minimum of x coordinate range
#' @param maxX Maximum of x coordinate range
#' @param minY Minimum of y coordinate range
#' @param maxY Maximum of y coordinate range
#' @param minT Minimum of t coordinate range
#' @param maxT Maximum of t coordinate range
#' @param L linear network
#'
#' @return A \code{stp} object
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @examples
#' \dontrun{
#'
#' h1 <- rstlpp(lambda = 500, seed = 2, L = chicagonet)
#'
#' inh <- rstlpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(4, 1.5),
#'         seed = 2, L = chicagonet)
#' }
#'
#'
rstlpp <- function(lambda = 500, nsim = 1, seed = NULL, verbose = FALSE,  par = NULL,
                  minX = 0, maxX = 1, minY = 0, maxY = 1, minT = 0, maxT = 1, L){


  if(is.numeric(lambda)){
    par <- log(lambda)
    lambda <- function(x,y,t, a) {exp(a[1])}
  }

  if(nsim != 1){pp0 <- list(l = nsim)}

  set.seed(seed)
  for(i in 1:nsim){
    if(verbose == T) progressreport(i, nsim)
    lam = lambda(x = 1, y = 1, t = 1, a = par)
    candn = rpois(1, lam)
    rand <- spatstat.random::datagen.runifpointOnLines(candn, spatstat.geom::as.psp(L))
    candx = rand$x
    candy = rand$y
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
      pp0[[i]] <- stp(cbind(lon, lat, t), L = L)
    } else {
      pp0 <- stp(cbind(lon, lat, t), L = L)
    }

  }

  return(pp0)
}





