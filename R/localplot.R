#' Plot the coefficients of a fitted local spatio-temporal Poisson process or local LGCP model
#'
#' The function plots the local estimates of a fitted local spatio-temporal Poisson process or local LGCP model
#'
#' @param x An object of class \code{locstppm} or \code{stlgcppm}
#' @param par Default to \code{TRUE}.
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{locstppm}, \link{stlgcppm}
#'
#'
#'
#' @examples
#' \dontrun{
#' 
#' # Local spatio-temporal Poisson process model
#' 
#' inh <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, 
#'              par = c(0.005, 5), seed = 2)
#' inh_local <- locstppm(inh, formula = ~ x)
#' 
#' localplot(inh_local)
#' 
#' # Local LGCP
#' catsub <- stp(greececatalog$df[1:200, ])
#' 
#' lgcp_loc <- stlgcppm(catsub, formula = ~ x, first = "local")
#'
#' localplot(lgcp_loc)
#' 
#'
#'}
#'
#'
#'
#'
#' @references
#' D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
#'
#'
localplot <- function(x, par = TRUE){
  if (!inherits(x, c("stlgcppm", "locstppm"))) stop("x should be either from class stlgcppm or locstppm")
  UseMethod("localplot")
}

