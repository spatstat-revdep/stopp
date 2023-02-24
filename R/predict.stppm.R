

#' Extract the the fitted intensity of a spatio-temporal Poisson process model
#'
#' @param object An object of class \code{stppm}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#'
#' @references
#' D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
#'
#' @seealso
#' \link{stppm}, \link{print.stppm}, \link{summary.stppm}, \link{plot.stppm}
#'
#' @return A vector containing the fitted intensity values
#'
#' @examples
#'
#' \dontrun{
#' ## Homogeneous
#' set.seed(2)
#' ph <- rstpp(lambda = 200, nsim = 1, seed = 2, verbose = T)
#' hom1 <- stppm(ph, formula = ~ 1)
#'
#' predict.stppm(hom1)
#'
#' ## Inhomogeneous
#' pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6),
#' nsim = 1, seed = 2, verbose = T)
#' inh1 <- stppm(pin, formula = ~ x)
#'
#' predict.stppm(inh1)
#'
#'}
#'
predict.stppm <- function(object, ...){
  if(!any(class(object) == "stppm")) stop("class(object) must be stppm")
  return(object$l)
}




