#' Extract the fitted intensity of a local spatio-temporal Poisson process model
#'
#' @param object An object of class \code{locstppm}
#' @param ... additional unused argument
#'
#' @return A vector containing the fitted intensity
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
#' \link{locstppm}, \link{print.locstppm}, \link{summary.locstppm}, \link{plot.locstppm}
#'
#'
#'
#' @examples
#'
#' \dontrun{
#' set.seed(2)
#' pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6),
#'              nsim = 1, seed = 2, verbose = T)
#' inh00_local <- locstppm(pin, formula = ~ 1)
#' inh01_local <- locstppm(pin, formula = ~ x)
#'
#' predict.locstppm(inh00_local)
#' predict.locstppm(inh01_local)
#'
#'}
#'
predict.locstppm <- function(object, ...){
  if(!any(class(object) == "locstppm")) stop("class(object) must be locstppm")
  return(object$l_local)

}
