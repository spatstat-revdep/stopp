#' Extract the fitted coefficients of a local spatio-temporal Poisson process model
#'
#' @param object An object of class \code{locstppm}
#' @param ... additional unused argument
#'
#' @return A list containing the global and local fitted coefficients
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
#' \donttest{
#' set.seed(2)
#' pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6),
#'              nsim = 1, seed = 2, verbose = TRUE)
#' inh00_local <- locstppm(pin, formula = ~ 1)
#' inh01_local <- locstppm(pin, formula = ~ x)
#'
#' coef(inh00_local)
#' coef(inh01_local)
#'
#'}
#'
coef.locstppm <- function(object, ...){
  if(!any(class(object) == "locstppm")) stop("class(object) must be locstppm")
  out <- list(
    globalpars = object$IntCoefs,
    localpars = summary(object$IntCoefs_local)
  )
  return(out)

}
