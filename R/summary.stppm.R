
#' Summary of a fitted spatio-temporal Poisson process model
#'
#'  The function summarises the main information of the fitted model.
#'
#' @param object An object of class \code{stppm}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stppm}, \link{print.stppm},
#' \link{plot.stppm}
#'
#'
#' @examples
#'
#' set.seed(2)
#' pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6))
#' inh1 <- stppm(pin, formula = ~ x)
#' 
#' summary(inh1)
#'
#'
#' @references
#' D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
#'
#'
#'
summary.stppm <- function(object, ...){

  if(length(object$IntCoefs) == 1){
    cat("Homogeneous Poisson process \n")
    cat("with Intensity: ")
    cat(as.numeric(round(as.numeric(exp(object$IntCoefs)), 5)))
    cat("\n\n")
  } else {
    cat("Inhomogeneous Poisson process \n")
    cat("with Trend: ")
    print(object$formula)
    cat("\n")
  }
  cat("Estimated coefficients: \n")
  print(round(object$IntCoefs, 3))
}

