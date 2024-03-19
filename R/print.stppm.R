
#' Print of a fitted spatio-temporal Poisson process model
#'
#'  The function prints the main information of the fitted model.
#'
#' @param x An object of class \code{stppm}
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
#'
#' @examples
#' \dontrun{
#' 
#' pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6),
#'              seed = 2)
#' inh1 <- stppm(pin, formula = ~ x)
#' 
#' inh1
#' 
#'}
#'
#'
#' @references
#' D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
#'
#'
#'
print.stppm <- function(x, ...){
  if (!inherits(x, c("stppm"))) stop("X should be from class stppm")
  
  if(length(x$IntCoefs) == 1){
    cat("Homogeneous Poisson process \n")
    cat("with Intensity: ")
    cat(as.numeric(round(as.numeric(exp(x$IntCoefs)), 5)))
    cat("\n\n")
  } else {
    cat("Inhomogeneous Poisson process \n")
    cat("with Trend: ")
    print(x$formula)
    cat("\n")
  }
  cat("Estimated coefficients: \n")
  print(round(x$IntCoefs, 3))
}


