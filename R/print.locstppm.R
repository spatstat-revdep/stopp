#' Print of a fitted local spatio-temporal Poisson process model
#'
#'  The function prints the main information of the distribution of the parameters
#'  of a fitted local spatio-temporal Poisson process model.
#'
#' @param x An object of class \code{locstppm}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{locstppm}, \link{summary.locstppm},
#' \link{plot.locstppm}
#'
#'
#'
#' @examples
#' \dontrun{
#' set.seed(2)
#' pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6),
#'              nsim = 1, seed = 2, verbose = T)
#' inh00_local <- locstppm(pin, formula = ~ 1)
#' inh01_local <- locstppm(pin, formula = ~ x)
#'
#' inh00_local
#' inh01_local
#'
#' }
#'
#'
#' @references
#' D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
#'
#'
#'
print.locstppm <- function(x, ...){
  if(!any(class(x) == "locstppm")) stop("class(x) must be locstppm")
  if(ncol(x$IntCoefs_local) == 1){
    cat("Homogeneous Poisson process \n")
    cat("with median Intensity: ")
    cat(as.numeric((as.numeric((median(x$IntCoefs_local[, 1])), 5))))
    cat("\n\n")
  } else {
    cat("Inhomogeneous Poisson process \n")
    cat("with Trend: ")
    print(x$formula)
    cat("\n")
  }
  cat("Summary of estimated coefficients \n")
  print(summary(x$IntCoefs_local))

}

