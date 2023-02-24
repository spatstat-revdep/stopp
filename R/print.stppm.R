
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
#' \link{plot.stppm}, \link{coef.stppm}
#'
#'
#'
#' @examples
#' \dontrun{
#' ## Homogeneous
#' set.seed(2)
#' ph <- rstpp(lambda = 200, nsim = 1, seed = 2, verbose = T)
#' hom1 <- stppm(ph, formula = ~ 1)
#'
#' hom1
#'
#' # Homogeneous Poisson process
#' # with Intensity: 202.093
#' #
#' # Estimated coefficients:
#' #   (Intercept)
#' # 5.309
#'
#' ## Inhomogeneous
#' pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6),
#' nsim = 1, seed = 2, verbose = T)
#' inh1 <- stppm(pin, formula = ~ x)
#'
#' inh1
#'
#' # Inhomogeneous Poisson process
#' # with Trend: ~x
#' #
#' # Estimated coefficients:
#' #   (Intercept)           x
#' # 2.180       5.783
#'}
#'
#'
#' @references
#' D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
#'
#'
#'
print.stppm <- function(x, ...){
  if(!any(class(x) == "stppm")) stop("class(x) must be stppm")
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


