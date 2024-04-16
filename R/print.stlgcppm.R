#' Print of a fitted LGCP model
#'
#'  The function prints the main information on the fitted model.
#'   In this case of local parameters (both first- and second-order),
#'   the summary function contains information on
#'  their distributions.
#'
#' @param x An object of class \code{stlgcppm}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo and Giada Adelfio
#'
#' @seealso
#' \link{stlgcppm}, \link{print.stlgcppm}, \link{localsummary},
#' \link{plot.stlgcppm}, \link{localplot}
#'
#'
#'
#' @examples
#' 
#' catsub <- stp(greececatalog$df[1:200, ])
#' 
#' lgcp1 <- stlgcppm(catsub)
#'
#' lgcp1
#' 
#'
#'
#' @references
#' D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
#'
#' Siino, M., Adelfio, G., and Mateu, J. (2018). Joint second-order parameter estimation for spatio-temporal log-Gaussian Cox processes. Stochastic environmental research and risk assessment, 32(12), 3525-3539.
#'
#'
print.stlgcppm <- function(x, ...){

  cat("\nJoint minimum contrast fit \n")
  cat("for a log-Gaussian Cox process with \n")
  if(inherits(x$IntCoefs, "numeric")){ cat("global ")} else { cat("local ")}
  cat("first-order intensity and \n")
  if(inherits(x$CovCoefs, "numeric")){ cat("global ")} else { cat("local ")}
  cat("second-order intensity \n")
  cat("--------------------------------------------------\n")
  if(inherits(x$IntCoefs, "numeric")){
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
    cat("Estimated coefficients of the first-order intensity: \n")
    print(round(x$IntCoefs, 3))
  } else {
    if(ncol(x$IntCoefs) == 1){
      cat("Homogeneous Poisson process \n")
      cat("with median Intensity: ")
      cat(as.numeric(round(as.numeric(exp(median(x$IntCoefs[, 1])), 5))))
      cat("\n\n")
    } else {
      cat("Inhomogeneous Poisson process \n")
      cat("with Trend: ")
      print(x$formula)
      cat("\n")
    }
    cat("Summary of estimated coefficients of the first-order intensity \n")
    print(summary(x$IntCoefs))
  }
  cat("--------------------------------------------------\n")
  cat("Covariance function: ")
  cat(switch(x$cov,
             "separable" = "separable \n\n",
             "gneiting" = "gneiting \n\n",
             "iaco-cesare" = "iaco-cesare \n\n"))
  if(inherits(x$CovCoefs, "numeric")){
    cat("Estimated coefficients of the second-order intensity: \n")
    print(round(x$CovCoefs, 3))
  } else {
    cat("Summary of estimated coefficients of the second-order intensity \n")
    print(summary(x$CovCoefs))
  }
  cat("--------------------------------------------------\n")
  cat("Model fitted in ")
  cat(x$time, "\n")
}


