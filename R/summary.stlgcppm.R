#' Summary of a fitted LGCP model
#'
#'  The function Summarises the main information on the fitted model.
#'  provided. In this case of local parameters (both first- and second-order),
#'   the summary function contains information on
#'  their distributions.
#'
#' @param object An object of class \code{stlgcppm}
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
#' \dontrun{
#' 
#' catsub <- stp(greececatalog$df[1:200, ])
#' 
#' lgcp1 <- stlgcppm(catsub)
#'
#' summary(lgcp1)
#' 
#'}
#'
#'
#' @references
#' D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
#'
#' Siino, M., Adelfio, G., and Mateu, J. (2018). Joint second-order parameter estimation for spatio-temporal log-Gaussian Cox processes. Stochastic environmental research and risk assessment, 32(12), 3525-3539.
#'
#'
summary.stlgcppm <- function(object, ...){
  if (!inherits(object, c("stlgcppm"))) stop("X should be from class stlgcppm")
  
  cat("\nJoint minimum contrast fit \n")
  cat("for a log-Gaussian Cox process with \n")
  if(inherits(object$IntCoefs, "numeric")){ cat("global ")} else { cat("local ")}
  cat("first-order intensity and \n")
  if(inherits(object$CovCoefs, "numeric")){ cat("global ")} else { cat("local ")}
  cat("second-order intensity \n")
  cat("--------------------------------------------------\n")
  if(inherits(object$IntCoefs, "numeric")){
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
    cat("Estimated coefficients of the first-order intensity: \n")
    print(round(object$IntCoefs, 3))
  } else {
    if(ncol(object$IntCoefs) == 1){
      cat("Homogeneous Poisson process \n")
      cat("with median Intensity: ")
      cat(as.numeric(round(as.numeric(exp(median(object$IntCoefs[, 1])), 5))))
      cat("\n\n")
    } else {
      cat("Inhomogeneous Poisson process \n")
      cat("with Trend: ")
      print(object$formula)
      cat("\n")
    }
    cat("Summary of estimated coefficients of the first-order intensity \n")
    print(summary(object$IntCoefs))
  }
  cat("--------------------------------------------------\n")
  cat("Covariance function: ")
  cat(switch(object$cov,
             "separable" = "separable \n\n",
             "gneiting" = "gneiting \n\n",
             "iaco-cesare" = "iaco-cesare \n\n"))
  if(inherits(object$CovCoefs, "numeric")){
    cat("Estimated coefficients of the second-order intensity: \n")
    print(round(object$CovCoefs, 3))
  } else {
    cat("Summary of estimated coefficients of the second-order intensity \n")
    print(summary(object$CovCoefs))
  }
  cat("--------------------------------------------------\n")
  cat("Model fitted in ")
  cat(object$time, "\n")
}

