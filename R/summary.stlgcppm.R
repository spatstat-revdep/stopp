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
#' \link{stlgcppm}, \link{print.stlgcppm}, \link{localsummary.stlgcppm},
#' \link{plot.stlgcppm}, \link{localplot.stlgcppm}
#'
#'
#'
#' @examples
#' \dontrun{
#' # Example with complex seismic point pattern
#' data("greececatalog")
#'
#' lgcp1 <- stlgcppm(greececatalog, formula = ~ 1, first = "global", second = "global")
#' lgcp2 <- stlgcppm(greececatalog, formula = ~ x, first = "local", second = "global")
#'
#' summary(lgcp1)
#'
#'
#' # Joint minimum contrast fit
#' # for a log-Gaussian Cox process with
#' # global first-order intensity and
#' # global second-order intensity
#' # --------------------------------------------------
#' #   Homogeneous Poisson process
#' # with Intensity: 0.0064
#' #
#' # Estimated coefficients of the first-order intensity:
#' #   (Intercept)
#' # -5.052
#' # --------------------------------------------------
#' #   Covariance function: separable
#' #
#' # Estimated coefficients of the second-order intensity:
#' #   sigma  alpha   beta
#' # 10.984  0.224 47.076
#' # --------------------------------------------------
#' #   Model fitted in 0.873 minutes
#'
#'
#' summary(lgcp2)
#'
#'
#' # Joint minimum contrast fit
#' # for a log-Gaussian Cox process with
#' # local first-order intensity and
#' # global second-order intensity
#' # --------------------------------------------------
#' #   Inhomogeneous Poisson process
#' # with Trend: ~x
#' #
#' # Summary of estimated coefficients of the first-order intensity
#' # (Intercept)           x
#' # Min.   :-6.282   Min.   :-0.96831
#' # 1st Qu.:-2.387   1st Qu.:-0.36685
#' # Median : 2.122   Median :-0.25871
#' # Mean   : 2.052   Mean   :-0.26309
#' # 3rd Qu.: 4.569   3rd Qu.:-0.07325
#' # Max.   :17.638   Max.   : 0.10269
#' # --------------------------------------------------
#' #   Covariance function: separable
#' #
#' # Estimated coefficients of the second-order intensity:
#' #   sigma  alpha   beta
#' # 2.612  0.001 36.415
#' # --------------------------------------------------
#' #   Model fitted in 3.503 minutes
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
  if(!any(class(object) == "stlgcppm")) stop("class(object) must be stlgcppm")
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
  cat(object$time)
}

