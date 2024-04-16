#' Print of the result of the permutation local test
#'
#' It prints the main information on the result of the local permutation test
#'  performed with \link{localtest}  on either a \code{stp} or \code{stlp} object:
#'  whether the local test was run on point patterns lying on a linear network or not;
#'  the number of points in the background \code{X} and alternative \code{Z} patterns;
#'  the number of points in \code{X} which exhibit local differences in the second-order
#'  structure with respect to \code{Z}, according to the performed test.
#'
#' @param x An object of class  \code{localtest}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{localtest}, \link{summary.localtest},
#' \link{plot.localtest}
#'
#' @examples
#'
#' \donttest{
#'
#' set.seed(2)
#' X <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)},
#'             par = c(.005, 5))
#' Z <- rstpp(lambda = 30)
#' 
#' test <- localtest(X, Z, method = "K", k = 3)
#' 
#' test
#'
#'}
#'
#'
#' @references
#' D’Angelo, N., Adelfio, G., and Mateu, J. (2021). Assessing local differences between the spatio-temporal second-order structure of two point patterns occurring on the same linear network. Spatial Statistics, 45, 100534.
#'
#' Siino, M., Rodríguez‐Cortés, F. J., Mateu, J. ,and Adelfio, G. (2018). Testing for local structure in spatiotemporal point pattern data. Environmetrics, 29(5-6), e2463.
#'
#'
print.localtest <- function(x, ...){

  if(inherits(x$Xsig, "stlp")){
    cat("Test for local differences between two \n")
    cat("spatio-temporal point patterns on a linear network \n")
    cat("--------------------------------------------------\n")
  } else {
    cat("Test for local differences between two \n")
    cat("spatio-temporal point patterns \n")
    cat("--------------------------------------\n")
  }
  cat(paste("Background pattern X:",nrow(x$X$df), " \n"))
  cat(paste("Alternative pattern Z:",nrow(x$Z$df), " \n  \n"))
  cat(paste(nrow(x$Xsig$df), "significant points at alpha =", x$alpha), "\n")
}

