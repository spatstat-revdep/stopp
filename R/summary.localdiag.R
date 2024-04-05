#' Summary of the diagnostics performed on a spatio-temporal point process model
#'
#' It summarises the main information on the result of the local diagnostics
#'  performed with \link{localdiag} on either a \code{stp} or \code{stlp} object:
#'  whether the local test was run on point patterns lying on a linear network or not;
#'  the number of points in the analysed spatio-temporal
#'  point pattern \code{X};
#'  the number of points of \code{X} which are identified as outlying by the
#'  previously performed local diagnostics.
#'
#' @param object A \code{localdiag} object
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo and Giada Adelfio
#'
#' @seealso
#' \link{infl}, \link{plot.localdiag}, \link{print.localdiag}
#'
#'
#' @examples
#' \dontrun{
#'
#' inh <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, 
#'              par = c(.3, 6), seed = 2)
#' 
#' mod1 <- stppm(inh, formula = ~ 1)
#' 
#' resmod1 <- localdiag(inh, mod1$l, p = .9)
#' 
#' summary(resmod1)
#' 
#' }
#'
#'
#'
#' @references
#' Adelfio, G., Siino, M., Mateu, J., and Rodríguez-Cortés, F. J. (2020). Some properties of local weighted second-order statistics for spatio-temporal point processes. Stochastic Environmental Research and Risk Assessment, 34(1), 149-168.
#'
#' D’Angelo, N., Adelfio, G.  and Mateu, J. (2022) Local inhomogeneous second-order characteristics for spatio-temporal point processes on linear networks. Stat Papers. https://doi.org/10.1007/s00362-022-01338-4
#'
#'
#'
#'
#'
summary.localdiag <- function(object, ...){
  if (!inherits(object, c("localdiag"))) stop("X should be from class localdiag")
  
  if(inherits(object$X, "stlp")){
    cat(paste("Points outlying from the", object$p , "percentile\n"))
    cat("of the analysed spatio-temporal point pattern on a linear network \n")
    cat("--------------------------------------------------\n")
  } else {
    cat(paste("Points outlying from the", object$p , "percentile\n"))
    cat("of the analysed spatio-temporal point pattern \n")
    cat("--------------------------------------\n")
  }
  cat(paste("Analysed pattern X:",nrow(object$X$df), "points \n"))
  cat(paste(length(object$ids), "outlying points"), "\n")
}
