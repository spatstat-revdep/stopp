#' Convert a \code{stpp} object to a \code{stp} object
#'
#' This function converts an object of class \code{stpp} to an object of class \code{stp}
#'
#' @param x An object of class \code{stpp}
#'
#' @return An object of class \code{stp}
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stp}, \link{print.stp}, \link{summary.stp},
#' \link{plot.stp}, \link{as.stpp}
#'
#' @examples
#'
#' \donttest{
#'
#' set.seed(12345)
#' rpp1 <- stpp::rpp(lambda = 200, replace = FALSE)
#'
#' class(as.stp(rpp1$xyt))
#'
#' }
#'
#' @references
#' Gabriel, E., Rowlingson, B. S., & Diggle, P. J. (2013). stpp: An R Package for Plotting, Simulating and Analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53(2), 1–29. https://doi.org/10.18637/jss.v053.i02
#'
as.stp <- function(x){

  if(!any(class(x) == "stpp")) stop("class(x) must be stpp")

  stp(cbind(x[, 1], x[, 2], x[, 3]))
}
