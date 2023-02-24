#' Convert a \code{stp} object to a \code{stpp} object
#'
#' This function converts an object of class \code{stp} to an object of class \code{stpp}
#'
#' @param x An object of class \code{stp}
#'
#' @return An object of class \code{stpp}
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stp}, \link{print.stp}, \link{summary.stp},
#' \link{plot.stp}, \link{as.stp}
#'
#' @examples
#'
#' \donttest{
#'
#' set.seed(12345)
#' rpp1 <- stpp::rpp(lambda = 200, replace = FALSE)
#' df0 <- cbind(rpp1$xyt[, 1], rpp1$xyt[, 2], rpp1$xyt[, 3])
#' stp1 <- stp(df0)
#'
#' class(as.stpp(stp1))
#'
#' }
#'
#' @references
#' Gabriel, E., Rowlingson, B. S., & Diggle, P. J. (2013). stpp: An R Package for Plotting, Simulating and Analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53(2), 1–29. https://doi.org/10.18637/jss.v053.i02
#'
#'
as.stpp <- function(x){

  if(!any(class(x) == "stp")) stop("class(x) must be stp")

  out <- cbind(x$df$x, x$df$y, x$df$t)
  colnames(out) <- c("x", "y", "t")
  class(out) <- "stpp"
  return(out)
}
