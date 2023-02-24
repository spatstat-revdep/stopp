#' Test whether an object is a spatio-temporal point pattern
#'
#' @param x Any object
#'
#' @return
#' \code{TRUE} is \code{x} is a spatio-temporal point pattern, otherwise \code{FALSE}
#'
#' @details
#' This function tests whether the object \code{x} is a spatio-temporal
#' point pattern object of class \code{stp}. See \link{stp} for details of this class
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso \link{stp}, \link{summary.stp}, \link{print.stp}, \link{plot.stp}
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' set.seed(12345)
#' rpp1 <- stpp::rpp(lambda = 200, replace = FALSE)
#' is.stp(rpp1)
#'
#' df0 <- cbind(rpp1$xyt[, 1], rpp1$xyt[, 2], rpp1$xyt[, 3])
#' stp1 <- stp(df0)
#'
#' is.stp(stp1)
#' }
is.stp <- function(x){
  if(any(class(x) == "stp")) print(TRUE) else print(FALSE)
}
