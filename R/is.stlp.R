
#' Test whether an object is a spatio-temporal point pattern on a linear network
#'
#' @param x Any object
#'
#' @return
#' \code{TRUE} is \code{x} is a spatio-temporal point pattern on a linear network,
#'  otherwise \code{FALSE}
#'
#' @details
#' This function tests whether the object \code{x} is a spatio-temporal
#' point pattern object of class \code{stlp}. See \link{stp} for details of this class
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso \link{stp}, \link{summary.stlp}, \link{print.stlp}, \link{plot.stlp}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(12345)
#' stlpp1 <- stlnpp::rpoistlpp(.2, a = 0, b = 5, L = stlnpp::easynet)
#' is.stlp(stlpp1)
#'
#' df0 <- cbind(stlpp1$data$x, stlpp1$data$y, stlpp1$data$t)
#' L0 <- stlpp1$domain
#' stlp1 <- stp(df0, L0)
#'
#' is.stlp(stlp1)
#' }
is.stlp <- function(x){
  if(any(class(x) == "stlp")) print(TRUE) else print(FALSE)
}
