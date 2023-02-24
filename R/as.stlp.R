#' Convert a \code{stlpp} object to a \code{stlp} object
#'
#' This function converts an object of class \code{stlpp} to an object of class \code{stlp}
#'
#' @param x An object of class \code{stlpp}
#'
#' @return An object of class \code{stlp}
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stp}, \link{summary.stlp}, \link{print.stlp},
#' \link{as.stlpp}, \link{plot.stlp}
#'
#' @examples
#'
#' \donttest{
#'
#' set.seed(12345)
#' stlpp1 <- stlnpp::rpoistlpp(.2, a = 0, b = 5, L = stlnpp::easynet)
#'
#' class(as.stlp(stlpp1))
#'
#' }
#'
#' @references
#' Moradi M, Cronie O, Mateu J (2020). stlnpp: Spatio-temporal analysis of point patterns on linear networks.
#'
as.stlp <- function(x){
  if(!any(class(x) == "stlpp")) stop("class(x) must be stlpp")

  stp(cbind(x$data$x, x$data$y, x$data$t), x$domain)
}
