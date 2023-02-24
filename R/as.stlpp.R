#' Convert a \code{stlp} object to a \code{stlpp} object
#'
#' This function converts an object of class \code{stlp} to an object of class \code{stlpp}
#'
#' @param x An object of class \code{stlp}
#'
#' @return An object of class \code{stlpp}
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stp}, \link{summary.stlp}, \link{print.stlp},
#' \link{as.stlp}, \link{plot.stlp}
#'
#' @examples
#'
#' \dontrun{
#'
#' set.seed(12345)
#' rpp1 <- stpp::rpp(lambda = 200, replace = FALSE)
#' df0 <- cbind(rpp1$xyt[, 1], rpp1$xyt[, 2], rpp1$xyt[, 3])
#' L0 <- spatstat.geom::domain(spatstat.data::chicago)
#' stlp1 <- stp(df0, L0)
#'
#' class(as.stlpp(stlp1))
#'
#' }
#'
#' @references
#' Moradi M, Cronie O, Mateu J (2020). stlnpp: Spatio-temporal analysis of point patterns on linear networks.
#'
#'
as.stlpp <- function(x){
  if(!any(class(x) == "stlp")) stop("class(x) must be stlp")

  stlnpp::as.stlpp(x$df$x, x$df$y, x$df$t, x$L)
}
