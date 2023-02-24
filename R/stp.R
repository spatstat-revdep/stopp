#' Create \code{stp} and \code{stlp} objects for point patterns storage
#'
#' This function creates a \code{stp} object as a dataframe with three columns:
#' \code{x}, \code{y}, and \code{t}.
#' If also the linear network \code{L}, of class linnet, is provided, a \code{stlp}
#' object is created instead.
#'
#' @seealso \link{summary.stp}, \link{print.stp}, \link{plot.stp}
#'
#' @param df A matrix with three columns, containing to two space and the
#' temporal coordinates
#' @param L Optional. The linear network of class \code{linnet}.
#'  If provided, the function returns a \code{stlp} object.
#'
#' @return An \code{stp} or \code{stlpp} object, depending on whether or not an object of class
#' linnet is provided for the \code{L} argument.
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stppm}, \link{AIC.stppm}, \link{BIC.stppm},
#' \link{print.stp}, \link{summary.stp}, \link{plot.stp},
#' \link{as.stp}, \link{as.stpp}, \link{print.stlp},
#' \link{summary.stlp}, \link{plot.stlp}, \link{as.stlp},  \link{as.stlpp}
#'
#' @examples
#'
#' \dontrun{
#'
#' # stp
#'
#' set.seed(12345)
#' rpp1 <- rpp(lambda = 200, replace = FALSE)
#' df0 <- cbind(rpp1$xyt[, 1], rpp1$xyt[, 2], rpp1$xyt[, 3])
#' stp1 <- stp(df0)
#'
#'
#'
#' #stlp
#'
#' set.seed(12345)
#' stlpp1 <- rpoistlpp(.2, a = 0, b = 5, L = easynet)
#' df0 <- cbind(stlpp1$data$x, stlpp1$data$y, stlpp1$data$t)
#' L0 <- stlpp1$domain
#' stlp1 <- stp(df0, L0)
#'
#' }
#'
#'
stp <- function(df, L){

  colnames(df) <- c("x", "y", "t")
  df <- data.frame(df)

  if(missing(L)){
    out <- list(df = df)
    class(out) <- "stp"
  } else {
    out <- list(df = df, L = L)
    class(out) <- "stlp"
  }

  return(out)

}





















