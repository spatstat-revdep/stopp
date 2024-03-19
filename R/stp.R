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
#' \link{stppm},
#' \link{print.stp}, \link{summary.stp}, \link{plot.stp},
#' \link{print.stlp},
#' \link{summary.stlp}, \link{plot.stlp}
#'
#' @examples
#'
#' \dontrun{
#'
#' set.seed(2)
#' df <- data.frame(runif(100), runif(100), runif(100))
#' 
#' stp1 <- stp(df)
#' 
#' set.seed(2)                       
#' df_net <- data.frame(runif(100, 0, 0.85), runif(100, 0, 0.85), runif(100))
#' 
#' stlp1 <- stp(df_net, L = chicagonet)
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
    suppressWarnings(lpp_new <- spatstat.linnet::lpp(df, L))
    df_new <- data.frame(x = lpp_new$data[, 1], y = lpp_new$data[, 2 ], t = df[, 3])
    out <- list(df = df_new, L = L)
    class(out) <- "stlp"
  }
  
  return(out)
  
}




















