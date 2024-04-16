#' Summary of a \code{stp} object
#'
#' It prints the summary statistics of the spatial and temporal coordinates
#' of  the spatio-temporal point pattern
#' stored in the \code{stp} object.
#'
#' @param object An object of class \code{stp}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stp}, \link{print.stp}, \link{plot.stp}
#'
#' @examples
#'
#'
#' set.seed(2)
#' df <- data.frame(cbind(runif(100), runif(100), runif(100)))
#' 
#' stp1 <- stp(df)
#' summary(stp1)
#'
#'
#'
summary.stp <- function(object, ...)
{
  summary(object$df)
}


