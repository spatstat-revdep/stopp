#' Summary of a \code{stlpm} object
#'
#' It prints the summary statistics of the spatio-temporal coordinates and the marks
#' of  the spatio-temporal point pattern
#' stored in the \code{stlpm} object.
#'
#' @param object An object of class \code{stlpm}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#' 
#' @examples
#'
#' set.seed(2)
#' df <- data.frame(x = runif(100, 0, 0.8), y = runif(100, 0, 0.8),
#'  t = runif(100), m = rpois(100, 15))
#'
#' stlpm1 <- stpm(df, L = chicagonet)
#'
#' summary(stlpm1)
#'
summary.stlpm <- function (object, ...) 
{
  summary(object$df)
}
