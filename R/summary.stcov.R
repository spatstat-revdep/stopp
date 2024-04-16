#' Summary of a \code{stcov} object
#'
#' It prints the summary statistics of the spatio-temporal coordinates and the 
#' covariates values
#' of  the spatio-temporal covariate
#' stored in the \code{stcov} object.
#'
#' @param object An object of class \code{stcov}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#' 
#' @examples
#'
#' set.seed(2)
#' df <- data.frame(runif(100), runif(100), runif(100), rpois(100, 15))
#'
#' cov <- stcov(df, interp = FALSE)
#' summary(cov)
#'
summary.stcov <- function (object, ...) 
{
  summary(object$df)
}
