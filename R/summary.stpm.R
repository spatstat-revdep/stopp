#' Summary of a \code{stpm} object
#'
#' It prints the summary statistics of the spatio-temporal coordinates and the marks
#' of  the spatio-temporal point pattern
#' stored in the \code{stpm} object.
#'
#' @param object An object of class \code{stpm}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#' 
#' @examples
#'
#' set.seed(2)
#' df <- data.frame(cbind(runif(100), runif(100), runif(100), rpois(100, 15),
#' rpois(100, 30)))
#'
#' stpm1 <- stpm(df)
#'
#' summary(stpm1)
#' 
#' ## Categorical marks
#' 
#' set.seed(2)
#' dfA <- data.frame(x = runif(100), y = runif(100), t = runif(100), 
#'                   m1 = rnorm(100), m2 = rep(c("C"), times = 100))
#' dfB <- data.frame(x = runif(50), y = runif(50), t = runif(50), 
#'                   m1 = rnorm(25), m2 = rep(c("D"), times = 50))
#' 
#' stpm2 <- stpm(rbind(dfA, dfB), names = c("continuous", "dichotomous"))
#' 
#' summary(stpm2)
#'
summary.stpm <- function (object, ...) 
{
  summary(object$df)
}
