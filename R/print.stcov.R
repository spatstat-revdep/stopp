#' Print a \code{stcov} object
#'
#' It prints the main information on the spatio-temporal covariate
#' stored in the \code{stcov} object: the number of points; the enclosing spatial window;
#' the temporal time period; information on the covariate values.
#'
#' @param x An object of class \code{stcov}
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
#' cov
#'
#'
print.stcov <- function (x, ...) 
{

  cat("Spatio-temporal covariate \n")
  {
    cat("observed at", paste0(nrow(x$df), " ", "sites"), "\n")
  }
  print(spatstat.geom::owin(range(x$df$x), range(x$df$y)), 
        prefix = "Enclosing window: ")
  cat(paste0("Time period: [", round(range(x$df$t)[1], 3), 
             ", ", round(range(x$df$t)[2], 3), "]"), "\n")
}
