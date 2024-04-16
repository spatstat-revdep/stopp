#' Print a \code{stlpm} object
#'
#' It prints the main information on the spatio-temporal point pattern
#' stored in the \code{stlpm} object: the number of points; the enclosing spatial window;
#' the temporal time period; information on marks.
#'
#' @param x An object of class \code{stlpm}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @examples
#'
#' set.seed(2)
#' df <- data.frame(x = runif(100, 0, 0.8), y = runif(100, 0, 0.8), t = runif(100), m = rpois(100, 15))
#'
#' stlpm1 <- stpm(df, L = chicagonet)
#'
#' stlpm1
#' 
#'
#'
print.stlpm <- function (x, ...) 
{

  cat("Spatio-temporal marked point pattern \n")
  {
    cat(paste0(nrow(x$df), " ", "points"), "\n")
  }
  print(spatstat.geom::owin(range(x$df$x), range(x$df$y)), 
        prefix = "Enclosing window: ")
  cat(paste0("Time period: [", round(range(x$df$t)[1], 3), 
             ", ", round(range(x$df$t)[2], 3), "]"), "\n")
  {
    cat(paste("with", ncol(x$df) - 3, "marks"), "\n")
  }
  print(x$L)
  cat(paste0("Time period: [",round(range(x$df$t)[1], 3),
             ", ", round(range(x$df$t)[2], 3),"]"),"\n")
}
