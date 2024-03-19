#' Print a \code{stpm} object
#'
#' It prints the main information on the spatio-temporal point pattern
#' stored in the \code{stpm} object: the number of points; the enclosing spatial window;
#' the temporal time period; information on marks.
#'
#' @param x An object of class \code{stpm}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @examples
#'
#' \dontrun{
#' set.seed(2)
#' df <- data.frame(cbind(runif(100), runif(100), runif(100), rpois(100, 15),
#' rpois(100, 30)))
#'
#' stpm1 <- stpm(df)
#'
#' summary(stpm1)
#' }
#'
#'
print.stpm <- function (x, ...) 
{
  if (!inherits(x, c("stpm"))) stop("X should be from class stpm")
  
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
}
