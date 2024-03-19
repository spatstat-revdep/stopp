#' Print a \code{stp} object
#'
#' It prints the main information on the spatio-temporal point pattern
#' stored in the \code{stp} object: the number of points; the enclosing spatial window;
#' the temporal time period.
#'
#' @param x An object of class \code{stp}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stp}, \link{summary.stp}, \link{plot.stp}
#'
#' @examples
#'
#' \dontrun{
#'
#' set.seed(2)
#' df <- data.frame(cbind(runif(100), runif(100), runif(100)))
#' 
#' stp1 <- stp(df)
#' stp1
#'
#' }
#'
#'
#'
#'
print.stp <- function(x, ...)
{
  if (!inherits(x, c("stp"))) stop("X should be from class stp")
  
  cat("Spatio-temporal point pattern \n")
  {cat(paste0(nrow(x$df)," ", "points"),"\n")}
  print(spatstat.geom::owin(range(x$df$x), range(x$df$y)), prefix = "Enclosing window: ")
  cat(paste0("Time period: [",round(range(x$df$t)[1], 3),
             ", ", round(range(x$df$t)[2], 3),"]"),"\n")
}
