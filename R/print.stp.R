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
#' \link{stp}, \link{summary.stp}, \link{plot.stp},
#' \link{as.stp}, \link{as.stpp}
#'
#' @examples
#'
#' \dontrun{
#'
#' set.seed(12345)
#' rpp1 <- rpp(lambda = 200, replace = FALSE)
#' df0 <- cbind(rpp1$xyt[, 1], rpp1$xyt[, 2], rpp1$xyt[, 3])
#' stp1 <- stp(df0)
#'
#' stp1
#'
#' # Spatio-temporal point pattern
#' # 208 points
#' # Enclosing window: rectangle = [0.0011366, 0.9933775] x [0.0155277,
#' #                                                          0.9960438] units
#' # Time period: [0.004, 0.997]
#'
#' }
#'
#'
#'
#'
print.stp <- function(x, ...)
{

  if(!any(class(x) == "stp")) stop("class(x) must be stp")

  cat("Spatio-temporal point pattern \n")
  {cat(paste0(nrow(x$df)," ", "points"),"\n")}
  print(spatstat.geom::owin(range(x$df$x), range(x$df$y)), prefix = "Enclosing window: ")
  cat(paste0("Time period: [",round(range(x$df$t)[1], 3),
             ", ", round(range(x$df$t)[2], 3),"]"),"\n")
}
