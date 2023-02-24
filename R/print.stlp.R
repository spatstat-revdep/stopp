#' Print a \code{stlp} object
#'
#' It prints the main information on the spatio-temporal point pattern on a linear
#' network stored in the \code{stlp} object: the number of points;
#' vertices and lines of the linear network; the enclosing spatial window;
#' the temporal time period.
#'
#' @param x An object of class \code{stlp}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @examples
#'
#' \dontrun{
#'
#' set.seed(12345)
#' stlpp1 <- rpoistlpp(.2, a = 0, b = 5, L = easynet)
#' df0 <- cbind(stlpp1$data$x, stlpp1$data$y, stlpp1$data$t)
#' L0 <- stlpp1$domain
#' stlp1 <- stp(df0, L0)
#'
#' stlp1
#'
#' # Spatio-temporal point pattern on a linear network
#' # 43 points
#' # Linear network with 19 vertices and 26 lines
#' # Enclosing window: rectangle = [-0.01, 5.1] x [-0.01, 5.1] units
#' # Time period: [0.043, 4.93]
#'
#' }
#'
#' @seealso
#' \link{stp}, \link{plot.stlp}, \link{summary.stlp},
#' \link{as.stlpp}, \link{as.stlp}
#'
#'
print.stlp <- function(x, ...)
{

  if(!any(class(x) == "stlp")) stop("class(x) must be stlp")

  cat("Spatio-temporal point pattern on a linear network \n")
  {cat(paste0(nrow(x$df)," ", "points"),"\n")}
  print(x$L)
  cat(paste0("Time period: [",round(range(x$df$t)[1], 3),
             ", ", round(range(x$df$t)[2], 3),"]"),"\n")
}
