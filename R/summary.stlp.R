#' Summary of a \code{stlp} object
#'
#' It prints the main information on the spatio-temporal point pattern on a linear
#' network stored in the \code{stlp} object: the number of points;
#' vertices and lines of the linear network; the enclosing spatial window;
#' the temporal time period.
#'
#' @param object An object of class \code{stlp}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stp}, \link{plot.stlp}, \link{print.stlp},
#' \link{as.stlpp}, \link{as.stlp}
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
#' summary(stlp1)
#'
#' }
#'
summary.stlp <- function(object, ...)
{
  if(!any(class(object) == "stlp")) stop("class(object) must be stlp")

  cat("Spatio-temporal point pattern on a linear network \n")
  {cat(paste0(nrow(object$df)," ", "points"),"\n")}
  print(object$L)
  cat(paste0("Time period: [",round(range(object$df$t)[1], 3),
             ", ", round(range(object$df$t)[2], 3),"]"),"\n")
}
