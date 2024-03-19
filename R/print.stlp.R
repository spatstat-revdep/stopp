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
#' set.seed(2)
#' df_net <- data.frame(cbind(runif(100, 0, 0.85), runif(100, 0, 0.85), runif(100)))
#'
#' stlp1 <- stp(df_net, L = chicagonet)
#' stlp1
#' }
#'
#' @seealso
#' \link{stp}, \link{plot.stlp}, \link{summary.stlp}
#'
#'
print.stlp <- function(x, ...)
{
  if (!inherits(x, c("stlp"))) stop("X should be from class stlp")
  
  cat("Spatio-temporal point pattern on a linear network \n")
  {cat(paste0(nrow(x$df)," ", "points"),"\n")}
  print(x$L)
  cat(paste0("Time period: [",round(range(x$df$t)[1], 3),
             ", ", round(range(x$df$t)[2], 3),"]"),"\n")
}
