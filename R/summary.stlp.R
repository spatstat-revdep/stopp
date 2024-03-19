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
#' \link{stp}, \link{plot.stlp}, \link{print.stlp}
#'
#' @examples
#'
#' \dontrun{
#'
#' set.seed(2)                       
#' df_net <- data.frame(cbind(runif(100, 0, 0.85), runif(100, 0, 0.85), runif(100)))
#'
#' stlp1 <- stp(df_net, L = chicagonet)
#' summary(stlp1)
#'
#' }
#'
summary.stlp <- function(object, ...)
{
  if (!inherits(object, c("stlp"))) stop("X should be from class stlp")
  
  cat("Spatio-temporal point pattern on a linear network \n")
  {cat(paste0(nrow(object$df)," ", "points"),"\n")}
  print(object$L)
  cat(paste0("Time period: [",round(range(object$df$t)[1], 3),
             ", ", round(range(object$df$t)[2], 3),"]"),"\n")
}
