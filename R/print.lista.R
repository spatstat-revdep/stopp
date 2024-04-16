#' Print a \code{lista} object
#'
#' It prints the main information on the local network summary statistics
#' stored in a \code{lista} object.
#'
#' @param x An object of class \code{lista}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @examples
#'
#' set.seed(2)
#' df_net <- data.frame(x = runif(25, 0, 0.85), y = runif(25, 0, 0.85), t = runif(25))
#' stlp1 <- stp(df_net, L = chicagonet)
#' lambda <- rep(diff(range(stlp1$df$x)) * diff(range(stlp1$df$y))
#'  * diff(range(stlp1$df$t)) / spatstat.geom::volume(stlp1$L),
#' nrow(stlp1$df))
#'
#' k <- localSTLKinhom(stlp1, lambda = lambda, normalize = TRUE)
#' 
#' k
#'
#'
print.lista <- function (x, ...) 
{

  cat("LISTA functions on a linear network \n")
  {
    cat("computed at", paste0(length(x), " ", "points"), "\n")
  }
}
