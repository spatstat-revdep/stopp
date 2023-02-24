#' Summary of a \code{stp} object
#'
#' It prints the summary statistics of the spatial and temporal coordinates
#' of  the spatio-temporal point pattern
#' stored in the \code{stp} object.
#'
#' @param object An object of class \code{stp}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stp}, \link{print.stp}, \link{plot.stp},
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
#' summary(stp1)
#'
#' #       x                  y                 t
#' # Min.   :0.001137   Min.   :0.01553   Min.   :0.004254
#' # 1st Qu.:0.290510   1st Qu.:0.27316   1st Qu.:0.266398
#' # Median :0.559250   Median :0.53679   Median :0.524548
#' # Mean   :0.536354   Mean   :0.52031   Mean   :0.513976
#' # 3rd Qu.:0.793498   3rd Qu.:0.78442   3rd Qu.:0.758390
#' # Max.   :0.993378   Max.   :0.99604   Max.   :0.996996
#'
#' }
#'
#'
summary.stp <- function(object, ...)
{
  if(!any(class(object) == "stp")) stop("class(object) must be stp")

  summary(object$df)
}


