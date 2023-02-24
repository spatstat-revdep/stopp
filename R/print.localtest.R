#' Print of the result of the permutation local test
#'
#' It prints the main information on the result of the local permutation test
#'  performed with \link{localtest}  on either a \code{stp} or \code{stlp} object:
#'  whether the local test was run on point patterns lying on a linear network or not;
#'  the number of points in the background \code{X} and alternative \code{Z} patterns;
#'  the number of points in \code{X} which exhibit local differences in the second-order
#'  structure with respect to \code{Z}, according to the performed test.
#'
#' @param x An object of class  \code{localtest}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{localtest}, \link{summary.localtest},
#' \link{plot.localtest}
#'
#' @examples
#'
#' \dontrun{
#'
#' # background pattern
#' set.seed(12345)
#' X <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(.05, 4),
#'            nsim = 1, seed = 2, verbose = T)
#'
#' # alternative pattern
#' set.seed(12345)
#' Z <- rstpp(lambda = 25, nsim = 1, seed = 2, verbose = T)
#'
#' # run the local test
#' test <- localtest(X, Z, method = "K", k = 9)
#'
#' test
#'
#' # Test for local differences between two
#' # spatio-temporal point patterns
#' # --------------------------------------
#' #   Backgound pattern X: 17
#' # Alternative pattern Z: 20
#' #
#' # 1 significant points at alpha = 0.05
#'
#'
#'# Linear networks
#'
#' # fix the linear network
#' L0 = domain(chicago)
#'
#' # background pattern
#' set.seed(12345)
#' X <- retastlp(cat = NULL, params = c(0.078915 / 1.95, 0.003696,  0.013362,  1.2,
#'                                         0.424466,  1.164793),
#'                  betacov = 0.5, m0 = 2.5, b = 1.0789, tmin = 0, t.lag = 200,
#'                  xmin = 600, xmax = 2200, ymin = 4000, ymax = 5300,
#'                  iprint = TRUE, covdiag = FALSE, covsim = FALSE, L = L0)
#'
#' # alternative pattern, on the same linear network
#' l <- 20 / (volume(domain(chicago)) * (200 - 25))
#' set.seed(12345)
#' stlppPOIS <- rpoistlpp(lambda = l, a = 25, b = 200, L = L0)
#' Z <- as.stlp(stlppPOIS)
#'
#' # run the local test
#' test <- localtest(X, Z, method = "K", k = 9)
#'
#' test
#'
#' # Test for local differences between two
#' # spatio-temporal point patterns on a linear network
#' # --------------------------------------------------
#' #   Backgound pattern X: 31
#' # Alternative pattern Z: 22
#' #
#' # 19 significant points at alpha = 0.05
#'
#'}
#'
#'
#' @references
#' D’Angelo, N., Adelfio, G., and Mateu, J. (2021). Assessing local differences between the spatio-temporal second-order structure of two point patterns occurring on the same linear network. Spatial Statistics, 45, 100534.
#'
#' Siino, M., Rodríguez‐Cortés, F. J., Mateu, J. ,and Adelfio, G. (2018). Testing for local structure in spatiotemporal point pattern data. Environmetrics, 29(5-6), e2463.
#'
#'
print.localtest <- function(x, ...){
  if(!any(class(x) == "localtest")) stop("class(x) must be localtest")
  if(inherits(x$Xsig, "stlp")){
    cat("Test for local differences between two \n")
    cat("spatio-temporal point patterns on a linear network \n")
    cat("--------------------------------------------------\n")
  } else {
    cat("Test for local differences between two \n")
    cat("spatio-temporal point patterns \n")
    cat("--------------------------------------\n")
  }
  cat(paste("Backgound pattern X:",nrow(x$X$df), " \n"))
  cat(paste("Alternative pattern Z:",nrow(x$Z$df), " \n  \n"))
  cat(paste(nrow(x$Xsig$df), "significant points at alpha =", x$alpha))
}

