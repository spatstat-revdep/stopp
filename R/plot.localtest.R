#' Plot of the result of the local permutation test
#'
#' @description
#' This function plots the result of the local permutation test performed with
#' \link{localtest} on either a \code{stp} or \code{stlp} object. It highlights the points of the background pattern \code{X},
#'  which exhibit local differences in the second-order
#'  structure with respect to \code{Z}, according to the previously performed test.
#'  The remaining points of \code{X} are also represented.
#'
#'  It also shows the underlying linear network, if the local test has been applied
#'  to point patterns occurring on the same linear network, that is, if \link{localtest}
#'  has been applied to a \code{stlp} object.
#'
#'
#' @param x An object of class  \code{localtest}
#' @param col Color of the significant points
#' @param cols Color of the linear network. If applicable.
#' @param col2 Color of the non-significant points
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{localtest},  \link{print.localtest},
#' \link{summary.localtest}
#'
#'
#' @examples
#'
#' \dontrun{
#'
#' # Euclidean
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
#'
#' plot(test)
#'
#' # Linear networks
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
#' plot(test)
#'
#'}
#'
#'
#' @references
#' D’Angelo, N., Adelfio, G., and Mateu, J. (2021). Assessing local differences between the spatio-temporal second-order structure of two point patterns occurring on the same linear network. Spatial Statistics, 45, 100534.
#'
#' Siino, M., Rodríguez‐Cortés, F. J., Mateu, J. ,and Adelfio, G. (2018). Testing for local structure in spatiotemporal point pattern data. Environmetrics, 29(5-6), e2463.
#'
plot.localtest <- function(x, col = "grey", cols = "lightgrey", col2 = "red", ...){
  if(!any(class(x) == "localtest")) stop("class(x) must be localtest")

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow = c(1, 3))

  plot3D::scatter3D(x$Xnosig$df$x, x$Xnosig$df$y, x$Xnosig$df$t,
                    theta = - 45, phi = 20, pch = 16,
                    ticktype = "detailed",
                    col = col,
                    xlab = "x", ylab = "y", zlab = "t",
                    main = c("Spatio-temporal point pattern"),
                    xlim = range(x$X$df$x), ylim = range(x$X$df$y), zlim = range(x$X$df$t))
  plot3D::scatter3D(x$Xsig$df$x, x$Xsig$df$y, x$Xsig$df$t,
                    theta = - 45, phi = 20, pch = 19,
                    ticktype = "detailed",
                    col = col2, add = T)

  if(inherits(x$X, "stlp")){
    plot(x$Xsig$L, col = cols, main = "Spatial point pattern")
    title(xlab = "x", ylab = "y")
    axis(1)
    axis(2)
    box()
  } else {
    plot(x$X$df$x, x$X$df$y, type = "n", xlab = "x", ylab = "y",
         main = "Spatial point pattern")
  }
  points(x$Xnosig$df$x, x$Xnosig$df$y, col = col, pch = 16)
  points(x$Xsig$df$x, x$Xsig$df$y, col = col2, pch = 19)

  plot(x$Xnosig$df$t, rep(1, length(x$Xnosig$df$t)),
       main = "Temporal point pattern", pch = 16,
       xlab = "t", ylab = "", yaxt = "n", col = col)
  points(x$Xsig$df$t, rep(1, length(x$Xsig$df$t)),
         col = col2, pch = 19)
}
