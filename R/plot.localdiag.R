#' Plot of the local diagnostics' result on a spatio-temporal point process model
#'
#' @description
#' This function plots the result of the local diagnostics performed with
#' \link{localdiag} on either a \code{stp} or \code{stlp} object.
#'  It highlights the  points of the analysed spatio-temporal point pattern \code{X}
#'  which are identified as outlying by the
#'  previously performed local diagnostics; the remaining points of \code{X}
#'  are also represented.
#'
#'  It also shows the underlying linear network, if the local diagnostics has been applied
#'  to point patterns occurring on the same linear network, that is, if \link{localdiag}
#'  has been applied to a \code{stlp} object.
#'
#' @param x A \code{localdiag} object
#' @param marg Default to \code{TRUE}. If \code{marg = F}, only the
#' spatio-temporal point pattern is plotted
#' @param col Color of the outlying points
#' @param col2 Color of the network (if applicable)
#' @param cols Color of the non-outlying points
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo and Giada Adelfio
#'
#' @seealso
#' \link{infl}, \link{print.localdiag}, \link{summary.localdiag}
#'
#'
#' @examples
#'
#' \dontrun{
#'
#' inh <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, 
#'              par = c(.3, 6), seed = 2)
#' 
#' mod1 <- stppm(inh, formula = ~ 1)
#' 
#' resmod1 <- localdiag(inh, mod1$l, p = .9)
#' 
#' plot(resmod1)
#' plot(resmod1, marg = FALSE)
#'
#' }
#'
#'
#' @references
#' Adelfio, G., Siino, M., Mateu, J., and Rodríguez-Cortés, F. J. (2020). Some properties of local weighted second-order statistics for spatio-temporal point processes. Stochastic Environmental Research and Risk Assessment, 34(1), 149-168.
#'
#' D’Angelo, N., Adelfio, G.  and Mateu, J. (2022) Local inhomogeneous second-order characteristics for spatio-temporal point processes on linear networks. Stat Papers. https://doi.org/10.1007/s00362-022-01338-4
#'
plot.localdiag <- function(x, marg = TRUE, col = "grey", col2 = "red",
                           cols = "lightgrey", ...){
  if (!inherits(x, c("localdiag"))) stop("x should be from class localdiag")
  

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  mask <- vector(length = nrow(x$X$df))
  mask[x$ids] <- "1"
  mask[- x$ids] <- "0"

  minx <- min(x$x2)
  maxx <- max(x$x2)
  breaksx <- seq(minx, maxx, length = 40)

  if(marg == T){
    par(mfrow = c(2, 2), mar = par("mar") - 1)
  } else {
    par(mfrow = c(1, 2))
  }

  par(mar = c(5, 4, 4, 2) + 0.1)

  hist(x$x2, breaks = breaksx, main =quote(Obj == chi[i]^2),
       xlab = " ")

  q <- quantile(x$x2, x$p)

  abline(v = q, col = 2)

  rug(x$x2)

  par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 1))

  plot3D::scatter3D(x$X$df$x[mask == "0"], x$X$df$y[mask == "0"], x$X$df$t[mask == "0"],
                    theta = - 45, phi = 20, pch = 20,
                    ticktype = "detailed",
                    col = col,
                    main = paste("Points with Obj >", round(q, 3)),
                    xlim = range(x$X$df$x), ylim = range(x$X$df$y), zlim = range(x$X$df$t))

  par(mar = c(5, 4, 4, 2) + 0.1)

  plot3D::scatter3D(x$X$df$x[mask == "1"], x$X$df$y[mask == "1"], x$X$df$t[mask == "1"],
                    add = TRUE, colkey = FALSE, pch = 20,
                    col = col2)

  if(marg == TRUE){
    if(inherits(x$X, "stlp")){
      plot(x$X$L, col = cols, main = "Spatial point pattern",
           xlab = "", ylab = "")
    } else {
      plot(x$X$df$x, x$X$df$y, type = "n",xlab = "", ylab = "",
           main = "Spatial point pattern")
    }
    title(xlab = "x", ylab = "y")
    axis(1)
    axis(2)
    box()

    points(x$X$df$x[mask == 0], x$X$df$y[mask == 0], col = col, pch = 20)
    points(x$X$df$x[mask == 1], x$X$df$y[mask == 1], col = col2, pch = 20)

    plot(x$X$df$t[mask == 0], rep(1, length(x$X$df$t[mask == 0])),
         main = "Temporal point pattern", pch = 20,
         xlab = "t", ylab = "", yaxt = "n", col = col)
    points(x$X$df$t[mask == 1], rep(1, length(x$X$df$t[mask == 1])),
           col = col2, pch = 20)
  }
}
