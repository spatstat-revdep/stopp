#' Plot the coefficients of a fitted local LGCP model
#'
#' The function plots the local estimates.
#' In the case of local covariance parameters, the function
#' displays the local estimates of the chosen covariance function.
#'
#' @param x An object of class \code{stlgcppm}
#' @param par Default to \code{TRUE}.
#'
#' @export
#'
#' @author Nicoletta D'Angelo and Giada Adelfio
#'
#' @seealso
#' \link{stlgcppm}, \link{print.stlgcppm}, \link{summary.stlgcppm},
#' \link{localsummary.stlgcppm}, \link{plot.stlgcppm}
#'
#'
#'
#' @examples
#' \dontrun{
#' catsub <- stp(greececatalog$df[1:100, ])
#' 
#' lgcp_loc <- stlgcppm(catsub, formula = ~ x, first = "local")
#'
#' localplot.stlgcppm(lgcp_loc)
#'
#'}
#'
#'
#'
#'
#' @references
#' D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
#'
#'
localplot.stlgcppm <- function(x, par = TRUE){
  if (!inherits(x, c("stlgcppm"))) stop("x should be from class stlgcppm")
  
  if(inherits(x$IntCoefs, "numeric") & inherits(x$CovCoefs, "numeric")){
    stop("No local parameters to plot")
  }
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if(!inherits(x$IntCoefs, "numeric")){
    nn <- length(names(x$IntCoefs))

    if(par == T){
      par(mfrow = c(1, nn))
    } else{
      par(mfrow = c(1, 1))
      par(ask = TRUE)
    }

    par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 0))
    plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                      theta = - 45, phi = 20,
                      col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                          range = range(x$IntCoefs$`(Intercept)`)),
                                 "stuff")$outputs,
                      ticktype = "detailed", pch = 19,
                      colvar = x$IntCoefs$`(Intercept)`,
                      xlab="x",ylab="y",zlab="t",
                      main = c("Intercept"))
    par(mar = c(5, 4, 4, 2) + 0.1)
    if(par != TRUE){
      par(ask = FALSE)
    }

    if(nn > 1){
      for(i in 2:nn){

        id <- x$IntCoefs[, i]

        if(par != TRUE){
          par(ask = TRUE)
        }
        par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, -1))
        plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                          theta = - 45, phi = 20,
                          col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                              range = range(id)),
                                     "stuff")$outputs,
                          ticktype = "detailed", pch = 19,
                          colvar = id,
                          xlab="x",ylab="y",zlab="t",
                          main = names(x$IntCoefs)[i])
        par(mar = c(5, 4, 4, 2) + 0.1)
        if(par != TRUE){
          par(ask = FALSE)
        }

      }
    }


  }


  if(!inherits(x$CovCoefs, "numeric")){

    if(ncol(x$CovCoefs) == 3){
      if(par == T){
        par(mfrow = c(1, 3))
      } else{
        par(mfrow = c(1, 1))
      }

      if(par != TRUE){
        par(ask = TRUE)
      }
      par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 1))
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$sigma)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$sigma,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(sigma))

      if(par != TRUE){
        par(ask = TRUE)
      }
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$alpha)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$alpha,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(alpha))
      if(par != TRUE){
        par(ask = TRUE)
      }
      par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 0))
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$beta)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$beta,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(beta))
      par(mar = c(5, 4, 4, 2) + 0.1)
      if(par != TRUE){
        par(ask = FALSE)
      }

    } else {
      if(par == TRUE){
        par(mfrow = c(2, 3))
      } else{
        par(mfrow = c(1, 1))
      }


      if(par != TRUE){
        par(ask = TRUE)
      }
      par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 1))
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$sigma)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$sigma,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(sigma))

      if(par != TRUE){
        par(ask = TRUE)
      }
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$alpha)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$alpha,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(alpha))
      par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 0))

      if(par != TRUE){
        par(ask = TRUE)
      }
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$beta)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$beta,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(beta))

      if(par != TRUE){
        par(ask = TRUE)
      }
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$gamma_s)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$gamma_s,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(gamma))

      if(par != TRUE){
        par(ask = TRUE)
      }
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$gamma_t)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$gamma_t,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(gamma[t]))
      par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 0))

      if(par != TRUE){
        par(ask = TRUE)
      }
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(x$CovCoefs$delta)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 19,
                        colvar = x$CovCoefs$delta,
                        xlab="x",ylab="y",zlab="t",
                        main = expression(delta))
      par(mar = c(5, 4, 4, 2) + 0.1)
      if(par != TRUE){
        par(ask = FALSE)
      }
    }
  }

}

