#' Plot the coefficients of a fitted local spatio-temporal Poisson process model
#'
#' The function plots the local estimates of a fitted local spatio-temporal Poisson process model.
#'
#' @param x An object of class \code{locstppm}
#' @param par Default to \code{TRUE}.
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{locstppm}, \link{print.locstppm}, \link{summary.locstppm},
#' \link{localsummary.locstppm}, \link{plot.locstppm}
#'
#'
#'
#' @examples
#' \dontrun{
#' inh <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, 
#'              par = c(0.005, 5), seed = 2)
#' inh_local <- locstppm(inh, formula = ~ x)
#' 
#' localplot.locstppm(inh_local)
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
localplot.locstppm <- function(x, par = TRUE){
  if (!inherits(x, c("locstppm"))) stop("x should be from class locstppm")
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  nn <- length(names(x$IntCoefs))

  if(par == TRUE){
    par(mfrow = c(1, nn))
  } else{
    par(mfrow = c(1, 1))
    par(ask = TRUE)
  }

  par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 0))
  plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                    theta = - 45, phi = 20,
                    col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                        range = range(x$IntCoefs_local[, 1])),
                               "stuff")$outputs,
                    ticktype = "detailed", pch = 19,
                    colvar = x$IntCoefs_local[, 1],
                    xlab="x",ylab="y",zlab="t",
                    main = c("Intercept"))
  par(mar = c(5, 4, 4, 2) + 0.1)
  if(par != TRUE){
    par(ask = FALSE)
  }

  if(nn > 1){
    for(i in 2:nn){

      id <- x$IntCoefs_local[, i]

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
                        main = names(x$IntCoefs_local)[i])
      par(mar = c(5, 4, 4, 2) + 0.1)
      if(par != TRUE){
        par(ask = FALSE)
      }

    }
  }

}

