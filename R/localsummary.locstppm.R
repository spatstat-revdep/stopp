#' Summary plots of the fitted coefficient of a local spatio-temporal Poisson process model
#'
#' The function breaks up the contribution of the local estimates
#' to the fitted intensity, by plotting the overall intensity and the
#' density kernel smoothing of some artificial intensities, obtained by
#' imputing the quartiles of the local parameters' distributions.
#'
#' @param x An object of class \code{locstppm}
#' @param scaler Optional. Controls the value for a scalar representation of the
#'  spatial scale of the data.
#'  Either a character string, \code{"silverman"} (default), \code{"IQR"},
#'   \code{"sd"}, or \code{"var"};
#'  or positive numeric value(s). See \link[sparr]{OS}.
#' @param do.points Add points to plot
#' @param print.bw It prints the estimated oversmoothing (\link[sparr]{OS}) bandwidth selector
#' @param zap Noise threshold factor (default to 0.00001). A numerical value greater than or equal to 1.
#'  If the range of pixel values is less than \code{zap * .Machine$double.eps},
#'   the image will be treated as constant. This avoids displaying images which
#'   should be constant but contain small numerical errors.
#' @param par Default to \code{TRUE}.
#'
#' @export
#'
#' @author Nicoletta D'Angelo and Giada Adelfio
#'
#' @seealso
#' \link{locstppm}, \link{print.locstppm}, \link{summary.locstppm},
#' \link{plot.locstppm}
#'
#'
#' @examples
#' \dontrun{
#' set.seed(2)
#' pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6),
#'              nsim = 1, seed = 2, verbose = T)
#' inh00_local <- locstppm(pin, formula = ~ 1)
#' inh01_local <- locstppm(pin, formula = ~ x)
#'
#' localsummary.locstppm(inh00_local)
#' localsummary.locstppm(inh01_local)
#'
#'}
#'
#'
#' @references
#' D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
#'
#' Davies, T.M. and Hazelton, M.L. (2010). Adaptive kernel estimation of spatial relative risk, Statistics in Medicine, 29(23) 2423-2437.
#'
#' Terrell, G.R. (1990). The maximal smoothing principle in density estimation, Journal of the American Statistical Association, 85, 470-477.
#'
localsummary.locstppm <- function(x,
                                  scaler = c("silverman", "IQR", "sd", "var"),
                                  do.points = TRUE,
                                  print.bw = FALSE,
                                  zap = 0.00001,
                                  par = TRUE){
  if(!any(class(x) == "locstppm")) stop("class(x) must be locstppm")

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  mod_Min <- mod_1stQu <- mod_Median <- mod_3rdQu <- mod_Max <- x$mod_global

  mod_Min$coefficients <- apply(x$IntCoefs_local, 2, summary)[1, ]
  mod_1stQu$coefficients <- apply(x$IntCoefs_local, 2, summary)[2, ]
  mod_Median$coefficients <- apply(x$IntCoefs_local, 2, summary)[3, ]
  mod_3rdQu$coefficients <- apply(x$IntCoefs_local, 2, summary)[5, ]
  mod_Max$coefficients <- apply(x$IntCoefs_local, 2, summary)[6, ]

  mark_int <- x$l_local
  l_Min <- predict(mod_Min, newdata = x$newdata)
  l_1stQu <- predict(mod_1stQu, newdata = x$newdata)
  l_Median <- predict(mod_Median, newdata = x$newdata)
  l_3rdQu <- predict(mod_3rdQu, newdata = x$newdata)
  l_Max <- predict(mod_Max, newdata = x$newdata)

  ppx_int <- spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = mark_int,
                                window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y)))
  ppx_Min <- spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = l_Min,
                                window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y)))
  ppx_1stQu <- spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = l_1stQu,
                                  window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y)))
  ppx_Median <- spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = l_Median,
                                   window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y)))
  ppx_3rdQu <- spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = l_3rdQu,
                                  window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y)))
  ppx_Max <- spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = l_Max,
                                window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y)))

  sig <- sparr::OS(unmark(ppx_int), scaler = scaler)

  s_int <- spatstat.explore::Smooth(ppx_int, sigma = sig)
  s_Min <- spatstat.explore::Smooth(ppx_Min, sigma = sig)
  s_1stQu <- spatstat.explore::Smooth(ppx_1stQu, sigma = sig)
  s_Median <- spatstat.explore::Smooth(ppx_Median, sigma = sig)
  s_3rdQu <- spatstat.explore::Smooth(ppx_3rdQu, sigma = sig)
  s_Max <- spatstat.explore::Smooth(ppx_Max, sigma = sig)

  g00 <- spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                  range = c(min(s_int, s_Min$v[1], s_1stQu$v[1], s_Median$v[1],
                                                s_3rdQu$v[1], s_Max$v[1]),
                                            max(s_int, s_Min$v[1], s_1stQu$v[1], s_Median$v[1],
                                                s_3rdQu$v[1], s_Max$v[1])))

  if(par == T){
    par(mfrow = c(2, 3))
  } else {
    par(mfrow = c(1, 1))
    par(ask = FALSE)
  }

  if(x$formula == "~1"){
    if(par != T){
      par(ask = TRUE)
    }
    par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 1))
    plot(s_int,
         col = g00,
         main = c("First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_int), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Min,
         col = g00,
         main = paste("Min. First-order Intensity in space \n lambda = ", round(s_Min$v[1], 2)))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Min), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_1stQu,
         col = g00,
         main = paste("1st Qu. First-order Intensity in space \n lambda = ",  round(s_1stQu$v[1], 2)))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_1stQu), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Median,
         col = g00,
         main = paste("Median First-order Intensity in space \n lambda = ",  round(s_Median$v[1], 2)))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Median), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_3rdQu,
         col = g00,
         main = paste("3rd Qu. First-order Intensity in space \n lambda = ",  round(s_3rdQu$v[1], 2)))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_3rdQu), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Max,
         col = g00,
         main = paste("Max. First-order Intensity in space \n lambda = ",  round(s_Max$v[1], 2)))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Max), add = T)}
    par(mar = c(5, 4, 4, 2) + 0.1)
    if(par != T){
      par(ask = FALSE)
    }
  } else {
    if(par != T){
      par(ask = TRUE)
    }
    par(mar = c(5, 4, 4, 2) + 0.1 - c(2, 1 , 1, 1))
    plot(s_int, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_int)),
         main = c("First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_int), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Min, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_Min)),
         main = c("Min. First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Min), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_1stQu, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_1stQu)),
         main = c("1st Qu. First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_1stQu), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Median, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_Median)),
         main = c("Median First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Median), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_3rdQu, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_3rdQu)),
         main = c("3rd Qu. First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_3rdQu), add = T)}
    if(par != T){
      par(ask = TRUE)
    }
    plot(s_Max, zap = zap,
         col = spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                        range = range(s_Max)),
         main = c("Max. First-order Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_Max), add = T)}
    par(mar = c(5, 4, 4, 2) + 0.1)
    if(par != T){
      par(ask = FALSE)
    }
  }


  if(print.bw == T){print(sig)}
}

