#' Global diagnostics of a spatio-temporal point process first-order intensity
#'
#' @description
#' This function performs global diagnostics of a model fitted for the
#' first-order intensity of a  spatio-temporal point pattern, by returning
#' the plots of the inhomogeneous K-function weighted by the
#' provided intensity to diagnose, its theoretical value,
#' and their difference.
#'
#' @details
#' If applied to a \code{stp} object, it resorts to  the
#' spatio-temporal inhomogeneous K-function (Gabriel and Diggle, 2009)
#'  documented by the function
#'  \link{STIKhat}  of the \code{stpp} package (Gabriel et al, 2013).
#'
#'  If applied to a \code{stlp} object, it uses the
#' spatio-temporal inhomogeneous K-function on a linear network (Moradi and Mateu, 2020)
#'  documented by the function
#'  \link{STLKinhom}  of the \code{stlnpp} package (Moradi et al., 2020).
#'
#'
#' @param X A \code{stp} object
#' @param intensity A vector of intensity values, of the same length as the number
#' of point in \code{X}
#' @param samescale Logical value. It indicates whether to plot the observed
#' and the theoretical K-function in the same or
#' different scale. Default to \code{TRUE}.
#'
#'
#' @return It plots three panels: the observed K-function, as returned by \link{STLKinhom};
#' the theoretical one; their difference. The function also prints the sum of
#'  squared differences between the observed and theoretical
#' K-function on the console.
#'
#'
#'
#' @export
#'
#' @author Nicoletta D'Angelo and Giada Adelfio
#'
#' @seealso
#' \link{infl}, \link{plot.localdiag}, \link{print.localdiag},
#' \link{summary.localdiag}, \link{localdiag},
#'
#'
#'
#'
#' @examples
#' \donttest{
#'
#' inh <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, 
#'                par = c(.3, 6), seed = 2)
#' 
#' mod1 <- stppm(inh, formula = ~ 1)
#' mod2 <- stppm(inh, formula = ~ x)
#' 
#' globaldiag(inh, mod1$l)
#' globaldiag(inh, mod2$l)
#' 
#' }
#'
#'
#'
#' @references
#' Adelfio, G., Siino, M., Mateu, J., and Rodríguez-Cortés, F. J. (2020). Some properties of local weighted second-order statistics for spatio-temporal point processes. Stochastic Environmental Research and Risk Assessment, 34(1), 149-168.
#'
#' D’Angelo, N., Adelfio, G.  and Mateu, J. (2022) Local inhomogeneous second-order characteristics for spatio-temporal point processes on linear networks. Stat Papers. https://doi.org/10.1007/s00362-022-01338-4
#'
#' Gabriel, E., and Diggle, P. J. (2009). Second‐order analysis of inhomogeneous spatio‐temporal point process data. Statistica Neerlandica, 63(1), 43-51.
#'
#' Gabriel, E., Rowlingson, B. S., & Diggle, P. J. (2013). stpp: An R Package for Plotting, Simulating and Analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53(2), 1–29. https://doi.org/10.18637/jss.v053.i02
#'
#' Moradi M, Cronie O, and Mateu J (2020). stlnpp: Spatio-temporal analysis of point patterns on linear networks.
#'
#' Moradi, M. M., and Mateu, J. (2020). First-and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics, 29(3), 432-443.
#'
#'
globaldiag <- function(X, intensity, samescale = TRUE){
  
  if (!inherits(X, c("stp", "stlp"))) stop("X should be either from class stp or stlp")

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if(any(class(X) == "stp")) {
    Kfunct <- stpp::STIKhat(xyt = as.stpp(X), lambda = intensity)

    est <- Kfunct$Khat
    theo <- Kfunct$Ktheo

    dist <- Kfunct$dist
    times <- Kfunct$times
  } else if(any(class(X) == "stlp")){
    Kfunct <- stlnpp::STLKinhom(X = as.stlpp(X), lambda = intensity)

    est <- Kfunct$Kinhom
    theo <- Kfunct$Ktheo

    dist <- Kfunct$r
    times <- Kfunct$t
  }

  diffK <- est - theo

  lims <- if(samescale){
    range(c(as.numeric(est), as.numeric(theo), as.numeric(diffK)))
  } else {
    NULL
  }

  par(mfrow = c(1, 3))

  fields::image.plot(dist, times, est, main = "Estimated", xlab = "r",
                     ylab = "h",
                     col = grDevices::hcl.colors(12, "YlOrRd", rev = TRUE),
                     zlim = lims,
                     legend.mar =  12, axes = FALSE)
  axis(1, at = seq(0, 1, l = length(dist)), labels = round(dist, 3))
  axis(2, at = seq(0, 1, l = length(times)), labels = round(times, 3))
  box()
  fields::image.plot(dist, times, theo, main = "Theoretical", xlab = "r",
                     ylab = "h",
                     col = grDevices::hcl.colors(12, "YlOrRd", rev = TRUE),
                     zlim = lims,
                     legend.mar =  12, axes = FALSE)
  axis(1, at = seq(0, 1, l = length(dist)), labels = round(dist, 3))
  axis(2, at = seq(0, 1, l = length(times)), labels = round(times, 3))
  box()

  fields::image.plot(dist, times, diffK, main = "Difference", xlab = "r",
                     ylab = "h",
                     col = grDevices::hcl.colors(12, "YlOrRd", rev = TRUE),
                     zlim = lims,
                     legend.mar =  15, axes = FALSE)
  axis(1, at = seq(0, 1, l = length(dist)), labels = round(dist, 3))
  axis(2, at = seq(0, 1, l = length(times)), labels = round(times, 3))
  box()

  paste("Sum of squared differences = ", round(sum(diffK ^ 2 / (est)), 3))

}
