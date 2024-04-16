#' Display outlying LISTA functions
#'
#' @description 
#' This function works on the objects of class \code{localdiag}, as returned by
#' \link{localdiag}, plotting the identified 'outlying'
#' LISTA functions. These correspond to the influential points in the fitting
#' of the model provided by \link{localdiag}
#'
#' @param x An object of class \code{localdiag}
#' @param id The id of the LISTA to display.
#' Default is set to the ids identified and stored in the \code{localdiag} object
#'
#' @export
#'
#' @author Nicoletta D'Angelo and Giada Adelfio
#'
#' @seealso
#' \link{localdiag}, \link{plot.localdiag},
#'  \link{print.localdiag}, \link{summary.localdiag}
#'
#'
#' @examples
#'
#' set.seed(2)
#' inh <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(.3, 6))
#' 
#' mod1 <- stppm(inh, formula = ~ 1)
#' 
#' resmod1 <- localdiag(inh, mod1$l, p = .9)
#' 
#' infl(resmod1)
#'
#'
#'
#' @references
#' Adelfio, G., Siino, M., Mateu, J., and Rodríguez-Cortés, F. J. (2020). Some properties of local weighted second-order statistics for spatio-temporal point processes. Stochastic Environmental Research and Risk Assessment, 34(1), 149-168.
#'
#' D’Angelo, N., Adelfio, G.  and Mateu, J. (2022) Local inhomogeneous second-order characteristics for spatio-temporal point processes on linear networks. Stat Papers. https://doi.org/10.1007/s00362-022-01338-4
#'
#'
infl <- function(x, id = NULL){

  if(!inherits(x,"localdiag")) stop("class(x) must be localdiag")
  
  if (!is.null(id)){
    if(!is.numeric(id)) {
      stop("id should be a numeric vector")
    } 
  }
  
  # normalize <- function(x){(x - min(x)) / ((max(x) - min(x)) + 0.0001) }

  if(is.null(id)){
    id <- x$ids
  }

  n <- length(id)

  if(n == 1) {
    oldpar <- par(mfrow = c(1, 1))
    on.exit(par(oldpar))
  } else {
    a <- as.numeric(max(spatstat.utils::primefactors(n)))
    b <- n / a
    oldpar <- par(mfrow = c(b, a))
    on.exit(par(oldpar))
  }

  if(inherits(x$X, "stlp")){
    
    for(i in id){
      fields::image.plot(
        # normalize(x$listas[[i]]$Kinhom)
        x$listas[[i]]$Kinhom
        , xlab = "r",
                      ylab = "h", #zlim = range(normalize(x$listas[[i]]$Kinhom)),
            main = paste("id", i),
            col = grDevices::hcl.colors(12, "YlOrRd",  rev = TRUE),
            axes = FALSE)
      axis(1, at = seq(0, 1, l = length(x$dist)), labels = round(x$dist, 3))
      axis(2, at = seq(0, 1, l = length(x$times)), labels = round(x$times, 3))
      box()
    }
  } else {
    for(i in id){
      fields::image.plot(
        # normalize(x$listas[, , i]),
        x$listas[, , i],
        xlab = "r",
                      ylab = "h", #zlim = range(normalize(x$listas[, , i])),
                      main = paste("id", i),
                      col = grDevices::hcl.colors(12, "YlOrRd", rev = TRUE),
                      axes = FALSE)
      
      axis(1, at = seq(0, 1, l = length(x$dist)), labels = round(x$dist, 3))
      axis(2, at = seq(0, 1, l = length(x$times)), labels = round(x$times, 3))
      box()
    }
  }

}
