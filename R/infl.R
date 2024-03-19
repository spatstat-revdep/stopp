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
#' \donttest{
#' inh <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, 
#'              par = c(.3, 6), seed = 2)
#' 
#' mod1 <- stppm(inh, formula = ~ 1)
#' 
#' resmod1 <- localdiag(inh, mod1$l, p = .9)
#' 
#' infl(resmod1)
#'
#' }
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
  
  normalize <- function(x){(x - min(x)) / ((max(x) - min(x)) + 0.0001) }

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
      graphics::image(normalize(x$listas[[i]]$Kinhom),
            main = paste("id", i),col = grDevices::hcl.colors(12, "YlOrRd",
                                                              rev = TRUE))
    }
  } else {
    for(i in id){
      graphics::image(normalize(x$listas[, , i]),
            main = paste("id", i), col = grDevices::hcl.colors(12, "YlOrRd",
                                                               rev = TRUE))
    }
  }

}
