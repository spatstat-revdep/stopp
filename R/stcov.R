#' Create \code{stcov} objects and interpolate spatio-temporal covariates on a regular grid
#'
#' @description This function interpolates the covariate values observed 
#' at some observed sites to a regular grid. The imput object
#' should be either a matrix or a dataframe with four columns:
#' \code{x}, \code{y}, \code{t}, and the covariate values, named as the
#'  covariate later called in the model formula (see \code{stppm}).
#'  The interpolation is performed through Inverse Distance Weighting (IDW).
#'  See the Details.
#'  
#' @details 
#' 
#' The function builds a regular grid with equispaced values along the three
#' coordinates and interpolates the covariate values at the new locations.
#' The interpolation at a point location \eqn{x_k} is performed
#' through the inverse-distance weighting smoothing procedure of the covariate 
#' values \eqn{Z(x_j)} at their sampling locations \eqn{j=1, \ldots, J}. 
#' In such a case, the smoothed value at location \eqn{x_k} is 
#' \deqn{Z(x_k) = \frac{\sum_j w_j Z(x_j)}{\sum_j w_j},}
#' where the weight \eqn{w_j} is the \eqn{j}-th element of the inverse \eqn{p}th powers
#' of distance,  
#' \deqn{\textbf{w}=\{w_j\}_{j=1}^J=\{\frac{1}{d(x_k-x_j)^p}\}_{j=1}^J,}
#' with \deqn{d(x_k-x_j) = ||x_k-x_j||} the Euclidean distance from \eqn{x_k}
#'  to \eqn{x_j}.
#'
#' @seealso \link{stppm}
#'
#' @param x A data.frame with four columns, containing the spatio-temporal coordinates
#' and the covariate values.
#' @param interp Logical value indicating whether to interpolate the covariate
#' on a regular grid. Default to \code{TRUE}.
#' @param nx Number of coordinates to generate for each dimention.
#' The default is \code{floor((mult * nrow(cov)) ^ (1/3))}.
#' @param mult The multiplicand of the number of points in the default for \code{nx}.
#' @param p Power of IDW distances.
#' @param names Factor string to name the covariate.
#'
#' @return A \code{stpm} object, to be imputed as list object in \code{stppm}.
#' 
#' @export
#'
#' @author Nicoletta D'Angelo
#' 
#' @examples
#'
#' \dontrun{
#'
#' set.seed(2)
#' df <- data.frame(runif(100), runif(100), runif(100), rpois(100, 15))
#'
#' cov <- stcov(df, interp = FALSE)
#'
#' }
#' 
stcov <- function(x, interp = TRUE, nx = NULL, mult = 1, p = 81, names = NULL){
  # cov must be a dataframe - x, y, t, cov
  if(!inherits(x, "data.frame")) stop("class(x) must be data.frame")
  if(!is.numeric(x[, 4])) stop("The covariate must be continuous")
  time1 <- Sys.time()
  if(interp){
    s.region <- splancs::sbox(cbind(x[, 1], x[, 2]), xfrac = 0.01, yfrac = 0.01)
    xr = range(x[, 3], na.rm = TRUE)
    xw = diff(xr)
    t.region <- c(xr[1] - 0.01 * xw, xr[2] + 0.01 * xw)
    if(is.null(nx)) ff <- floor((mult * nrow(x)) ^ (1/3))
    x0 <- y0 <- t0 <- 1:ff
    scale_to_range <- function(x, new_min, new_max) {
      ((x - min(x)) / (max(x) - min(x))) * (new_max - new_min) + new_min
    }
    x0 <- scale_to_range(x0, s.region[1, 1], s.region[2, 1])
    y0 <- scale_to_range(y0, s.region[1, 2], s.region[3, 2])
    t0 <- scale_to_range(t0, t.region[1], t.region[2])
    df0 <- expand.grid(x0, y0, t0)
    colnames(df0) <- c("xx", "xy", "xt")
    colnames(x) <- c("x", "y", "t", names(x[4]))
    
    gu <- interp3D(df0, x, p = p, d = 3, parallel = FALSE)
    
    df <- cbind(df0, gu)
    if(is.null(names)){
      colnames(df) <- c("x", "y", "t", "cov") 
    } else {
      colnames(df) <- c("x", "y", "t", names) 
    }
    out <- list(df = df)
  } else {
    df <- data.frame(df)
    if(is.null(names)){
      colnames(df) <- c("x", "y", "t", "cov") 
    } else {
      colnames(df) <- c("x", "y", "t", names) 
    }
      out <- list(df = df)
  }
  
  time2 <- Sys.time()
  if(interp) print(paste0(round(as.numeric(difftime(time1 = time2, 
                                         time2 = time1, units = "mins")), 3),
               " minutes"))
  class(out) <- "stcov"
  return(out)
}
