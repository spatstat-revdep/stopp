#' Fit a local Poisson process model to a spatio-temporal point pattern
#'
#' @description
#' This function fits a Poisson process model to an observed spatio-temporal
#' point pattern stored in a \code{stp} object, that is, a Poisson model with
#' a set of parameters \eqn{\theta_i} for each point \eqn{i}.
#'
#' @details
#' 
#' We assume that the template model is a Poisson process, with a parametric
#' intensity or rate function \eqn{\lambda(\textbf{u}, t; \theta_i)}  with space
#' and time locations \eqn{\textbf{u} \in
#' W,  t \in T} and parameters \eqn{\theta_i \in \Theta.}
#'
#' Estimation is performed through the fitting of a \code{glm} using a localized
#' version of the quadrature scheme by Berman and Turner (1992), firstly introduced
#' in the purely spatial context by Baddeley (2017), and in the spatio-temporal
#' framework by D'Angelo et al. (2023).
#' 
#' @param X A \code{stp} object
#' @param formula An object of class \code{"formula"}: a symbolic description of
#'  the model to be fitted.
#' The current version only supports formulas depending on the spatial and temporal coordinates:
#' \code{x}, \code{y}, \code{t}.
#' @param verbose Default to \code{TRUE}
#' @param mult The multiplicand of the number of data points,
#'  for setting the number of dummy
#' points to generate for the quadrature scheme
#' @param seed The seed used for the simulation of the dummy points. Default to
#'  \code{NULL}.
#' @param hs Character string indicating whether to select fixed or variable bandwidths
#' for the kernel weights to be used in the log-likelihood.
#' In any of those cases, the well-supported rule-of-thumb for choosing the
#' bandwidth of a Gaussian kernel density estimator is employed.
#' If \code{hs = "global"} (default), a fixed bandwidth is selected.
#' If \code{hs = "local"}, an individual bandwidth is selected for each point in the
#' pattern  \code{X}.
#' @param npx0 Number of lags for the space grid period for variable bandwidths kernel
#' @param npt0 Number of lags for the time period for variable bandwidths kernel
#'
#' @return An object of class \code{locstppm}. A list of
#' \describe{
#' \item{\code{IntCoefs}}{The fitted global coefficients}
#' \item{\code{IntCoefs_local}}{The fitted local coefficients}
#' \item{\code{X}}{The \code{stp} object provided as input}
#' \item{\code{nX}}{The number of points in \code{X}}
#' \item{\code{I}}{Vector indicating which points are dummy or data}
#' \item{\code{y_resp}}{The response variable of the model fitted to the quadrature scheme}
#' \item{\code{formula}}{The formula provided as input}
#' \item{\code{l}}{Fitted intensity through the global parameters}
#' \item{\code{l_local}}{Fitted intensity through the local parameters}
#' \item{\code{mod_global}}{The \code{glm} object of the model fitted to the quadrature scheme}
#' \item{\code{newdata}}{The data used to fit the model, without the dummy points}
#' \item{\code{time}}{Time elapsed to fit the model, in minutes}
#' }
#'
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stppm}
#'
#'
#' @examples
#'
#' set.seed(2)
#' inh <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, 
#'              par = c(0.005, 5))
#' inh_local <- locstppm(inh, formula = ~ x)
#' 
#'
#'
#' @references
#' Baddeley, A. (2017). Local composite likelihood for spatial point processes. Spatial Statistics, 22, 261-295.
#'
#' D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
#'
#'
locstppm <- function(X, formula, verbose = TRUE, mult = 4, seed = NULL,
                     hs = c("global", "local"), npx0 = 10, npt0 = 10){
  
  if (!inherits(X, c("stp"))) stop("X should be from class stp")
  
  time1 <- Sys.time()
  
  if (!is.numeric(mult)) {
    stop("mult should be a numeric value")
  } else {
    if(mult <= 0) {
      stop("mult should be mult > 0")
    }
  } 
  
  if (!is.numeric(npx0)) {
    stop("npx0 should be a numeric value")
  } else {
    if(npx0 <= 0) {
      stop("npx0 should be npx0 > 0")
    }
  } 
  
  if (!is.numeric(npt0)) {
    stop("npt0 should be a numeric value")
  } else {
    if(npt0 <= 0) {
      stop("npt0 should be npt0 > 0")
    }
  } 
  
  hs <- match.arg(hs)
  nX <- nrow(X$df)

  x <- X$df$x
  y <- X$df$y
  t <- X$df$t

  s.region <- splancs::sbox(cbind(x, y), xfrac = 0.01, yfrac = 0.01)
  xr = range(t, na.rm = TRUE)
  xw = diff(xr)
  t.region <- c(xr[1] - 0.01 * xw, xr[2] + 0.01 * xw)

  HomLambda <- nX

  rho <- mult * HomLambda

  set.seed(seed)
  dummy_points <- rstpp(lambda = rho, nsim = 1, verbose = F,
                        minX = s.region[1, 1], maxX = s.region[2, 1],
                        minY = s.region[1, 2], maxY = s.region[3, 2],
                        minT = t.region[1], maxT = t.region[2])$df

  quad_p <- rbind(X$df, dummy_points)

  z <- c(rep(1, nX), rep(0, dim(dummy_points)[1]))

  xx   <- quad_p[, 1]
  xy   <- quad_p[, 2]
  xt   <- quad_p[, 3]
  win <- spatstat.geom::box3(xrange = range(xx), yrange = range(xy), zrange = range(xt))

  ncube <- .default.ncube(quad_p)
  length(ncube) == 1
  ncube <- rep.int(ncube, 3)
  nx <- ncube[1]
  ny <- ncube[2]
  nt <- ncube[3]

  nxyt <- nx * ny * nt
  cubevolume <-  spatstat.geom::volume(win) / nxyt
  volumes <- rep.int(cubevolume, nxyt)

  id <- .grid.index(xx, xy, xt,
                    win$xrange, win$yrange, win$zrange, nx, ny, nt)$index

  w <- .counting.weights(id, volumes)

  y_resp <- z / w

  dati.modello <- cbind(y_resp, w, quad_p[, 1], quad_p[, 2], quad_p[, 3])
  colnames(dati.modello) <- c("y_resp", "w", "x", "y", "t")

  dati.modello <- as.data.frame(dati.modello)

  suppressWarnings(mod_global <- try(glm(as.formula(paste("y_resp", paste(formula, collapse = " "), sep = " ")),
                        weights = w, family = poisson, data = dati.modello), silent = T))
  res_global <- mod_global$coefficients
  pred_global <- exp(predict(mod_global, newdata = dati.modello[1:nX, ]))

  nU <- dim(quad_p)[1]

  h_x <- MASS::bandwidth.nrd(x)
  h_y <- MASS::bandwidth.nrd(y)
  h_t <- MASS::bandwidth.nrd(t)

    if(hs == "local"){
      h_x <- h_y <- kde2d.new.var(x, y, gx = x, gy = y, np = npx0)$hx
      h_t <- density.new.var(t, gx = t, np = npt0)$h
    }

    localwt <- matrix(NA, nrow = nX, ncol = nU)

    if(verbose) cat(paste("\n", "Computing Kernel Densities to the", nX,
                          "points", "\n", "\n"))
    for(j in 1:nX) {
      if(verbose) spatstat.geom::progressreport(j, nX)
      localwt[j, ] <- switch(hs,
                             "local" = dnorm(xx - x[j],
                                             sd = h_x[j]) * dnorm(xy - y[j],
                                                                  sd = h_y[j]) * dnorm(xt - t[j],
                                                                                       sd = h_t[j]),
                             "global" = dnorm(xx - x[j],
                                              sd = h_x) * dnorm(xy - y[j],
                                                                sd = h_y) * dnorm(xt - t[j],
                                                                                  sd = h_t))
    }

    a_s <- localwt * w

    res_local <- matrix(NA, nrow = nX, ncol = length(mod_global$coefficients))
    pred_local <- vector(length = nX)
    if(verbose) cat(paste("\n", "Fitting local model to the", nX,
                          "points", "\n", "\n"))
    for(i in 1:nX){
      spatstat.geom::progressreport(i, nX)
      suppressWarnings(mod_local <- try(glm(as.formula(paste("y_resp", paste(formula, collapse = " "), sep = " ")),
                                            weights = a_s[i, ], family = poisson, data = dati.modello), silent = T))
      res_local[i, ] <- mod_local$coefficients
      pred_local[i] <- exp(predict(mod_local, newdata = dati.modello[i, ]))
    }

    res_local <- as.data.frame(res_local)
    colnames(res_local) <- names(res_global)
    
  time2 <- Sys.time()

    list.obj <- list(IntCoefs = res_global,
                     IntCoefs_local = res_local,
                     X = X,
                     nX = nX,
                     I = z,
                     y_resp = y_resp,
                     formula = formula,
                     l = as.vector(pred_global),
                     l_local = as.vector(pred_local),
                     mod_global = mod_global,
                     newdata = dati.modello[1:nX, ],
                     time = paste0(round(as.numeric(difftime(time1 = time2, time2 = time1, units = "mins")), 3), " minutes"))

  class(list.obj) <- "locstppm"
  return(list.obj)

}















