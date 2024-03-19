#' Fit a log-Gaussian Cox process model to a spatio-temporal point pattern
#'
#' @description
#'
#' This function estimates a log-Gaussian Cox process (LGCP), following the
#' **joint minimum contrast** procedure introduced in Siino et al. (2018)
#' .
#'
#' Three covariances are available: separable exponential,
#' Gneiting, and De Iaco-Cesare.
#'
#' If the \code{first} and \code{second} arguments are set to \code{local}, a local
#'  log-Gaussian
#' Cox process is fitted by means of the ** locally weighted minimum contrast**
#'  procedure proposed in
#' D'Angelo et al. (2023).
#'
#' @details
#'
#' Following the inhomogeneous specification in Diggle et al. (2013),
#' we consider LGCPs with intensity
#' \deqn{
#' \Lambda(\textbf{u},t)=\lambda(\textbf{u},t)\exp(S(\textbf{u},t)).
#' }
#'
#'
#' @param X A \code{stp} object
#' @param formula An object of class \code{formula}: a symbolic description of the first-order intensity to be fitted.
#' The current version only supports formulas depending on the spatial and temporal coordinates:
#' \code{x}, \code{y}, \code{t}. Default to \code{formula = ~ 1} which provides an homogeneous
#' first-order intensity.
#' @param verbose Default to TRUE
#' @param cov Covariance function to be fitted for the second-order intensity function.
#' Default to \code{separable}. Other options are \code{gneiting} and \code{iaco-cesare}".
#' @param first Character string indicating whether to fit a first-order intensity function
#'  with global or local parameters:
#'  either \code{global} (default) or \code{local}.
#' @param second Character string indicating whether to fit a second-order intensity function
#'  with global or local parameters:
#'  either \code{global} (default) or \code{local}.
#' @param mult The multiplicand of the number of data points,
#'  for setting the number of dummy
#' points to generate for the quadrature scheme
#' @param hs Character string indicating whether to select fixed or variable bandwidths
#' for the kernel weights to be used in the log-likelihood.
#' In any of those cases, the well-supported rule-of-thumb for choosing the
#' bandwidth of a Gaussian kernel density estimator is employed.
#' If \code{hs = "global"} (default), a fixed bandwidth is selected.
#' If \code{hs = "local"}, an individual bandwidth is selected for each point in the
#' pattern  \code{X}.
#' @param npx0 A positive integer representing the spatial distance to np-th closest event.
#' Used in the computation of the local bandwidth.
#' Suitable values are in the range from 10 (default) to 100.
#' @param npt0 A positive integer representing the temporal distance to np-th closest event.
#' Used in the computation of the local bandwidth.
#' Suitable values are in the range from 10 (default) to 100.
#' @param itnmax Maximum number of iterations to run in the optimization procedure
#' for the estimation of the second-order intensity parameters.
#' @param min_vals Minimum values of the optimization procedure for the minimum contrast.
#' @param max_vals Maximum values of the optimization procedure for the minimum contrast.
#'
#' @return A list of the class \code{stlgcppm}, containing
#' \describe{
#' \item{\code{IntCoefs}}{The fitted coefficients of the first-order intensity function}
#' \item{\code{CovCoefs}}{The fitted coefficients of the second-order intensity function}
#' \item{\code{X}}{The stp object provided as input}
#' \item{\code{formula}}{The formula provided as input}
#' \item{\code{cov}}{A string with the chosen covariance type}
#' \item{\code{l}}{Fitted first-order intensity}
#' \item{\code{mu}}{Mean function of the random intensity}
#' \item{\code{mod_global}}{The glm object of the model fitted to the quadrature scheme
#' for the first-order intensity parameters estimation}
#' \item{\code{newdata}}{The data used to fit the model, without the dummy points}
#' \item{\code{time}}{Time elapsed to fit the model, in minutes}
#' }
#'
#'
#'
#' @export
#'
#' @author Nicoletta D'Angelo, Giada Adelfio, and Marianna Siino
#'
#' @seealso
#' \link{print.stlgcppm}, \link{summary.stlgcppm}, \link{localsummary.stlgcppm},
#'  \link{plot.stlgcppm}, \link{localplot.stlgcppm}
#'
#'
#' @examples
#'
#' \dontrun{
#' 
#' catsub <- stp(greececatalog$df[1:100, ])
#' 
#' lgcp1 <- stlgcppm(catsub)
#' 
#'}
#'
#'
#' @references
#'
#' Baddeley, A. (2017). Local composite likelihood for spatial point processes. Spatial Statistics, 22, 261-295.
#'
#' D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
#'
#' Diggle, P. J., Moraga, P., Rowlingson, B., and Taylor, B. M. (2013). Spatial and spatio-temporal log-gaussian cox processes: extending the geostatistical paradigm. Statistical Science, 28(4):542–563.
#'
#' Gabriel, E., Rowlingson, B. S., and Diggle, P. J. (2013). stpp: An R Package for Plotting, Simulating and Analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53(2), 1–29. https://doi.org/10.18637/jss.v053.i02
#'
#' Siino, M., Adelfio, G., and Mateu, J. (2018). Joint second-order parameter estimation for spatio-temporal log-Gaussian Cox processes. Stochastic environmental research and risk assessment, 32(12), 3525-3539.
#'
stlgcppm <- function(X, formula = ~ 1, verbose = TRUE,
                     cov = c("separable", "gneiting", "iaco-cesare"),
                     first = c("global", "local"), second = c("global", "local"), mult = 4,
                     hs = c("global", "local"), npx0 = 10, npt0 = 10, itnmax = 100,
                     min_vals = NULL, max_vals = NULL){
  
  if (!inherits(X, c("stp"))) stop("X should be from class stp")
  
  time1 <- Sys.time()
  
  cov <- match.arg(cov)
  first <- match.arg(first)
  second <- match.arg(second)
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
  
  dummy_points <- rstpp(lambda = rho, nsim = 1, seed = 2, verbose = F,
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
  xx <- quad_p[, 1]
  xy <- quad_p[, 2]
  xt <- quad_p[, 3]
  
  h_x <- MASS::bandwidth.nrd(x)
  h_y <- MASS::bandwidth.nrd(y)
  h_t <- MASS::bandwidth.nrd(t)
  
  if(first == "local"){
    
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
      progressreport(i, nX)
      suppressWarnings(mod_local <- try(glm(as.formula(paste("y_resp", paste(formula, collapse = " "), sep = " ")),
                                            weights = a_s[i, ], family = poisson, data = dati.modello), silent = T))
      res_local[i, ] <- mod_local$coefficients
      pred_local[i] <- exp(predict(mod_local, newdata = dati.modello[i, ]))
    }
    
  }
  
  # covariance
  
  tlim <- max(t)
  cube <- (range(x) - min(x))[2]
  
  lgrid_s <- 25
  rsup_t <- tlim / 4
  
  dt <- dist(t)
  bw_vector_t <- KernSmooth::dpik(dt, kernel = "box", range.x = c(min(dt), max(dt)))
  
  rmin_vector_t2 <- (bw_vector_t * ((1 + lgrid_s) / lgrid_s))
  
  vals_t <- cbind(rmin_vector_t2, bw_vector_t)
  
  rsup_s <- cube / 4
  
  ds <- dist(cbind(x, y))
  bw2_vector <- KernSmooth::dpik(ds, kernel = "epanech", range.x = c(min(ds), max(ds)))
  
  rmin3_vector <- bw2_vector * ((1 + lgrid_s) / lgrid_s)
  
  vals_s <- cbind(rmin3_vector, bw2_vector)
  
  s.region <- matrix(c(range(x)[1] - 0.05,
                       range(x)[2] + 0.05,
                       range(x)[2] + 0.05,
                       range(x)[1] - 0.05,
                       range(y)[1] - 0.05,
                       range(y)[1] - 0.05,
                       range(y)[2] + 0.05,
                       range(y)[2] + 0.05), ncol = 2)
  
  t.region <- c(min(t) - 1, tlim + 1)
  
  
  u <- seq(vals_s[1], rsup_s, len = lgrid_s)
  v <- seq(vals_t[1], rsup_t, len = lgrid_s)
  us0 <- vals_s[2]
  vt0 <- vals_t[2]
  
  g0 <- if(first == "global"){
    stpp::LISTAhat(xyt = as.stpp(X), lambda = pred_global, s.region = s.region, t.region = t.region,
                   dist = u, times = v, ks = "epanech", hs = us0, kt = "box",
                   ht = vt0)
  } else {
    stpp::LISTAhat(xyt = as.stpp(X), lambda = pred_local, s.region = s.region, t.region = t.region,
                   dist = u, times = v, ks = "epanech", hs = us0, kt = "box",
                   ht = vt0)
  }
  
  max_dist <- max(dist(s.region))
  max_dist_t <- max(diff(t.region))
  
  nonpar_g_st_finite <- is.finite(g0$list.LISTA)
  g0$list.LISTA[!nonpar_g_st_finite] <- 0
  
  if(cov == "separable"){
    
    if(is.null(min_vals)) min_vals <- c(0.01, 0.001, 0.01)#c(0.01, 0.001, 20)
    if(is.null(max_vals)) max_vals <-  c(50, 10, 1500)
    
    start_par <- c(log(nX) / 2, max_dist / 100, max_dist_t / 100)
    
    if(second == "global"){
      
      MINCON.EXP_EXP <- c(NA, NA, NA)
      
      MINCON.EXP_EXP <- optimx::optimx(par = start_par,
                                       fn = g.sep_st_exp_exp2,
                                       useq = g0$dist,
                                       vseq = g0$times,
                                       ghat = g0$list.LISTA,
                                       transform = NULL,
                                       power = 1, method = c("nlm"),
                                       lower = min_vals,
                                       upper = max_vals,
                                       itnmax = itnmax)
      
      res <- as.numeric(as.vector(MINCON.EXP_EXP[1 : 3]))
      
    } else {
      
      res <- matrix(NA, nrow = nX, ncol = 3)
      
      if(verbose) cat(paste("\n", "Fitting local Separable covariance to the", nX,
                            "points", "\n", "\n"))
      
      for(id in 1:nX){
        
        spatstat.geom::progressreport(id, nX)
        
        MINCON.EXP_EXP <- c(NA, NA, NA)
        
        wi <- dnorm(x - x[id], sd = h_x) * dnorm(y - y[id], sd = h_y) * dnorm(t - t[id], sd = h_t)
        
        numer0 <- sapply(1:nX, function(x) g0$list.LISTA[, , x] * wi[x], simplify = "array")
        
        numer <- apply(numer0, 1:2, sum)
        
        denom0 <- sapply(1:nX, function(x) nonpar_g_st_finite[, , x] * wi[x], simplify = "array")
        
        denom <- apply(denom0, 1:2, sum)
        
        avg_lista <- numer / denom
        
        MINCON.EXP_EXP <- optimx::optimx(par = start_par,
                                         fn = g.sep_st_exp_exp2,
                                         useq = g0$dist,
                                         vseq = g0$times,
                                         ghat = avg_lista,
                                         transform = NULL,
                                         power = 1, method = c("nlm"),
                                         lower = min_vals,
                                         upper = max_vals,
                                         itnmax = itnmax)
        
        
        res[id, ] <- as.numeric(as.vector(MINCON.EXP_EXP[1 : 3]))
      }
      
    }
    
  }
  
  if(cov == "gneiting") {
    
    if(is.null(min_vals)) min_vals <- c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
    if(is.null(max_vals)) max_vals <- c(nX, max_dist * 250, max_dist_t * 250, 2, 2, 2)
    
    start_par <- c(log(nX) / 2, max_dist, max_dist_t, 1, 1, 1)
    
    if(second == "global"){
      
      MINCON_GN <- c(NA, NA, NA, NA, NA, NA)
      
      MINCON_GN <- optimx::optimx(par = start_par,
                                  fn = g_st,
                                  useq = g0$dist,
                                  vseq = g0$times,
                                  ghat = g0$list.LISTA,
                                  transform = NULL,
                                  power = 1,
                                  method = c("nlm"),
                                  lower = min_vals,
                                  upper = max_vals,
                                  itnmax = itnmax)
      
      res <- as.numeric(as.vector(MINCON_GN[1 : 6]))
      
    } else {
      
      res <- matrix(NA, nrow = nX, ncol = 6)
      
      if(verbose) cat(paste("\n", "Fitting local Gneiting covariance to the", nX,
                            "points", "\n", "\n"))
      
      for(id in 1:nX){
        
        spatstat.geom::progressreport(id, nX)
        
        MINCON_GN <- c(NA, NA, NA, NA, NA, NA)
        
        wi <- dnorm(x - x[id], sd = h_x) * dnorm(y - y[id], sd = h_y) * dnorm(t - t[id], sd = h_t)
        
        numer0 <- sapply(1:nX, function(x) g0$list.LISTA[, , x] * wi[x], simplify = "array")
        
        numer <- apply(numer0, 1:2, sum)
        
        denom0 <- sapply(1:nX, function(x) nonpar_g_st_finite[, , x] * wi[x], simplify = "array")
        
        denom <- apply(denom0, 1:2, sum)
        
        avg_lista <- numer / denom
        
        MINCON_GN <- optimx::optimx(par = start_par,
                                    fn = g_st,
                                    useq = g0$dist,
                                    vseq = g0$times,
                                    ghat = avg_lista,
                                    transform = NULL,
                                    power = 1,
                                    method = c("nlm"),
                                    lower = min_vals,
                                    upper = max_vals,
                                    itnmax = itnmax)
        
        res[id, ] <- as.numeric(as.vector(MINCON_GN[1 : 6]))
      }
    }
  }
  
  if(cov == "iaco-cesare") {
    
    if(is.null(min_vals)) min_vals <- c(0.01, 0.01, 0.01, 0.01, 0.01, 1.5)
    if(is.null(max_vals)) max_vals <- c(nX, max_dist * 250, max_dist_t * 250, 2, 2, 15)
    
    start_par <- c(log(nX) / 2, max_dist, max_dist_t, 1, 1, 8)
    
    if(second == "global"){
      
      MINCON_IACO <- c(NA, NA, NA, NA, NA, NA)
      
      MINCON_IACO <- optimx::optimx(par = start_par,
                                    fn = g_st_iaco,
                                    useq = g0$dist,
                                    vseq = g0$times,
                                    ghat = g0$list.LISTA,
                                    transform = NULL,
                                    power = 1,
                                    method = c("nlm"),
                                    lower = min_vals,
                                    upper = max_vals,
                                    itnmax = itnmax)
      
      
      res <- as.numeric(as.vector(MINCON_IACO[1 : 6]))
      
    } else {
      
      res <- matrix(NA, nrow = nX, ncol = 6)
      
      if(verbose) cat(paste("\n", "Fitting local Iaco-Cesare covariance to the", nX,
                            "points", "\n", "\n"))
      
      for(id in 1:nX){
        
        spatstat.geom::progressreport(id, nX)
        
        MINCON_IACO <- c(NA, NA, NA, NA, NA, NA)
        
        wi <- dnorm(x - x[id], sd = h_x) * dnorm(y - y[id], sd = h_y) * dnorm(t - t[id], sd = h_t)
        
        numer0 <- sapply(1:nX, function(x) g0$list.LISTA[, , x] * wi[x], simplify = "array")
        
        numer <- apply(numer0, 1:2, sum)
        
        denom0 <- sapply(1:nX, function(x) nonpar_g_st_finite[, , x] * wi[x], simplify = "array")
        
        denom <- apply(denom0, 1:2, sum)
        
        avg_lista <- numer / denom
        
        MINCON_IACO <- optimx::optimx(par = start_par,
                                      fn = g_st_iaco,
                                      useq = g0$dist,
                                      vseq = g0$times,
                                      ghat = avg_lista,
                                      transform = NULL,
                                      power = 1,
                                      method = c("nlm"),
                                      lower = min_vals,
                                      upper = max_vals,
                                      itnmax = itnmax)
        
        
        res[id, ] <- as.numeric(as.vector(MINCON_IACO[1 : 6]))
      }
    }
  }
  
  time2 <- Sys.time()
  
  if(second == "local"){
    
    colnames(res) <- switch(cov,
                            "separable" = c("sigma", "alpha", "beta"),
                            "gneiting" = c("sigma", "alpha", "beta", "gamma_s", "gamma_t", "delta"),
                            "iaco-cesare" = c("sigma", "alpha", "beta", "gamma_s", "gamma_t", "delta"))
    
    
    res <- as.data.frame(res)
    
  } else {
    
    names(res) <- switch(cov,
                         "separable" = c("sigma", "alpha", "beta"),
                         "gneiting" = c("sigma", "alpha", "beta", "gamma_s", "gamma_t", "delta"),
                         "iaco-cesare" = c("sigma", "alpha", "beta", "gamma_s", "gamma_t", "delta"))
    
  }
  
  if(inherits(res, "numeric")){
    int2 <- rep(res[1], nX) / 2
  } else {
    int2 <- res$sigma / 2
  }
  
  if(first == "local"){
    
    names(res_global) <- names(mod_global$coefficients)
    res_local <- data.frame(res_local)
    colnames(res_local) <- names(mod_global$coefficients)
    
    
    if(hs == "global"){
      bw <- c(round(h_x, 3), round(h_y, 3), round(h_t, 3))
      names(bw) <- c("h_x", "h_y", "h_t")
    } else {
      bw <- cbind(round(h_x, 3), round(h_y, 3), round(h_t, 3))
      colnames(bw) <- c("h_x", "h_y", "h_t")
    }
    
    quad_p <- rbind(X$df)
    quad_p <- as.data.frame(quad_p)
    
    list.obj <- list(IntCoefs = res_local,
                     CovCoefs = res,
                     X = X,
                     formula = formula,
                     cov = cov,
                     l = as.vector(pred_local),
                     mu = as.vector(pred_local - int2),
                     mod_global = mod_global,
                     newdata = dati.modello[1:nX, ],
                     time = paste0(round(as.numeric(difftime(time1 = time2, time2 = time1, units = "mins")), 3), " minutes"))
  } else {
    
    
    list.obj <- list(IntCoefs = res_global,
                     CovCoefs = res,
                     X = X,
                     formula = formula,
                     cov = cov,
                     l = as.vector(pred_global),
                     mu = as.vector(pred_global - int2),
                     mod_global = mod_global,
                     newdata = dati.modello[1:nX, ],
                     time = paste0(round(as.numeric(difftime(time1 = time2, time2 = time1, units = "mins")), 3), " minutes"))
    
  }
  
  
  
  
  class(list.obj) <- "stlgcppm"
  return(list.obj)
  
}













