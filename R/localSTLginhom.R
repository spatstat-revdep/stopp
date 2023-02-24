#' Local inhomogeneous Spatio-temporal pair correlation functions on a linear network
#'
#' @description
#' The functions \code{localSTLKinhom} and \code{localSTLginhom} implement the
#'inhomogeneous LISTA functions proposed in D'Angelo et al. (2022).
#'
#' The homogeneous K-function and pair correlation functions, in
#' D'Angelo et al. (2021), can be obtained easily with \code{localSTLKinhom} and
#' \code{localSTLginhom}, by imputing a lambda vector of constant intensity
#' values, the same for each point.
#'
#' These local functions are the local counterparts of \link{STLKinhom} and
#' \link{STLginhom} by
#'  Moradi and Mateu (2020), available in the \code{stlnpp} package (Moradi et al., 2020).
#' Basically, we avoid summing up all
#' the points as in the global statistics counterparts, and we denote the
#' individual contribution to the global statistics with the index \eqn{i}.
#'
#'
#' @details
#'
#' The *local spatio-temporal inhomogeneous*
#'   K-function for the i-th event \eqn{(\boldsymbol{u}_i,t_i)} on a linear network
#' is \deqn{\hat{K}^i_{L,I}(r,h)=\frac{1}{ \vert L  \vert  \vert T  \vert}\sum_{(\boldsymbol{u}_i,t_i)\ne (\boldsymbol{v},s)} \frac{I\{ d_L(\boldsymbol{u}_i,\boldsymbol{v})<r,\vert t_i-s\vert <h\} }{\hat{\lambda}(\boldsymbol{u}_i,t_i)\hat{\lambda}(\boldsymbol{v},s)M((\boldsymbol{u}_i,t_i),d_L(\boldsymbol{u}_i,\boldsymbol{v}),\vert t_i-s\vert )},}
#'  and the corresponding *local pair correlation function* (pcf)
#'   \deqn{\hat{g}^i_{L,I}(r,h)=\frac{1}{ \vert L  \vert  \vert T  \vert}\sum_{(\boldsymbol{u}_i,t_i)\ne (\boldsymbol{v},s)} \frac{\kappa( d_L(\boldsymbol{u}_i,\boldsymbol{v})-r)\kappa(\vert t_i-s\vert -h) }{\hat{\lambda}(\boldsymbol{u}_i,t_i)\hat{\lambda}(\boldsymbol{v},s)M((\boldsymbol{u}_i,t_i),d_L(\boldsymbol{u}_i,\boldsymbol{v}),\vert t_i-s\vert )},}
#' with
#'   \deqn{D(X) = \frac{n-1}{ \vert L  \vert  \vert T  \vert}\sum_{i=1}^n\sum_{i \ne j}}
#'   normalization factor. This leads to the unbieased estimators \eqn{\frac{1}{D(X)}\hat{K}^i_{L,I}(r,h)} and
#' \eqn{\frac{1}{D(X)}\hat{g}^i_{L,I}(r,h)} (D'Angelo et al., 2022) .
#'
#' The homogeneous versions (D'Angelo et al., 2021) can be obtained by weighting the second-order
#' summary statistics (either K or pcf) by a constant intensity
#' \eqn{\hat{\lambda}=n/( \vert L \vert \vert T \vert)}, giving
#'   \deqn{\hat{K}_L^i(r,h)=\frac{1}{\hat{\lambda}^{2} \vert L  \vert  \vert T  \vert}\sum_{(\boldsymbol{u}_i,t_i)\ne (\boldsymbol{v},s)} \frac{I\{ d_L(\boldsymbol{u}_i,\boldsymbol{v})<r,\vert t_i-s\vert <h\} }{M((\boldsymbol{u}_i,t_i),d_L(\boldsymbol{u}_i,\boldsymbol{v}),\vert t_i-s\vert )},}
#' and
#'  \deqn{\hat{g}_L^i(r,h)=\frac{1}{\hat{\lambda}^{2} \vert L  \vert  \vert T  \vert}\sum_{(\boldsymbol{u}_i,t_i)\ne (\boldsymbol{v},s)} \frac{\kappa( d_L(\boldsymbol{u}_i,\boldsymbol{v})-r)\kappa(\vert t_i-s\vert -h) }{M((\boldsymbol{u}_i,t_i),d_L(\boldsymbol{u}_i,\boldsymbol{v}),\vert t_i-s\vert )}.}
#'
#'
#' @param x A realisation of a spatio-temporal point processes on a linear network in \code{stlp} format
#' @param lambda values of estimated intensity.
#' @param normalize normalization factor to be considered.
#' @param r values of argument r where pair correlation function will be evaluated. optional.
#' @param t values of argument t where pair correlation function will be evaluated. optional.
#' @param nxy pixel array dimensions. optional.
#' @return A list of objects of class \code{sumstlpp} (Moradi and Mateu, 2020).
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso \link{localSTLginhom}, \link[stlnpp]{STLKinhom}, \link[stlnpp]{STLginhom}
#'
#' @examples
#'
#'\dontrun{
#'set.seed(10)
#'X <- stlnpp::rpoistlpp(.2, a = 0, b = 5, L = stlnpp::easynet)
#'X
#'lambda <- density(X, at = "points")
#'x <- as.stlp(X)
#'g <- localSTLginhom(x, lambda = lambda, normalize = TRUE)
#'#select an individual point
#'j = 1
#'g[[j]]
#'#plot the lista function
#'inhom <- list(x = g[[j]]$r, y = g[[j]]$t, z = g[[j]]$ginhom)
#'
#'fields::image.plot(inhom, main = "ginhom", col = hcl.colors(12, "YlOrRd", rev = FALSE),
#' xlab = "Spatial distance", ylab = "Temporal distance")
#'
#'}
#'
#'
#'
#' @references
#' D'Angelo, N., Adelfio, G., and Mateu, J. (2021). Assessing local differences between the spatio-temporal second-order structure of two point patterns occurring on the same linear network. Spatial Statistics, 45, 100534.
#'
#' D'Angelo, N., Adelfio, G.  and Mateu, J. (2022). Local inhomogeneous second-order characteristics for spatio-temporal point processes on linear networks. Stat Papers. https://doi.org/10.1007/s00362-022-01338-4
#'
#' Moradi M, Cronie O, and Mateu J (2020). stlnpp: Spatio-temporal analysis of point patterns on linear networks.
#'
#' Moradi, M. M., and Mateu, J. (2020). First-and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics, 29(3), 432-443.
#'

localSTLginhom = function(x, lambda, normalize = FALSE, r = NULL, t = NULL, nxy = 10)
{
  if (!inherits(x, "stlp"))
    stop("x should be from class stlp")

  X <- stopp::as.stlpp(x)

  Y <- stlnpp::as.lpp.stlpp(X)
  l <- spatstat.geom::domain(Y)
  tleng <- summary(l)$totlength
  n <- spatstat.geom::npoints(Y)
  a <- X$time[1]
  b <- X$time[2]
  trange <- b - a
  timev <- X$data$t
  sdist <- spatstat.linnet::pairdist.lpp(Y)
  tdist <- as.matrix(dist(timev))
  toler <- spatstat.linnet::default.linnet.tolerance(l)
  ml <- matrix(1, n, n)
  for (j in 1:n) {
    ml[-j, j] <- spatstat.linnet::countends(l, Y[-j], sdist[-j, j], toler = toler)
  }
  mtplus <- matrix(timev, n, n, byrow = T) + tdist
  mtminus <- matrix(timev, n, n, byrow = T) - tdist
  mtedge <- (mtplus <= b) + (mtminus >= a)
  diag(mtedge) <- 1
  lamden <- outer(lambda, lambda, FUN = "*")
  edgetl <- mtedge * ml * lamden
  maxs <- 0.7 * max(sdist[!is.infinite(sdist)])
  maxt <- 0.7 * (trange/2)
  if (is.null(r))
    r <- seq((maxs/nxy),maxs,by=(maxs-(maxs/nxy))/(nxy-1))
  if (is.null(t))
    t <- seq((maxt/nxy),maxt,by=(maxt-(maxt/nxy))/(nxy-1))

  g_local = list()
  for(k in 1:n){

    g <- matrix(NA, nrow = nxy, ncol = nxy)
    no <- sdist == 0&tdist == 0|sdist == Inf|sdist > maxs|tdist > maxt
    bwl <- stats::bw.nrd0(as.numeric(sdist[!no]))
    bwt <- stats::bw.nrd0(as.numeric(tdist[!no]))

    for (i in 1:length(r)) {
      for (j in 1:length(t)) {

        outl <- spatstat.explore::dkernel(as.numeric(sdist[!no] - r[i]), sd = bwl)
        outt <- spatstat.explore::dkernel(as.numeric(tdist[!no] - t[j]), sd = bwt)
        g1 <- outl * outt / (edgetl[!no])
        no2 <- no
        no2[no] <- NA
        no2[!no] <- g1
        no2 <- no2[, k]
        g[i, j] <- sum(no2[!is.na(no2) & !is.infinite(no2)])

      }
    }
    g_local[[k]] <- g
  }

  if (normalize) {
    revrho <- outer(1 / lambda, 1 / lambda, FUN = "*")
    appx <- (n - 1)*(tleng * trange)/(sum(revrho[lower.tri(revrho,diag = FALSE)]) * 2)
    gval <- lapply(g_local, "*" ,appx)
  }
  else {
    gval <- lapply(g_local, FUN= function(K) (n-1) * K / (trange * tleng) )

  }

  g_theo <- list()
  for(k in 1:n){
    g_theo[[k]] <- matrix(rep(1, length(t) * length(r)), ncol = nxy)
  }

  gout <- list()
  for(k in 1:n){
    gout[[k]] <- list(ginhom = gval[[k]], gtheo = g_theo[[k]], r = r, t = t)
  }
  class(gout[[k]]) <- c("sumstlpp")
  attr(gout[[k]], "nxy") <- nxy
  return(gout)
}
