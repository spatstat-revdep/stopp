
.counting.weights <- function(id, volumes) {
  id <- as.integer(id)
  fid <- factor(id, levels = seq_along(volumes))
  counts <- table(fid)
  w <- volumes[id] / counts[id]
  w <- as.vector(w)
  names(w) <- NULL
  return(w)
}




.default.ncube <- function(X){
  guess.ngrid <- floor((splancs::npts(X) / 2) ^ (1 / 3))
  max(5, guess.ngrid)
}



density.new.var  <-function(x,np,gx){
  nx  <-  length(x)
  n =length(gx)
  h<-numeric(nx)
  mat.d=as.matrix(dist((x),upper=TRUE,diag=TRUE))
  diag(mat.d)=999999
  if(missing(np))  np=trunc(nx/2)


  for(i in 1:ncol(mat.d)){
    o=order(mat.d[,i])
    pp=o[np]
    h[i]=mat.d[pp,i]
  }
  h   <-  h/4

  ax  <-  outer(gx,x,"-")/h
  z   <-  rowSums(matrix(dnorm(ax),n,nx)/(nx*h))
  return(list(x=gx,z=z,h=h))
}


g_st <- function(param, useq, vseq, ghat, transform, power){
  uv <- expand.grid(useq, vseq)
  g.parametric <- exp(
    param[1] / ((uv[, 2] / param[3]) ^ param[5] + 1) ^ (param[6] / param[5]) *
      exp( - (uv[, 1] / param[2]) ^ param[4] / ((uv[, 2] / param[3]) ^ param[5] + 1) ^ (param[6] / param[5] * param[4] / 2))
  )
  ghat2 <- as.vector(ghat)
  return(sum((((ghat2)) ^ power - (g.parametric) ^ power) ^ 2))
}

g_st_iaco <-function (param, useq, vseq,ghat, transform, power)
{ uv<-expand.grid(useq,vseq)
g.parametric <- exp((as.numeric(param[1])*
                       (1+(uv[,1]/as.numeric(param[2]))^as.numeric(param[4])+
                          (uv[,2]/as.numeric(param[3]))^as.numeric(param[5]))^(-as.numeric(param[6]))
))
ghat2<-as.vector(ghat)
return(sum((((ghat2))^power - (g.parametric)^power)^2))
}



.grid1.index <- function(x, xrange, nx) {
  i <- ceiling(nx * (x - xrange[1]) / diff(xrange))
  i <- pmax.int(1, i)
  i <- pmin.int(i, nx)
  i
}


.grid.index <- function(x, y, t, xrange, yrange, trange, nx, ny, nt) {

  ix <- .grid1.index(x, xrange, nx)
  iy <- .grid1.index(y, yrange, ny)
  it <- .grid1.index(t, trange, nt)

  return(list(ix = ix, iy = iy, it = it, index = as.integer((iy - 1) * nx + ix + (it - 1) * nx * ny)))
}



g.sep_st_exp_exp2 <- function(param, useq, vseq, ghat, transform, power){
  uv <- expand.grid(useq, vseq)
  g.parametric <- exp(param[1] * exp( - uv[, 2] / param[3]) * exp( - uv[, 1] / param[2]))
  ghat2 <- as.vector(ghat)
  return(sum((ghat2 ^ power - g.parametric ^ power) ^ 2))
}



kde2d.new.var <- function(x, y, gx = x, gy = y, np, var.bin = TRUE, y.var = TRUE, hx, hy){

  nx  <-  length(x)
  if(length(y) != nx)   stop("Data vectors must be the same length")
  if(var.bin){
    h <- numeric()
    if(y.var){
      mat.d <- as.matrix(dist(cbind(x, y), upper = TRUE, diag = TRUE))
    }
    if(!y.var){
      mat.d <- as.matrix(dist(cbind(x), upper = TRUE, diag = TRUE))
    }
    diag(mat.d) <- 999999
    if(missing(np))  np <- trunc(nx / 2)#
    np <- trunc(np)
    if(np <= 1) np <- 1
    for(i in 1:ncol(mat.d)){
      o <- order(mat.d[, i])
      pp <- o[np]
      h[i] <- mat.d[pp, i]
    }
    h <- h / 4
    hx <- h
    hy <- h
    rm(mat.d)
  } else{ np <- 0
  if(missing(hx) || missing(hy)){
    h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))
    h <- h / 4 / 1.06 * 0.9
    hx <- rep(h[1], nx)
    hy <- rep(h[2], nx)
  }
  }
  n <- length(gx)
  ax <- (outer(gx, x, "-") / hx)
  ay <- (outer(gy, y, "-") / hy)
  z <- rowSums(t(t(matrix(dnorm(ax), n, nx) * matrix(dnorm(ay), n, nx)) / (nx * hx * hy)))
  return(list(x = gx, y = gy, z = z, hx = hx, hy = hy, np = np))

}



norm2_etas <- function(n = 1, d = 1, q = 2, x0, y0){
  theta   = runif(n, max = 2 * pi)
  U       = runif(n)
  R       = sqrt(d * (U ^ (1 / (1 - q)) - 1))
  return(cbind(x0 + R * cos(theta), y0 + R * sin(theta)))
}


permutest.stlp <- function(perm, Q_n, Xi, L0, method){
  n_perm <- length(perm)
  n_points <- sum(perm == 2) + 1
  A <- sample(perm, size = n_perm, replace = FALSE, prob = NULL)
  cod <- which(A == 2)
  Q_perm <- Q_n[cod, 1:3]
  Q <- stlnpp::as.stlpp(c(Xi$x, Q_perm$x),
                c(Xi$y, Q_perm$y),
                c(Xi$t, Q_perm$t),
                L = L0)

  Q_dens <- rep(n_points / (spatstat.geom::volume(L0) * (Q$time[2] - Q$time[1])), n_points)
  lista_0 <- if (method == "K"){
    STLKinhom_i(Q, lambda = Q_dens)
  } else {
    STLginhom_i(Q, lambda = Q_dens)
  }
  lista_0
}



permutest.stp <- function(perm, Q_n, Xi, method){

  n_perm <- length(perm)
  n_points <- sum(perm == 2) + 1
  A <- sample(perm, size = n_perm, replace = FALSE, prob = NULL)
  cod <- which(A == 2)
  Q_perm <- Q_n[cod, 1:3]

  Q <- cbind(c(Xi$x, Q_perm$x),
             c(Xi$y, Q_perm$y),
             c(Xi$t, Q_perm$t))
  class(Q) <- "stpp"

  lista_0 <- if (method == "K"){
    KLISTAhat(Q)$list.KLISTA
  } else {
    LISTAhat(Q)$list.LISTA
  }
  lista_0
}




STLginhom_i <- function(X, lambda, normalize = FALSE, r = NULL, t = NULL, nxy = 10)
{
  if (!inherits(X, "stlpp"))
    stop("X should be from class stlpp")

  Y <- stlnpp::as.lpp.stlpp(X)
  l <- spatstat.geom::domain(Y)
  tleng <- summary(l)$totlength
  n <- spatstat.geom::npoints(Y)
  a <- X$time[1]
  b <- X$time[2]
  trange <- b - a
  timev <- X$data$t
  sdist <- spatstat.geom::pairdist(Y)
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

  g_local <- list()
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
    appx <- (n-1) * (tleng * trange) / (sum(revrho[lower.tri(revrho,diag = FALSE)]) * 2)
    gval <- lapply(g_local, "*" , appx)
  }
  else {
    g_local <- lapply(g_local, FUN= function(K) (n-1) * K / (trange * tleng) )
  }

  arr <- array( unlist(g_local) , c(nxy, nxy, n) )

  return(arr)
}




STLKinhom_i <- function(X, lambda = lambda, normalize = FALSE, r = NULL, t = NULL,
                        nxy = 10)
{
  if (!inherits(X, "stlpp"))
    stop("X should be from class stlpp")

  Y <- as.lpp.stlpp(X)
  l <- domain(Y)
  tleng <- summary(l)$totlength
  n <- npoints(X)
  a <- X$time[1]
  b <- X$time[2]
  trange <- b - a
  timev <- X$data$t
  sdist <- pairdist.lpp(Y)
  tdist <- as.matrix(dist(timev))
  toler <- default.linnet.tolerance(l)
  ml <- matrix(1, n, n)
  for (j in 1:n) {
    ml[-j, j] <- countends(l, Y[-j], sdist[-j, j], toler = toler)
  }
  mtplus <- matrix(timev, n, n, byrow = T) + tdist
  mtminus <- matrix(timev, n, n, byrow = T) - tdist
  mtedge <- (mtplus <= b) + (mtminus >= a)
  diag(mtedge) <- 1
  lamden <- outer(lambda, lambda, FUN = "*")
  diag(lamden) <- 1
  edgetl <- mtedge* ml * lamden

  maxs <- 0.7 * max(sdist[!is.infinite(sdist)])
  maxt <- 0.7 * (trange/2)
  if (is.null(r))
    r <- seq((maxs/nxy), maxs, by = (maxs - (maxs/nxy))/(nxy - 1))
  if (is.null(t))
    t <- seq((maxt/nxy), maxt, by = (maxt - (maxt/nxy))/(nxy - 1))

  K_local <- list()
  for(k in 1:n){
    K <- matrix(NA, nrow = nxy, ncol = nxy)
    for (i in 1:length(r)) {
      for (j in 1:length(t)) {
        out <- (sdist <= r[i]) * (tdist <= t[j])
        diag(out) <- 0
        kout <- out / edgetl
        kout <- kout[, k]
        K[i, j] <- sum(kout[is.finite(kout)])
      }
    }
    K_local[[k]] <- K
  }

  if (normalize) {
    revrho <- outer(1 / lambda, 1 / lambda, FUN = "*")
    appx <- (n - 1) * (tleng * trange) / (sum(revrho[lower.tri(revrho, diag = FALSE)]) * 2)
    K_local <- lapply(K_local, "*" , appx)#

  }
  else {
    K_local <- lapply(K_local, FUN= function(K) (n - 1) * K / (trange * tleng) )
  }

  arr <- array(unlist(K_local), c(nxy, nxy, n))

  return(arr)
}




































