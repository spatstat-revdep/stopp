#' Simulation of a spatio-temporal ETAS (Epidemic Type Aftershock Sequence) model
#'
#'
#'
#' @description
#' This function simulates a spatio-temporal ETAS
#' (Epidemic Type Aftershock Sequence) process as a \code{stpm} object.
#'
#' It follows the generating scheme for simulating a pattern from an
#' Epidemic Type Aftershocks-Sequences (ETAS) process
#' (Ogata and Katsura 1988) with conditional intensity function (CIF) as in
#' Adelfio and Chiodi (2020), adapted for the space location of events
#' to be constrained.
#'
#' See the 'Details' section.
#'
#'
#' @details
#' The CIF of an ETAS
#' process as in Adelfio and Chiodi (2020) can be written as \deqn{
#'   \lambda_{\theta}(t,\textbf{u}|\mathcal{H}_t)=\mu f(\textbf{u})+\sum_{t_j<t} \frac{\kappa_0 \exp(\eta_j)}{(t-t_j+c)^p} \{    (\textbf{u}-\textbf{u}_j)^2+d   \}^{-q}  ,
#'} where
#'
#' \eqn{\mathcal{H}_t} is the past history of the process up to time
#' \eqn{t}
#'
#' \eqn{\mu} is the large-scale general intensity
#'
#' \eqn{f(\textbf{u})} is
#' the spatial density
#'
#' \eqn{\eta_j=\boldsymbol{\beta}' \textbf{Z}_j} is a linear predictor
#'
#' \eqn{\textbf{Z}_j} the external known covariate vector, including the
#' magnitude
#'
#' \eqn{\boldsymbol{\theta}= (\mu, \kappa_0, c, p, d, q, \boldsymbol{\beta})}
#' are the parameters to be estimated
#'
#' \eqn{\kappa_0} is a
#' normalising constant
#'
#' \eqn{c} and \eqn{p} are characteristic parameters of the
#' seismic activity of the given region,
#'
#'  and \eqn{d} and \eqn{q} are two parameters
#' related to the spatial influence of the mainshock
#'
#' In the usual ETAS
#' model for seismic analyses, the only external covariate represents the magnitude,
#' \eqn{\boldsymbol{\beta}=\alpha}, as
#' \eqn{\eta_j = \boldsymbol{\beta}' \textbf{Z}_j = \alpha (m_j-m_0)}, where
#' \eqn{m_j} is the magnitude of the \eqn{j^{th}} event and \eqn{m_0} the threshold
#' magnitude, that is, the lower bound for which earthquakes with higher
#' values of magnitude are surely recorded in the catalogue.
#'
#'
#' @param pars A vector of parameters of the ETAS model to be simulated.
#'  See the 'Details' section.
#' @param betacov Numerical array. Parameters of the ETAS model covariates.
#' @param m0 Parameter for the background general intensity of the ETAS model.
#' In the common seismic analyses it represents the threshold
#' magnitude.
#' @param b 1.0789
#' @param t.lag 200
#' @param tmin Minimum value of time.
#' @param xmin Minimum of x coordinate range
#' @param xmax Maximum of x coordinate range
#' @param ymin Minimum of y coordinate range
#' @param ymax Maximum of y coordinate range
#' @param covsim Default \code{FALSE}
#' @param all.marks Logical value indicating whether to store
#' all the simulation information as marks in the \code{stpm} object.
#' If \code{FALSE} (default option) only the magnitude is returned.
#' @param seed Seed to set, if ones wished to reproduce the analyses
#' @return A \code{stpm} object
#' @export
#'
#' @author Nicoletta D'Angelo and Marcello Chiodi
#'
#'
#' @examples
#'
#'\dontrun{
#' X <- rETASp(pars = c(0.1293688525, 0.003696, 0.013362, 1.2,0.424466,  1.164793),
#'          betacov = 0.5, 
#'          xmin = 600, xmax = 2200, ymin = 4000, ymax = 5300,
#'          seed = 95)
#' 
#' plot(X)
#' 
#'}
#'
#' @references
#' Adelfio, G., and Chiodi, M. (2021). Including covariates in a space-time point process with application to seismicity. Statistical Methods & Applications, 30(3), 947-971.
#'
#' Ogata, Y., and Katsura, K. (1988). Likelihood analysis of spatial inhomogeneity for marked point patterns. Annals of the Institute of Statistical Mathematics, 40(1), 29-39.
#'
rETASp <-function(pars=NULL,
                           betacov=0.39,m0=2.5,b=1.0789,tmin=0,t.lag=200,
                  xmin=0,xmax=1,ymin=0,ymax=1,
                  covsim=FALSE, all.marks = FALSE,
                  seed = NULL){
  set.seed(seed)
  if(is.null(pars)) stop("Please provide some parameters")
  
  cat=NULL

  mu	= pars[1]
  k0  	= pars[2]
  c   	= pars[3]
  p   	= pars[4]
  gamma   =0
  d       = pars[5]
  q       = pars[6]
  ncov    = length(betacov)
  beta    =log(10)*b
  mm      =(ymax-ymin)/(xmax-xmin)
  qq      =ymin-mm*xmin
  cat.pois=NULL
  cat.sim =NULL
  cat.new =NULL
  tmax    =tmin+t.lag
  ak      =k0*c^(1-p)/(p-1)
  sk      =(pi*d^(1-q))/(q-1)
  muback  =mu*(tmax-tmin)
  n0      =rpois(1,muback)
  nstart  =n0

  if(nstart>0){

    lgen       =array(0,nstart)
    ind     =array(0,nstart)
    father  =array(0,nstart)
    
    tback   =runif(n0,tmin,tmax)
    xback   =runif(n0,xmin,xmax)
    yback   =runif(n0,ymin,ymax)
    mback   =m0+rexp(n0,rate=beta)
    zback   =array(0,n0)
    if(ncov==1){
      cov2back=array(0,n0)
    }
    else
    {
      if(covsim)
      {
        cov2back=rnorm(n0)^2
      }
      else
      {
        cov2back=abs(yback-mm*xback-qq)/(sqrt(1+mm*mm))
      }
    }


    cat.new=cbind(tback,xback,yback,mback,zback,cov2back)
    cat.pois=cat.new

    cat.new             =cbind(cat.new,lgen,ind,father)
    colnames(cat.new)   =c("time","long","lat","magn1","z","cov2","lgen","ind","father")
    cat.new             =as.data.frame(cat.new)
    cat.new=cat.new[order(cat.new$time),]
    i=0
    while(!prod(cat.new$ind)){
      i   =i+1
      cat.new$ind[i]  =TRUE
      pred            =(cat.new$magn1[i]-m0)*betacov[1]
      if(ncov>1) pred=pred+cat.new$cov2[i]*betacov[2]
      nexpected       =ak*sk*exp(pred)
      ni              =rpois(1,nexpected)
      if(ni>0){
        xy              =norm2_etas(ni,d,q,cat.new$long[i],cat.new$lat[i])
        t             =c*runif(ni)^(-1/(p-1))-c+cat.new$time[i]
        ind.txy=(t>tmin)&(t<tmax)&(xy[,1]>xmin)&(xy[,1]<xmax)&(xy[,2]>ymin)&(xy[,2]<ymax)
        ntxy=sum(ind.txy)
        if(ntxy>0){
          m1      =m0+rexp(ntxy,rate=beta)
          z1      =array(0,ntxy)
          cov2    =array(0,ntxy)
          lgen    =array(1,ntxy)+cat.new$lgen[i]
          ind     =array(0,ntxy)
          father  =array(i,ntxy)
          if (ncov>1)
          {

            if(covsim)
            {
              cov2=rnorm(ntxy)^2
            }
            else
            {
              xo      =xy[ind.txy,1]
              yo      =xy[ind.txy,2]
              cov2    =abs(yo-mm*xo-qq)/(sqrt(1+mm*mm))
            }


          }
          cat.son =cbind(t[ind.txy],xy[ind.txy,1],xy[ind.txy,2],m1,z1,cov2,lgen,ind,father)
          colnames(cat.son)   =c("time","long","lat","magn1","z","cov2","lgen","ind","father")
          cat.son =as.data.frame(cat.son)
          #               print(cat.son)
          cat.new =rbind(cat.new,cat.son)
          cat.new=cat.new[order(cat.new$time),]
          colnames(cat.new)   = colnames(cat.son) 

        }
      }
    }
  }
  nson=nrow(cat.new)-n0
  nson = nrow(cat.new) - n0
  cat.new[, 1:3] <- cat.new[c(2, 3, 1)]
  if(all.marks){
    return(stpm(cat.new, names = colnames(cat.new)[-c(1:3)]))
  } else {
    return(stpm(cat.new[, 1:4], names = "Magnitude"))
  }
}




