#' Test of local structure for spatio-temporal point processes
#'
#' @description
#' This function performs the permutation test of the local structure for
#' spatio-temporal point pattern data, proposed in Siino et al. (2018), as well as for
#' spatio-temporal point pattern data occurring on the same linear network,
#'  following
#' D'Angelo et al. (2021).
#'
#' @param X Background spatio-temporal point pattern. Usually, the most clustered
#' between \code{X} and \code{Z}. Must be either a \code{stp} or \code{stlp} object.
#' @param Z Other spatio-temporal point pattern. Must also be of the same class as \code{X}.
#' @param method Character string indicating which version of LISTA function to use:
#' either
#'  \code{"K"} or \code{"g"}.
#' If \code{"K"}, the local spatio-temporal K-function is used to run the test.
#' If \code{"g"}, the local spatio-temporal pair correlation function is used.
#' @param k Number of permutations
#' @param alpha Significance level
#' @param verbose If TRUE (default) the progress of the test is printed
#'
#'
#' @details
#' The test detects local differences between  \eqn{\textbf{x}} and \eqn{\textbf{z}}
#' occurring on the same space-time region.
#' This procedure was firstly introduced in Moraga and Montes (2011), extended in
#' the spatio-temporal context by Siino at al, (2018). Finally,
#'  test has been made suitable also for spatio-temporal point patterns
#' with spatial domain coinciding with a linear network by D'Angelo et al. (2021).
#'
#'
#' In general, for each point \eqn{(\textbf{u},t)} in the spatio-temporal observed
#' point pattern \eqn{\textbf{x}}, we test
#' \deqn{
#' \begin{cases}
#' \mathcal{H}_{0}: & \text{no difference in the second-order local  structure of }  (\textbf{u},t) \quad   \text{   w.r.t  } \quad \{  \{  \textbf{x} \setminus    (\textbf{u},t)   \} \cup  \textbf{z} \}\\
#' \mathcal{H}_{1}: &  \text{significant difference in the second-order local  }   \text{structure of} (\textbf{u},t) \quad \text{   w.r.t  } \quad \{  \{  \textbf{x} \setminus    (\textbf{u},t)   \} \cup  \textbf{z} \}
#' \end{cases}}
#'
#' The sketch of the test is as follows:
#'
#'
#' 1. Set \eqn{k} as the number of permutations
#'
#' 2. For each point  \eqn{(\textbf{u}_i,t_i) \in \textbf{x}, i = 1, \ldots, n}:
#'
#' 2.1. Estimate the LISTA function  \eqn{\hat{L}^{(i)}(r,h)}
#'   and Compute the local  deviation test
#'   \eqn{T^i=\int_{0}^{t_0} \int_{0}^{r_0} \Big(
#'     \hat{L}^{(i)}(r,h)- \hat{L}^{-(i)}_{H_0}(r,h)
#'     \Big)^2 \text{d}r \text{d}h,
#'   }
#'     where  \eqn{\hat{L}^{-(i)}_{H_0}(r,h)}
#'      is the LISTA function for the \eqn{i^{th}} point,
#'       averaged over the \eqn{j=1,\dots,k} permutations
#'
#' 2.2 Compute a \eqn{p}-value as
#'       \eqn{p^i=\sum_{j=1}^{k}  \textbf{1}(T^{i,j}_{H_0} \geq T^i)/k}
#'
#'
#'    The test ends providing a vector \eqn{p} of  \eqn{p}- values, one for each point
#'    in \eqn{\textbf{x}}.
#'
#'  If the test is  performed for spatio-temporal point patterns as in
#'    Siino et al. (2018), that is, on an object of class \code{stp}, the LISTA
#'    functions \eqn{\hat{L}^{(i)}} employed are the local functions of
#'    Adelfio et al. (2020), documented in
#'  \link{KLISTAhat} and \link{LISTAhat} of the \code{stpp} package (Gabriel et al, 2013).
#'
#'  If the function is applied to a \code{stlp} object, that is, on two spatio-temporal
#'  point patterns observed on the same linear network \code{L}, the LISTA function
#' \eqn{\hat{L}^{(i)}} used are the ones proposed in D'Angelo et al. (2021), documented
#' in  \link{localSTLKinhom} and \link{localSTLginhom}.
#'
#' Details on the performance of the test are found in Siino et al. (2018) and
#' D'Angelo et al. (2021), for Euclidean and network spaces, respectively.
#'
#'
#' @return A list of class \code{localtest}, containing
#' \describe{
#' \item{\code{p}}{A vector of p-values, one for each of the points in \code{X}}
#' \item{\code{X}}{The background spatio-temporal point pattern given in input}
#' \item{\code{Z}}{The alternative spatio-temporal point pattern given in input}
#' \item{\code{alpha}}{The threshold given in input}
#' \item{\code{Xsig}}{A \code{stp} object storing the resulting significant points}
#' \item{\code{Xnosig}}{A \code{stp} object storing the resulting non-significant points}
#' \item{\code{id}}{The ids of the resulting significan points}
#' }
#'
#'
#' @export
#'
#' @author Nicoletta D'Angelo and Marianna Siino
#'
#' @seealso
#' \link{print.localtest}, \link{summary.localtest}, \link{plot.localtest}
#'
#'
#' @examples
#'
#'
#'
#'\dontrun{
#' # background pattern
#' set.seed(12345)
#' X <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(.05, 4),
#'            nsim = 1, seed = 2, verbose = TRUE)
#'
#' # alternative pattern
#' set.seed(12345)
#' Z <- rstpp(lambda = 25, nsim = 1, seed = 2, verbose = T)
#'
#'
#' # run the local test
#' test <- localtest(X, Z, method = "K", k = 9)
#'
#' ## Example on a linear network
#'
#' # fix the linear network
#' L0 <- spatstat.geom::domain(spatstat.data::chicago)
#'
#' # background pattern
#' set.seed(12345)
#' X <- retastlp(cat = NULL, params = c(0.078915 / 1.95, 0.003696,  0.013362,  1.2,
#'                                         0.424466,  1.164793),
#'                  betacov = 0.5, m0 = 2.5, b = 1.0789, tmin = 0, t.lag = 200,
#'                  xmin = 600, xmax = 2200, ymin = 4000, ymax = 5300,
#'                  iprint = TRUE, covdiag = FALSE, covsim = FALSE, L = L0)
#'
#' # alternative pattern, on the same linear network
#' l <- 20 / (spatstat.geom::volume(L0) * (200 - 25))
#' set.seed(12345)
#' stlppPOIS <- stlnpp::rpoistlpp(lambda = l, a = 25, b = 200, L = L0)
#' Z <- as.stlp(stlppPOIS)
#'
#' # run the local test
#' test <- localtest(X, Z, method = "K", k = 9)
#'
#'
#'
#'}
#'
#'
#' @references
#' Adelfio, G., Siino, M., Mateu, J., and Rodríguez-Cortés, F. J. (2020). Some properties of local weighted second-order statistics for spatio-temporal point processes. Stochastic Environmental Research and Risk Assessment, 34(1), 149-168.
#'
#' D’Angelo, N., Adelfio, G., and Mateu, J. (2021). Assessing local differences between the spatio-temporal second-order structure of two point patterns occurring on the same linear network. Spatial Statistics, 45, 100534.
#'
#' Gabriel, E., Rowlingson, B. S., and Diggle, P. J. (2013). stpp: An R Package for Plotting, Simulating and Analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53(2), 1–29. https://doi.org/10.18637/jss.v053.i02
#'
#' Moraga, P. and Montes, F. (2011). Detection of spatial disease clusters with lisa functions. Statistics in Medicine, 30(10):1057–1071
#'
#' Siino, M., Rodríguez‐Cortés, F. J., Mateu, J. ,and Adelfio, G. (2018). Testing for local structure in spatiotemporal point pattern data. Environmetrics, 29(5-6), e2463.
#'
localtest <- function(X, Z, method = "K", k, alpha = 0.05, verbose = TRUE){

  nX <- nrow(X$df)
  nZ <- nrow(Z$df)
  perm <- c(rep(1, nZ), rep(2, nX - 1))
  p <- double(nX)

  if ((inherits(X,"stlp")|inherits(Z,"stlp"))) L0 <- X$L

  res_test <- double(nX)
  res_test_d <- matrix(NA, k, nX)

  if ((inherits(X,"stp")|inherits(Z,"stp"))){

    X_lista <- if (method == "K"){
      stpp::KLISTAhat(as.stpp(X))$list.KLISTA
    } else if (method == "g"){
      stpp::LISTAhat(as.stpp(X))$list.LISTA
    } else {stop(" 'method' argument must be either \"g\" or \"K\ ")}

  } else if ((inherits(X,"stlp")|inherits(Z,"stlp"))){
    X_lista <- if (method == "K"){
      STLKinhom_i(as.stlpp(X), lambda = rep(nX / (spatstat.geom::volume(L0) * (range(X$df$t)[2] - range(X$df$t)[1])), nX))
    } else if (method == "g"){
      STLginhom_i(as.stlpp(X), lambda = rep(nX / (spatstat.geom::volume(L0) * (range(X$df$t)[2] - range(X$df$t)[1])), nX))
    } else {stop(" 'method' argument must be either \"g\" or \"K\ ")}

  }

  lista_0 <- array(NA, dim = c(dim(X_lista)[1:3], k))

  for (i in 1:nX){

    print(i)

    a <- cbind(Z$df, 1)
    b <- cbind(X$df[-i, ], 2)
    colnames(a) <- colnames(b) <- c("x", "y", "t", "id")

    Q_n <- data.frame(rbind(a, b))
    Xi <- X$df[i, ]

    if(verbose) pb <- utils::txtProgressBar(min = 0, max = k, style = 3)
    lista_0 <- sapply(1:k, function(q) {
      if(verbose)  utils::setTxtProgressBar(pb, q)
      if ((inherits(X,"stp")|inherits(Z,"stp"))){
        permutest.stp(perm = perm, Q_n = Q_n, Xi = Xi, method = method)
      } else if ((inherits(X,"stlp")|inherits(Z,"stlp"))){
        permutest.stlp(perm = perm, Q_n = Q_n, Xi = Xi, L0 = L0, method = method)
      }
    }, simplify = "array")
    if(verbose) close(pb)
  }

  mlistas <- apply(lista_0, 1:3, mean, na.rm = TRUE)

  for (i in 1:nX){
    lista_a <- X_lista[, , i]
    lista_b <- mlistas[, , i]
    t2 <- sum((lista_a - lista_b) ^ 2)
    res_test[i] <- t2
    for (j in 1:k){
      lista_a <- lista_0[, , i, j]
      t2_0 <- sum((lista_a - lista_b) ^ 2)
      res_test_d[j,i] <- t2_0
    }
    p[i] <- mean(res_test_d[,i] >= res_test[i])
  }

  id <- which(p <= alpha)

  if ((inherits(X,"stp")|inherits(Z,"stp"))){
    Xsig <- stp(X$df[id, ])
    Xnosig <- stp(X$df[- id, ])
  } else if ((inherits(X,"stlp")|inherits(Z,"stlp"))){
    Xsig <- stp(X$df[id, ], L = L0)
    Xnosig <- stp(X$df[- id, ], L = L0)
  }


  out <- list(p = p, X = X, Z = Z, alpha = alpha, Xsig = Xsig, Xnosig = Xnosig,
              id = id)
  class(out) <- "localtest"
  return(out)
}








