#' Summary of the diagnostics performed on a spatio-temporal point process model
#'
#' It summarises the main information on the result of the local diagnostics
#'  performed with \link{localdiag} on either a \code{stp} or \code{stlp} object:
#'  whether the local test was run on point patterns lying on a linear network or not;
#'  the number of points in the analysed spatio-temporal
#'  point pattern \code{X};
#'  the number of points of \code{X} which are identified as outlying by the
#'  previously performed local diagnostics.
#'
#' @param object A \code{localdiag} object
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo and Giada Adelfio
#'
#' @seealso
#' \link{infl}, \link{plot.localdiag}, \link{print.localdiag}
#'
#'
#' @examples
#' \dontrun{
#'
#' #load data
#' set.seed(12345)
#' id <- sample(1:nrow(catalog.withcov), 200)
#' cat <- catalog.withcov[id, ]
#' stp1 <- stp(cat[, 5:3])
#'
#' #fit two competitor models
#' # and extract the fitted spatio-temporal intensity
#'
#' lETAS <- etasclass(cat.orig = cat, magn.threshold = 2.5, magn.threshold.back = 3.9,
#' mu = 0.3, k0 = 0.02, c = 0.015, p = 1.01, gamma = 0, d = 1,
#' q = 1.5, params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE,
#'                         TRUE), formula1 = "time ~  magnitude- 1",
#'                         declustering = TRUE,
#'                         thinning = FALSE, flp = TRUE, ndeclust = 15, onlytime = FALSE,
#'                         is.backconstant = FALSE, sectoday = FALSE, usenlm = TRUE,
#'                         compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36)$l
#'
#' lPOIS <- etasclass(cat.orig = cat, magn.threshold = 2.5, magn.threshold.back = 3.9,
#' mu = 0.3, k0 = 0.02, c = 0.015, p = 1.01, gamma = 0, d = 1,
#' q = 1.5, params.ind = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
#'                         FALSE), formula1 = "time ~  magnitude- 1",
#'                         declustering = TRUE,
#'                         thinning = FALSE, flp = TRUE, ndeclust = 15, onlytime = FALSE,
#'                         is.backconstant = FALSE, sectoday = FALSE, usenlm = TRUE,
#'                         compsqm = TRUE, epsmax = 1e-04, iterlim = 100, ntheta = 36)$l
#'
#' # let's identify the outlying points at a .9 percentile
#'
#' resETAS <- localdiag(stp1, lETAS, p = .9)
#' resPOIS <- localdiag(stp1, lPOIS, p = .9)
#'
#' summary(resETAS)
#'
#' # Points outlying from the 0.9 percentile
#' # of the anaysed spatio-temporal point pattern
#' # --------------------------------------
#' #   Analysed pattern X: 200 points
#' # 20 outlying points
#'
#'
#' summary(resPOIS)
#'
#' # Points outlying from the 0.9 percentile
#' # of the anaysed spatio-temporal point pattern
#' # --------------------------------------
#' #   Analysed pattern X: 200 points
#' # 20 outlying points
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
#'
#'
#'
#'
summary.localdiag <- function(object, ...){
  if(!any(class(object) == "localdiag")) stop("class(object) must be localdiag")
  if(inherits(object$X, "stlp")){
    cat(paste("Points outlying from the", object$p , "percentile\n"))
    cat("of the anaysed spatio-temporal point pattern on a linear network \n")
    cat("--------------------------------------------------\n")
  } else {
    cat(paste("Points outlying from the", object$p , "percentile\n"))
    cat("of the anaysed spatio-temporal point pattern \n")
    cat("--------------------------------------\n")
  }
  cat(paste("Analysed pattern X:",nrow(object$X$df), "points \n"))
  cat(paste(length(object$ids), "outlying points"))
}
