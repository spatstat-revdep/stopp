#' Create \code{stpm} and \code{stlpm} objects for marked point patterns storage
#'
#' This function creates a \code{stpm} object as a dataframe with 3 + m columns:
#' \code{x}, \code{y}, \code{t}, and \code{m} columns to store different marks.
#' If also the linear network \code{L}, of class linnet, is provided, a \code{stlp}
#' object is created instead.
#' 
#' @param df A matrix with three columns + m marks
#' @param names Factor string to name the marks columns.
#' @param L Optional. The linear network of class \code{linnet}.
#' If provided, the function returns a \code{stlpm} object.
#'
#' @return An \code{stpm} or \code{stlppm} object, depending on whether or not an object of class
#' linnet is provided for the \code{L} argument.
#' @export
#'
#' @author Nicoletta D'Angelo
#' 
#' @examples
#'
#'
#' set.seed(2)
#' df <- data.frame(cbind(runif(100), runif(100), runif(100), rpois(100, 15),
#' rpois(100, 30)))
#'
#' stpm1 <- stpm(df)
#'
#' ## Categorical marks
#' 
#' set.seed(2)
#' dfA <- data.frame(x = runif(100), y = runif(100), t = runif(100), 
#'                   m1 = rnorm(100), m2 = rep(c("C"), times = 100))
#' dfB <- data.frame(x = runif(50), y = runif(50), t = runif(50), 
#'                   m1 = rnorm(25), m2 = rep(c("D"), times = 50))
#' 
#' stpm2 <- stpm(rbind(dfA, dfB), names = c("continuous", "dichotomous"))
#' 
#' ## Linear network
#' 
#' set.seed(2)
#' dfL <- data.frame(cbind(runif(100, 0, 0.85), runif(100, 0, 0.85), runif(100), 
#'                        rpois(100, 15)))
#' 
#' stlpm1 <- stpm(dfL, L = chicagonet)
#' 
#'
#'
stpm <- function(df, names = NULL, L) {
  if (!inherits(df, c("matrix", "data.frame"))){
    stop("df should be either a matrix or a data.frame")
  } 
  
  if (!inherits(df, c("data.frame"))) df <- data.frame(df)
  if(is.null(names)){
    nm <- ncol(df) - 3
    colnames(df) <- c("x", "y", "t", paste0("m", 1:nm)) 
  } else {
    colnames(df) <- c("x", "y", "t", names) 
  }
  if(missing(L)){
    out <- list(df = df)
    class(out) <- "stpm"
  } else {
    
    if (!inherits(L, c("linnet"))){
      stop("L should be a linnet object")
    } 
    
    suppressWarnings(lpp_new <- spatstat.linnet::lpp(df, L))
    df_new <- data.frame(x = lpp_new$data[, 1], y = lpp_new$data[, 2 ])
    df[, 1:2] <- df_new
    out <- list(df = df, L = L)
    class(out) <- "stlpm"
  }
  return(out)
}
