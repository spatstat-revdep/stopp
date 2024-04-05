
#' Summary of a fitted fitted separable spatio-temporal Poisson process model on a linear network
#'
#'  The function summarises the main information of the fitted model.
#'
#' @param object An object of class \code{sepstlppm}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{sepstlppm}
#'
#'
#' @examples
#' \dontrun{
#' crimesub <- stpm(valenciacrimes$df[101:200, ],
#'            names = colnames(valenciacrimes$df)[-c(1:3)],
#' L = valencianet)
#' 
#' mod1 <- sepstlppm(crimesub, spaceformula = ~x ,
#'                   timeformula = ~ day)
#'                   
#' summary(mod1)
#'}
#'
#'
#'
summary.sepstlppm <- function(object, ...){
  if(!inherits(object, "sepstlppm")) stop("class(object) must be sepstlppm")
  
  
  cat("Fitted separable spatio-temporal Poisson process model \n")
  cat("on a linear network \n \n")
  cat("with spatial estimates: \n \n")
  print(summary(object$spacemod$fit)$coefs.SE.CI)
  cat("\n")
  cat("and temporal estimates: \n \n")
  print(summary(object$timemod)$coefficients)
  cat("\n")
 
}

