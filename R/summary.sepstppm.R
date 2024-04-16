
#' Summary of a fitted separable spatio-temporal Poisson process model 
#'
#'  The function summarises the main information of the fitted model.
#'
#' @param object An object of class \code{sepstppm}
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
#'
#' crimesub <- stpm(valenciacrimes$df[101:200, ],
#'            names = colnames(valenciacrimes$df)[-c(1:3)])
#' 
#' mod1 <- sepstppm(crimesub, spaceformula = ~x ,
#'                   timeformula = ~ day)
#'                   
#'summary(mod1)
#'
#'
#'
summary.sepstppm <- function(object, ...){

  cat("Fitted separable spatio-temporal Poisson process model \n")
  cat("with spatial estimates: \n \n")
  print(summary(object$spacemod)$coefs.SE.CI)
  cat("\n")
  cat("and temporal estimates: \n \n")
  print(summary(object$timemod)$coefficients)
    cat("\n")

}

