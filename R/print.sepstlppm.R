
#' Print of a fitted separable spatio-temporal Poisson process model on a linear network
#'
#'  The function prints the main information of the fitted model.
#'
#' @param x An object of class \code{sepstlppm }
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
#' mod1
#' 
#'}
#'
#'
#'
print.sepstlppm <- function(x, ...){
  if (!inherits(x, c("sepstlppm"))) stop("X should be from class sepstlppm")
  
    cat("Fitted separable spatio-temporal Poisson process model \n")
    cat("on a linear network \n")
}


