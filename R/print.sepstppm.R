
#' Print of a fitted separable spatio-temporal Poisson process model
#'
#'  The function prints the main information of the fitted model.
#'
#' @param x An object of class \code{sepstppm}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{sepstppm}
#'
#'
#'
#' @examples
#' \dontrun{
#' crimesub <- stpm(valenciacrimes$df[101:200, ],
#'            names = colnames(valenciacrimes$df)[-c(1:3)])
#' 
#' mod1 <- sepstppm(crimesub, spaceformula = ~x ,
#'                   timeformula = ~ day)
#' mod1
#'}
#'
#'
#' 
#'
print.sepstppm  <- function(x, ...){
  if (!inherits(x, c("sepstppm"))) stop("X should be from class sepstppm")
  cat("Fitted separable spatio-temporal Poisson process model \n")
}


