#' Fit a separable spatio-temporal Poisson process model on a linear network
#'
#' @param x A dataframe
#' @param spaceformula A formula for the spatial component. See \link{lppm} for details
#' @param timeformula A formula for the temporal component. It fits a log-linear model with the \link{glm} function
#' @param L A linear network of class \code{linnet}
#'
#' @return An object of class \code{sepstlppm}
#' @export
#'
#' @examples
#' \dontrun{
#' mod1 <- sepstlppm(valenciacrimes[1:2500, ], spaceformula = ~x,
#' timeformula = ~ crime_hour + week_day, L = valencianet)
#'
#'
#' }
sepstlppm <- function(x, spaceformula, timeformula, L){

  if (!inherits(x, "data.frame")) stop("x should a dataframe")

  ot <- x$t

  n <- nrow(x)

  X <- spatstat.linnet::lpp(spatstat.geom::ppp(x$x, x$y,
                            window = spatstat.geom::owin(range(x$x), range(x$y))), L = L)
  spacemod <- spatstat.linnet::lppm(as.formula(paste("X", paste(spaceformula, collapse = " "), sep = " ")))
  spaceint <- predict(spacemod, locations = X)
  spaceint_plot <- predict(spacemod)

  stab.expanse1 <- table(x[, all.vars(timeformula)]) |> as.data.frame()
  if(length(all.vars(timeformula))  == 1){
    colnames(stab.expanse1)[1] <- attr(terms(timeformula), "term.labels")
  }
  timemod <- glm(as.formula(paste("Freq", paste(timeformula, collapse = " "), sep = " ")),
                 data = stab.expanse1, family = poisson)

  newdata0 <- data.frame(apply(x[, all.vars(timeformula), drop = F], 2, as.factor))
  timeint <- exp(predict(timemod, newdata = newdata0))

  stint <-  spaceint * timeint

  stint <-  as.numeric(stint / sum(stint) * n)


  out = list(l = stint, x = x, spaceformula = spaceformula, timeformula = timeformula,
             spacemod = spacemod, timemod = timemod, spaceint_plot = spaceint_plot,
             spaceint = spaceint, timeint = timeint)
  class(out) <- "sepstlppm"
  return(out)

}



