#' Fit a separable spatio-temporal Poisson process model
#'
#' @param x A dataframe
#' @param spaceformula A formula for the spatial component. See \link{ppm} for details
#' @param timeformula A formula for the temporal component. It fits a log-linear model with the \link{glm} function
#'
#' @return An object of class \code{sepstppm}
#' @export
#'
#' @examples
#' \dontrun{
#' df1 <- valenciacrimes[valenciacrimes$x < 210000 & valenciacrimes$x > 206000
#' & valenciacrimes$y < 4377000 & valenciacrimes$y > 4373000, ]
#'
#' mod1 <- sepstppm(df1, spaceformula = ~x * y,
#'                 timeformula = ~ crime_hour + week_day)
#' }
#'
#'
#'
sepstppm <- function(x, spaceformula, timeformula){

  if (!inherits(x, "data.frame")) stop("x should a dataframe")

  ot <- x$t

  n <- nrow(x)

  X <- spatstat.geom::ppp(x$x, x$y,
                          spatstat.geom::owin(range(x$x), range(x$y)))
  spacemod <- spatstat.model::ppm(as.formula(paste("X", paste(spaceformula, collapse = " "), sep = " ")))
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
             spacemod = spacemod, timemod = timemod, spaceint = spaceint, timeint = timeint)
  class(out) <- "sepstppm"
  return(out)

}
