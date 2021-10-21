#' Compute standard error of a sample mean
#'
#' Given a distribution, returns its standard error.
#'
#' @param x A vector of values.
#' @return The standard error as a numeric value.
#' @export
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
