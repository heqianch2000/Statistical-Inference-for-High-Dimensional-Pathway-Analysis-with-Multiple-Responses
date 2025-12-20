#' sign function
#'
#' @param x a numeric value
#' @returns 1 if positive and -1 if negative
tsign <- function(x){
  ifelse(x >= 0, 1, -1)
}
