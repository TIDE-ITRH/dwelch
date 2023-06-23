#' Calculate Sinc Function
#'
#' @param x value
#'
#' @return value of sinc(x)
#' @export
#'
#' @examples
sinc <- function(x) {
    value <- sin(x) / x
    value[1] <- 1
    value
}