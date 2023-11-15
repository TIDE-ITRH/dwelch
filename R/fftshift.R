#' Shift FFT output to be negative/positive
#'
#' @param x vector
#'
#' @return A shifted vector
#' @export
#'
#' @examples
fftshift <- function(x) {
    len <- length(x)
    is_even <- len %% 2 == 0
    half <- floor(len / 2)

    if (is_even) {
        x <- c(x[(half + 1):len], x[1:half])
    } else {
        x <- c(x[(half + 2):len], x[1:(half + 1)])
    }

    x
}