#' Calculate AR spectrum
#'
#' @inheritParams pwelch
#' @param ff frequencies at which to compute spectrum
#' @param phis autoregressive parameters
#' @param sd standard deviation of the noise
#'
#' @return A vector of the AR spectrum corresponding the ff frequencies
#' @export
#'
#' @examples
ar_spectrum <- function(ff, phis, sd, delta = 1) {
    d <- 1
    d_conj <- 1

    ff <- ff * delta

    for (j in seq_along(phis)) {
        d <- d - phis[j] * exp(-2i * j * pi * ff)
        d_conj <- d_conj - phis[j] * exp(2i * j * pi * ff)
    }

    sd^2 / Re(d * d_conj)

}