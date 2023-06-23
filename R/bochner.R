#' Compute Bochner's Theorem
#'
#' @inheritParams pwelch
#' @param acf acf of the time-series
#' @param bias boolean value to bias as per Sykulski (2019)
#'
#' @return A vector of the positive spectral density,
#' excluding zero and Nyquist frequencies
#' @export
#'
#' @examples
bochner <- function(acf, delta = 1, bias = FALSE) {
    n <- length(acf)
    nfreq <- get_nfreq(n)

    if (bias) {
        acf <- (1 - (0:nfreq) / n) * acf
    }

    acf <- c(acf[1] / 2, acf[2:(n - 1)], acf[n] / 2)

    2 * delta * Re(stats::fft(acf))[1:nfreq + 1]

}