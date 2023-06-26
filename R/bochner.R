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
bochner <- function(acf, delta = 1, h = NULL) {
    n <- length(acf)
    nfreq <- get_nfreq(n)

    if (is.null(h)) {
        acf <- (1 - (0:(n - 1)) / n) * acf
    } else {
        h <- h / sqrt(sum(h^2)) # normalise h
        h_conv <- stats::convolve(h, h, type = "open")[n:(2 * n - 1)]
        acf <- h_conv * acf
    }

    acf <- c(acf[1] / 2, acf[2:(n - 1)], acf[n] / 2)

    2 * delta * Re(stats::fft(acf))[2:(nfreq + 1)]

}