#' Calculate Welch's Estimator
#'
#' @param ts time-series as a vector
#' @param m number of blocks
#' @param l length of blocks
#' @param s shift factor
#' @param delta sampling interval
#' @param window optional window to apply to data. Defaults to rectangular
#'
#' @return A tibble containing Welch's estimate
#' @export
#'
#' @examples
pwelch <- function(ts, m, l, s, delta = 1, window = NULL) {

    n <- (m - 1) * s + l
    nfreq <- get_nfreq(l)

    if (length(ts) != n) {
        warning("n does not equal length of ts")
        stop() # Change to zero-pad and continue
    }

    if (is.null(window)) {
        window <- rep(1, l)
    }

    ff <- (1:(nfreq)) / (l * delta)
    pxx <- matrix(nrow = nfreq, ncol = m)

    window <- window / sqrt(sum(window^2))

    for (i in 1:m){
        start <- (i - 1) * s + 1
        end <- start + l - 1
        ts_tmp <- ts[start:end] * window
        pxx[, i] <- (Mod(stats::fft(ts_tmp))^2)[2:(nfreq + 1)]
    }

    welch_estimate <- rowMeans(pxx)

    dplyr::tibble(
        ff = ff,
        welch = welch_estimate
    )
}

# Notes for future development
# Add some tests and warning messages
# Change so it zero-pads the last one
# Add functionality for percentage overlap