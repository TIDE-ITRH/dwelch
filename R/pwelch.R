#' Calculate Welch's Estimator
#'
#' @param ts time-series as a vector
#' @param m number of blocks
#' @param l length of blocks
#' @param s shift factor
#' @param delta sampling interval
#'
#' @return A tibble containing Welch's estimate
#' @export
#'
#' @examples
pwelch <- function(ts, m, l, s, delta = 1) {

    is_even <- l %% 2 == 0
    n <- (m - 1) * s + l

    if (length(ts) != n) {
        warning("n does not equal length of ts")
        stop() # Change to zero-pad and continue
    }

    if (is_even) {
        n_spec <- floor(l / 2) - 1
    } else {
        n_spec <- floor(l / 2)
    }

    ff <- (1:(n_spec)) / (l * delta)
    pxx <- matrix(nrow = n_spec, ncol = m)

    for (i in 1:m){
        start <- (i - 1) * s + 1
        end <- start + l - 1
        ts_tmp <- ts[start:end]
        pxx[, i] <- 1 / l * (Mod(stats::fft(ts_tmp))^2)[2:(n_spec + 1)]
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