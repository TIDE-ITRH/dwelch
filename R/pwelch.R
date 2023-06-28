#' Calculate Welch's Estimator
#'
#' @param ts time-series as a vector
#' @param m number of blocks
#' @param l length of blocks
#' @param s shift factor, alternatively specified by the overlap parameter
#' @param delta sampling interval
#' @param h optional window to apply to data. Defaults to rectangular
#' @param overlap optional overlap parameter specified between 0 and 1
#'
#' @return A tibble containing Welch's estimate
#' @export
#'
#' @examples
pwelch <- function(ts, m, l, s = NULL, delta = 1, h = NULL, overlap = NULL) {

    if (is.null(s) && is.null(overlap)) {
        stop("either s or overlap must be specified")
    }

    if ((!is.null(s)) && (!is.null(overlap))) {
        warning("both s and overlap are specified. Taking s value")
    }

    if (is.null(s) && !is.null(overlap)) {
        if (!0 <= overlap && overlap <= 1) {
            stop("overlap must be between 0 and 1")
        }
        s <- round(overlap * l, 0)
    }

    n <- (m - 1) * s + l
    nfreq <- get_nfreq(l)

    if (length(ts) != n) {
        warning("n does not equal length of ts")
        stop() # Change to zero-pad and continue
    }

    if (is.null(h)) {
        h <- rep(1, l)
    }

    ff <- (1:(nfreq)) / (l * delta)
    pxx <- matrix(nrow = nfreq, ncol = m)

    h <- h / sqrt(sum(h^2))

    for (i in 1:m){
        start <- (i - 1) * s + 1
        end <- start + l - 1
        ts_tmp <- ts[start:end] * h
        pxx[, i] <- (Mod(stats::fft(ts_tmp))^2)[2:(nfreq + 1)]
    }

    pwelch <- rowMeans(pxx)

    dplyr::tibble(
        ff = ff,
        pwelch = pwelch
    )
}

# Notes for future development
# Add some tests and warning messages
# Change so it zero-pads the last one