#' Whittle Likelihood
#'
#' @param psd Power spectral density
#' @param periodogram Periodogram
#' @param return_log Return the log likelihood
#'
#' @return The value of the Whittle likelihood
#' @export
#'
#' @examples
whittle <- function(psd, periodogram, return_log = TRUE) {
    wll <- sum(log(psd) + periodogram / psd)

    if (return_log) {
        wll
    } else {
        exp(wll)
    }
}