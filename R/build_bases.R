#' Calculate biased bases (slow)
#' 
#' @inheritParams pwelch
#' @param k number of bases
#'
#' @return A l x k matrix of each of the bases
#' @export
#'
#' @examples
build_bases <- function(l, k, h = NULL) {

    is_even <- l %% 2 == 0
    nfreq <- get_nfreq(l)

    if (is_even) {
        width <- 0.5 / (k + 1)
    } else {
        width <- 0.5 / (k + 0.5)
    }

    centres <- seq(width, k * width, width)

    bases <- matrix(nrow = nfreq, ncol = k)

    tt <- 0:(l - 1)

    for (i in 1:k) {
        acf_tmp <- 2 * width * sinc(pi * tt * width) *
            cos(2 * pi * centres[i] * tt)

        bases[, i] <- bochner(acf_tmp, h = h)
    }

    colnames(bases) <- seq(1, k, 1)

    bases

}