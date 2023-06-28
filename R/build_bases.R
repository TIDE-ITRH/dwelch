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

    nfreq <- dwelch::get_nfreq(l)

    centres <- dwelch::get_centres(l, k)$centres
    width <- dwelch::get_centres(l, k)$width

    bases <- matrix(nrow = nfreq, ncol = k)

    tt <- 0:(l - 1)

    for (i in 1:k) {
        acf_tmp <- 2 * width * dwelch::sinc(pi * tt * width) *
            cos(2 * pi * centres[i] * tt)

        bases[, i] <- dwelch::bochner(acf_tmp, h = h)
    }

    colnames(bases) <- seq(1, k, 1)

    bases

}