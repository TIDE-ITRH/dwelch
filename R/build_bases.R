#' Calculate biased bases (slow)
#'
#' @inheritParams pwelch
#' @param k number of bases
#'
#' @return A l x k matrix of each of the bases
#' @export
#'
#' @examples
build_bases <- function(
    l, k = NULL, h = NULL, delta = 1, centres = NULL, 
    widths = NULL, lowers = NULL, uppers = NULL
) {

    nfreq <- dwelch::get_nfreq(l)

    if (length(k) == 1){
        
        bases_type <- "even"
        centres <- dwelch::get_centres(l, k, delta)$centres
        widths <- rep(dwelch::get_centres(l, k, delta)$width, k)

    } else if ((!is_null(centres) && !is_null(widths))) {
        
        bases_type <- "centred"

    } else if ((!is_null(lowers) && !is_null(uppers))) {
        
        bases_type <- "bounded"
        widths <- uppers-lowers
        centres <- lowers + widths/2

    }

    k <- length(centres)
    bases <- matrix(nrow = nfreq, ncol = k)

    tt <- 0:(l - 1)

    for (i in 1:k) {
        acf_tmp <- 2 * widths[i] * dwelch::sinc(pi * tt * widths[i]) *
            cos(2 * pi * centres[i] * tt)

        bases[, i] <- dwelch::bochner(acf_tmp, h = h)
    }

    colnames(bases) <- seq(1, k, 1)

    bases

}