#' Calculate locations of dwelch centres
#'
#' @inheritParams pwelch
#' @inheritParams build_bases
#'
#' @return A list containing the centres and the bin width
#' @export
#'
#' @examples
get_centres <- function(l, k, delta = 1) {
    is_even <- l %% 2 == 0

    if (is_even) {
        width <- 0.5 / (k + 1)
    } else {
        width <- 0.5 / (k + 0.5)
    }

    centres <- seq(width, k * width, width) / delta

    return(list(centres = centres, width = width))

}