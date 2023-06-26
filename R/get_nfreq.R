#' Calculate number of frequencies exluding zero and Nyquist
#'
#' @param l length of blocks
#'
#' @return An integer of number of frequencies
#' @export
#'
#' @examples
get_nfreq <- function(l) {

    is_even <- l %% 2 == 0

    if (is_even) {
        floor(l / 2) - 1
    } else {
        floor(l / 2)
    }

}