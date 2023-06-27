#' Calculate the Debiased Welch Estimator (slow)
#'
#' @inheritParams pwelch
#' @inheritParams build_bases
#' @param pwelch output tibble from function pwelch
#'
#' @return A tibble containing the debiased Welch estimate
#' @export
#'
#' @examples
dwelch <- function(pwelch, l, k, delta = 1, h = NULL) {

    W <-  diag(1 / (pwelch$pwelch^2))
    pw <- matrix(pwelch$pwelch)

    B <- build_bases(l, k, h)

    dwelch <- solve(t(B) %*% W %*% B) %*% t(B) %*% W %*% pw

    dplyr::tibble(
        ff = get_centres(l, k, delta = delta)$centres,
        dwelch = as.numeric(dwelch)
    )

}

# I think we could probably figure out l from figuring out
# odd/even and backing out. Probably not worth it for now.
# Add in functionality through an option on how to compute