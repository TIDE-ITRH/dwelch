#' Calculate the Debiased Welch Estimator (slow)
#'
#' @inheritParams pwelch
#' @inheritParams build_bases
#'
#' @return A tibble containing the debiased Welch estimate
#' @export
#'
#' @examples
dwelch <- function(ts, m, l, s, k, delta = 1, h = NULL) {

    pw_tbl <- dwelch::pwelch(ts, m, l, s, delta, h)

    W <-  diag(1 / (pw_tbl$pwelch^2)) # nolint: object_name_linter.
    pw <- matrix(pw_tbl$pwelch)

    B <- build_bases(l, k, h) # nolint: object_name_linter.

    dwelch <- solve(t(B) %*% W %*% B) %*% t(B) %*% W %*% pw

    dplyr::tibble(
        ff = get_centres(l, k, delta = delta)$centres,
        dwelch = as.numeric(dwelch)
    )

}

# I think we could probably figure out l from figuring out
# odd/even and backing out. Probably not worth it for now.
# Add in functionality through an option on how to compute
# dwelch should run it's own pwelch??