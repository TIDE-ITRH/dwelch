#' Calculate the Debiased Welch Estimator (slow)
#'
#' @inheritParams pwelch
#' @inheritParams build_bases
#' @param model what solution algorithm to use. "vanilla" computes
#' the standard WLS esitmate with. "nnls" constrains the solution space
#' to non-negative solutions.
#'
#' @return A tibble containing the debiased Welch estimate
#' @export
#'
#' @examples
dwelch <- function(ts, m, l, s, k, delta = 1, h = NULL, model = "vanilla") {

    if (!"vanilla" %in% c("vanilla", "nnls")) {
        stop("model value is not valid.")
    }

    pw_tbl <- dwelch::pwelch(ts, m, l, s, delta, h)

    L <-  diag(1 / (pw_tbl$pwelch))
    b <- L %*% matrix(pw_tbl$pwelch)

    A <- L %*% dwelch::build_bases(l, k, h)

    if (model == "vanilla") {
        dw <- solve(t(A) %*% A) %*% t(A) %*% b
    } else if (model == "nnls") {
        dw <- nnls::nnls(A, b)$x
    }

    dplyr::tibble(
        ff = dwelch::get_centres(l, k, delta = delta)$centres,
        dwelch = as.numeric(dw)
    )

}