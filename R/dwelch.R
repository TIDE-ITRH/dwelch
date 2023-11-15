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
dwelch <- function(
    ts, m, l, s, k = NULL, delta = 1, h = NULL, model = c("vanilla", "nnls"), 
    centres = NULL, widths = NULL, lowers = NULL, uppers = NULL
) {
    
    model <- match.arg(model)

    pw_tbl <- dwelch::pwelch(ts, m, l, s, delta, h)

    L <-  diag(1 / (pw_tbl$pwelch))
    b <- L %*% matrix(pw_tbl$pwelch)

    if (length(k) == 1){

        A <- L %*% dwelch::build_bases(l, k, h, delta = delta)
        centres <- dwelch::get_centres(l, k, delta = delta)$centres

    } else if ((!is_null(centres) && !is_null(widths))) {
        
        A <- L %*% dwelch::build_bases(
            l, k, h, delta = delta, centres = centres, widths = widths
        )

    } else if ((!is_null(lowers) && !is_null(uppers))) {
        
        A <- L %*% dwelch::build_bases(
            l, k, h, delta = delta, lowers = lowers, uppers = uppers
        )
        widths <- uppers-lowers
        centres <- lowers + widths/2

    } else {

        stop("I need a working basis description")

    }

    if (model == "vanilla") {
        dw <- solve(t(A) %*% A) %*% t(A) %*% b
    } else if (model == "nnls") {
        dw <- nnls::nnls(A, b)$x
    }

    dplyr::tibble(
        ff = centres,
        dwelch = as.numeric(dw)
    )

}