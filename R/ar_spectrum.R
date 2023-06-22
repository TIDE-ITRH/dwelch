
ar_spectrum <- function(phis, sd, ff) {
    d <- 1
    d_conj <- 1

    for (j in seq_along(phis)) {
        d <- d - phis[j] * exp(-2i * j * pi * ff)
        d_conj <- d_conj - phis[j] * exp(2i * j * pi * ff)
    }

    sd^2 / Re(d * d_conj)

}