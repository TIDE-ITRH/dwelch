
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dwelch

<!-- badges: start -->
<!-- badges: end -->

<tt> dwelch </tt> provides code to calculate the debiased Welch
estimator developed in LINK.

## Installation

You can install the development version of <tt> dwelch </tt> from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("astfalckl/dwelch")
```

## AR(4) Example

We demonstrate the basic functionality of <tt> dwelch </tt> with the
classic AR(4) problem. First, import some basic packages.

``` r
library(dwelch)
library(ggplot2)
library(tidyverse)
library(gsignal)

theme_set(theme_bw())
```

Next, we set our parameters and generate an AR process.

``` r
m <- 32 # Number of blocks
l <- 1024 # Block length
s <- ceiling(l / 2) # Shift factor
n <- (m - 1) * s + l # Total data length

phis <- c(2.7607, -3.8106, 2.6535, -0.9238)
sd <- 1

sampled_ar <- stats::arima.sim(list(ar = phis), n, n.start = 1000, sd = sd)
```

Define our data taper, <tt> h </tt>, and calculate Welch’s estimate of
the AR process. Here we have assumed a rectangular taper, but
alternative tapers can be used. See <tt> gsignal </tt> for a fairly
comprehensive list; for example, <tt> gsignal::hamming(l) </tt>.

``` r
h <- rep(1, l)

welch_estimate <- dwelch::pwelch(sampled_ar, m, l, s, window = h)

welch_estimate %>%
    mutate(ar = ar_spectrum(welch_estimate$ff, phis, sd)) %>%
    ggplot() +
        geom_line(aes(x = ff, y = ar), colour = "#e17c1d") +
        geom_line(aes(x = ff, y = welch)) +
        scale_x_continuous(
            "frequency [Hz]",
            limits = c(0, 0.5),
            expand = c(0, 0)
        ) +
        scale_y_continuous(
            "log spectral density",
            limits = c(1e-3, 1e5),
            expand = c(0, 0),
            trans = "log10"
        )
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />
