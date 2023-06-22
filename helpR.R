
setwd("/Users/00077063/Documents/PostDoc/20220903_bandlimited_gps/dwelch")

library(devtools)
library(ggplot2)
library(tidyverse)
library(gsignal)
library(dwelch)

theme_set(theme_bw())

m <- 34
l <- 1111
s <- ceiling(l / 2)
n <- (m - 1) * s + l

phis <- c(2.7607, -3.8106, 2.6535, -0.9238)
sd <- 1e-3
sampled_ar <- stats::arima.sim(list(ar = phis), n, n.start = 1000)

# h <- gsignal::hamming(l)
h <- rep(1, l)

welch_estimate <- pwelch(sampled_ar, m, l, s, window = h)

ggplot(welch_estimate) +
    geom_line(aes(x = ff, y = welch)) +
    scale_x_continuous(
        limits = c(0, 0.5),
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        limits = c(1e-2, 1e5),
        expand = c(0, 0),
        trans = "log10"
    )

ar_spectrum(phis, noise, welch_estimate$ff)


# Calculate theoretical AR spectrum
# SAMPLE AR
# Calculate fft
# Do again but with a taper
# Plot