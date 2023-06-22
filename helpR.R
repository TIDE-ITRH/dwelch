
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
sampled_ar <- stats::arima.sim(list(ar = phis), n, n.start = 1000)

pwelch(sampled_ar, m, l, s) %>%
    ggplot() +
        geom_line(aes(x = ff, y = welch)) +
        scale_x_continuous(limits = c(0, 0.5)) +
        scale_y_continuous(limits = c(0, 25000))











((stats::fft(ts_tmp) %>% abs())^2)[1:floor(l/2)]
((stats::fft(ts_tmp) %>% abs())^2)[2:(l/2)]
((stats::fft(ts_tmp) %>% abs())^2)[l:((l/2)+2)]


((stats::fft(ts_tmp) %>% abs())^2)[floor(l / 2)]

tibble(y = ((stats::fft(ts_tmp) %>% abs())^2), x = (0:(l-1))/l) %>%
    ggplot() +
        geom_line(aes(x = x, y = y))

0:floor(l / 2)

for (i in 1:m){
    start <- (i - 1) * n_star + 1
    end <- start + l - 1
    ts_tmp <- sampled_ar$ar[start:end]
    periodogram <- ((stats::fft(ts_tmp) %>% abs())^2)[1:floor(l / 2)]
}

# First even
# Then odd




# sampled_ar <- tibble(
#     time = seq(1, n),
#     ar = stats::arima.sim(list(ar = phis), n, n.start = 1000)
# )

# ggplot(sampled_ar) +
#     geom_line(aes(x = time, y = ar))

# NOTE Okay, write these out into a table. average, get into a pwelch function,
# save and write it out. Move on.


tibble(
    ff = (1:floor(l / 2)) / l,
    I = periodogram
) %>%
    ggplot() +
        geom_line(aes(x = ff, y = I)) +
        scale_y_log10()

# }

window <- rep(1 / sqrt(l), l)

pwelch(sampled_ar$ar, window, 0.5)




# Calculate theoretical AR spectrum
# SAMPLE AR
# Calculate fft
# Do again but with a taper
# Plot