library(truncdist)
library(VGAM)

l <- 1; u <- 5
rate_back <- 3
rate_sig <- 1 #or 
shape_sig <- 0.01
rate_gb <- 2
eta <- 0.25

back <- function(t) dtrunc(t, spec = 'exp', a = l, b = u, rate = rate_back)
sig <- function(t) dtrunc(t, spec = 'exp', a = l, b = u, rate = rate_sig)
gb <- function(t) dtrunc(t, spec = 'exp', a = l, b = u, rate = rate_gb)
mix <- function(t) (1-eta)*back(t) + eta*sig(t)


curve(back, l, u, col = 'brown', lwd = 2)
curve(sig, l, u, col = 'blue', lwd = 2, add = TRUE)
curve(gb, l, u, col = 'red', lwd = 2, add = TRUE)
curve(mix, l, u, col = 'black', lwd = 2, add = TRUE)

# delta:
integrate(function(t) (sig(t)/gb(t) - 1)*back(t), l, u)



back <- function(t) dtrunc(t, spec = 'exp', a = l, b = u, rate = rate_back)
sig <- function(t) dtrunc(t, spec = 'pareto', a = l, b = u, shape = shape_sig, scale = l)
gb <- function(t) dtrunc(t, spec = 'exp', a = l, b = u, rate = rate_gb)
mix <- function(t) (1-eta)*back(t) + eta*sig(t)

curve(back, l, u, col = 'brown', lwd = 2)
curve(sig, l, u, col = 'blue', lwd = 2, add = TRUE)
curve(gb, l, u, col = 'red', lwd = 2, add = TRUE)
curve(mix, l, u, col = 'black', lwd = 2, add = TRUE)

# delta:
integrate(function(t) (sig(t)/gb(t) - 1)*back(t), l, u)
