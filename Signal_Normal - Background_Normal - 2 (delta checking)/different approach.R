rm(list = ls())
library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)

l <- 1; u <- 5

mean_sig <- 3.5; sd_sig <- 0.25
mean_back <- -0.5; sd_back <- 3.25
# rate_gb <- 0.25
# rate_gb <- 0.15
rate_gb <- 0.2
eta_true <- 0.05
fs_prop <- 0.02
# fs_prop <- 0

fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(dtrunc(x, spec = 'norm', a = l, b = u,
                mean = mean, sd = sd))
}

Fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(ptrunc(x, spec = 'norm', a = l, b = u,
                mean = mean, sd = sd))
}

f <- function(x)
{
  return(eta_true*dtrunc(x, spec = 'norm', a = l, b = u,
                         mean = mean_sig, sd = sd_sig) + (1-eta_true)*dtrunc(x, spec = 'norm', a = l, b = u,
                                                                             mean = mean_back,
                                                                             sd = sd_back))
}


set.seed(12345)
n <- 5e3

back <- rtrunc(n, spec = 'norm',
               mean = mean_back, sd = sd_back,
               a = l, b = u)
sig <- rtrunc(n, spec = 'norm',
              mean = mean_sig, sd = sd_sig,
              a = l, b = u)
u_mask <- rbinom(n, size = 1, prob = eta_true)

obs <- ifelse(u_mask, sig, back)

gb <- function(x)
{
  # fs_val <- fs(x, mean = mean_sig, sd = sd_sig)
  fs_val <- fs(x, mean = mean_sig, sd = 0.5)
  gb_val <- dtrunc(x, spec = 'exp', rate = rate_gb, a = l, b = u)
  return(fs_prop*fs_val + (1-fs_prop)*gb_val)
}

fb <- function(x)
{
  dtrunc(x, spec = 'norm',
         mean = mean_back, sd = sd_back,
         a = l, b = u)
}

fb_by_gb <- function(x) fb(x)/gb(x)


calc_norm_gb <- function(fun, ...)
{
  integrate(function(t) gb(t)*fun(t,...)^2, l, u)$value |> sqrt()
}

norm_S <- calc_norm_gb(function(t) fs(t, mean= mean_sig, sd = sd_sig)/gb(t) - 1)

S1 <- function(x, ...)
{
  f_sig <- fs(x, ...)
  g_b <- gb(x)
  return((f_sig/g_b -1)/norm_S)
}

# checking if delta < 0:
integrate(function(t) S1(t, mean = mean_sig,
                         sd = sd_sig)*dtrunc(t, mean = mean_back,
                                             sd = sd_back, spec = 'norm',
                                             a = l, b = u),
          lower = l, upper = u)$value

find_intersection <- function(t)
{
  return(gb(t) - fs(t))
}

s1 <- uniroot(find_intersection, interval = c(l, mean_sig))$root
s2 <- uniroot(find_intersection, interval = c(mean_sig, u))$root



curve(gb, col = 'red',l,u, lwd = 2, ylim = c(0, 1))
curve(fs, col = 'blue',l,u, lwd = 2, add = TRUE)
curve(fb, col = 'brown',l,u, lwd = 2, add = TRUE)
curve(f, col = 'skyblue',l,u, lwd = 2, add = TRUE)
abline(v = c(s1, s2), col = 'black', lwd = 2, lty = 2)
legend('topleft', legend = c('gb', 'fs', 'fb','f','v'),
       col = c('red', 'blue', 'brown', 'skyblue'),
       lwd = 2, bty = 'n')

curve(fb_by_gb, s1, s2)

curve(fb_by_gb, l, u)
abline(v = c(s1, s2), col = 'black', lwd = 2, lty = 2)

xs <- seq(l,u,0.01)
fb_by_gb_vec <- sapply(xs, fb_by_gb)

min(fb_by_gb_vec)
max(fb_by_gb_vec)


inf_val <- min(sapply(xs[xs<s1 | xs> s2], fb_by_gb))
sup_val <- max(sapply(xs[xs>s1 & xs<s2], fb_by_gb))
