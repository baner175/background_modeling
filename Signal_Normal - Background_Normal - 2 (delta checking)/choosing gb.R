rm(list = ls())
library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)

l <- 1; u <- 5

mean_sig <- 3.5; sd_sig <- 0.25
# mean_sig <- 3.5; sd_sig <- 0.15
mean_back <- -0.5; sd_back <- 3.25
# rate_gb <- 0.18
rate_gb <- 0.3
eta_true <- 0.05
fs_prop <- 0.05

eps <- 0.001

fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(dtrunc(x, spec = 'norm', a = l, b = u,
                mean = mean, sd = sd))
}

fb <- function(x)
{
  dtrunc(x, spec = 'norm',
         mean = mean_back, sd = sd_back,
         a = l, b = u)
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

find_d <- function(d)
{
  pl <- Fs(mean_sig-d, mean = mean_sig, sd = sd_sig)
  pu <- Fs(mean_sig+d, mean = mean_sig, sd = sd_sig)
  return(pu-pl-1+eps)
}

sol <- uniroot(find_d, lower = 0, upper = min(mean_sig - l,u - mean_sig))

r <- sol$root

d <- seq(0, min(mean_sig - l,u - mean_sig),0.001)
vals <- sapply(d, find_d)
plot(x = d, y = vals, type = 'l')
abline(v = r, h=0, col = 'red', lty = 2)

M_lower <- mean_sig - r
M_upper <- mean_sig + r

# mean_in_gb <- M_lower; sd_in_gb <- sd_sig

mean_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig

gb <- function(x)
{
  fs_val <- fs(x, mean = mean_in_gb, sd = sd_in_gb)
  # fs_val <- dtrunc(x, location = mean_in_gb, scale = sd_in_gb,
                    # spec = 'cauchy', l,u)
  gb_val <- dtrunc(x, spec = 'exp', rate = rate_gb, a = l, b = u)
  return(fs_prop*fs_val + (1-fs_prop)*gb_val)
}

Gb <- function(x)
{
  Fs_val <- Fs(x, mean = mean_in_gb, sd = sd_in_gb)
  # Fs_val <- ptrunc(x, location = mean_in_gb, scale = sd_in_gb,
          # spec = 'cauchy', l,u)
  Qb_val <- ptrunc(x, spec = 'exp', rate = rate_gb, a = l, b = u)
  return(fs_prop*Fs_val + (1-fs_prop)*Qb_val)
}

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



# kde <- kdensity::kdensity(obs)

emp_cdf <- ecdf(obs)

curve(gb, l, u, col = 'red', lwd = 2, ylim = c(0,0.5))
# curve(kde, l, u, col = 'pink',add = TRUE, lwd = 2)
curve(f, l, u, col = 'blue', add = TRUE, lwd = 2, lty = 2)
curve(fb, l, u, col = 'brown', add = TRUE, lwd = 2, lty = 2)

legend('topright', legend = c('gb', 'true f', 'true fb'),
       col = c('red', 'blue', 'brown'), lwd = 2)
abline(v = c(M_lower, M_upper))

curve(Gb, l, u, lwd = 2, col = 'red')
curve(emp_cdf,l,u, col = 'blue', add = TRUE, lwd = 2)

legend('bottomright', legend = c('Gb','empirical F'),
       col = c('red', 'blue'), lwd = 2)

abline(v = c(M_lower))


# checking if delta < 0:
integrate(function(t) S1(t, mean = mean_sig,
                         sd = sd_sig)*fb(t),
          lower = l, upper = u)$value

fb_by_gb <- function(t) fb(t)/gb(t)

xs <- seq(l,u,0.001)
sup_in_M <- max(sapply(xs[xs<M_upper & xs>M_lower], fb_by_gb))
sup_out_M <- max(sapply(xs[xs>M_upper | xs<M_lower], fb_by_gb))
curve(fb_by_gb, l, u)
abline(v = c(M_lower, M_upper), h = sup_in_M, lty = 2)

(1-eps)*sup_in_M 
eps*sup_out_M


# testing for signal:
S1_vec <- sapply(obs, function(t) S1(t, mean = mean_sig,
                                     sd = sd_sig))
(theta <- mean(S1_vec))

theta/norm_S

(se_theta <- sqrt((mean(S1_vec^2) - theta^2)/n))
(t_stat <- theta/se_theta)
t_stat>qnorm(0.95)


