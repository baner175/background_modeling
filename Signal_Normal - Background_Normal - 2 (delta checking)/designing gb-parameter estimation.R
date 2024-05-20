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

M_lower <- mean_sig - r
M_upper <- mean_sig + r

# mean_in_gb <- M_lower; sd_in_gb <- sd_sig

mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;


x_bar <- mean(obs)
x_sq_bar <- mean(obs^2)
fs_star <- function(x) {
  fs_val_1 <- fs(x, mean = mean1_in_gb, sd = sd_in_gb)
  fs_val_2 <- fs(x, mean = mean2_in_gb, sd = sd_in_gb)
  return((fs_val_1+fs_val_2)/2)
}
mu1_fs_star <- integrate(function(t) fs_star(t)*t, l, u)$value
mu2_fs_star <- integrate(function(t) fs_star(t)*(t^2), l, u)$value

# solve_beta <- function(beta){
#   mu1_qb <- 1/beta + 
#     (l*exp(-beta*l) - u*exp(-beta*u))/(exp(-beta*l) - exp(-beta*u))
#   mu2_qb <- ((l^2)*exp(-beta*l) - (u^2)*exp(-beta*u))/(exp(-beta*l) - exp(-beta*u)) + 
#     (2/beta)*mu1_qb
#   
#   val1 <- (x_bar - mu1_fs_star)/(mu1_qb - mu1_fs_star)
#   val2 <- (x_sq_bar - mu2_fs_star)/(mu2_qb - mu2_fs_star)
#   return(val1 - val2)
# }

solve_beta <- function(beta){
  mu1_qb <- ((1-beta)/(2-beta))*(u^(2-beta) - l^(2-beta))/(u^(1-beta) - l^(1-beta))
  mu2_qb <- ((1-beta)/(3-beta))*(u^(3-beta) - l^(3-beta))/(u^(1-beta) - l^(1-beta))

  val1 <- (x_bar - mu1_fs_star)/(mu1_qb - mu1_fs_star)
  val2 <- (x_sq_bar - mu2_fs_star)/(mu2_qb - mu2_fs_star)
  return(val1 - val2)
}

curve(solve_beta, 0.01, 15)
abline(h = 0, col = 'red')

# gb <- function(x)
# {
#   fs_val_1 <- fs(x, mean = mean1_in_gb, sd = sd_in_gb)
#   fs_val_2 <- fs(x, mean = mean2_in_gb, sd = sd_in_gb)
#   gb_val <- dtrunc(x, spec = 'exp', rate = rate_gb, a = l, b = u)
#   return((fs_prop)*fs_val_1 + (fs_prop)*fs_val_2 + (1-2*fs_prop)*gb_val)
# }
# 
# 
# 
# calc_norm_gb <- function(fun, ...)
# {
#   integrate(function(t) gb(t)*fun(t,...)^2, l, u)$value |> sqrt()
# }
# 
# norm_S <- calc_norm_gb(function(t) fs(t, mean= mean_sig, sd = sd_sig)/gb(t) - 1)
# 
# S1 <- function(x, ...)
# {
#   f_sig <- fs(x, ...)
#   g_b <- gb(x)
#   return((f_sig/g_b -1)/norm_S)
# }
# 
# 
# 
# kde <- kdensity::kdensity(obs)
# 
# curve(gb, l, u, col = 'red', lwd = 2, ylim = c(0,0.5))
# curve(kde, l, u, col = 'pink',add = TRUE, lwd = 2)
# curve(f, l, u, col = 'blue', add = TRUE, lwd = 2, lty = 2)
# curve(fb, l, u, col = 'brown', add = TRUE, lwd = 2, lty = 2)
# 
# legend('topright', legend = c('gb', 'true f', 'true fb', 'kde'),
#        col = c('red', 'blue', 'brown', 'pink'), lwd = 2)
# abline(v = c(M_lower, M_upper))
# 
# abline(v = c(M_lower))
# 
# 
# # checking if delta < 0:
# integrate(function(t) S1(t, mean = mean_sig,
#                          sd = sd_sig)*fb(t),
#           lower = l, upper = u)$value
# 
# 
