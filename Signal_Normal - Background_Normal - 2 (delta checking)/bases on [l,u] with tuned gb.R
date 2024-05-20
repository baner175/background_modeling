rm(list = ls())
library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)

l <- 1; u <- 5

mean_sig <- 3.5; sd_sig <- 0.25
mean_back <- -0.5; sd_back <- 3.25
rate_gb <- 0.3
eta_true <- 0.05
fs_prop <- 0.05
eps <- 1e-3

# SIGNAL DENSITY:
fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(dtrunc(x, spec = 'norm', a = l, b = u,
                mean = mean, sd = sd))
}

# SIGNAL CDF:
Fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(ptrunc(x, spec = 'norm', a = l, b = u,
                mean = mean, sd = sd))
}

# Figuring out (mu_s -d, mu_s + d) that covers an area of 1-eps:
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
# plot(x = d, y = vals, type = 'l')
# abline(v = r, h=0, col = 'red', lty = 2)

M_lower <- mean_sig - r
M_upper <- mean_sig + r

# PROPOSAL BACKGROUND DENSITY:
mean_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
gb <- function(x)
{
  fs_val <- fs(x, mean = mean_in_gb, sd = sd_in_gb)
  qb_val <- dtrunc(x, spec = 'exp', rate = rate_gb, a = l, b = u)
  return(fs_prop*fs_val + (1-fs_prop)*qb_val)
}


# DATA GENERATION
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


# Area under gb
integrate(gb,l,u)

# Function to calculate norm under Gb
calc_norm_gb <- function(fun, ...)
{
  integrate(function(t) gb(t)*fun(t,...)^2, l, u)$value |> sqrt()
}

# Convarting fs into S1:
norm_S <- calc_norm_gb(function(t) fs(t, mean= mean_sig, sd = sd_sig)/gb(t) - 1)

S1 <- function(x, ...)
{
  f_sig <- fs(x, ...)
  g_b <- gb(x)
  return((f_sig/g_b -1)/norm_S)
}

# Basis functions:
T1_inner_S1 <- integrate(function(t) t*S1(t,
                                          mean = mean_sig, sd = sd_sig)*gb(t), l, u)$value
T1_inner_1 <- integrate(function(t) t*gb(t), l, u)$value
T1 <- function(x)
{
  return(x - T1_inner_S1*S1(x, mean = mean_sig, sd = sd_sig) - T1_inner_1)
}

T1_norm <- calc_norm_gb(T1)
T1_normed <- function(x) T1(x)/T1_norm

T2_inner_S1 <- integrate(function(t) (t^2)*S1(t,
                                              mean = mean_sig, sd = sd_sig)*gb(t), l, u)$value
T2_inner_T1 <- integrate(function(t) (t^2)*T1_normed(t)*gb(t), l, u)$value
T2_inner_1 <- integrate(function(t) (t^2)*gb(t), l, u)$value

T2 <- function(x)
{
  return(x^2 - T2_inner_S1*S1(x, mean = mean_sig, sd = sd_sig) -
           T2_inner_T1*T1_normed(x) - T2_inner_1)
}

T2_norm <- calc_norm_gb(T2)
T2_normed <- function(x) T2(x)/T2_norm


T3_inner_S1 <- integrate(function(t) (t^3)*S1(t,
                                              mean = mean_sig, sd = sd_sig)*gb(t), l, u)$value
T3_inner_T1 <- integrate(function(t) (t^3)*T1_normed(t)*gb(t), l, u)$value
T3_inner_T2 <- integrate(function(t) (t^3)*T2_normed(t)*gb(t), l, u)$value
T3_inner_1 <- integrate(function(t) (t^3)*gb(t), l, u)$value

T3 <- function(x)
{
  return(x^3 - T3_inner_S1*S1(x, mean = mean_sig, sd = sd_sig) -
           T3_inner_T1*T1_normed(x) - T3_inner_T2*T2_normed(x) - 
           T3_inner_1)
}

T3_norm <- calc_norm_gb(T3)
T3_normed <- function(x) T3(x)/T3_norm


T4_inner_S1 <- integrate(function(t) (t^4)*S1(t,
                                              mean = mean_sig, sd = sd_sig)*gb(t), l, u)$value
T4_inner_T1 <- integrate(function(t) (t^4)*T1_normed(t)*gb(t), l, u)$value
T4_inner_T2 <- integrate(function(t) (t^4)*T2_normed(t)*gb(t), l, u)$value
T4_inner_T3 <- integrate(function(t) (t^4)*T3_normed(t)*gb(t), l, u)$value
T4_inner_1 <- integrate(function(t) (t^4)*gb(t), l, u)$value

T4 <- function(x)
{
  return(x^4 - T4_inner_S1*S1(x, mean = mean_sig, sd = sd_sig) -
           T4_inner_T1*T1_normed(x) - T4_inner_T2*T2_normed(x) - 
           T4_inner_T3*T3_normed(x) - T4_inner_1)
}

T4_norm <- calc_norm_gb(T4)
T4_normed <- function(x) T4(x)/T4_norm


T5_inner_S1 <- integrate(function(t) (t^5)*S1(t,
                                              mean = mean_sig,
                                              sd = sd_sig)*gb(t), l, u)$value
T5_inner_T1 <- integrate(function(t) (t^5)*T1_normed(t)*gb(t), l, u)$value
T5_inner_T2 <- integrate(function(t) (t^5)*T2_normed(t)*gb(t), l, u)$value
T5_inner_T3 <- integrate(function(t) (t^5)*T3_normed(t)*gb(t), l, u)$value
T5_inner_T4 <- integrate(function(t) (t^5)*T4_normed(t)*gb(t), l, u)$value
T5_inner_1 <- integrate(function(t) (t^5)*gb(t), l, u)$value

T5 <- function(x)
{
  return(x^5 - T5_inner_S1*S1(x, mean = mean_sig, sd = sd_sig) -
           T5_inner_T1*T1_normed(x) - T5_inner_T2*T2_normed(x) - 
           T5_inner_T3*T3_normed(x) - T5_inner_T4*T4_normed(x) - 
           T5_inner_1)
}

T5_norm <- calc_norm_gb(T5)
T5_normed <- function(x) T5(x)/T5_norm

# The basis:
T_basis <- c(T1_normed, T2_normed, T3_normed, T4_normed, T5_normed)

full_basis <- c(S1, T_basis)

# model:
mod <- function(x, theta, tau)
{
  tau_len <- length(tau)
  Tvec <- c()
  for(j in 1:tau_len) {
    Tvec <- c(Tvec,
              T_basis[[j]](x))
  }
  mix <- 1 + theta*S1(x, mean = mean_sig, sd = sd_sig) + as.numeric(crossprod(tau,Tvec))
  mix <- mix*gb(x)
  eta <- theta/norm_S
  beta <- tau/(1-eta)
  back <- 1 + as.numeric(crossprod(beta,Tvec))
  back <- back*gb(x)
  return(list(mix = mix, back = back))
}

#log-likelihood:
log_lik <- function(data, theta, tau)
{
  y <- sapply(data, function(t) mod(t, theta = theta, tau = tau)$mix)
  return(sum(log(y)))
}
