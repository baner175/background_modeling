rm(list = ls())
library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)
library(latex2exp)
library(Hmisc)

options(digits = 22)

mean_sig <- 125; sd_sig <- 3
eps <- 1e-3
u <- 160; l <- 110
dat <- read.csv('Data/VBF_Cat2.csv', header = FALSE)

ni <- dat[,2]
N <- sum(ni)
k <- nrow(dat)
bins <- seq(l, u, length.out = k+1)
xi <- (bins[1:k] + bins[2:(k+1)])/2

# DEFINING SIGNAL DENSITY AND CALCULATING SIGNAL REGION

fs <- function(x) dtrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)
Fs <- function(x) ptrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)

# FIGURING OUT (mu_s -d, mu_s + d):
find_d <- function(d)
{
  pl <- Fs(mean_sig-d)
  pu <- Fs(mean_sig+d)
  return(pu-pl-1+eps)
}

sol <- uniroot(find_d, lower = 0, upper = min(mean_sig - l,u - mean_sig))

r <- sol$root

(M_lower <- mean_sig - r)
(M_upper <- mean_sig + r)

round(integrate(fs, M_lower, M_upper)$value,5) == 1-eps


# DEFINING q_b AND g_b

# generating pseudo unbinned data
set.seed(123456)
obs <- c()
for(i in 1:(length(bins))-1)
{
  obs <- c(obs, runif(ni[i], bins[i], bins[i+1]))
}

# defining gb:
bkg_loc <- 91.2
bkg_scale <- 2.49/2

qb <- function(x)
{
  dtrunc(x, a = l, b= u, spec = 'cauchy', 
         location = bkg_loc, 
         scale = bkg_scale)
}

# PROPOSAL BACKGROUND DENSITY:
mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;

gb_test <- function(x, fs_prop = 0)
{
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    spec = 'norm', a = l, b = u)
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     spec = 'norm', a = l, b = u)
  qb_val <- dtrunc(x, a = l, b= u, spec = 'cauchy', 
                   location = bkg_loc, 
                   scale = bkg_scale)
  return(fs_prop*fs_val1 + fs_prop*fs_val2 + (1-2*fs_prop)*qb_val)
}

gb <- function(x) gb_test(x, fs_prop = 0.05)
gb <- Vectorize(gb)

Gb <- function(x, ...)
{
  integrate(function(t) gb(t, ...), l, x)$value
}

Gb <- Vectorize(Gb)

# Function to calculate norm under Gb
calc_norm_gb <- function(fun, ...)
{
  integrate(function(t) gb(t)*fun(t,...)^2, l, u)$value |> sqrt()
}

# Convarting fs into S1:
norm_S <- integrate(function(t) gb(t)*(fs(t)/gb(t)-1)^2, l, u)$value|> sqrt()

S1 <- function(x)
{
  f_sig <- fs(x)
  g_b <- gb(x)
  return((f_sig/g_b -1)/norm_S)
}

# Basis functions:
T1_inner_S1 <- integrate(function(t) t*S1(t)*gb(t), l, u)$value
T1_inner_1 <- integrate(function(t) t*gb(t), l, u)$value
T1 <- function(x)
{
  return(x - T1_inner_S1*S1(x) - T1_inner_1)
}
T1_norm <- calc_norm_gb(T1)
T1_normed <- function(x) T1(x)/T1_norm

T2_inner_S1 <- integrate(function(t) (t^2)*S1(t)*gb(t), l, u)$value
T2_inner_T1 <- integrate(function(t) (t^2)*T1_normed(t)*gb(t), l, u)$value
T2_inner_1 <- integrate(function(t) (t^2)*gb(t), l, u)$value

T2 <- function(x)
{
  return(x^2 - T2_inner_S1*S1(x) -
           T2_inner_T1*T1_normed(x) - T2_inner_1)
}

T2_norm <- calc_norm_gb(T2)
T2_normed <- function(x) T2(x)/T2_norm

T3_inner_S1 <- integrate(function(t) (t^3)*S1(t)*gb(t), l, u)$value
T3_inner_T2 <- integrate(function(t) (t^3)*T2_normed(t)*gb(t), l, u)$value
T3_inner_T1 <- integrate(function(t) (t^3)*T1_normed(t)*gb(t), l, u)$value
T3_inner_1 <- integrate(function(t) (t^3)*gb(t), l, u)$value

T3 <- function(x)
{
  return(x^3 - T3_inner_S1*S1(x) - T3_inner_T2*T2_normed(x) -
           T3_inner_T1*T1_normed(x) - T3_inner_1)
}

T3_norm <- calc_norm_gb(T3)
T3_normed <- function(x) T3(x)/T3_norm


# T4_inner_S1 <- integrate(function(t) (t^4)*S1(t)*gb(t), l, u)$value
# T4_inner_T3 <- integrate(function(t) (t^4)*T3_normed(t)*gb(t), l, u)$value
# T4_inner_T2 <- integrate(function(t) (t^4)*T2_normed(t)*gb(t), l, u)$value
# T4_inner_T1 <- integrate(function(t) (t^4)*T1_normed(t)*gb(t), l, u)$value
# T4_inner_1 <- integrate(function(t) (t^4)*gb(t), l, u)$value
# 
# T4 <- function(x)
# {
#   return(x^4 - T4_inner_S1*S1(x) - T4_inner_T3*T3_normed(x) -
#            T4_inner_T2*T2_normed(x) -
#            T4_inner_T1*T1_normed(x) - T4_inner_1)
# }
# 
# T4_norm <- calc_norm_gb(T4)
# T4_normed <- function(x) T4(x)/T4_norm


# T_basis <- c(T1_normed, T2_normed, T3_normed, T4_normed)

curve(S1, l, u, lwd = 2, col = 'skyblue', ylim = c(-1,3))
curve(T1_normed, l, u, lwd = 2, col = 'red', add = TRUE)
curve(T2_normed, l, u, lwd = 2, col = 'blue', add = TRUE)
curve(T3_normed, l, u, lwd = 2, col = 'green', add = TRUE)
# curve(T4_normed, l, u, lwd = 2, col = 'orange', add = TRUE)

T_basis <- c(T1_normed, T2_normed, T3_normed)

T_basis <- sapply(T_basis, Vectorize)
tau <- sapply(T_basis, function(f) sapply(obs, f) |> mean())

# fm_null <- function(x) gb(x)*(1 + tau[1]*T1_normed(x) + tau[2]*T2_normed(x) +
#                                 tau[3]*T3_normed(x) + tau[4]*T4_normed(x))

fm_null <- function(x)
{
  T_vec <- sapply(T_basis, function(f) f(x))
  return(gb(x)*(1 + crossprod(tau, T_vec)[1,1]))
}

fm_null <- Vectorize(fm_null)

# sampmling from fm_null using acceptance rejection sampling

## sampling from Gb
u_obs <- Gb(obs)
extrap <- approxExtrap(obs, fm_null(obs)/gb(obs), xout = c(l, obs, u))
dG_GF <- approxfun(extrap$x,extrap$y, rule=2) # d(Gb(x),; Gb, Fm_hat)
extrap <- approxExtrap(u_obs, fm_null(obs)/gb(obs), xout = c(0, u_obs, 1))
du_GF <- approxfun(extrap$x,extrap$y, rule=2)
M <- max(sapply(u_obs, du_GF))
n_sim <- 1e5
fs_prop <- 0.05
xFs1 <- rtrunc(n_sim, spec = 'norm', a = l, b = u,
              mean = mean1_in_gb, sd = sd_in_gb)
xFs2 <- rtrunc(n_sim, spec = 'norm', a = l, b = u,
              mean = mean2_in_gb, sd = sd_in_gb)
vG_sim1 <- runif(n_sim)
xG_bump <- ifelse(vG_sim1<0.5, xFs1, xFs2)
xQ <- rtrunc(n_sim, a = l, b= u, spec = 'cauchy', 
             location = bkg_loc, scale = bkg_scale)
vG_sim2 <- runif(n_sim)
xG <- ifelse(vG_sim2<1-2*fs_prop, xQ, xG_bump) # sample from Gb

## sampling from fm_null
dG_GF.xG <- sapply(xG, dG_GF)
v_sim <- runif(n_sim)
xF <- xG[v_sim*M<dG_GF.xG]

hist(xF, probability = TRUE, breaks = 50)
curve(fm_null, add = TRUE, lwd = 2, col = 'brown')
curve(gb, add = TRUE, lwd = 2, col = 'red', lty = 2)
