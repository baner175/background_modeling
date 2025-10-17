rm(list = ls())

library(VGAM)
library(truncdist)
library(latex2exp)
library(knitr)
library(kableExtra)
library(ggplot2)

real_l <- 1; real_u <- 35
l <- log(real_l); u <- log(real_u)

mean_sig <- 3.5; sd_sig <- sqrt(0.01*3.5^2)
eps <- 1e-3
mu_in_qb <- -1; sigma_factor_in_qb <- 2

# SIGNAL DENSITY:
fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(dtrunc(exp(x), spec = 'norm', a = real_l, b = real_u,
                mean = mean, sd = sd)*exp(x))
}

# SIGNAL CDF:
Fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(ptrunc(exp(x), spec = 'norm', a = real_l, b = real_u,
                mean = mean, sd = sd))
}

phys_data <- read.table('Data_ex1.txt', header = TRUE)
bkg_data <- read.table('cal.txt', header = TRUE)
y <- log(bkg_data$x)
x <- log(phys_data$x)

N <- length(x); M <- length(y)
k <- 1e2
bin_ends <- seq(l, u, length.out = k+1)
xi <- (bin_ends[1:k] + bin_ends[2:(k+1)])/2

ni <- sapply(1:k, function(i){
  sum((x > bin_ends[i])&(x <= bin_ends[i+1]))
})
mi <- sapply(1:k, function(i){
  sum((y > bin_ends[i])&(y <= bin_ends[i+1]))
})

gb <- function(t) dunif(t, min = l, max = u)

norm_S <- integrate(function(t){
  fs <- dtrunc(exp(t), spec = 'norm', a = real_l, b = real_u,
               mean = mean_sig, sd = sd_sig)*exp(t)
  gb <- gb(t)
  return(((fs/gb - 1)^2)*gb)
}, l, u)$value |> sqrt()

S2_vec <- sapply(xi, function(t){
  fs <- dtrunc(exp(t), spec = 'norm', a = real_l, b = real_u,
               mean = mean_sig, sd = sd_sig)*exp(t)
  gb <- gb(t)
  S_val <- (fs/gb - 1)
  return(S_val/(norm_S^2))
})

theta0_hat <- sum(S2_vec*ni)/N; delta0_hat <- sum(S2_vec*mi)/M
eta_hat <- (theta0_hat-delta0_hat)/(1-delta0_hat)

test_num <- sqrt(M*N)*eta_hat


d_theta_T <- 1/(1-delta0_hat)
d_delta_T <- (theta0_hat - 1)/((1-delta0_hat)^2)

var_S2_F_hat <- sum((S2_vec^2)*ni)/N - (theta0_hat^2)
var_S2_Fb_hat <- sum((S2_vec^2)*mi)/M - (delta0_hat^2)

denom1 <- M*var_S2_F_hat*(d_theta_T^2)
denom2 <- N*var_S2_Fb_hat*(d_delta_T^2)

test_denom <- sqrt(denom1 + denom2)

test_stat <- test_num/test_denom
p_val <- pnorm(test_stat, lower.tail = FALSE)

std_err <- test_denom/(sqrt(M*N))
ci_95 <- eta_hat + c(-1,1)*qnorm(0.975)*std_err

plot(y = mi, x = xi,
     pch = 16, col = 'grey', xlab = 'log(x)',
     main = TeX('Estimated proposal bkg $q_b(x; \\hat{beta})$ on bkg data'))

curve(M*gb(x)*(u-l)/(k+1), col = alpha('blue', 0.6), add = TRUE, lwd = 2.2,
      lty = 1)

plot(y = ni, x = xi,
     pch = 16, col = 'grey', xlab = 'log(x)',
     main = TeX('Estimated proposal bkg $q_b(x; \\hat{beta})$ on physics data'))

curve(N*gb(x)*(u-l)/(k+1), col = alpha('blue', 0.6), add = TRUE, lwd = 2.2,
      lty = 1)
text(x = 2.5, y = 60, 
     TeX(sprintf(
       paste0('$\\hat{\\eta} = %.4f$, $p$-value: ', as.character(round(p_val, 6)), '; $\\sigma$-signif.: %.1f'),
       eta_hat, qnorm(p_val, lower.tail = FALSE))
     ))
text(x = 2.5, y = 50, 
     TeX(sprintf(
       paste0('95 \\%% CI for $\\hat{\\eta}$: [ %.4f, %.4f]'),
       ci_95[1], ci_95[2])
     ))
