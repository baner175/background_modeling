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
n <- length(x)
m <- length(y)

qb_y_model <- function(beta){
  qb_i <- dtrunc(y, spec = 'norm',
                 mean = mu_in_qb,
                 sd = sqrt(sigma_factor_in_qb*beta),
                 a = l, b = u)
  return(-sum(log(qb_i)))
}

beta_hat <- nlminb(start = 0.01,
                   objective = qb_y_model,
                   upper = Inf, lower = 0)$par
h <- function(t) exp(-((t-mu_in_qb)^2)/(2*sigma_factor_in_qb*beta_hat))
qb <- function(t) dtrunc(t, spec = 'norm',
                         mean = mu_in_qb,
                         sd = sqrt(sigma_factor_in_qb*beta_hat),
                         a = l, b = u)
d_log_h <- function(t) ((t-mu_in_qb)^2)/(2*sigma_factor_in_qb*(beta_hat^2))
d2_log_h <- function(t) (-(t-mu_in_qb)^2)/(sigma_factor_in_qb*(beta_hat^3))
E_qb_d_log_h <- integrate(function(t) d_log_h(t)*qb(t), l, u)$value

d2_log_qb_int_1 <- integrate(function(y){
  d2_log_h <- d2_log_h(y)
  d_log_h <- d_log_h(y)
  d2_h_by_h <- d2_log_h + d_log_h^2
  qb <- qb(y)
  return(d2_h_by_h*qb)
}, l, u)$value
d2_log_qb_int_2 <- integrate(function(y){
  d_log_h <- d_log_h(y)
  qb <- qb(y)
  return(d_log_h*qb)
}, l, u)$value

d_log_qb <- function(t) d_log_h(t) - E_qb_d_log_h

d2_log_qb <- function(t){
  d2_log_h <- d2_log_h(t)
  int_val <- d2_log_qb_int_1 - (d2_log_qb_int_2^2)
  return(d2_log_h - int_val)
}

norm_S <- integrate(function(t){
  fs <- dtrunc(exp(t), spec = 'norm', a = real_l, b = real_u,
               mean = mean_sig, sd = sd_sig)*exp(t)
  qb <- qb(t)
  return(((fs/qb - 1)^2)*qb)
}, l, u)$value |> sqrt()


S2_phys_vec <- sapply(x, function(t){
  fs <- dtrunc(exp(t), spec = 'norm', a = real_l, b = real_u,
               mean = mean_sig, sd = sd_sig)*exp(t)
  qb <- qb(t)
  S_val <- (fs/qb - 1)
  return(S_val/(norm_S^2))
})

S2_bkg_vec <- sapply(y, function(t){
  fs <- dtrunc(exp(t), spec = 'norm', a = real_l, b = real_u,
               mean = mean_sig, sd = sd_sig)*exp(t)
  qb <- qb(t)
  S_val <- (fs/qb - 1)
  return(S_val/(norm_S^2))
})

theta0_hat <- mean(S2_phys_vec); delta0_hat <- mean(S2_bkg_vec)
eta_hat <- (theta0_hat-delta0_hat)/(1-delta0_hat)

test_num <- sqrt(m*n)*eta_hat

d_normS_sq <- -integrate(function(t){
  fs <- dtrunc(exp(t), spec = 'norm', a = real_l, b = real_u,
               mean = mean_sig, sd = sd_sig)*exp(t)
  qb <- qb(t)
  d_log_qb <- d_log_qb(t)
  return((fs^2)*d_log_qb/qb)
}, l, u)$value

d_S2 <- function(t){
  fs <- dtrunc(exp(t), spec = 'norm', a = real_l, b = real_u,
               mean = mean_sig, sd = sd_sig)*exp(t)
  qb <- qb(t)
  d_log_qb <- d_log_qb(t)
  
  return(-((norm_S^2)*(fs/qb)*qb*d_log_qb + (fs/qb-1)*d_normS_sq)/(norm_S^4))
}

d_log_qb_yi <- sapply(y, d_log_qb)

d2_log_qb_yi <- sapply(y, d2_log_qb)

d_theta0 <- sapply(x, d_S2) |> mean()
d_delta0 <- sapply(y, d_S2) |> mean()
d_theta_T <- 1/(1-delta0_hat)
d_delta_T <- (theta0_hat - 1)/((1-delta0_hat)^2)
cov_term <- mean(d_log_qb_yi*S2_bkg_vec)
V_hat <- mean((d_log_qb_yi)^2)
J_hat <- -mean(d2_log_qb_yi)

var_S2_F_hat <- mean(S2_phys_vec^2) - (theta0_hat^2)
var_S2_Fb_hat <- mean(S2_bkg_vec^2) - (delta0_hat^2)

denom1 <- m*(d_theta_T^2)*var_S2_F_hat
denom2 <- n*(d_delta_T^2)*var_S2_Fb_hat
denom3 <- n*(V_hat/(J_hat^2))*((d_theta_T*d_theta0 + d_delta_T*d_delta0)^2)
denom4 <- 2*(n/J_hat)*d_delta_T*cov_term*(d_theta_T*d_theta0 + d_delta_T*d_delta0)

test_denom <- sqrt(denom1 + denom2 + denom3 + denom4)

test_stat <- test_num/test_denom
p_val <- pnorm(test_stat, lower.tail = FALSE)

sig_hat <- test_denom/(sqrt(m*n))
std_err <- sig_hat*sqrt((m+n)/(m*n))
ci_95 <- eta_hat + c(-1,1)*qnorm(0.975)*std_err

# hist(y, probability = TRUE, breaks = 50,
#      col = 'white', xlab = 'log(x)',
#      main = TeX('Estimated proposal bkg $g_{\\hat{beta}}(x)$ on bkg data'),
#      cex.main = 2,
#      cex.lab = 2)

# curve(qb, col = alpha('blue', 0.6), add = TRUE, lwd = 3,
#       lty = 1)

hist(x, probability = TRUE, breaks = 50,
     col = 'white', xlab = 'log(x)',
     main = TeX('$g_{\\beta}(x) \\propto \\exp\\left(
                -\\frac{(x+1)^2}{4\\beta}
                \\right)$'),
     cex.main = 1.5,
     cex.lab = 2)

curve(qb, col = alpha('blue', 0.6), add = TRUE, lwd = 3,
      lty = 1)
text(x = 2.5, y = 0.8, 
     TeX(sprintf(
       paste0('$\\hat{\\eta} = %.4f$, $p$-value: ', as.character(round(p_val, 6)), '; $\\sigma$-signif.: %.1f'),
       eta_hat, qnorm(p_val, lower.tail = FALSE))
     ),
     cex = 1.5)
text(x = 2.5, y = 0.65, 
     TeX(sprintf(
       paste0('95 \\%% CI for $\\hat{\\eta}$: [ %.4f, %.4f]'),
       ci_95[1], ci_95[2])
     ),
     cex = 1.5)
