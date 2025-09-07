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
  qb_i <- dtrunc(y, spec = 'exp', rate = beta, a = l, b = u)
  return(-sum(log(qb_i)))
}

beta_hat <- nlminb(start = 0.01,
                   objective = qb_y_model,
                   upper = Inf, lower = -Inf)$par

norm_S <- integrate(function(t){
  fs <- dtrunc(exp(t), spec = 'norm', a = real_l, b = real_u,
               mean = mean_sig, sd = sd_sig)*exp(t)
  qb <- dtrunc(t, spec = 'exp', rate = beta_hat, a = l, b = u)
  return(((fs/qb - 1)^2)*qb)
}, l, u)$value |> sqrt()

S2_phys_vec <- sapply(x, function(t){
  fs <- dtrunc(exp(t), spec = 'norm', a = real_l, b = real_u,
               mean = mean_sig, sd = sd_sig)*exp(t)
  qb <- dtrunc(t, spec = 'exp', rate = beta_hat, a = l, b = u)
  S_val <- (fs/qb - 1)
  return(S_val/(norm_S^2))
})

S2_bkg_vec <- sapply(y, function(t){
  fs <- dtrunc(exp(t), spec = 'norm', a = real_l, b = real_u,
               mean = mean_sig, sd = sd_sig)*exp(t)
  qb <- dtrunc(t, spec = 'exp', rate = beta_hat, a = l, b = u)
  S_val <- (fs/qb - 1)
  return(S_val/(norm_S^2))
})

theta0_hat <- mean(S2_phys_vec); delta0_hat <- mean(S2_bkg_vec)
eta_hat <- (theta0_hat-delta0_hat)/(1-delta0_hat)

test_num <- sqrt(m*n)*eta_hat

d_normS_sq <- -integrate(function(t){
  fs <- dtrunc(exp(t), spec = 'norm', a = real_l, b = real_u,
               mean = mean_sig, sd = sd_sig)*exp(t)
  qb <- dtrunc(t, spec = 'exp', rate = beta_hat, a = l, b = u)
  d_log_qb <- (1/beta_hat) - t - 
    (u*exp(-beta_hat*u) - l*exp(-beta_hat*l))/(exp(-beta_hat*l) - exp(-beta_hat*u))
  return((fs^2)*d_log_qb/qb)
}, l, u)$value

d_S2 <- function(t){
    fs <- dtrunc(exp(t), spec = 'norm', a = real_l, b = real_u,
                 mean = mean_sig, sd = sd_sig)*exp(t)
    qb <- dtrunc(t, spec = 'exp', rate = beta_hat, a = l, b = u)
    d_log_qb <- (1/beta_hat) - t - 
      (u*exp(-beta_hat*u) - l*exp(-beta_hat*l))/(exp(-beta_hat*l) - exp(-beta_hat*u))
    
    return(-((norm_S^2)*(fs/qb)*qb*d_log_qb + (fs/qb-1)*d_normS_sq)/(norm_S^4))
}

d_log_qb_yi <- sapply(y, function(t){
  return(
    (1/beta_hat) - t - 
      (u*exp(-beta_hat*u) - l*exp(-beta_hat*l))/(exp(-beta_hat*l) - exp(-beta_hat*u)) 
  )
})

d2_log_qb_y <- -(1/(beta_hat^2)) - (((l^2)*exp(-beta_hat*l) - (u^2)*exp(-beta_hat*u))/(exp(-beta_hat*l) - exp(-beta_hat*u)) - 
                                      ((u*exp(-beta_hat*u) - l*exp(-beta_hat*l))/(exp(-beta_hat*l) - exp(-beta_hat*u)))^2)


d_theta0 <- sapply(x, d_S2) |> mean()
d_delta0 <- sapply(y, d_S2) |> mean()
d_theta_T <- 1/(1-delta0_hat)
d_delta_T <- (theta0_hat - 1)/((1-delta0_hat)^2)
cov_term <- mean(d_log_qb_yi*S2_bkg_vec)
V_hat <- mean((d_log_qb_yi)^2)
J_hat <- -d2_log_qb_y

var_S2_F_hat <- mean(S2_phys_vec^2) - (theta0_hat^2)
var_S2_Fb_hat <- mean(S2_bkg_vec^2) - (delta0_hat^2)

denom1 <- m*(d_theta_T^2)*var_S2_F_hat
denom2 <- n*(d_delta_T^2)*var_S2_Fb_hat
denom3 <- n*(V_hat/(J_hat^2))*((d_theta_T*d_theta0 + d_delta_T*d_delta0)^2)
denom4 <- 2*(n/J_hat)*d_delta_T*cov_term*(d_theta_T*d_theta0 + d_delta_T*d_delta0)

test_denom <- sqrt(denom1 + denom2 + denom3 + denom4)

test_stat <- test_num/test_denom
p_val <- pnorm(test_stat, lower.tail = FALSE)

qb <- function(y) dtrunc(y, spec = 'exp', rate = beta_hat, a = l, b = u)

hist(y, probability = TRUE, breaks = 50,
     col = 'white', xlab = 'log(x)',
     main = TeX('Estimated proposal bkg $q_b(x; \\hat{beta})$ on bkg data'))

curve(qb, col = alpha('blue', 0.6), add = TRUE, lwd = 2.2,
      lty = 1)

hist(x, probability = TRUE, breaks = 50,
     col = 'white', xlab = 'log(x)',
     main = TeX('Estimated proposal bkg $q_b(x; \\hat{beta})$ on physics data'))

curve(qb, col = alpha('blue', 0.6), add = TRUE, lwd = 2.2,
      lty = 1)
text(x = 2.5, y = 0.8, 
     TeX(sprintf(
       paste0('$\\hat{\\eta} = %.4f$, $p$-value: ', as.character(round(p_val, 6)), '; $\\sigma$-signif.: %.1f'),
       eta_hat, qnorm(p_val, lower.tail = FALSE))
       ))
