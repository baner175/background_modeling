rm(list = ls())
source('bases on [l,u] - gb with 2 bumps.R')

n_rep <- 5e3 # do it for 1e4
n_samp <- 1e3
mean_back <- 0.5; sd_back <- 2.5; eta_true <- 0.03

theta_vec <- c()
t_stat_vec <- c()
set.seed(12345)
for(i in 1:n_rep){
  obs_new <- rtrunc(n_samp, spec = 'norm', a = l, b = u,
                    mean = mean_back, sd = sd_back)
  S1_vec <- sapply(obs_new, function(t) S1(t, mean = mean_sig,
                                           sd = sd_sig))
  theta <- mean(S1_vec)
  se_theta <- sqrt((mean(S1_vec^2) - theta^2)/n)
  theta_vec <- c(theta_vec, theta)
  t_stat_vec <- c(t_stat_vec, theta/se_theta)
  message(sprintf('Iteration: %d/%d',i,n_rep))
}

mean(t_stat_vec>qnorm(0.95)) #came out 0.0522

fb <- function(t) dtrunc(t, spec = 'norm', a = l, b = u,
                         mean = mean_back, sd = sd_back)
f_mix <- function(t) eta_true*fs(t) + (1-eta_true)*fb(t)

(delta <- integrate(function(t) fb(t)*S1(t), l, u)$value)
sd_theta <- (integrate(function(t) fb(t)*S1(t)^2, l, u)$value - delta^2)/sqrt(n_samp)
hist(theta_vec, probability = TRUE, breaks = 50)
true_dist <- function(t) dnorm(t, mean = delta, sd =sd_theta)
curve(true_dist, add = TRUE, col = 'blue')
abline(v = delta, col = 'red', lty = 2)
