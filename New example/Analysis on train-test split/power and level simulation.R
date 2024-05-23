rm(list = ls())
source('bases on [l,u] - gb with 2 bumps.R')
library(doSNOW)
library(foreach)
library(iterators)

n_rep <- 1e4 # do it for 1e4
n_samp <- 1e3

# UNKNOWN DENSITIES:
mean_back <- 0.5; sd_back <- 2.5; eta_true <- 0.03
fb <- function(t) dtrunc(t, spec = 'norm', a = l, b = u,
                         mean = mean_back, sd = sd_back)
f_mix <- function(t) eta_true*fs(t) + (1-eta_true)*fb(t)


eta_seq <- seq(0, 0.06, 0.03)
theta_vec <- c()
t_stat_vec <- c()
set.seed(12345)

numCores <- detectCores()
cl <- makeSOCKcluster(numCores)
registerDoSNOW(cl)
pb <- txtProgressBar(max=n_rep*length(eta_seq), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- options(progress = progress)

time_elapsed <- system.time({
result_eta <- foreach(j = icount(length(eta_seq)),
                      .packages = c('truncdist', 'tcltk')) %:%
          foreach(i = icount(n_rep), .combine = rbind,
                      .packages = c('truncdist', 'tcltk'),
                      .options.snow = opts) %dopar%
      {
        sig <- rtrunc(n_samp, spec = 'norm', a = l, b = u,
                      mean = mean_sig, sd = sd_sig)
        back <- rtrunc(n_samp, spec = 'norm', a = l, b = u,
                       mean = mean_back, sd = sd_back)
        u_mask <- runif(n_samp)
        obs_new <- ifelse(u_mask<eta_seq[j], sig, back)
        
        S1_vec <- sapply(obs_new, function(t) S1(t, mean = mean_sig,
                                                 sd = sd_sig))
        theta <- mean(S1_vec)
        se_theta <- sqrt((mean(S1_vec^2) - theta^2)/n_samp)
        t_stat <- theta/se_theta
        c(t_stat, theta, se_theta)
      }
})
stopCluster(cl)

print(time_elapsed) #took 23 minutes

mean(result_eta[[1]][,1]>qnorm(0.95)) # 0.0078
mean(result_eta[[2]][,1]>qnorm(0.95)) # 0.5896
mean(result_eta[[3]][,1]>qnorm(0.95)) # 0.9947

theta_vec_null <- result_eta[[1]][,2]
(delta <- integrate(function(t) fb(t)*S1(t), l, u)$value)
sd_theta <- (integrate(function(t) fb(t)*S1(t)^2, l, u)$value - delta^2)/sqrt(n_samp)
hist(theta_vec_null, probability = TRUE, breaks = 50)
true_dist <- function(t) dnorm(t, mean = delta, sd =sd_theta)
curve(true_dist, add = TRUE, col = 'blue')
abline(v = delta, col = 'red', lty = 2)
