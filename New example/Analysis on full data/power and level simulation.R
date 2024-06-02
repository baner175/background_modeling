rm(list = ls())
source('bases on [l,u] - gb with 2 bumps.R')
library(foreach)
library(doSNOW)
library(iterators)
library(parallel)

n_rep <- 1e4 # do it for 1e4
n_samp <- 5e3

# UNKNOWN DENSITIES:
mean_back <- 0.5; sd_back <- 2.5; eta_true <- 0.03
fb <- function(t) dtrunc(t, spec = 'norm', a = l, b = u,
                         mean = mean_back, sd = sd_back)
f_mix <- function(t) eta_true*fs(t) + (1-eta_true)*fb(t)


eta_seq <- seq(0, 0.09, 0.01)
theta_vec <- c()
t_stat_vec <- c()
seed <- 12345
set.seed(seed)

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

print(time_elapsed) #took 17 minutes

# mean(result_eta[[1]][,1]>qnorm(0.95)) # 0.0077
# mean(result_eta[[2]][,1]>qnorm(0.95)) # 0.5904
# mean(result_eta[[3]][,1]>qnorm(0.95)) # 0.9939

(filename <- paste0('nrep_',n_rep, '__nsamp_', n_samp, '__seed_',seed,'.csv'))

levels <- sapply(result_eta, function(x) mean(x[,1] > qnorm(0.95)))
res_df <- data.frame(eta = eta_seq, levels = levels)
write.csv(res_df, file = file.path('simulation_results/', filename), row.names = FALSE)
message('simulation completed...')

# theta_vec_null <- result_eta[[1]][,2]
# (delta <- integrate(function(t) fb(t)*S1(t), l, u)$value)
# sd_theta <- (integrate(function(t) fb(t)*S1(t)^2, l, u)$value - delta^2)/sqrt(n_samp)
# hist(theta_vec_null, probability = TRUE, breaks = 50)
# true_dist <- function(t) dnorm(t, mean = delta, sd =sd_theta)
# curve(true_dist, add = TRUE, col = 'blue')
# abline(v = delta, col = 'red', lty = 2)
