rm(list = ls())

library(truncdist)
library(VGAM)
library(doSNOW)
library(parallel)
library(foreach)
l <- 1; u <- 2
################################################################
################ SIGNAL AND SIGNAL REGION ######################
################################################################

#parameters for the signal
mean_sig <- 1.28
sd_sig <- 0.02

eta_true <- 0

# signal density
fs <- function(x, mean = mean_sig) dtrunc(x, a = l, b = u, spec = 'norm', mean = mean, sd = sd_sig)

#parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5

# true bkg density
fb_true <- function(x) dtrunc(x, a = l, b = u, spec = 'gamma',
                              rate = bkg_rate, shape = bkg_shape)
# true mixture density
f <- function(x) eta_true*fs(x)+(1-eta_true)*fb_true(x)

beta0 <- 3.87
gb <- function(x, beta = beta0){
  dtrunc(x, spec = 'pareto', a = l, b = u,
         scale = l, shape = beta)
}

S <- function(x, beta = beta0) fs(x)/gb(x, beta) - 1
norm_S <- function(beta = beta0) integrate(function(x) (S(x, beta)^2)*gb(x,beta),
                                   l, u)$value |> sqrt()
S2 <- function(x, beta = beta0) (fs(x)/gb(x, beta) - 1)/norm_S(beta)^2

B <- 1e4
n_samp <- 5e2
set.seed(12345)
seeds <- sample.int(.Machine$integer.max, B)

cl <- makeCluster(8)
registerDoSNOW(cl)
pb <- txtProgressBar(max = B, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

start_time <- Sys.time()
test_stat_eta <- foreach(i = 1:B, .combine = c,
                         .packages = c('truncdist', 'VGAM'),
                         .options.snow = opts) %dopar%
  {
    set.seed(seeds[i])
    # bkg-only sample:
    bkg_samp <- rtrunc(n_samp, a = l, b = u, spec = 'gamma',
                       rate = bkg_rate, shape = bkg_shape)
    
    # physics-sample:
    s_samp <- rtrunc(n_samp, a = l, b = u, spec = 'norm',
                     mean = mean_sig, sd = sd_sig)
    b_samp <- rtrunc(n_samp, a = l, b = u, spec = 'gamma',
                     rate = bkg_rate, shape = bkg_shape)
    u_mask <- runif(n_samp)
    phys_samp <- ifelse(u_mask <= eta_true, s_samp, b_samp)
    
    
    S2_phys_vec <- sapply(phys_samp, 
                         function(x) S2(x))
    S2_bkg_vec <- sapply(bkg_samp, 
                        function(x) S2(x))
    theta_0_hat <- mean(S2_phys_vec)
    delta_0_hat <- mean(S2_bkg_vec)
    
    var_F_S <- mean(S2_phys_vec^2) - theta_0_hat^2
    var_Fb_S <- mean(S2_bkg_vec^2) - delta_0_hat^2
    
    eta_hat <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat)
    sigma_eta_hat <- sqrt(
      (var_F_S*(1-delta_0_hat)^2 + var_Fb_S*(1-theta_0_hat)^2)/(1-delta_0_hat)^4
    )
    
    sqrt(n_samp)*(eta_hat-eta_true)/sigma_eta_hat
  }
close(pb)

Sys.time() - start_time
stopCluster(cl)

file_name <- paste0('Results/beta_',beta0,'_known','_B_',B,'_n_samp_',n_samp,'_eta_',eta_true,'.csv')

write.csv(data.frame(test_stat = test_stat_eta),
          file_name, row.names = FALSE)