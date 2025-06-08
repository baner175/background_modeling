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

gb <- function(x, beta){
  dtrunc(x, spec = 'pareto', a = l, b = u,
         scale = l, shape = beta)
}

gb_bkg_model <- function(beta, data){
  gb_i <- sapply(data, function(x){
    dtrunc(x, spec = 'pareto', a = l, b = u,
           scale = l, shape = beta)
  })
  return(-sum(log(gb_i)))
}

d_log_gb <- function(beta, x){
  1/beta - log(x) - (log(u)*u^(-beta) - log(l)*l^(-beta))/(l^(-beta) - u^(-beta))
}
d2_log_gb <- function(beta){
  -1/(beta^2) - 
    ((log(l)^2)*l^(-beta) - (log(u)^2)*u^(-beta))/(l^(-beta) - u^(-beta)) -
    ((log(u)*u^(-beta) - log(l)*l^(-beta))/(l^(-beta) - u^(-beta)))^2
}

S <- function(x, beta) fs(x)/gb(x, beta) - 1
norm_S <- function(beta) integrate(function(x) (S(x, beta)^2)*gb(x,beta),
                                   l, u)$value |> sqrt()
S2 <- function(x, beta) (fs(x)/gb(x, beta) - 1)/norm_S(beta)^2

d_S2 <- function(x, beta){
  nrm_S <- norm_S(beta)
  d_S <- (-fs(x)/gb(x, beta))*d_log_gb(beta, x)
  d_normS2 <- -integrate(function(t) (fs(t)^2/gb(t,beta))*d_log_gb(t,beta),
                         l, u)$value
  S_val <- fs(x)/gb(x, beta) - 1
  
  return(((nrm_S^2)*d_S - S_val*d_normS2)/(nrm_S^4))
}

bkg_to_phys_ratio <- 2
B <- 1e4
n_phys <- 5e3
n_bkg <- n_phys*bkg_to_phys_ratio
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
    bkg_samp <- rtrunc(n_bkg, a = l, b = u, spec = 'gamma',
                       rate = bkg_rate, shape = bkg_shape)
    
    # physics-sample:
    s_samp <- rtrunc(n_phys, a = l, b = u, spec = 'norm',
                     mean = mean_sig, sd = sd_sig)
    b_samp <- rtrunc(n_phys, a = l, b = u, spec = 'gamma',
                     rate = bkg_rate, shape = bkg_shape)
    u_mask <- runif(n_phys)
    phys_samp <- ifelse(u_mask <= eta_true, s_samp, b_samp)
    
    opt <- tryCatch(
      nlminb(start = 0.01,
             objective = gb_bkg_model,
             lower = 0, upper = Inf,
             data = bkg_samp),
      error = function(e) return(NULL)
    )
    
    if (is.null(opt)) return(NA_real_)  # skip this iteration
    
    beta_hat <- opt$par
    
    S2_phys_vec <- sapply(phys_samp, 
                          function(x) S2(x, beta = beta_hat))
    S2_bkg_vec <- sapply(bkg_samp, 
                         function(x) S2(x, beta = beta_hat))
    theta_0_hat <- mean(S2_phys_vec)
    delta_0_hat <- mean(S2_bkg_vec)
    J1_hat <- -d2_log_gb(beta_hat)
    V1_hat <- sapply(bkg_samp,
                     function(x) d_log_gb(beta = beta_hat, x)^2) |> mean()
    
    eta_hat <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat)
    
    d_theta_hat <- sapply(phys_samp, function(x) d_S2(x, beta = beta_hat)) |> mean()
    d_delta_hat <- sapply(bkg_samp, function(x) d_S2(x, beta = beta_hat)) |> mean()
    
    test_num <- sqrt(n_phys*n_bkg)*(eta_hat - eta_true)
    test_denom <- sqrt(
      (n_bkg/(1-delta_0_hat)^2) * (mean(S2_phys_vec^2) - theta_0_hat^2) +
        n_phys*((theta_0_hat - 1)^2/(1-delta_0_hat)^4) * (mean(S2_bkg_vec^2) - delta_0_hat^2) + 
        n_phys*V1_hat/(J1_hat^2) * (d_theta_hat/(1-delta_0_hat) + ((theta_0_hat-1)/((1-delta_0_hat)^2))*d_delta_hat)^2 + 
        2 * (n_phys/J1_hat) * (theta_0_hat-1)/((1-delta_0_hat)^2) * 
        (d_theta_hat/(1-delta_0_hat) + ((theta_0_hat-1)/((1-delta_0_hat)^2))*d_delta_hat) * 
        sapply(bkg_samp, function(x) S2(x, beta = beta_hat)*d_log_gb(beta_hat, x)) |> mean()
    )
    
    test_num/test_denom
  }
close(pb)

Sys.time() - start_time
stopCluster(cl)

file_name <- paste0('Results/unbinned_test_eta_w_bkg__',
                    'beta_estimated_',
                    'B(',B,')_',
                    'n_phys(',n_phys,')_',
                    'bkg_to_phys(',bkg_to_phys_ratio,')_',
                    'eta(',eta_true,')','.csv')

write.csv(data.frame(test_stat = test_stat_eta),
          file_name, row.names = FALSE)
