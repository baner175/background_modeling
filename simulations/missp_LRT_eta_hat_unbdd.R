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

neg_KL_div <- function(eta){
  -integrate(function(t) log(eta*dtrunc(t, spec = 'norm', a = l, b = u,
                                        mean = mean_sig, sd = sd_sig)+
                               (1-eta)*dtrunc(t, spec = 'pareto', a = l, b = u,
                                              scale = l, shape = beta0))*
               (eta_true*dtrunc(t, spec = 'norm', a = l, b = u,
                                mean = mean_sig, sd = sd_sig) + 
                  (1-eta_true)*dtrunc(t, a = l, b = u, spec = 'gamma',
                                      rate = bkg_rate, shape = bkg_shape)),
             l, u)$value
}

eta_star <- nlminb(start = 0.01,
                   objective = neg_KL_div,
                   lower = 0, upper = 1)$par

J1 <- integrate(function(t){
  val <- (fs(t) - gb(t))/(eta_star*fs(t)+(1-eta_star)*gb(t))
  return(val^2)
}, l, u)$value

neg_ll <- function(eta, data){
  fi <- sapply(data, function(t){
    eta*dtrunc(t, spec = 'norm', a = l, b = u,
               mean = mean_sig, sd = sd_sig) + 
      (1-eta)*dtrunc(t, spec = 'pareto', a = l, b = u,
                     scale = l, shape = beta0)
  })
  return(-sum(log(fi)))
}

B <- 1e4
n_samp <- 1e3
set.seed(12345)
seeds <- sample.int(.Machine$integer.max, B)
cl <- makeCluster(8)
registerDoSNOW(cl)
pb <- txtProgressBar(max = B, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

start_time <- Sys.time()
test_stat_LRT <- foreach(i = 1:B, .combine = rbind,
                         .packages = c('truncdist', 'VGAM'),
                         .options.snow = opts) %dopar%
  {
    set.seed(seeds[i])
    
    # physics-sample:
    s_samp <- rtrunc(n_samp, a = l, b = u, spec = 'norm',
                     mean = mean_sig, sd = sd_sig)
    b_samp <- rtrunc(n_samp, a = l, b = u, spec = 'gamma',
                     rate = bkg_rate, shape = bkg_shape)
    u_mask <- runif(n_samp)
    phys_samp <- ifelse(u_mask <= eta_true, s_samp, b_samp)
    
    ll_res <- nlminb(start = 0.01,
                     objective = neg_ll,
                     lower = -Inf, upper = Inf,
                     data = phys_samp)
    eta_hat <- ll_res$par
    ll1 <- -ll_res$objective
    ll_star <- -neg_ll(eta = eta_star, data = phys_samp)
    ll_0 <- -neg_ll(eta = 0, data = phys_samp)
    test_stat_star <- -2*(ll_star-ll1)
    test_stat_0 <- -2*(ll_0-ll1)
    
    c(test_stat_star, test_stat_0)
  }
close(pb)

Sys.time() - start_time
stopCluster(cl)

df <- data.frame('test_stat_star' = test_stat_LRT[,1],
                 'test_stat_0' = test_stat_LRT[,2])

if(eta_true > 0)
{
  df <- df[,1]
}

file_name <- paste0('Results/LRT', 
                    '_B_', B, 
                    '_n_samp_', n_samp,
                    '_eta_', eta_true,
                    '.csv')
write.csv(df, file_name,
          row.names = FALSE)